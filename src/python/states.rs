use anyhow::Result;
use crate::numerics::float::{Float, Float3};
use crate::transport::PhotonState;
use pyo3::prelude::*;
use pyo3::exceptions::{PyKeyError, PyTypeError, PyValueError};
use pyo3::types::{PyDict, PyType};
use self::ErrorKind::{KeyError, TypeError, ValueError};
use super::numpy::{Dtype, PyArray, ShapeArg};


// ===============================================================================================
// C representation of a photon state.
// ===============================================================================================
#[repr(C)]
#[derive(Clone, Copy)]
pub(crate) struct CState {
    pub energy: Float,
    pub position: [Float; 3],
    pub direction: [Float; 3],
    pub length: Float,
    pub weight: Float,
}

impl From<CState> for PhotonState {
    fn from(state: CState) -> Self {
        Self {
            energy: state.energy,
            position: state.position.into(),
            direction: state.direction.into(),
            length: state.length,
            weight: state.weight,
        }
    }
}

impl From<PhotonState> for CState {
    fn from(state: PhotonState) -> Self {
        Self {
            energy: state.energy,
            position: state.position.into(),
            direction: state.direction.into(),
            length: state.length,
            weight: state.weight
        }
    }
}


// ===============================================================================================
// Utility function for creating a numpy array of photon states.
// ===============================================================================================

#[pyfunction]
#[pyo3(signature=(shape=None, **kwargs))]
pub fn states(
    py: Python,
    shape: Option<ShapeArg>,
    kwargs: Option<&Bound<PyDict>>
    ) -> Result<PyObject> {
    let shape: Vec<usize> = match shape {
        None => vec![0],
        Some(shape) => shape.into(),
    };
    let array = PyArray::<CState>::zeros(py, &shape)?
        .into_any();
    let mut has_direction = false;
    let mut has_energy = false;
    let mut has_weight = false;
    if let Some(kwargs) = kwargs {
        for (key, value) in kwargs.iter() {
            {
                let key: String = key.extract()?;
                match key.as_str() {
                    "direction" => { has_direction = true; },
                    "energy" => { has_energy = true; },
                    "weight" => { has_weight = true; },
                    _ => {},
                }
            }
            array.set_item(key, value)?;
        }
    }
    if !has_direction {
        array.set_item("direction", (0.0, 0.0, 1.0))?;
    }
    if !has_energy {
        array.set_item("energy", 1.0)?;
    }
    if !has_weight {
        array.set_item("weight", 1.0)?;
    }
    Ok(array.unbind())
}


// ===============================================================================================
// Utility struct for managing an array of photon states.
// ===============================================================================================

pub struct States<'a> {
    energy: &'a PyArray<Float>,
    position: &'a PyArray<Float>,
    direction: &'a PyArray<Float>,
    length: Option<&'a PyArray<Float>>,
    pid: Option<&'a PyArray<i32>>,
    weight: Option<&'a PyArray<Float>>,
    size: usize,
}

impl<'a> States<'a> {
    pub fn new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Self> {
        let energy = extract(elements, "energy")?;
        let position = extract(elements, "position")?;
        let direction = extract(elements, "direction")?;
        let length = maybe_extract(elements, "length")?;
        let pid = maybe_extract(elements, "pid")?;
        let weight = maybe_extract(elements, "weight")?;
        if *position.shape().last().unwrap_or(&0) != 3 {
            let why = format!("expected a shape '[..,3]' array, found '{:?}'", position.shape());
            let err = Error::new(TypeError, "position", why).into_err();
            return Err(err);
        }
        if *direction.shape().last().unwrap_or(&0) != 3 {
            let why = format!("expected a shape '[..,3]' array, found '{:?}'", direction.shape());
            let err = Error::new(TypeError, "direction", why).into_err();
            return Err(err);
        }
        let size = energy.size();
        let others = [
            position.size() / 3,
            direction.size() / 3,
            length.map(|a| a.size()).unwrap_or(size),
            pid.map(|a| a.size()).unwrap_or(size),
            weight.map(|a| a.size()).unwrap_or(size),
        ];
        if others.iter().any(|x| *x != size) {
            let err = Error::new(TypeError, "states", "differing arrays sizes").into_err();
            return Err(err);
        }
        let states = Self { energy, position, direction, length, pid, weight, size };
        Ok(states)
    }

    pub fn get(&self, index: usize) -> PyResult<PhotonState> {
        let length = match self.length {
            None => 0.0,
            Some(length) => length.get(index)?,
        };
        if let Some(pid) = self.pid {
            let pid = pid.get(index)?;
            if pid != 22 {
                let why = format!("expected '22', found '{}'", pid);
                let err = Error::new(ValueError, "pid", why).into_err();
                return Err(err)
            }
        };
        let weight = match self.weight {
            None => 1.0,
            Some(weight) => weight.get(index)?,
        };
        let state = PhotonState {
            energy: self.energy.get(index)?,
            position: [
                self.position.get(3 * index)?,
                self.position.get(3 * index + 1)?,
                self.position.get(3 * index + 2)?,
            ].into(),
            direction: [
                self.direction.get(3 * index)?,
                self.direction.get(3 * index + 1)?,
                self.direction.get(3 * index + 2)?,
            ].into(),
            length,
            weight,
        };
        Ok(state)
    }

    pub fn has_weights(&self) -> bool {
        self.weight.is_some()
    }

    pub fn set(&self, index: usize, state: &PhotonState) -> PyResult<()> {
        self.energy.set(index, state.energy)?;
        self.position.set(3 * index, state.position.0)?;
        self.position.set(3 * index + 1, state.position.1)?;
        self.position.set(3 * index + 2, state.position.2)?;
        self.direction.set(3 * index, state.direction.0)?;
        self.direction.set(3 * index + 1, state.direction.1)?;
        self.direction.set(3 * index + 2, state.direction.2)?;
        if let Some(length) = self.length {
            length.set(index, state.length)?;
        }
        if let Some(weight) = self.weight {
            weight.set(index, state.weight)?;
        }
        Ok(())
    }

    pub fn shape(&self) -> Vec<usize> {
        self.energy.shape()
    }

    pub fn size(&self) -> usize {
        self.size
    }
}

fn extract<'a, 'py, T>(elements: &'a Bound<'py, PyAny>, key: &str) -> PyResult<&'a PyArray<T>>
where
    'py: 'a,
    T: Dtype,
{
    let py = elements.py();
    elements.get_item(key)
        .map_err(|err| Error::new(KeyError, "states", err.value_bound(py).to_string()).into_err())?
        .extract()
        .map_err(|err| Error::new(TypeError, key, err.value_bound(py).to_string()).into_err())
}

fn maybe_extract<'a, 'py, T>(
    elements: &'a Bound<'py, PyAny>,
    key: &str
) -> PyResult<Option<&'a PyArray<T>>>
where
    'py: 'a,
    T: Dtype,
{
    let py = elements.py();
    match extract(elements, key) {
        Ok(property) => Ok(Some(property)),
        Err(err) => if err.get_type_bound(py).is(&PyType::new_bound::<PyKeyError>(py)) {
            Ok(None)
        } else {
            Err(err)
        },
    }
}

struct Error<'a> {
    kind: ErrorKind,
    what: &'a str,
    why: String
}

enum ErrorKind {
    KeyError,
    TypeError,
    ValueError,
}

impl<'a> Error<'a> {
    fn new<T: Into<String>>(kind: ErrorKind, what: &'a str, why: T) -> Self {
        Self { kind, what, why: why.into() }
    }

    fn into_err(self) -> PyErr {
        let msg = format!("bad {} ({})", self.what, self.why);
        match self.kind {
            ErrorKind::KeyError => PyKeyError::new_err(msg),
            ErrorKind::TypeError => PyTypeError::new_err(msg),
            ErrorKind::ValueError => PyValueError::new_err(msg),
        }
    }
}


// ===============================================================================================
// Utility struct for managing an array of coordinates.
// ===============================================================================================

pub struct Coordinates<'a> {
    position: &'a PyArray<Float>,
    direction: Option<&'a PyArray<Float>>,
    size: usize,
}

impl<'a> Coordinates<'a> {
    pub fn new_with_direction<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Self> {
        let coordinates = Self::new(elements)?;
        if coordinates.direction.is_none() {
            let err = Error::new(
                TypeError,
                "position",
                "missing 'direction'".to_owned()
            ).into_err();
            return Err(err);
        }
        Ok(coordinates)
    }

    pub fn new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Self> {
        let position = extract(elements, "position")?;
        let direction = maybe_extract(elements, "direction")?;
        if *position.shape().last().unwrap_or(&0) != 3 {
            let why = format!("expected a shape '[..,3]' array, found '{:?}'", position.shape());
            let err = Error::new(TypeError, "position", why).into_err();
            return Err(err);
        }
        if let Some(direction) = direction {
            if *direction.shape().last().unwrap_or(&0) != 3 {
                let why = format!(
                    "expected a shape '[..,3]' array, found '{:?}'",
                    direction.shape()
                );
                let err = Error::new(TypeError, "direction", why).into_err();
                return Err(err);
            }
        }
        let size = position.size() / 3;
        if let Some(direction) = direction {
            if position.size() != direction.size() {
                let err = Error::new(TypeError, "states", "differing arrays sizes").into_err();
                return Err(err);
            }
        }
        let coordinates = Self { position, direction, size };
        Ok(coordinates)
    }

    pub fn get_direction(&self, index: usize) -> PyResult<Float3> {
        let direction = match self.direction.as_ref() {
            None => Float3::new(0.0, 0.0, 1.0),
            Some(direction) => Float3::new(
                direction.get(3 * index)?,
                direction.get(3 * index + 1)?,
                direction.get(3 * index + 2)?,
            ),
        };
        Ok(direction)
    }

    pub fn get_position(&self, index: usize) -> PyResult<Float3> {
        let position = Float3::new(
                self.position.get(3 * index)?,
                self.position.get(3 * index + 1)?,
                self.position.get(3 * index + 2)?,
            );
        Ok(position)
    }

    pub fn shape(&self) -> Vec<usize> {
        let mut shape = self.position.shape();
        shape.pop();
        shape
    }

    pub fn size(&self) -> usize {
        self.size
    }
}
