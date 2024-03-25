use anyhow::{anyhow, Result};
use crate::numerics::consts::PI;
use crate::numerics::{Float, Float3, Float3x3, FloatRng};
use crate::pretty_enumerate;
use crate::transport::{BoxShape, GeometryShape, PhotonState, SphereShape, TransportBoundary,
                       TransportMode};
use enum_iterator::{all, Sequence};
use pyo3::prelude::*;
use std::fmt;
use super::numpy::{FloatOrFloat3, PyArray};
use super::rand::PyRandomStream;
use super::transport::{CState, PyTransportEngine};


// ===============================================================================================
// Unresolved boundary.
// ===============================================================================================

#[derive(FromPyObject)]
pub enum PyTransportBoundary<'py> {
    Box(PyRef<'py, PyBoxShape>),
    Sector(usize),
    Sphere(PyRef<'py, PySphereShape>),
}

impl<'py> From<PyTransportBoundary<'py>> for TransportBoundary {
    fn from(boundary: PyTransportBoundary) -> Self {
        match boundary {
            PyTransportBoundary::Box(shape) => Self::Box(shape.0),
            PyTransportBoundary::Sector(index) => Self::Sector(index),
            PyTransportBoundary::Sphere(shape) => Self::Sphere(shape.0),
        }
    }
}

impl IntoPy<PyObject> for TransportBoundary {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            Self::None => py.None(),
            Self::Box(shape) => {
                let shape = PyBoxShape(shape);
                shape.into_py(py)
            },
            Self::Sector(index) => index.into_py(py),
            Self::Sphere(shape) => {
                let shape = PySphereShape(shape);
                shape.into_py(py)
            },
        }
    }
}


// ===============================================================================================
// Python wrapper for a box shape.
// ===============================================================================================

#[pyclass(name = "BoxShape", module = "goupil")]
pub struct PyBoxShape(BoxShape);

#[pymethods]
impl PyBoxShape {
    #[new]
    fn new(
        center: Option<Float3>,
        size: Option<FloatOrFloat3>,
        rotation: Option<Float3x3>
    ) -> Self {
        let center = center.unwrap_or(Float3::zero());
        let size = match size {
            None => Float3::splat(1.0),
            Some(size) => match size {
                FloatOrFloat3::Float(size) => Float3::splat(size),
                FloatOrFloat3::Float3(size) => size,
            }
        };
        let shape = BoxShape { center, size, rotation };
        Self(shape)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    #[getter]
    fn get_center(&self) -> Float3 {
        self.0.center
    }

    #[getter]
    fn get_rotation(&self) -> Option<Float3x3> {
        self.0.rotation
    }

    #[getter]
    fn get_size(&self) -> Float3 {
        self.0.size
    }

    fn distance(&self, states: &PyArray<CState>) -> Result<PyObject> {
        self.0.distance_v(states)
    }

    fn inside(&self, states: &PyArray<CState>) -> Result<PyObject> {
        self.0.inside_v(states)
    }

    fn randomise(
        &self,
        states: &PyArray<CState>,
        engine: Option<&PyTransportEngine>,
        rng: Option<Py<PyRandomStream>>,
        side: Option<&str>,
        direction: Option<&str>,
        weight: Option<bool>,
    ) -> Result<()> {
        let generator = BoxGenerator::new(&self.0);
        generator.generate_on_v(states, engine, rng, side, direction, weight)
    }
}

struct BoxGenerator<'a> {
    shape: &'a BoxShape,
    cumulative_surface: [Float; 6],
}

impl<'a> BoxGenerator<'a> {
    fn new(shape: &'a BoxShape) -> Self {
        let cumulative_surface = {
            let mut tmp = [
                shape.size.1 * shape.size.2,
                shape.size.1 * shape.size.2,
                shape.size.2 * shape.size.0,
                shape.size.2 * shape.size.0,
                shape.size.0 * shape.size.1,
                shape.size.0 * shape.size.1,
            ];
            for i in 1..6 { tmp[i] += tmp[i - 1] };
            tmp
        };
        Self { shape, cumulative_surface }
    }
}

impl<'a> Generate for BoxGenerator<'a> {
    fn generate_on<R: FloatRng>(
        &self,
        rng: &mut R,
        side: Side,
        direction: Direction,
        weight: bool,
        state: &mut PhotonState
    ) {
        // Randomise the face according to respective surfaces.
        let (axis, dir, total_surface) = {
            let n = self.cumulative_surface.len();
            let r = self.cumulative_surface[n - 1] * rng.uniform01();
            let face = {
                let mut i = 0;
                loop {
                    if r <= self.cumulative_surface[i] { break i }
                    else if i == n - 2 { break n - 1 }
                    else { i += 1 }
                }
            };
            let axis = face / 2;
            let dir = if (face % 2) == 0 { -1.0 } else { 1.0 };
            (axis, dir, self.cumulative_surface[n - 1])
        };

        // Randomise the position.
        let position: Float3 = {
            let mut position: [Float; 3] = [ 0.0, 0.0, 0.0];
            let size: [Float; 3] = self.shape.size.into();
            let eps = 1E-04; // For numeric safety.
            let sgn = match side {
                Side::Inside => -1.0,
                Side::Outside => 1.0,
            };
            position[axis] = dir * ((0.5 * size[axis] + eps * sgn).max(0.0));
            for i in 0..2 {
                let ii = (axis + 1 + i) % 3;
                position[ii] = size[ii] * (rng.uniform01() - 0.5);
            }
            position.into()
        };
        match self.shape.rotation.as_ref() {
            None => { state.position = position + self.shape.center },
            Some(rotation) => { state.position = rotation * position + self.shape.center },
        }

        // Randomise the direction using a cosine distribution.
        let xi = rng.uniform01();
        let cos_theta = xi.sqrt();
        let sin_theta = (1.0 - xi).sqrt();
        let phi = 2.0 * PI * rng.uniform01();

        let sgn = match direction {
            Direction::Ingoing => -1.0,
            Direction::Outgoing => 1.0,
        };
        let direction: Float3 = {
            let mut u: [Float; 3] = [ 0.0, 0.0, 0.0 ];
            u[(axis + 1) % 3] = sgn * dir * sin_theta * phi.cos();
            u[(axis + 2) % 3] = sgn * dir * sin_theta * phi.sin();
            u[axis] = sgn * dir * cos_theta;
            u.into()
        };
        match self.shape.rotation.as_ref() {
            None => { state.direction = direction },
            Some(rotation) => { state.direction = rotation * direction },
        }

        // Weight according to 1 / PDF.
        if weight {
            state.weight *= total_surface * PI;
        }
    }
}

impl VectorisedOperations for BoxShape {}


// ===============================================================================================
// Python wrapper for a spherical shape.
// ===============================================================================================

#[pyclass(name = "SphereShape", module = "goupil")]
pub struct PySphereShape(SphereShape);

#[pymethods]
impl PySphereShape {
    #[new]
    fn new(center: Option<Float3>, radius: Option<Float>) -> Self {
        let center = center.unwrap_or(Float3::zero());
        let radius = radius.unwrap_or(1.0);
        let shape = SphereShape { center, radius };
        Self(shape)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    #[getter]
    fn get_center(&self) -> Float3 {
        self.0.center
    }

    #[getter]
    fn get_radius(&self) -> Float {
        self.0.radius
    }

    fn distance(&self, states: &PyArray<CState>) -> Result<PyObject> {
        self.0.distance_v(states)
    }

    fn inside(&self, states: &PyArray<CState>) -> Result<PyObject> {
        self.0.inside_v(states)
    }

    fn randomise(
        &self,
        states: &PyArray<CState>,
        engine: Option<&PyTransportEngine>,
        rng: Option<Py<PyRandomStream>>,
        side: Option<&str>,
        direction: Option<&str>,
        weight: Option<bool>,
    ) -> Result<()> {
        self.0.generate_on_v(states, engine, rng, side, direction, weight)
    }
}

impl Generate for SphereShape {
    fn generate_on<R: FloatRng>(
        &self,
        rng: &mut R,
        side: Side,
        direction: Direction,
        weight: bool,
        state: &mut PhotonState
    ) {
        // Generate position.
        let cos_theta = rng.uniform01();
        let sin_theta = {
            let s2 = 1.0 - cos_theta * cos_theta;
            if s2 > 0.0 { s2.sqrt() } else { 0.0 }
        };
        let phi = 2.0 * PI * rng.uniform01();
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        let ur = Float3::new(cos_phi * sin_theta, sin_phi * sin_theta, cos_theta);
        let ut = Float3::new(cos_phi * cos_theta, sin_phi * cos_theta, -sin_theta);
        let up = Float3::new(-sin_phi, cos_phi, 0.0);
        let eps = 1E-04; // For numeric safety.
        let sgn = match side {
            Side::Inside => -1.0,
            Side::Outside => 1.0,
        };
        state.position = (self.radius + eps * sgn).max(0.0) * ur + self.center;

        // Generate direction.
        let xi = rng.uniform01();
        let cos_theta = xi.sqrt();
        let sin_theta = (1.0 - xi).sqrt();
        let phi = 2.0 * PI * rng.uniform01();
        let sgn = match direction {
            Direction::Ingoing => -1.0,
            Direction::Outgoing => 1.0,
        };
        state.direction =
            (sgn * phi.cos() * sin_theta) * ut +
            (sgn * phi.sin() * sin_theta) * up +
            (sgn * cos_theta) * ur;

        // Weight according to 1 / PDF.
        if weight {
            let tmp = 2.0 * PI * self.radius;
            state.weight *= tmp * tmp;
        }
    }
}

impl VectorisedOperations for SphereShape {}


// ===============================================================================================
// Vectorised geometry operations for shapes.
// ===============================================================================================

trait VectorisedOperations: GeometryShape {
    fn distance_v(&self, states: &PyArray<CState>) -> Result<PyObject> {
        let py = states.py();
        let result = PyArray::<Float>::empty(py, &states.shape())?;
        let n = states.size();
        for i in 0..n {
            let state: PhotonState = states.get(i)?.into();
            let distance = self
                .distance(state.position, state.direction)
                .unwrap_or(Float::NAN);
            result.set(i, distance)?;
        }
        let result: &PyAny = result;
        Ok(result.into_py(py))
    }

    fn inside_v(&self, states: &PyArray<CState>) -> Result<PyObject> {
        let py = states.py();
        let result = PyArray::<bool>::empty(py, &states.shape())?;
        let n = states.size();
        for i in 0..n {
            let state: PhotonState = states.get(i)?.into();
            let inside = self.inside(state.position);
            result.set(i, inside)?;
        }
        let result: &PyAny = result;
        Ok(result.into_py(py))
    }
}


// ===============================================================================================
// Generator interface.
// ===============================================================================================

trait Generate {
    fn generate_on<R: FloatRng>(
        &self,
        rng: &mut R,
        side: Side,
        direction: Direction,
        weight: bool,
        state: &mut PhotonState
    );

    fn generate_on_v<'p>(
        &self,
        states: &PyArray<CState>,
        engine: Option<&PyTransportEngine>,
        rng: Option<Py<PyRandomStream>>,
        side: Option<&str>,
        direction: Option<&str>,
        weight: Option<bool>,
    ) -> Result<()> {
        let py = states.py();
        let default_rng: Py<PyRandomStream>;
        let rng = match rng.as_ref() {
            None => match engine.as_ref() {
                None => {
                    default_rng = Py::new(py, PyRandomStream::new(None)?)?;
                    &default_rng
                },
                Some(engine) => &engine.random,
            },
            Some(rng) => rng,
        };
        let rng: &mut PyRandomStream = &mut rng.borrow_mut(py);
        let side = match side {
            None => match engine.as_ref() {
                None => Side::Inside,
                Some(engine) => match engine.settings.borrow(py).inner.mode {
                    TransportMode::Forward => Side::Inside,
                    TransportMode::Backward => Side::Outside,
                },
            },
            Some(side) => {
                let side: Side = side.try_into()?;
                side
            },
        };
        let direction = match direction {
            None => Direction::Ingoing,
            Some(direction) => {
                let direction: Direction = direction.try_into()?;
                direction
            },
        };
        let weight = match weight {
            None => match engine.as_ref() {
                None => false,
                Some(engine) => match engine.settings.borrow(py).inner.mode {
                    TransportMode::Forward => false,
                    TransportMode::Backward => true,
                },
            },
            Some(weight) => weight,
        };
        let n = states.size();
        for i in 0..n {
            let mut state: PhotonState = states.get(i)?.into();
            self.generate_on(rng, side, direction, weight, &mut state);
            states.set(i, state.into())?;
        }
        Ok(())
    }
}


// ===============================================================================================
// Generator direction and side flags.
// ===============================================================================================

#[derive(Clone, Copy, Sequence)]
enum Direction {
    Ingoing,
    Outgoing,
}

impl Direction {
    const INGOING: &str = "Ingoing";
    const OUTGOING: &str = "Outgoing";

    fn pretty_variants() -> String {
        let variants: Vec<_> = all::<Self>()
            .map(|e| format!("'{}'", e))
            .collect();
        pretty_enumerate(&variants)
    }
}

impl fmt::Display for Direction {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s: &str = (*self).into();
        write!(f, "{}", s)
    }
}

impl From<Direction> for &str {
    fn from(value: Direction) -> Self {
        match value {
            Direction::Ingoing => Direction::INGOING,
            Direction::Outgoing => Direction::OUTGOING,
        }
    }
}

impl TryFrom<&str> for Direction {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            Self::INGOING => Ok(Self::Ingoing),
            Self::OUTGOING => Ok(Self::Outgoing),
            _ => Err(anyhow!(
                "bad direction (expected {}, found '{}')",
                Self::pretty_variants(),
                value,
            )),
        }
    }
}

#[derive(Clone, Copy, Sequence)]
enum Side {
    Inside,
    Outside,
}

impl Side {
    const INSIDE: &str = "Inside";
    const OUTSIDE: &str = "Outside";

    fn pretty_variants() -> String {
        let variants: Vec<_> = all::<Self>()
            .map(|e| format!("'{}'", e))
            .collect();
        pretty_enumerate(&variants)
    }
}

impl fmt::Display for Side {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s: &str = (*self).into();
        write!(f, "{}", s)
    }
}

impl From<Side> for &str {
    fn from(value: Side) -> Self {
        match value {
            Side::Inside => Side::INSIDE,
            Side::Outside => Side::OUTSIDE,
        }
    }
}

impl TryFrom<&str> for Side {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            Self::INSIDE => Ok(Self::Inside),
            Self::OUTSIDE => Ok(Self::Outside),
            _ => Err(anyhow!(
                "bad side (expected {}, found '{}')",
                Self::pretty_variants(),
                value,
            )),
        }
    }
}
