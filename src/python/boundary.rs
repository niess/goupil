use anyhow::Result;
use crate::numerics::{Float, Float3, Float3x3};
use crate::transport::{BoxShape, GeometryShape, PhotonState, SphereShape, TransportBoundary};
use pyo3::prelude::*;
use super::numpy::{FloatOrFloat3, PyArray};
use super::transport::CState;


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
