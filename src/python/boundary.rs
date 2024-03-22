use anyhow::Result;
use crate::numerics::{Float, Float3};
use crate::transport::{PhotonState, SphereShape, TransportBoundary};
use pyo3::prelude::*;
use super::numpy::PyArray;
use super::transport::CState;


// ===============================================================================================
// Unresolved boundary.
// ===============================================================================================

#[derive(FromPyObject)]
pub enum PyTransportBoundary<'py> {
    Sector(usize),
    Sphere(PyRef<'py, PySphereShape>),
}

impl<'py> From<PyTransportBoundary<'py>> for TransportBoundary {
    fn from(boundary: PyTransportBoundary) -> Self {
        match boundary {
            PyTransportBoundary::Sector(index) => Self::Sector(index),
            PyTransportBoundary::Sphere(sphere) => Self::Sphere(sphere.0),
        }
    }
}

impl IntoPy<PyObject> for TransportBoundary {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            Self::None => py.None(),
            Self::Sector(index) => index.into_py(py),
            Self::Sphere(sphere) => {
                let sphere = PySphereShape(sphere);
                sphere.into_py(py)
            },
        }
    }
}


// ===============================================================================================
// Python wrapper for a spherical shape.
// ===============================================================================================

#[pyclass(name = "SphereShape", module = "goupil")]
pub struct PySphereShape(SphereShape);

#[pymethods]
impl PySphereShape {
    #[new]
    fn new(center: Option<Float3>, radius: Option<Float>) -> Self {
        let center = center.unwrap_or(Float3::new(0.0, 0.0, 0.0));
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
        let py = states.py();
        let result = PyArray::<Float>::empty(py, &states.shape())?;
        let n = states.size();
        for i in 0..n {
            let state: PhotonState = states.get(i)?.into();
            let distance = self.0
                .distance(state.position, state.direction)
                .unwrap_or(Float::NAN);
            result.set(i, distance)?;
        }
        let result: &PyAny = result;
        Ok(result.into_py(py))
    }

    fn inside(&self, states: &PyArray<CState>) -> Result<PyObject> {
        let py = states.py();
        let result = PyArray::<bool>::empty(py, &states.shape())?;
        let n = states.size();
        for i in 0..n {
            let state: PhotonState = states.get(i)?.into();
            let inside = self.0.inside(state.position);
            result.set(i, inside)?;
        }
        let result: &PyAny = result;
        Ok(result.into_py(py))
    }
}
