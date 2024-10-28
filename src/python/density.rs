use anyhow::Result;
use crate::numerics::float::{Float, Float3};
use crate::transport::density::DensityModel;
use pyo3::prelude::*;
use super::macros::value_error;
use super::numpy::{ArrayOrFloat3, PyArray, PyArrayMethods};


// ===============================================================================================
// Python wrapper for a simple geometry object.
// ===============================================================================================

#[pyclass(name = "DensityGradient", module = "goupil")]
pub struct PyDensityGradient (pub DensityModel);

#[pymethods]
impl PyDensityGradient {
    const DEFAULT_DIRECTION: Float3 = Float3(0.0, 0.0, -1.0);
    const DEFAULT_ORIGIN: Float3 = Float3(0.0, 0.0, 0.0);

    #[new]
    #[pyo3(signature = (density, scale, /, *, direction=None, origin=None))]
    fn new(
        density: Float,
        scale: Float,
        direction: Option<Float3>,
        origin: Option<Float3>,
    ) -> Result<Self> {
        let direction = direction.unwrap_or(Self::DEFAULT_DIRECTION);
        let origin = origin.unwrap_or(Self::DEFAULT_ORIGIN);
        let gradient = DensityModel::gradient(density, origin, scale, direction)?;
        Ok(Self(gradient))
    }

    fn __repr__(&self) -> String {
        match &self.0 {
            DensityModel::Gradient { rho0, origin, lambda, direction } => {
                if *origin != Self::DEFAULT_ORIGIN {
                    format!("DensityGradient({}, {}, {}, {})", rho0, lambda, direction, origin)
                } else if *direction != Self::DEFAULT_DIRECTION {
                    format!("DensityGradient({}, {}, {})", rho0, lambda, direction)
                } else {
                    format!("DensityGradient({}, {})", rho0, lambda)
                }
            },
            _ => unreachable!(),
        }
    }

    #[getter]
    fn get_density(&self) -> Float {
        match &self.0 {
            DensityModel::Gradient { rho0, .. } => *rho0,
            _ => unreachable!(),
        }
    }

    #[getter]
    fn get_direction(&self) -> Float3 {
        match &self.0 {
            DensityModel::Gradient { direction, .. } => *direction,
            _ => unreachable!(),
        }
    }

    #[getter]
    fn get_origin(&self) -> Float3 {
        match &self.0 {
            DensityModel::Gradient { origin, .. } => *origin,
            _ => unreachable!(),
        }
    }

    #[getter]
    fn get_scale(&self) -> Float {
        match &self.0 {
            DensityModel::Gradient { lambda, .. } => *lambda,
            _ => unreachable!(),
        }
    }

    fn __call__(&self, py: Python, position: ArrayOrFloat3) -> Result<PyObject> {
        let result: PyObject = match position {
            ArrayOrFloat3::Array(position) => {
                let shape = position.shape();
                let n = shape.len();
                if (n < 1) || (shape[n - 1] != 3) {
                    let shape: PyObject = shape.into_py(py);
                    value_error!("bad shape (expected [3] or [..., 3], found {})", shape);
                }
                let result = PyArray::<Float>::empty(py, &shape[0..(n-1)])?;
                let m = result.size();
                for i in 0..m {
                    let r = Float3::new(
                        position.get(3 * i)?,
                        position.get(3 * i + 1)?,
                        position.get(3 * i + 2)?,
                    );
                    let v = self.0.value(r);
                    result.set(i, v)?;
                }
                result.into_py(py)
            },
            ArrayOrFloat3::Float3(position) => self.0.value(position).into_py(py),
        };
        Ok(result)
    }
}


// ===============================================================================================
// Conversion between DensityModel and Pyobject.
// ===============================================================================================

#[derive(FromPyObject)]
pub(crate) enum DensityArg<'py> {
    Gradient(PyRef<'py, PyDensityGradient>),
    Uniform(Float),
}

impl<'py> FromPyObject<'py> for DensityModel {
    fn extract_bound(ob: &Bound<'py, PyAny>) -> PyResult<Self> {
        let density: DensityArg = ob.extract()?;
        let density = match density {
            DensityArg::Gradient(gradient) => gradient.0.clone(),
            DensityArg::Uniform(density) => DensityModel::uniform(density)?,
        };
        Ok(density)
    }
}

impl<'py> IntoPy<PyObject> for DensityModel {
    fn into_py(self, py: Python) -> PyObject {
        match &self {
            Self::Gradient { .. } => PyDensityGradient(self).into_py(py),
            Self::Uniform(density) => density.into_py(py),
        }
    }
}

impl ToPyObject for DensityModel {
    fn to_object(&self, py: Python) -> PyObject {
        self.into_py(py)
    }
}
