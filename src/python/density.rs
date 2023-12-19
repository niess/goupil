use anyhow::Result;
use crate::numerics::float::{Float, Float3};
use crate::transport::density::DensityModel;
use pyo3::prelude::*;


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
            _ => unreachable!()
        }
    }
}


// ===============================================================================================
// Conversion between DensityModel and Pyobject.
// ===============================================================================================

#[derive(FromPyObject)]
enum DensityArg<'py> {
    Gradient(PyRef<'py, PyDensityGradient>),
    Uniform(Float),
}

impl<'source> FromPyObject<'source> for DensityModel {
    fn extract(ob: &'source PyAny) -> PyResult<Self> {
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
