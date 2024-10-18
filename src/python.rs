use anyhow::Result;
use crate::numerics::{Float, Float3, Float3x3};
use pyo3::prelude::*;
use pyo3::exceptions::PyKeyboardInterrupt;
use pyo3::ffi;
use pyo3::sync::GILOnceCell;
use self::boundary::{PyBoxShape, PySphereShape};
use self::density::PyDensityGradient;
use self::elements::{elements as elements_fun, PyAtomicElement};
use self::geometry::{
    PyExternalGeometry,
    PyGeometrySector,
    PySimpleGeometry,
    PyStratifiedGeometry,
    PyTopographyMap,
    PyTopographySurface,
};
use self::materials::{
    PyCrossSection,
    PyDistributionFunction,
    PyElectronicStructure,
    PyFormFactor,
    PyInverseDistribution,
    PyMaterialDefinition,
    PyMaterialRecord,
    PyMaterialRegistry
};
use self::rand::PyRandomStream;
use process_path::get_dylib_path;
use self::process::{PyAbsorptionProcess, PyComptonProcess, PyRayleighProcess};
use self::spectrum::PyDiscreteSpectrum;
use self::states::{states as states_fun};
use self::transport::{PyTransportEngine, PyTransportSettings, PyTransportStatus};
use std::path::PathBuf;

mod boundary;
mod density;
mod elements;
mod export;
mod geometry;
mod materials;
mod namespace;
mod numpy;
mod rand;
mod process;
mod spectrum;
mod states;
mod transport;


//================================================================================================
// Check for a keyboard interrupt.
//================================================================================================

pub fn ctrlc_catched() -> Result<()> {
    if unsafe { ffi::PyErr_CheckSignals() } == -1 {
        Err(PyKeyboardInterrupt::new_err("").into())
    } else {
        Ok(())
    }
}


//================================================================================================
// Module prefix.
//================================================================================================

fn prefix<'py>(py: Python<'py>) -> Result<&'py PathBuf> {
    static PREFIX: GILOnceCell<PathBuf> = GILOnceCell::new();
    let prefix = PREFIX.get_or_try_init(py, || {
        let path = get_dylib_path();
        if let Some(mut path) = path {
            if path.pop() {
                return Ok::<_, PyErr>(path);
            }
        }
        Ok(".".into())
    })?;
    Ok(prefix)
}


//================================================================================================
// Helper macro(s) for bailing Python exceptions.
//================================================================================================

mod macros {
    macro_rules! key_error {
        ($($tts:tt)*) => {
            return Err(pyo3::exceptions::PyKeyError::new_err(format!($($tts)*)).into())
        }
    }
    pub(crate) use key_error;

    macro_rules! not_implemented_error {
        ($($tts:tt)*) => {
            return Err(pyo3::exceptions::PyNotImplementedError::new_err(format!($($tts)*)).into())
        }
    }
    pub(crate) use not_implemented_error;

    macro_rules! type_error {
        ($($tts:tt)*) => {
            return Err(pyo3::exceptions::PyTypeError::new_err(format!($($tts)*)).into())
        }
    }
    pub(crate) use type_error;

    macro_rules! value_error {
        ($($tts:tt)*) => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!($($tts)*)).into())
        }
    }
    pub(crate) use value_error;
}


//================================================================================================
// Implement from Python for Float3 and Float3x3.
//================================================================================================

impl<'py> FromPyObject<'py> for Float3 {
    fn extract_bound(ob: &Bound<'py, PyAny>) -> PyResult<Self> {
        let v: [Float; 3] = ob.extract()?;
        Ok(v.into())
    }
}

impl<'py> FromPyObject<'py> for Float3x3 {
    fn extract_bound(ob: &Bound<'py, PyAny>) -> PyResult<Self> {
        let v: [[Float; 3]; 3] = ob.extract()?;
        Ok(v.into())
    }
}


// ===============================================================================================
// Goupil Python 3 module.
// ===============================================================================================

#[pymodule]
fn goupil(module: &Bound<PyModule>) -> PyResult<()> {

    // Initialise Numpy array interface.
    let py = module.py();
    numpy::initialise(py)?;

    // Register attributes.
    module.add("PREFIX", prefix(py)?)?;

    // Register class object(s).
    module.add_class::<PyAbsorptionProcess>()?;
    module.add_class::<PyAtomicElement>()?;
    module.add_class::<PyBoxShape>()?;
    module.add_class::<PyComptonProcess>()?;
    module.add_class::<PyCrossSection>()?;
    module.add_class::<PyDensityGradient>()?;
    module.add_class::<PyDiscreteSpectrum>()?;
    module.add_class::<PyDistributionFunction>()?;
    module.add_class::<PyElectronicStructure>()?;
    module.add_class::<PyExternalGeometry>()?;
    module.add_class::<PyFormFactor>()?;
    module.add_class::<PyGeometrySector>()?;
    module.add_class::<PyInverseDistribution>()?;
    module.add_class::<PyMaterialDefinition>()?;
    module.add_class::<PyMaterialRecord>()?;
    module.add_class::<PyMaterialRegistry>()?;
    module.add_class::<PySimpleGeometry>()?;
    module.add_class::<PyStratifiedGeometry>()?;
    module.add_class::<PyRandomStream>()?;
    module.add_class::<PyRayleighProcess>()?;
    module.add_class::<PySphereShape>()?;
    module.add_class::<PyTopographyMap>()?;
    module.add_class::<PyTopographySurface>()?;
    module.add_class::<PyTransportEngine>()?;
    module.add_class::<PyTransportSettings>()?;
    module.add_class::<PyTransportStatus>()?;

    // Register function(s).
    module.add_function(wrap_pyfunction!(elements_fun, module)?)?;
    module.add_function(wrap_pyfunction!(states_fun, module)?)?;

    Ok(())
}
