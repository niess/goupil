use anyhow::Result;
use crate::numerics::float::Float;
use crate::physics::elements::AtomicElement;
use pyo3::prelude::*;
use pyo3::class::basic::CompareOp;
use pyo3::exceptions::PyNotImplementedError;
use pyo3::types::{PyBytes, PyTuple};
use rmp_serde::{Deserializer, Serializer};
use serde::{Deserialize, Serialize};
use super::materials::PyElectronicStructure;


// ===============================================================================================
// Python wrapper for an atomic element
// ===============================================================================================

#[pyclass(name = "AtomicElement", module = "goupil")]
pub struct PyAtomicElement (pub(crate) &'static AtomicElement);

#[pymethods]
impl PyAtomicElement {
    #[allow(non_snake_case)]
    #[new]
    #[pyo3(signature = (arg=None, /))]
    fn new(arg: Option<AtomArg>) -> Result<Self> {
        let element = match arg.as_ref() {
            None => AtomicElement::none(),
            Some(arg) => match arg {
                AtomArg::Z(Z) => AtomicElement::from_Z(*Z)?,
                AtomArg::Symbol(symbol) => AtomicElement::from_symbol(symbol)?,
            },
        };
        Ok(Self(element))
    }

    #[allow(non_snake_case)]
    #[getter]
    fn get_A(&self) -> Float {
        self.0.A
    }

    #[getter]
    fn get_name(&self) -> &str {
        &self.0.name
    }

    #[getter]
    fn get_symbol(&self) -> &str {
        &self.0.symbol
    }

    #[allow(non_snake_case)]
    #[getter]
    fn get_Z(&self) -> i32 {
        self.0.Z
    }

    fn __getstate__<'py>(&self, py: Python<'py>) -> Result<Bound<'py, PyBytes>> {
        let mut buffer = Vec::new();
        self.0.serialize(&mut Serializer::new(&mut buffer))?;
        Ok(PyBytes::new_bound(py, &buffer))
    }

    fn __setstate__(&mut self, state: &Bound<PyBytes>) -> Result<()> {
        self.0 = Deserialize::deserialize(&mut Deserializer::new(state.as_bytes()))?;
        Ok(())
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool> {
        match op {
            CompareOp::Eq => Ok(self.0 == other.0),
            CompareOp::Ne => Ok(self.0 != other.0),
            _ => Err(PyNotImplementedError::new_err("").into()),
        }
    }

    fn __str__(&self) -> &str {
        &self.0.symbol
    }

    fn __repr__(&self) -> &str {
        &self.0.symbol
    }

    fn electrons(&self, py: Python) -> Result<PyObject> {
        let electrons = self.0.electrons()?;
        let electrons = PyElectronicStructure::new(electrons.clone(), false)?;
        Ok(electrons.into_py(py))
    }
}

#[derive(FromPyObject)]
enum AtomArg {
    Symbol(String),
    Z(i32),
}

#[allow(non_snake_case)]
#[pyfunction]
#[pyo3(signature = (*args,))]
pub fn elements(args: &Bound<PyTuple>) -> Result<PyObject> {
    let py = args.py();
    let args: Vec<AtomArg> = args.extract()?;
    let mut elements = Vec::<PyObject>::with_capacity(args.len());
    for arg in args.iter() {
        match arg {
            AtomArg::Symbol(symbols) => {
                let symbols = symbols
                    .split(",")
                    .map(|s| s.trim());
                for symbol in symbols {
                    let element = PyAtomicElement(AtomicElement::from_symbol(symbol)?);
                    elements.push(element.into_py(py));
                }
            },
            AtomArg::Z(z) => {
                let element = PyAtomicElement(AtomicElement::from_Z(*z)?);
                elements.push(element.into_py(py));
            },
        };
    }
    let result = match elements.len() {
        0 => py.None(),
        1 => elements.pop().unwrap(),
        _ => PyTuple::new_bound(py, elements).into_any().unbind(),
    };
    Ok(result)
}
