use anyhow::{anyhow, Result};
use crate::numerics::{
    float::Float,
    grids::Grid,
};
use crate::physics::elements::AtomicElement;
use crate::physics::materials::{
    electronic::ElectronicStructure,
    MaterialDefinition,
    MaterialRecord,
    MaterialRegistry,
};
use crate::physics::process::absorption::{
    AbsorptionMode::{self, Discrete},
    table::AbsorptionCrossSection,
};
use crate::physics::process::compton::{
    ComptonModel::{self, ScatteringFunction},
    table::{ComptonCrossSection, ComptonCDF, ComptonInverseCDF},
    ComptonMethod::{self, RejectionSampling},
    ComptonMode::{self, Adjoint, Direct, Inverse},
};
use crate::physics::process::rayleigh::table::{RayleighCrossSection, RayleighFormFactor};
use crate::transport::{
    TransportMode::{self, Backward, Forward},
    TransportSettings,
};
use pyo3::{
    prelude::*,
    exceptions::PyKeyError,
    gc::PyVisit,
    PyTraverseError,
    types::{PyBytes, PyTuple},
};
use rmp_serde::{Deserializer, Serializer};
use serde::{Deserialize, Serialize};
use super::{
    elements::PyAtomicElement,
    macros::value_error,
    numpy::{ArrayOrFloat, PyArray, PyArrayMethods, PyArrayFlags},
    prefix,
    transport::PyTransportSettings,
};
use std::borrow::Cow;
use std::collections::HashMap;
use std::ops::Deref;
use std::path::PathBuf;


// ===============================================================================================
// Python wrapper for a material definition
// ===============================================================================================

#[pyclass(name = "MaterialDefinition", module = "goupil")]
#[derive(Clone)]
pub struct PyMaterialDefinition (pub MaterialDefinition);

enum Element<'py> {
    Object(Bound<'py, PyAtomicElement>),
    Symbol(String),
}

impl<'py> FromPyObject<'py> for Element<'py> {
    fn extract_bound(ob: &Bound<'py, PyAny>) -> PyResult<Self> {
        let result: PyResult<Bound<PyAtomicElement>> = ob.extract();
        if let PyResult::Ok(result) = result {
            return PyResult::Ok(Self::Object(result));
        }
        let result: PyResult<String> = ob.extract();
        match result {
            PyResult::Ok(symbol) => {
                let _unused = AtomicElement::from_symbol(symbol.as_str())?;
                PyResult::Ok(Self::Symbol(symbol))
            },
            PyResult::Err(err) => PyResult::Err(err),
        }
    }
}

impl<'py> Element<'py> {
    fn get(&self) -> Result<&'static AtomicElement> {
        let element = match self {
            Element::Object(element) => element.borrow().0,
            Element::Symbol(symbol) => AtomicElement::from_symbol(symbol)?,
        };
        Ok(element)
    }

    fn weight(&self, value: Float) -> Result<(Float, &'static AtomicElement)> {
        let element = self.get()?;
        Ok((value, element))
    }
}

#[derive(FromPyObject)]
enum Material<'py> {
    Formula(String),
    Object(PyRef<'py, PyMaterialDefinition>),
}

impl<'py> Material<'py> {
    fn get<'a>(&'a self) -> Result<Cow<'a, MaterialDefinition>> {
        let material = match self {
            Material::Formula(material) => {
                Cow::Owned(MaterialDefinition::from_formula(material.as_str())?)
            },
            Material::Object(material) => Cow::Borrowed(&material.deref().0),
        };
        Ok(material)
    }

    fn weight(&self, value: Float) -> Result<(Float, Cow<MaterialDefinition>)> {
        let material = self.get()?;
        Ok((value, material))
    }
}

#[derive(FromPyObject)]
enum PyMassComposition<'py> {
    Atomic(Vec<(Float, Element<'py>)>),
    Material(Vec<(Float, Material<'py>)>),
}

type PyMoleComposition<'py> = Vec<(Float, Element<'py>)>;

#[pymethods]
impl PyMaterialDefinition {
    #[new]
    #[pyo3(signature = (name=None, /, *, mass_composition=None, mole_composition=None))]
    fn new(
        name: Option<&str>,
        mass_composition: Option<PyMassComposition>,
        mole_composition: Option<PyMoleComposition>,
    ) -> Result<Self> {
        let definition = match name {
            None => {
                if !mass_composition.is_none() || !mole_composition.is_none() {
                    value_error!("bad material name (expected a string value, found None)")
                }
                MaterialDefinition::default()
            },
            Some(name) => match mass_composition {
                None => match mole_composition {
                    None => {
                        // Try to interpret name as a chemical formula.
                        MaterialDefinition::from_formula(name)?
                    },
                    Some(composition) => {
                        let composition: Result<Vec<_>> = composition
                            .iter()
                            .map(|(weight, element)| element.weight(*weight))
                            .collect();
                        let composition = composition?;
                        MaterialDefinition::from_mole(name, &composition)
                    },
                },
                Some(composition) => {
                    if let Some(_) = mole_composition {
                        value_error!(
                            "bad composition for '{}' (expected one of 'mass_composition' or \
                                'mole_composition', found both of them)",
                            name,
                        )
                    }
                    match composition {
                        PyMassComposition::Atomic(composition) => {
                            let composition: Result<Vec<_>> = composition
                                .iter()
                                .map(|(weight, element)| element.weight(*weight))
                                .collect();
                            let composition = composition?;
                            MaterialDefinition::from_mass(name, &composition)
                        },
                        PyMassComposition::Material(composition) => {
                            let composition: Result<Vec<_>> = composition
                                .iter()
                                .map(|(weight, material)| material.weight(*weight))
                                .collect();
                            let composition = composition?;
                            let references: Vec<_> = composition
                                .iter()
                                .map(|(weight, material)| (*weight, material.as_ref()))
                                .collect();
                            MaterialDefinition::from_others(name, &references)
                        },
                    }
                },
            }
        };
        Ok(Self(definition))
    }

    // Implementation of equality test.
    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    // Implementation of pickling protocol.
    pub fn __setstate__(&mut self, state: &Bound<PyBytes>) -> Result<()> {
        self.0 = Deserialize::deserialize(&mut Deserializer::new(state.as_bytes()))?;
        Ok(())
    }

    fn __getstate__<'py>(&self, py: Python<'py>) -> Result<Bound<'py, PyBytes>> {
        let mut buffer = Vec::new();
        self.0.serialize(&mut Serializer::new(&mut buffer))?;
        Ok(PyBytes::new_bound(py, &buffer))
    }

    fn __repr__(&self) -> &str {
        self.0.name()
    }

    // Public interface, as getters.
    #[getter]
    fn get_mass(&self) -> Float {
        self.0.mass()
    }

    #[getter]
    fn get_mass_composition<'py>(&self, py: Python<'py>) -> Bound<'py, PyTuple> {
        let composition: Vec<_> = self.0
            .mass_composition()
            .iter()
            .map(|(weight, element)| (*weight, PyAtomicElement(*element).into_py(py)))
            .collect();
        PyTuple::new_bound(py, composition)
    }

    #[getter]
    fn get_mole_composition<'py>(&self, py: Python<'py>) -> Bound<'py, PyTuple> {
        let composition: Vec<_> = self.0
            .mole_composition()
            .iter()
            .map(|(weight, element)| (*weight, PyAtomicElement(*element).into_py(py)))
            .collect();
        PyTuple::new_bound(py, composition)
    }

    #[getter]
    fn get_name(&self) -> &str {
        self.0.name()
    }

    fn electrons(&self, py: Python) -> Result<PyObject> {
        let electrons = self.0.compute_electrons()?;
        let electrons = PyElectronicStructure::new(electrons, false)?;
        Ok(electrons.into_py(py))
    }
}


// ===============================================================================================
// Arguments alike material definitions.
// ===============================================================================================

#[derive(FromPyObject)]
pub enum MaterialLike<'py> {
    Definition(PyRef<'py, PyMaterialDefinition>),
    Formula(String),
    Record(PyRefMut<'py, PyMaterialRecord>),
}

// Public interface.
impl<'py> MaterialLike<'py> {
    pub fn get_electrons(&self) -> Result<Cow<'py, ElectronicStructure>> {
        // Fetch reference to any computed structure.
        if let Self::Record(record) = self {
            let py = record.py();
            if let Some(electrons) = record.get(py)?.electrons() {
                return Ok(Cow::Borrowed(electrons))
            }
        }

        // Build a new instance.
        let definition = self.unpack()?;
        let electrons = definition.compute_electrons()?;
        Ok(Cow::Owned(electrons))
    }

    pub fn pack(self, py: Python) -> Result<PyObject> {
        let result: PyObject = match self {
            Self::Definition(definition) => definition.into_py(py),
            Self::Formula(formula) => {
                let definition = PyMaterialDefinition::new(Some(formula.as_str()), None, None)?;
                let definition = Py::new(py, definition)?;
                definition.into_py(py)
            },
            Self::Record(mut record) => record.get_definition(py)?.into_py(py),
        };
        Ok(result)
    }

    pub fn absorption_cross_section(&self, py: Python) -> Result<Cow<AbsorptionCrossSection>> {
        // Fetch reference to any computed table.
        if let Self::Record(record) = self {
            let py = record.py();
            if let Some(table) = record.get(py)?.absorption_cross_section() {
                return Ok(Cow::Borrowed(table))
            }
        }

        // Build a new table.
        let definition = self.unpack()?;
        let registry = Self::new_registry(py, &definition)?;
        let mut composition = Vec::<(Float, &AbsorptionCrossSection)>::default();
        for (weight, element) in definition.mole_composition().iter() {
            let cross_section = match registry.absorption.get(element) {
                None => value_error!(
                    "missing scattering cross-section for '{}'",
                    element.symbol,
                ),
                Some(table) => table,
            };
            composition.push((*weight, cross_section));
        }
        let table = AbsorptionCrossSection::from_others(&composition).unwrap();
        Ok(Cow::Owned(table))
    }

    pub fn rayleigh_cross_section(&self, py: Python) -> Result<Cow<RayleighCrossSection>> {
        // Fetch reference to any computed table.
        if let Self::Record(record) = self {
            let py = record.py();
            if let Some(table) = record.get(py)?.rayleigh_cross_section() {
                return Ok(Cow::Borrowed(table))
            }
        }

        // Build a new table.
        let definition = self.unpack()?;
        let registry = Self::new_registry(py, &definition)?;
        let mut composition = Vec::<(Float, &RayleighCrossSection)>::default();
        for (weight, element) in definition.mole_composition().iter() {
            let cross_section = match registry.scattering_cs.get(element) {
                None => value_error!(
                    "missing scattering cross-section for '{}'",
                    element.symbol,
                ),
                Some(table) => table,
            };
            composition.push((*weight, cross_section));
        }
        let table = RayleighCrossSection::from_others(&composition).unwrap();
        Ok(Cow::Owned(table))
    }

    pub fn rayleigh_form_factor(&self, py: Python) -> Result<Cow<RayleighFormFactor>> {
        // Fetch reference to any computed table.
        if let Self::Record(record) = self {
            let py = record.py();
            if let Some(table) = record.get(py)?.rayleigh_form_factor() {
                return Ok(Cow::Borrowed(table))
            }
        }

        // Build a new table.
        let definition = self.unpack()?;
        let registry = Self::new_registry(py, &definition)?;
        let mut composition = Vec::<(Float, &RayleighFormFactor)>::default();
        for (weight, element) in definition.mole_composition().iter() {
            let form_factor = match registry.scattering_ff.get(element) {
                None => value_error!(
                    "missing scattering form-factor for '{}'",
                    element.symbol,
                ),
                Some(table) => table,
            };
            composition.push((*weight, form_factor));
        }
        let table = RayleighFormFactor::from_others(&composition).unwrap();
        Ok(Cow::Owned(table))
    }

    pub fn unpack(&self) -> Result<Cow<MaterialDefinition>> {
        let result = match self {
            Self::Definition(definition) => Cow::Borrowed(&definition.0),
            Self::Formula(formula) => {
                let definition = MaterialDefinition::from_formula(formula)?;
                Cow::Owned(definition)
            },
            Self::Record(record) => match &record.proxy {
                RecordProxy::Borrowed { name, registry } => {
                    let py = record.py();
                    let registry: PyRef<PyMaterialRegistry> = registry.extract(py)?;
                    let record = registry.inner.get(name)?;
                    Cow::Owned(record.definition.clone())
                },
                RecordProxy::Owned(record) => Cow::Borrowed(&record.definition),
            },
        };
        Ok(result)
    }
}

// Private interface.
impl<'py> MaterialLike<'py> {
    pub(crate) fn new_registry(
        py: Python,
        definition: &MaterialDefinition
    ) -> Result<MaterialRegistry> {
        let mut registry = MaterialRegistry::default();
        let path = PyMaterialRegistry::default_elements_path(py)?;
        registry.load_elements(&path)?;
        registry.add(&definition)?;
        Ok(registry)
    }
}


// ===============================================================================================
// Python wrapper for a material registry.
// ===============================================================================================

#[pyclass(name = "MaterialRegistry", module = "goupil")]
pub struct PyMaterialRegistry {
    pub inner: MaterialRegistry,
    proxies: HashMap<String, Py<PyMaterialRecord>>,
}

#[pymethods]
impl PyMaterialRegistry {
    #[new]
    #[pyo3(signature = (*args))]
    pub fn new(args: Vec<MaterialLike>) -> Result<Self> {
        let mut registry = MaterialRegistry::default();
        for definition in args.iter() {
            let definition = definition.unpack()?;
            registry.add(&definition)?;
        }
        Ok(Self{
            inner: registry,
            proxies: HashMap::<String, Py<PyMaterialRecord>>::default(),
        })
    }

    // Implementation of GC protocol.
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        for record in self.proxies.values() {
            visit.call(record)?
        }
        Ok(())
    }

    fn __clear__(&mut self) {
        self.proxies.clear();
    }

    // Implementation of mapping protocol.
    fn __delitem__(&mut self, py: Python, key: &str) -> Result<()> {
        let record = self.inner.remove(key)
            .map_err(|e| PyKeyError::new_err(e.to_string()))?;
        let proxy = self.proxies.remove(key);
        if let Some(proxy) = proxy {
            Self::into_owned(py, proxy, record);
        }
        Ok(())
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__(
        slf: Py<PyMaterialRegistry>,
        py: Python,
        key: &str
    ) -> Result<Py<PyMaterialRecord>> {
        let mut registry = slf.bind(py).borrow_mut();
        registry.inner.get(key)
            .map_err(|e| PyKeyError::new_err(e.to_string()))?;
        if let Some(record) = registry.proxies.get(key) {
            Ok(record.clone())
        } else {
            let proxy = RecordProxy::Borrowed {
                name: key.to_string(),
                registry: slf.clone(),
            };
            let record = Py::new(py, PyMaterialRecord::new(proxy))?;
            registry.proxies.insert(key.to_string(), record.clone());
            Ok(record)
        }
    }

    fn __repr__(&self) -> String {
        let items = self.inner.keys().join(", ");
        format!("{{{}}}", items)
    }

    // Implementation of pickling protocol.
    pub fn __setstate__(&mut self, py: Python, state: &Bound<PyBytes>) -> Result<()> {
        // Detach pending record(s).
        for (k, v) in self.proxies.drain() {
            let record = self.inner.remove(&k).unwrap();
            Self::into_owned(py, v, record);
        }

        // Update inner registry.
        self.inner = Deserialize::deserialize(&mut Deserializer::new(state.as_bytes()))?;
        Ok(())
    }

    fn __getstate__<'py>(&self, py: Python<'py>) -> Result<Bound<'py, PyBytes>> {
        let mut buffer = Vec::new();
        self.inner.serialize(&mut Serializer::new(&mut buffer))?;
        Ok(PyBytes::new_bound(py, &buffer))
    }

    // Direct public interface.
    fn add(&mut self, definition: MaterialLike) -> Result<()> {
        let definition = definition.unpack()?;
        self.inner.add(&definition)?;
        Ok(())
    }

    #[pyo3(signature = (settings=None, /, *, shape=None, precision=None, mode=None, absorption=None, compton_method=None, compton_model=None, compton_mode=None, constraint=None, energy_max=None, energy_min=None))]
    fn compute(
        &mut self,
        py: Python<'_>,
        settings: Option<&PyTransportSettings>,
        shape: Option<PyObject>,
        precision: Option<Float>,
        // Override arguments.
        mode: Option<&str>,
        absorption: Option<&str>,
        compton_method: Option<&str>,
        compton_model: Option<&str>,
        compton_mode: Option<&str>,
        constraint: Option<bool>,
        energy_max: Option<Float>,
        energy_min: Option<Float>,
    ) -> Result<()> {
        let (length, width) = {
            // Parse shape, trying various patterns.
            match shape {
                None => (None, None),
                Some(shape) => match shape.extract::<(usize, usize)>(py) {
                    Ok(v) => (Some(v.0), Some(v.1)),
                    _ => match shape.extract::<[usize; 2]>(py) {
                        Ok(v) => (Some(v[0]), Some(v[1])),
                        _ => match shape.extract::<usize>(py) {
                            Ok(v) => (Some(v), None),
                            _ => value_error!(
                                "bad shape (expected an integer or a size 2 sequence, found {:})",
                                shape.to_string(),
                            ),
                        },
                    },
                },
            }
        };

        let mut config = TransportSettings::default();

        let mode = mode
            .map(|s| TransportMode::try_from(s))
            .or(settings.map(|settings| Ok(settings.inner.mode)));
        let mode = match mode {
            None => None,
            Some(mode) => {
                let mode = mode?;
                config.mode = mode;
                match mode {
                    Backward => config.compton_mode = ComptonMode::Adjoint,
                    Forward => config.compton_mode = ComptonMode::Direct,
                }
                Some(mode)
            },
        };

        config.absorption = absorption
            .map(|s| AbsorptionMode::try_from(s))
            .or(settings.map(|settings| Ok(settings.inner.absorption)))
            .unwrap_or(Ok(Discrete))?;

        config.compton_model = compton_model
            .map(|s| ComptonModel::try_from(s))
            .or(settings.map(|settings| Ok(settings.inner.compton_model)))
            .unwrap_or(Ok(ScatteringFunction))?;

        config.energy_max = energy_max
            .or(settings.and_then(|settings| settings.inner.energy_max));

        config.energy_min = energy_min
            .or(settings.and_then(|settings| settings.inner.energy_min));

        config.compton_method = compton_method
            .map(|s| ComptonMethod::try_from(s))
            .or(settings.map(|settings| Ok(settings.inner.compton_method)))
            .unwrap_or(Ok(RejectionSampling))?;

        config.compton_mode = compton_mode
            .map(|s| ComptonMode::try_from(s))
            .or(settings.map(|settings| Ok(settings.inner.compton_mode)))
            .unwrap_or(Ok(Direct))?;

        match &config.compton_mode {
            Adjoint | Inverse => match mode {
                None => config.mode = Backward,
                Some(mode) => if let Forward = mode {
                    value_error!(
                        "bad transport mode for compton mode '{}' (expected '{}', found '{}')",
                        config.compton_mode,
                        Backward,
                        mode
                    )
                },
            }
            Direct => match mode {
                None => config.mode = Forward,
                Some(mode) => if let Backward = mode {
                    value_error!(
                        "bad transport mode for compton mode '{}' (expected '{}', found '{}')",
                        config.compton_mode,
                        Forward,
                        mode
                    )
                },
            },
            ComptonMode::None => (),
        }

        config.constraint = match config.compton_mode {
            Direct => None,
            _ => match constraint {
                None => settings.and_then(|settings| settings.inner.constraint ),
                Some(constraint) => if constraint { Some(1.0) } else { None },
            },
        };

        if !self.inner.atomic_data_loaded() {
            // Load default atomic data.
            let mut path = prefix(py)?.clone();
            path.push(PyMaterialRegistry::ELEMENTS_DATA);
            self.inner.load_elements(&path)?;
        }

        self.inner.compute(
            &config,
            length,
            width,
            precision,
        )
    }

    #[pyo3(signature = (path=None, /))]
    fn load_elements(&mut self, py: Python, path: Option<String>) -> Result<()> {
        let path = match path {
            None => Self::default_elements_path(py)?,
            Some(path) => path.into(),
        };
        self.inner.load_elements(&path)?;
        Ok(())
    }
}

impl PyMaterialRegistry {
    pub(crate) const ELEMENTS_DATA: &str = "data/elements";

    fn default_elements_path(py: Python) -> Result<PathBuf> {
        let mut path = prefix(py)?.clone();
        path.push(Self::ELEMENTS_DATA);
        Ok(path)
    }

    // Transforms a borrowed material record to an owned one.
    fn into_owned(py: Python, proxy: Py<PyMaterialRecord>, record: MaterialRecord) {
        if proxy.get_refcnt(py) >= 2 {
            let mut proxy = proxy.bind(py).borrow_mut();
            proxy.proxy = RecordProxy::Owned(record);
        }
    }
}

impl Drop for PyMaterialRegistry {
    fn drop(&mut self) {
        Python::with_gil(|py| {
            for (k, v) in self.proxies.drain() {
                let record = self.inner.remove(&k).unwrap();
                Self::into_owned(py, v, record);
            }
        })
    }
}


// ===============================================================================================
// Python wrapper for a material record.
// ===============================================================================================

#[pyclass(name = "MaterialRecord", module = "goupil")]
pub struct PyMaterialRecord {
    proxy: RecordProxy,
    definition: Option<Py<PyMaterialDefinition>>,
    electrons: Option<Py<PyElectronicStructure>>,
}

pub enum RecordProxy {
    Borrowed { name: String, registry: Py<PyMaterialRegistry> },
    Owned(MaterialRecord),
}

impl PyMaterialRecord {
    // Returns a reference to the underlying record, with lifetime bounded by the GIL.
    // Note that the current implementation uses unsafe pointer to (owned) PyObject.
    pub(crate) fn get<'py>(&self, py: Python<'py>) -> Result<&'py MaterialRecord> {
        let ptr = match &self.proxy {
            RecordProxy::Borrowed {name, registry} => registry
                .clone()
                .bind(py)
                .borrow()
                .inner
                .get(name)? as *const MaterialRecord,
            RecordProxy::Owned(record) => record as *const MaterialRecord,
        };
        unsafe { ptr.as_ref() }
            .ok_or_else(|| anyhow!("null pointer"))
    }

    fn new(proxy: RecordProxy) -> Self {
        Self {
            proxy,
            definition: None,
            electrons: None,
        }
    }
}

#[pymethods]
impl PyMaterialRecord {
    // Implementation of GC protocol.
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Some(definition) = self.definition.as_ref() {
            visit.call(definition)?
        }
        if let Some(electrons) = self.electrons.as_ref() {
            visit.call(electrons)?
        }
        Ok(())
    }

    fn __clear__(&mut self) {
        if self.definition.is_some() {
            self.definition = None;
        }
        if self.electrons.is_some() {
            self.electrons = None;
        }
    }

    fn __repr__(&self) -> &str {
        match &self.proxy {
            RecordProxy::Borrowed {name, ..} => name.as_str(),
            RecordProxy::Owned(record) => record.definition.name(),
        }
    }

    #[getter]
    fn get_definition(&mut self, py: Python) -> Result<Py<PyMaterialDefinition>> {
        match &self.definition {
            None => {
                let definition = PyMaterialDefinition(self.get(py)?.definition().clone());
                let definition = Py::new(py, definition)?;
                self.definition = Some(definition.clone());
                Ok(definition)
            },
            Some(definition) => Ok(definition.clone()),
        }
    }

    #[getter]
    fn get_electrons(&mut self, py: Python) -> Result<PyObject> {
        match &self.electrons {
            None => {
                let record = self.get(py)?;
                let electrons = match record.electrons() {
                    None => record.definition().compute_electrons()?,
                    Some(electrons) => electrons.clone(),
                };
                let electrons = PyElectronicStructure::new(electrons.clone(), false)?;
                let electrons = Py::new(py, electrons)?;
                self.electrons = Some(electrons.clone());
                let object: PyObject = electrons.into_py(py);
                Ok(object)
            },
            Some(electrons) => Ok(electrons.clone().into_py(py)),
        }
    }

    fn absorption_cross_section(this: &Bound<Self>) -> Result<PyObject> {
        PyCrossSection::new_absorption(this)
    }

    #[pyo3(signature = (*, model=None, mode=None))]
    fn compton_cdf(
        this: &Bound<Self>,
        py: Python,
        model: Option<&str>,
        mode: Option<&str>,
    ) -> Result<PyObject> {
        let model = model
            .map(|model| ComptonModel::try_from(model))
            .unwrap_or(Ok(ScatteringFunction))?;
        let mode = mode
            .map(|mode| ComptonMode::try_from(mode))
            .unwrap_or(Ok(Direct))?;
        PyDistributionFunction::new(py, this, model, mode)
    }

    #[pyo3(signature = (*, model=None, mode=None))]
    fn compton_cross_section(
        this: &Bound<Self>,
        model: Option<&str>,
        mode: Option<&str>,
    ) -> Result<PyObject> {
        let model = model
            .map(|model| ComptonModel::try_from(model))
            .unwrap_or(Ok(ScatteringFunction))?;
        let mode = mode
            .map(|mode| ComptonMode::try_from(mode))
            .unwrap_or(Ok(Direct))?;
        PyCrossSection::new_compton(this, model, mode)
    }

    #[pyo3(signature = (*, model=None, mode=None))]
    fn compton_inverse_cdf(
        this: &Bound<Self>,
        model: Option<&str>,
        mode: Option<&str>,
    ) -> Result<PyObject> {
        let model = model
            .map(|model| ComptonModel::try_from(model))
            .unwrap_or(Ok(ScatteringFunction))?;
        let mode = mode
            .map(|mode| ComptonMode::try_from(mode))
            .unwrap_or(Ok(Direct))?;
        PyInverseDistribution::new(this, model, mode)
    }

    #[pyo3(signature = (energy_in, energy_out, *, model=None, mode=None))]
    fn compton_weight(
        &mut self,
        py: Python,
        energy_in: Float,
        energy_out: Float,
        model: Option<&str>,
        mode: Option<&str>,
    ) -> Result<Float> {
        let model = model
            .map(|model| ComptonModel::try_from(model))
            .unwrap_or(Ok(ScatteringFunction))?;
        let mode = mode
            .map(|mode| ComptonMode::try_from(mode))
            .unwrap_or(Ok(Adjoint))?;
        self.get(py)?.compton_weight(model, mode, energy_in, energy_out)
    }

    fn rayleigh_cross_section(this: &Bound<Self>) -> Result<PyObject> {
        PyCrossSection::new_rayleigh(this)
    }

    fn rayleigh_form_factor(this: &Bound<Self>) -> Result<PyObject> {
        PyFormFactor::new(this)
    }
}


// ===============================================================================================
// Python wrapper for an ElectronicStructure object.
// ===============================================================================================

#[pyclass(name = "ElectronicStructure", module = "goupil")]
pub struct PyElectronicStructure {
    electrons: ElectronicStructure,
    writable: bool,
    shells: Option<PyObject>,
}

impl PyElectronicStructure {
    pub(crate) fn new(electrons: ElectronicStructure, writable: bool) -> Result<Self> {
        Ok(Self {
            electrons,
            writable,
            shells: None,
        })
    }
}

#[pymethods]
impl PyElectronicStructure {
    // Implementation of GC protocol.
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Some(shells) = self.shells.as_ref() {
            visit.call(shells)?;
        }
        Ok(())
    }

    fn __clear__(&mut self) {
        if self.shells.is_some() {
            self.shells = None;
        }
    }

    #[getter]
    fn get_charge(&self) -> Float {
        self.electrons.charge()
    }

    #[getter]
    fn get_shells(slf: &Bound<Self>) -> Result<PyObject> {
        let py = slf.py();
        let mut obj = slf.borrow_mut();
        if obj.shells.is_none() {
            let flags = if obj.writable {
                PyArrayFlags::ReadWrite
            } else {
                PyArrayFlags::ReadOnly
            };
            let shells = PyArray::from_data(
                py,
                &obj.electrons,
                slf,
                flags,
                None,
            )?;
            obj.shells = Some(shells.into_any().unbind());
        }
        let shells = obj.shells
            .as_ref()
            .unwrap()
            .clone_ref(py);
        Ok(shells)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.electrons == other.electrons
    }
}


// ===============================================================================================
// Process wrapper.
// ===============================================================================================

#[derive(Clone, Copy)]
enum Process {
    Absorption,
    Compton(ComptonModel, ComptonMode),
    Rayleigh,
}

impl Into<String> for Process {
    fn into(self) -> String {
        match self {
            Self::Absorption => "Absorption".to_string(),
            Self::Compton(model, mode) => format!("Compton::{}::{}", model, mode),
            Self::Rayleigh => "Rayleigh".to_string(),
        }
    }
}


// ===============================================================================================
// Python wrapper for a CrossSection object.
// ===============================================================================================

#[pyclass(name = "CrossSection", module="goupil")]
pub struct PyCrossSection {
    #[pyo3(get)]
    energies: PyObject,
    #[pyo3(get)]
    material: PyObject,
    #[pyo3(get)]
    values: PyObject,

    process: Process
}

impl PyCrossSection {
    fn new_absorption(record: &Bound<PyMaterialRecord>) -> Result<PyObject> {
        let py = record.py();
        let (energies, values) = {
            match Self::table_absorption(py, &record.borrow())? {
                None => return Ok(py.None()),
                Some(table) => {
                    let energies = readonly1(table.energies.as_ref(), record)?;
                    let values = readonly1(table.values.as_ref(), record)?;
                    (energies, values)
                },
            }
        };
        let material: PyObject = record.into_py(py);
        let process = Process::Absorption;
        let result = Self { energies, material, values, process };
        Ok(result.into_py(py))
    }

    fn new_compton(
        record: &Bound<PyMaterialRecord>,
        model: ComptonModel,
        mode: ComptonMode,
    ) -> Result<PyObject> {
        let py = record.py();
        let (energies, values) = {
            match Self::table_compton(py, &record.borrow(), model, mode)? {
                None => return Ok(py.None()),
                Some(table) => {
                    let energies = readonly1(table.energies.as_ref(), record)?;
                    let values = readonly1(table.values.as_ref(), record)?;
                    (energies, values)
                },
            }
        };
        let material: PyObject = record.into_py(py);
        let process = Process::Compton(model, mode);
        let result = Self { energies, material, values, process };
        Ok(result.into_py(py))
    }

    fn new_rayleigh(record: &Bound<PyMaterialRecord>) -> Result<PyObject> {
        let py = record.py();
        let (energies, values) = {
            match Self::table_rayleigh(py, &record.borrow())? {
                None => return Ok(py.None()),
                Some(table) => {
                    let energies = readonly1(table.energies.as_ref(), record)?;
                    let values = readonly1(table.values.as_ref(), record)?;
                    (energies, values)
                },
            }
        };
        let material: PyObject = record.into_py(py);
        let process = Process::Rayleigh;
        let result = Self { energies, material, values, process };
        Ok(result.into_py(py))
    }

    fn table_absorption<'py>(
        py: Python<'py>,
        record: &PyMaterialRecord
    ) -> Result<Option<&'py AbsorptionCrossSection>> {
        Ok(record.get(py)?.table.absorption.as_ref())
    }

    fn table_compton<'py>(
        py: Python<'py>,
        record: &PyMaterialRecord,
        model: ComptonModel,
        mode: ComptonMode
    ) -> Result<Option<&'py ComptonCrossSection>> {
        Ok(record.get(py)?.table.compton.get(model).get(mode).cross_section.as_ref())
    }

    fn table_rayleigh<'py>(
        py: Python<'py>,
        record: &PyMaterialRecord,
    ) -> Result<Option<&'py RayleighCrossSection>> {
        Ok(record.get(py)?.table.rayleigh.cross_section.as_ref())
    }
}

#[pymethods]
impl PyCrossSection {
    // Implementation of GC protocol.
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.energies)?;
        visit.call(&self.material)?;
        visit.call(&self.values)?;
        Ok(())
    }

    #[getter]
    fn get_process(&self) -> String {
        self.process.into()
    }

    fn __call__(
        &self,
        py: Python,
        energy: ArrayOrFloat
    ) -> Result<PyObject> {
        let record: PyRef<PyMaterialRecord> = self.material.extract(py)?;
        let compute_cross_section = |energy: Float| -> Result<Float> {
            let value = match self.process {
                Process::Absorption => {
                    let table = Self::table_absorption(py, &record)?.unwrap();
                    table.interpolate(energy)
                },
                Process::Compton(model, mode) => {
                    let table = Self::table_compton(py, &record, model, mode)?.unwrap();
                    table.interpolate(energy)
                },
                Process::Rayleigh => {
                    let table = Self::table_rayleigh(py, &record)?.unwrap();
                    table.interpolate(energy)
                }
            };
            Ok(value)
        };
        let value: PyObject = match energy {
            ArrayOrFloat::Array(energy) => {
                let value = PyArray::<Float>::empty(py, &energy.shape())?;
                let n = energy.size();
                for i in 0..n {
                    let v = compute_cross_section(energy.get(i)?)?;
                    value.set(i, v)?;
                }
                value.into_py(py)
            },
            ArrayOrFloat::Float(energy) => {
                let value = compute_cross_section(energy)?;
                value.into_py(py)
            },
        };
        Ok(value)
    }
}


// ===============================================================================================
// Python wrapper for a ComptonCDF object.
// ===============================================================================================

#[pyclass(name = "DistributionFunction", module="goupil")]
pub struct PyDistributionFunction {
    #[pyo3(get)]
    energies_in: PyObject,
    #[pyo3(get)]
    material: PyObject,
    #[pyo3(get)]
    values: PyObject,
    #[pyo3(get)]
    x: PyObject,

    model: ComptonModel,
    mode: ComptonMode,
}

impl PyDistributionFunction {
    fn new(
        py: Python,
        record: &Bound<PyMaterialRecord>,
        model: ComptonModel,
        mode: ComptonMode,
    ) -> Result<PyObject> {
        let (energies_in, x, values) = {
            match Self::table(py, &record.borrow(), model, mode)? {
                None => return Ok(py.None()),
                Some(table) => {
                    let energies_in = readonly1(table.energies_in.as_ref(), record)?;
                    let x = copy1(py, table.x.len(), table.x.iter())?;
                    let values = readonly2(table.values.as_ref(), table.shape(), record)?;
                    (energies_in, x, values)
                },
            }
        };
        let material: PyObject = record.into_py(py);
        let result = Self { energies_in, material, values, x, model, mode };
        Ok(result.into_py(py))
    }

    fn table<'py>(
        py: Python<'py>,
        record: &PyMaterialRecord,
        model: ComptonModel,
        mode: ComptonMode
    ) -> Result<Option<&'py ComptonCDF>> {
        Ok(record.get(py)?.table.compton.get(model).get(mode).cdf.as_ref())
    }
}

#[pymethods]
impl PyDistributionFunction {
    // Implementation of GC protocol.
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.energies_in)?;
        visit.call(&self.material)?;
        visit.call(&self.values)?;
        visit.call(&self.x)?;
        Ok(())
    }

    #[getter]
    fn get_process(&self) -> String {
        Process::Compton(self.model, self.mode).into()
    }

    fn __call__(
        &self,
        py: Python,
        energy_in: Float,
        energy_out: ArrayOrFloat
    ) -> Result<PyObject> {
        let record: PyRef<PyMaterialRecord> = self.material.extract(py)?;
        let table = Self::table(py, &record, self.model, self.mode)?.unwrap();
        let result: PyObject = match energy_out {
            ArrayOrFloat::Array(energy) => {
                let result = PyArray::<Float>::empty(py, &energy.shape())?;
                let n = energy.size();
                for i in 0..n {
                    let v = table.interpolate(
                        energy_in,
                        energy.get(i)?,
                    );
                    result.set(i, v)?;
                }
                result.into_py(py)
            },
            ArrayOrFloat::Float(energy) => table.interpolate(energy_in, energy).into_py(py),
        };
        Ok(result)
    }

    fn energies_out(
        &self,
        py: Python,
        i: usize,
    ) -> Result<PyObject> {
        let record: PyRef<PyMaterialRecord> = self.material.extract(py)?;
        let table = Self::table(py, &record, self.model, self.mode)?.unwrap();
        let (_, m) = table.shape();
        let energies = (0..m).map(|j| table.energy_out(i, j));
        let array = copy1(py, m, energies)?;
        Ok(array)
    }
}


// ===============================================================================================
// Python wrapper for a Compton InverseCDF object.
// ===============================================================================================

#[pyclass(name = "InverseDistribution", module="goupil")]
pub struct PyInverseDistribution {
    #[pyo3(get)]
    cdf: PyObject,
    #[pyo3(get)]
    energies: PyObject,
    #[pyo3(get)]
    material: PyObject,
    #[pyo3(get)]
    values: PyObject,
    #[pyo3(get)]
    weights: PyObject,

    model: ComptonModel,
    mode: ComptonMode,
}

impl PyInverseDistribution {
    fn new(
        record: &Bound<PyMaterialRecord>,
        model: ComptonModel,
        mode: ComptonMode,
    ) -> Result<PyObject> {
        let py = record.py();
        let (energies, cdf, values, weights) = {
            match Self::table(py, &record.borrow(), model, mode)? {
                None => return Ok(py.None()),
                Some(table) => {
                    let energies = readonly1(table.energies.as_ref(), record)?;
                    let cdf = copy1(py, table.cdf.len(), table.cdf.iter())?;
                    let values = readonly2(table.values.as_ref(), table.shape(), record)?;
                    let weights = match &table.weights {
                        None => py.None(),
                        Some(table) => readonly2(table.as_ref(), table.shape(), record)?,
                    };
                    (energies, cdf, values, weights)
                },
            }
        };
        let material: PyObject = record.into_py(py);
        let result = Self { cdf, energies, material, values, weights, model, mode };
        Ok(result.into_py(py))
    }

    fn table<'py>(
        py: Python<'py>,
        record: &PyMaterialRecord,
        model: ComptonModel,
        mode: ComptonMode
    ) -> Result<Option<&'py ComptonInverseCDF>> {
        Ok(record.get(py)?.table.compton.get(model).get(mode).inverse_cdf.as_ref())
    }
}

#[pymethods]
impl PyInverseDistribution {
    // Implementation of GC protocol.
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.cdf)?;
        visit.call(&self.energies)?;
        visit.call(&self.material)?;
        visit.call(&self.values)?;
        visit.call(&self.weights)?;
        Ok(())
    }

    #[getter]
    fn get_process(&self) -> String {
        Process::Compton(self.model, self.mode).into()
    }

    fn __call__(
        &self,
        py: Python,
        energy: Float,
        cdf: ArrayOrFloat
    ) -> Result<PyObject> {
        let record: PyRef<PyMaterialRecord> = self.material.extract(py)?;
        let table = Self::table(py, &record, self.model, self.mode)?.unwrap();
        let result = match cdf {
            ArrayOrFloat::Array(cdf) => match &table.weights {
                None => {
                    let values = PyArray::<Float>::empty(py, &cdf.shape())?;
                    let n = cdf.size();
                    for i in 0..n {
                        let v = table.interpolate(energy, cdf.get(i)?);
                        values.set(i, v.0)?;
                    }
                    values.into_py(py)
                },
                Some(_) => {
                    let values = PyArray::<Float>::empty(py, &cdf.shape())?;
                    let weights = PyArray::<Float>::empty(py, &cdf.shape())?;
                    let n = cdf.size();
                    for i in 0..n {
                        let v = table.interpolate(energy, cdf.get(i)?);
                        values.set(i, v.0)?;
                        weights.set(i, v.1)?;
                    }
                    let values: PyObject = values.into_py(py);
                    let weights: PyObject = weights.into_py(py);
                    (values, weights).into_py(py)
                },
            },
            ArrayOrFloat::Float(cdf) => {
                let v = table.interpolate(energy, cdf);
                match &table.weights {
                    None => v.0.into_py(py),
                    Some(_) => v.into_py(py),
                }
            },
        };
        Ok(result)
    }
}


// ===============================================================================================
// Python wrapper for a Rayleigh FormFactor object.
// ===============================================================================================

#[pyclass(name = "FormFactor", module = "goupil")]
pub struct PyFormFactor {
    #[pyo3(get)]
    material: PyObject,
    #[pyo3(get)]
    momenta: PyObject,
    #[pyo3(get)]
    values: PyObject,
}

impl PyFormFactor {
    fn new(record: &Bound<PyMaterialRecord>) -> Result<PyObject> {
        let py = record.py();
        let (momenta, values) = {
            match Self::table(py, &record.borrow())? {
                None => return Ok(py.None()),
                Some(table) => {
                    let momenta = readonly1(table.momenta.as_ref(), record)?;
                    let values = readonly1(table.values.as_ref(), record)?;
                    (momenta, values)
                },
            }
        };
        let material: PyObject = record.into_py(py);
        let result = Self { material, momenta, values };
        Ok(result.into_py(py))
    }

    fn table<'py>(
        py: Python<'py>,
        record: &PyMaterialRecord,
    ) -> Result<Option<&'py RayleighFormFactor>> {
        Ok(record.get(py)?.table.rayleigh.form_factor.as_ref())
    }
}

#[pymethods]
impl PyFormFactor {
    // Implementation of GC protocol.
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.material)?;
        visit.call(&self.momenta)?;
        visit.call(&self.values)?;
        Ok(())
    }

    #[getter]
    fn get_process(&self) -> String {
        "Rayleigh".to_string()
    }

    fn __call__(
        &self,
        py: Python,
        momentum: ArrayOrFloat
    ) -> Result<PyObject> {
        let record: PyRef<PyMaterialRecord> = self.material.extract(py)?;
        let table = Self::table(py, &record)?.unwrap();
        let result = match momentum {
            ArrayOrFloat::Array(momentum) => {
                let result = PyArray::<Float>::empty(py, &momentum.shape())?;
                let n = momentum.size();
                for i in 0..n {
                    let v = table.interpolate(momentum.get(i)?);
                    result.set(i, v)?;
                }
                result.into_py(py)
            },
            ArrayOrFloat::Float(momentum) => table.interpolate(momentum).into_py(py),
        };
        Ok(result)
    }
}


// ===============================================================================================
// Some routines for wrapping Float data as numpy arrays.
// ===============================================================================================

fn readonly1(data: &[Float], owner: &Bound<PyAny>) -> Result<PyObject> {
    let array = PyArray::from_data(
        owner.py(),
        data,
        owner,
        PyArrayFlags::ReadOnly,
        None
    )?;
    Ok(array.into_any().unbind())
}

fn readonly2(data: &[Float], shape: (usize, usize), owner: &Bound<PyAny>) -> Result<PyObject> {
    let shape: [usize; 2] = shape.into();
    let array = PyArray::from_data(
        owner.py(),
        data,
        owner,
        PyArrayFlags::ReadOnly,
        Some(&shape),
    )?;
    Ok(array.into_any().unbind())
}

fn copy1<I>(py: Python, n: usize, iter: I) -> Result<PyObject>
where
    I: Iterator<Item=Float>,
{
    let array = PyArray::<Float>::from_iter(py, &[n], iter)?;
    array.readonly();
    Ok(array.into_any().unbind())
}
