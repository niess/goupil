use anyhow::{anyhow, bail, Error, Result};
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
    sync::GILOnceCell,
    types::{PyBytes, PyTuple},
};
use regex::Regex;
use rmp_serde::{Deserializer, Serializer};
use serde::{Deserialize, Serialize};
use super::{
    elements::PyAtomicElement,
    macros::value_error,
    numpy::{PyArray, PyArrayFlags},
    prefix,
    transport::PyTransportSettings,
};
use std::collections::HashMap;


// ===============================================================================================
// Python wrapper for a material definition
// ===============================================================================================

#[pyclass(name = "MaterialDefinition", module = "goupil")]
pub struct PyMaterialDefinition (pub MaterialDefinition);

#[derive(FromPyObject)]
enum Element<'py> {
    Object(PyRef<'py, PyAtomicElement>),
    Symbol(&'py str),
}

impl<'py> Element<'py> {
    fn get(&self) -> Result<&'static AtomicElement> {
        let element = match self {
            Element::Object(element) => element.0,
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
enum PyMassComposition<'py> {
    Atomic(Vec<(Float, Element<'py>)>),
    Material(Vec<(Float, PyRef<'py, PyMaterialDefinition>)>),
}

type PyMoleComposition<'py> = Vec<(Float, Element<'py>)>;

#[pymethods]
impl PyMaterialDefinition {
    #[new]
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
                        let re = Regex::new(r"([A-Z][a-z]?)([0-9]*)")?;
                        let mut composition = Vec::<(Float, &AtomicElement)>::default();
                        for captures in re.captures_iter(name) {
                            let symbol = captures.get(1).unwrap().as_str();
                            let element = AtomicElement::from_symbol(symbol)?;
                            let weight = captures.get(2).unwrap().as_str();
                            let weight: Float = if weight.len() == 0 {
                                1.0
                            } else {
                                weight.parse::<Float>()?
                            };
                            composition.push((weight, element));
                        }
                        MaterialDefinition::from_mole(name, &composition)
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
                            let composition: Vec<_> = composition
                                .iter()
                                .map(|(weight, material)| (*weight, &material.0))
                                .collect();
                            MaterialDefinition::from_others(name, &composition)
                        },
                    }
                },
            }
        };
        Ok(Self(definition))
    }

    // Implementation of pickling protocol.
    pub fn __setstate__(&mut self, state: &PyBytes) -> Result<()> {
        self.0 = Deserialize::deserialize(&mut Deserializer::new(state.as_bytes()))?;
        Ok(())
    }

    fn __getstate__<'py>(&self, py: Python<'py>) -> Result<&'py PyBytes> {
        let mut buffer = Vec::new();
        self.0.serialize(&mut Serializer::new(&mut buffer))?;
        Ok(PyBytes::new(py, &buffer))
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
    fn get_mass_composition<'p>(&self, py: Python<'p>) -> &'p PyTuple {
        let composition: Vec<_> = self.0
            .mass_composition()
            .iter()
            .map(|(weight, element)| (*weight, PyAtomicElement(*element).into_py(py)))
            .collect();
        PyTuple::new(py, composition)
    }

    #[getter]
    fn get_mole_composition<'p>(&self, py: Python<'p>) -> &'p PyTuple {
        let composition: Vec<_> = self.0
            .mole_composition()
            .iter()
            .map(|(weight, element)| (*weight, PyAtomicElement(*element).into_py(py)))
            .collect();
        PyTuple::new(py, composition)
    }

    #[getter]
    fn get_name(&self) -> &str {
        self.0.name()
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
    pub fn new(definitions: Option<Vec<PyRef<PyMaterialDefinition>>>) -> Result<Self> {
        let mut registry = MaterialRegistry::default();
        if let Some(definitions) = definitions {
            for definition in definitions.iter() {
                registry.add(&definition.0)?;
            }
        };
        Ok(Self{
            inner: registry,
            proxies: HashMap::<String, Py<PyMaterialRecord>>::default(),
        })
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
        let mut registry = slf.as_ref(py).borrow_mut();
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

    // Implementation of pickling protocol.
    pub fn __setstate__(&mut self, py: Python, state: &PyBytes) -> Result<()> {
        // Detach pending record(s).
        for (k, v) in self.proxies.drain() {
            let record = self.inner.remove(&k).unwrap();
            Self::into_owned(py, v, record);
        }

        // Update inner registry.
        self.inner = Deserialize::deserialize(&mut Deserializer::new(state.as_bytes()))?;
        Ok(())
    }

    fn __getstate__<'py>(&self, py: Python<'py>) -> Result<&'py PyBytes> {
        let mut buffer = Vec::new();
        self.inner.serialize(&mut Serializer::new(&mut buffer))?;
        Ok(PyBytes::new(py, &buffer))
    }

    // Direct public interface.
    fn add(&mut self, definition: &PyMaterialDefinition) -> Result<()> {
        self.inner.add(&definition.0)?;
        Ok(())
    }

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
        constraint: Option<PyConstraint>,
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
            .or(settings.map(|settings| Ok(settings.0.mode)));
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
            .or(settings.map(|settings| Ok(settings.0.absorption)))
            .unwrap_or(Ok(Discrete))?;

        config.compton_model = compton_model
            .map(|s| ComptonModel::try_from(s))
            .or(settings.map(|settings| Ok(settings.0.compton_model)))
            .unwrap_or(Ok(ScatteringFunction))?;

        config.energy_max = energy_max
            .or(settings.and_then(|settings| settings.0.energy_max));

        config.energy_min = energy_min
            .or(settings.and_then(|settings| settings.0.energy_min));

        config.compton_method = compton_method
            .map(|s| ComptonMethod::try_from(s))
            .or(settings.map(|settings| Ok(settings.0.compton_method)))
            .unwrap_or(Ok(RejectionSampling))?;

        config.compton_mode = compton_mode
            .map(|s| ComptonMode::try_from(s))
            .or(settings.map(|settings| Ok(settings.0.compton_mode)))
            .unwrap_or(Ok(Direct))?;

        match &config.compton_mode {
            Adjoint | Inverse => match mode {
                None => config.mode = Backward,
                Some(mode) => if let Forward = mode {
                    bail!(
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
                    bail!(
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
                None => settings.and_then(|settings| settings.0.constraint),
                Some(constraint) => match constraint {
                    PyConstraint::Bool(b) => if b { Some(1.0) } else { None },
                    PyConstraint::Float(f) => Some(f),
                }
            },
        };

        self.inner.compute(
            &config,
            length,
            width,
            precision,
        )
    }

    fn load_elements(&mut self, py: Python, path: Option<String>) -> Result<()> {
        let path = match path {
            None => {
                let mut path = prefix(py)?.clone();
                path.push(Self::ELEMENTS_DATA);
                path
            },
            Some(path) => path.into(),
        };
        self.inner.load_elements(&path)?;
        Ok(())
    }
}

#[derive(Clone, FromPyObject)]
pub enum PyConstraint {
    Bool(bool),
    Float(Float),
}

impl PyMaterialRegistry {
    pub(crate) const ELEMENTS_DATA: &str = "data/elements";

    // Transforms a borrowed material record to an owned one.
    fn into_owned(py: Python, proxy: Py<PyMaterialRecord>, record: MaterialRecord) {
        if proxy.get_refcnt(py) >= 2 {
            let mut proxy = proxy.as_ref(py).borrow_mut();
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
                .into_ref(py)
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
                let object: PyObject = match self.get(py)?.electrons() {
                    None => py.None(),
                    Some(electrons) => {
                        let electrons = PyElectronicStructure::new(electrons.clone(), false)?;
                        let electrons = Py::new(py, electrons)?;
                        self.electrons = Some(electrons.clone());
                        electrons.into_py(py)
                    },
                };
                Ok(object)
            },
            Some(electrons) => Ok(electrons.clone().into_py(py)),
        }
    }

    fn absorption_cross_section(this: &PyCell<Self>, py: Python) -> Result<PyObject> {
        PyAbsorptionCrossSection::new(py, this)
    }

    fn compton_cdf(
        this: &PyCell<Self>,
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
        PyComptonCDF::new(py, this, model, mode)
    }

    fn compton_cross_section(
        this: &PyCell<Self>,
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
        PyComptonCrossSection::new(py, this, model, mode)
    }

    fn compton_inverse_cdf(
        this: &PyCell<Self>,
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
        PyComptonInverseCDF::new(py, this, model, mode)
    }

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

    fn rayleigh_cross_section(
        this: &PyCell<Self>,
        py: Python,
    ) -> Result<PyObject> {
        PyRayleighCrossSection::new(py, this)
    }

    fn rayleigh_form_factor(
        this: &PyCell<Self>,
        py: Python,
    ) -> Result<PyObject> {
        PyRayleighFormFactor::new(py, this)
    }
}


// ===============================================================================================
// Python wrapper for an ElectronicStructure object.
// ===============================================================================================

#[pyclass(name = "ElectronicStructure", module = "goupil")]
pub struct PyElectronicStructure {
    electrons: ElectronicStructure,
    writable: bool,
    shells: GILOnceCell<PyObject>,
}

impl PyElectronicStructure {
    fn new(electrons: ElectronicStructure, writable: bool) -> Result<Self> {
        Ok(Self {
            electrons,
            writable,
            shells: GILOnceCell::new(),
        })
    }
}

#[pymethods]
impl PyElectronicStructure {
    #[getter]
    fn get_charge(&self) -> Float {
        self.electrons.charge()
    }

    #[getter]
    fn get_shells(slf: &PyCell<Self>, py: Python) -> Result<PyObject> {
        let obj = slf.borrow();
        let shells = obj.shells.get_or_try_init(py, || {
            let flags = if obj.writable {
                PyArrayFlags::ReadWrite
            } else {
                PyArrayFlags::ReadOnly
            };
            let shells: &PyAny = PyArray::from_data(
                py,
                &obj.electrons,
                slf,
                flags,
                None,
            )?;
            Ok::<Py<PyAny>, Error>(shells.into())
        })?;
        Ok(shells.clone_ref(py))
    }

    #[staticmethod]
    fn from_others(composition: Vec<(Float, PyRef<Self>)>) -> Result<Self> {
        let composition: Vec<_> = composition
            .iter()
            .map(|(weight, electrons)| (*weight, &electrons.electrons))
            .collect();
        let electrons = ElectronicStructure::from_others(&composition);
        Ok(Self::new(electrons, true)?)
    }
}


// ===============================================================================================
// Python wrapper for an AbsorptionCrossSection object.
// ===============================================================================================

#[pyclass(name = "AbsorptionCrossSection", module = "goupil")]
pub struct PyAbsorptionCrossSection {
    #[pyo3(get)]
    record: PyObject,
    #[pyo3(get)]
    energies: PyObject,
    #[pyo3(get)]
    values: PyObject,
}

impl PyAbsorptionCrossSection {
    fn new(
        py: Python,
        record: &PyCell<PyMaterialRecord>,
    ) -> Result<PyObject> {
        let (energies, values) = {
            match Self::table(py, &record.borrow())? {
                None => return Ok(py.None()),
                Some(table) => {
                    let energies = readonly1(table.energies.as_ref(), record)?;
                    let values = readonly1(table.values.as_ref(), record)?;
                    (energies, values)
                },
            }
        };
        let record: PyObject = record.into_py(py);
        let result = Self { record, energies, values };
        Ok(result.into_py(py))
    }

    fn table<'py>(
        py: Python<'py>,
        record: &PyMaterialRecord
    ) -> Result<Option<&'py AbsorptionCrossSection>> {
        Ok(record.get(py)?.table.absorption.as_ref())
    }
}

#[pymethods]
impl PyAbsorptionCrossSection {
    fn __call__(
        &self,
        py: Python,
        energy: Float
    ) -> Result<Float> {
        let record: PyRef<PyMaterialRecord> = self.record.extract(py)?;
        let table = Self::table(py, &record)?.unwrap();
        Ok(table.interpolate(energy))
    }
}


// ===============================================================================================
// Python wrapper for a ComptonCrossSection object.
// ===============================================================================================

#[pyclass(name = "ComptonCrossSection", module="goupil")]
pub struct PyComptonCrossSection {
    #[pyo3(get)]
    record: PyObject,
    #[pyo3(get)]
    energies: PyObject,
    #[pyo3(get)]
    values: PyObject,

    model: ComptonModel,
    mode: ComptonMode,
}

impl PyComptonCrossSection {
    fn new(
        py: Python,
        record: &PyCell<PyMaterialRecord>,
        model: ComptonModel,
        mode: ComptonMode,
    ) -> Result<PyObject> {
        let (energies, values) = {
            match Self::table(py, &record.borrow(), model, mode)? {
                None => return Ok(py.None()),
                Some(table) => {
                    let energies = readonly1(table.energies.as_ref(), record)?;
                    let values = readonly1(table.values.as_ref(), record)?;
                    (energies, values)
                },
            }
        };
        let record: PyObject = record.into_py(py);
        let result = Self { record, energies, values, model, mode };
        Ok(result.into_py(py))
    }

    fn table<'py>(
        py: Python<'py>,
        record: &PyMaterialRecord,
        model: ComptonModel,
        mode: ComptonMode
    ) -> Result<Option<&'py ComptonCrossSection>> {
        Ok(record.get(py)?.table.compton.get(model).get(mode).cross_section.as_ref())
    }
}

#[pymethods]
impl PyComptonCrossSection {
    fn __call__(
        &self,
        py: Python,
        energy: Float
    ) -> Result<Float> {
        let record: PyRef<PyMaterialRecord> = self.record.extract(py)?;
        let table = Self::table(py, &record, self.model, self.mode)?.unwrap();
        Ok(table.interpolate(energy))
    }

    #[getter]
    fn mode(&self) -> &str {
        self.mode.into()
    }

    #[getter]
    fn model(&self) -> &str {
        self.model.into()
    }
}


// ===============================================================================================
// Python wrapper for a ComptonCDF object.
// ===============================================================================================

#[pyclass(name = "ComptonCDF", module="goupil")]
pub struct PyComptonCDF {
    #[pyo3(get)]
    record: PyObject,
    #[pyo3(get)]
    energies_in: PyObject,
    #[pyo3(get)]
    x: PyObject,
    #[pyo3(get)]
    values: PyObject,

    model: ComptonModel,
    mode: ComptonMode,
}

impl PyComptonCDF {
    fn new(
        py: Python,
        record: &PyCell<PyMaterialRecord>,
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
        let record: PyObject = record.into_py(py);
        let result = Self { record, energies_in, x, values, model, mode };
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
impl PyComptonCDF {
    fn __call__(
        &self,
        py: Python,
        energy_in: Float,
        energy_out: Float
    ) -> Result<Float> {
        let record: PyRef<PyMaterialRecord> = self.record.extract(py)?;
        let table = Self::table(py, &record, self.model, self.mode)?.unwrap();
        Ok(table.interpolate(energy_in, energy_out))
    }

    fn energies_out(
        &self,
        py: Python,
        i: usize,
    ) -> Result<PyObject> {
        let record: PyRef<PyMaterialRecord> = self.record.extract(py)?;
        let table = Self::table(py, &record, self.model, self.mode)?.unwrap();
        let (_, m) = table.shape();
        let energies = (0..m).map(|j| table.energy_out(i, j));
        let array = copy1(py, m, energies)?;
        Ok(array)
    }

    #[getter]
    fn mode(&self) -> &str {
        self.mode.into()
    }

    #[getter]
    fn model(&self) -> &str {
        self.model.into()
    }
}


// ===============================================================================================
// Python wrapper for a ComptonInverseCDF object.
// ===============================================================================================

#[pyclass(name = "ComptonInverseCDF", module="goupil")]
pub struct PyComptonInverseCDF {
    #[pyo3(get)]
    record: PyObject,
    #[pyo3(get)]
    energies: PyObject,
    #[pyo3(get)]
    cdf: PyObject,
    #[pyo3(get)]
    values: PyObject,
    #[pyo3(get)]
    weights: PyObject,

    model: ComptonModel,
    mode: ComptonMode,
}

impl PyComptonInverseCDF {
    fn new(
        py: Python,
        record: &PyCell<PyMaterialRecord>,
        model: ComptonModel,
        mode: ComptonMode,
    ) -> Result<PyObject> {
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
        let record: PyObject = record.into_py(py);
        let result = Self { record, energies, cdf, values, weights, model, mode };
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
impl PyComptonInverseCDF {
    fn __call__(
        &self,
        py: Python,
        energy: Float,
        cdf: Float
    ) -> Result<PyObject> {
        let record: PyRef<PyMaterialRecord> = self.record.extract(py)?;
        let table = Self::table(py, &record, self.model, self.mode)?.unwrap();
        let result = table.interpolate(energy, cdf);
        let result = match &table.weights {
            None => result.0.into_py(py),
            Some(_) => result.into_py(py),
        };
        Ok(result)
    }

    #[getter]
    fn mode(&self) -> &str {
        self.mode.into()
    }

    #[getter]
    fn model(&self) -> &str {
        self.model.into()
    }
}


// ===============================================================================================
// Python wrapper for a RayleighCrossSection object.
// ===============================================================================================

#[pyclass(name = "RayleighCrossSection", module = "goupil")]
pub struct PyRayleighCrossSection {
    #[pyo3(get)]
    record: PyObject,
    #[pyo3(get)]
    energies: PyObject,
    #[pyo3(get)]
    values: PyObject,
}

impl PyRayleighCrossSection {
    fn new(
        py: Python,
        record: &PyCell<PyMaterialRecord>,
    ) -> Result<PyObject> {
        let (energies, values) = {
            match Self::table(py, &record.borrow())? {
                None => return Ok(py.None()),
                Some(table) => {
                    let energies = readonly1(table.energies.as_ref(), record)?;
                    let values = readonly1(table.values.as_ref(), record)?;
                    (energies, values)
                },
            }
        };
        let record: PyObject = record.into_py(py);
        let result = Self { record, energies, values };
        Ok(result.into_py(py))
    }

    fn table<'py>(
        py: Python<'py>,
        record: &PyMaterialRecord,
    ) -> Result<Option<&'py RayleighCrossSection>> {
        Ok(record.get(py)?.table.rayleigh.cross_section.as_ref())
    }
}

#[pymethods]
impl PyRayleighCrossSection {
    fn __call__(
        &self,
        py: Python,
        energy: Float
    ) -> Result<Float> {
        let record: PyRef<PyMaterialRecord> = self.record.extract(py)?;
        let table = Self::table(py, &record)?.unwrap();
        Ok(table.interpolate(energy))
    }
}


// ===============================================================================================
// Python wrapper for a RayleighFormFactor object.
// ===============================================================================================

#[pyclass(name = "RayleighFormFactor", module = "goupil")]
pub struct PyRayleighFormFactor {
    #[pyo3(get)]
    record: PyObject,
    #[pyo3(get)]
    momenta: PyObject,
    #[pyo3(get)]
    values: PyObject,
}

impl PyRayleighFormFactor {
    fn new(
        py: Python,
        record: &PyCell<PyMaterialRecord>,
    ) -> Result<PyObject> {
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
        let record: PyObject = record.into_py(py);
        let result = Self { record, momenta, values };
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
impl PyRayleighFormFactor {
    fn __call__(
        &self,
        py: Python,
        momentum: Float // XXX Vectorize these functions.
    ) -> Result<Float> {
        let record: PyRef<PyMaterialRecord> = self.record.extract(py)?;
        let table = Self::table(py, &record)?.unwrap();
        Ok(table.interpolate(momentum))
    }
}


// ===============================================================================================
// Some routines for wrapping Float data as numpy arrays.
// ===============================================================================================

fn readonly1(data: &[Float], owner: &PyAny) -> Result<PyObject> {
    let array: &PyAny = PyArray::from_data(
        owner.py(),
        data,
        owner,
        PyArrayFlags::ReadOnly,
        None
    )?;
    Ok(array.into())
}

fn readonly2(data: &[Float], shape: (usize, usize), owner: &PyAny) -> Result<PyObject> {
    let shape: [usize; 2] = shape.into();
    let array: &PyAny = PyArray::from_data(
        owner.py(),
        data,
        owner,
        PyArrayFlags::ReadOnly,
        Some(&shape),
    )?;
    Ok(array.into())
}

fn copy1<I>(py: Python, n: usize, iter: I) -> Result<PyObject>
where
    I: Iterator<Item=Float>,
{
    let array = PyArray::<Float>::from_iter(py, &[n], iter)?;
    array.readonly();
    let array: &PyAny = array;
    Ok(array.into())
}
