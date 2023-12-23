use anyhow::{bail, Result};
use crate::numerics::float::Float;
use crate::physics::materials::MaterialRegistry;
use crate::physics::process::{
    absorption::AbsorptionMode,
    compton::{ComptonModel, ComptonMethod, ComptonMode::{self, Adjoint, Direct, Inverse}},
    rayleigh::RayleighMode,
};
use crate::transport::{
    agent::{TransportAgent, TransportBoundary, TransportStatus},
    geometry::{ExternalTracer, GeometryDefinition, GeometryTracer, SimpleTracer},
    PhotonState,
    TransportMode::{self, Backward, Forward},
    TransportSettings,
};
use pyo3::{
    prelude::*,
    types::{PyBytes, PyDict, PyString},
};
use rmp_serde::{Deserializer, Serializer};
use serde::{Deserialize, Serialize};
use super::{
    ctrlc_catched,
    geometry::{PyExternalGeometry, PyGeometryDefinition},
    macros::type_error,
    materials::PyMaterialRegistry,
    numpy::{Dtype, PyArray, PyScalar},
    rand::PyRandomStream,
    prefix,
};


// ===============================================================================================
// Python wrapper for a Goupil Monte Carlo engine.
// ===============================================================================================

#[pyclass(name = "TransportSettings", module = "goupil")]
pub(crate) struct PyTransportSettings (pub TransportSettings);

// Convert from an optional string.
macro_rules! from_optstr {
    ($type:ty, $var:expr, $value:expr) => {
        $var = match $value {
            None => <$type>::None,
            Some(s) => <$type>::try_from(s)?,
        };
    }
}

// Convert to an optional string.
macro_rules! to_optstr {
    ($type:ty, $var:expr) => {
        match $var {
            <$type>::None => None,
            _ => Some($var.into()),
        }
    }
}

#[pymethods]
impl PyTransportSettings {
    #[new]
    fn new() -> Self {
        Self(TransportSettings::default())
    }

    #[getter]
    fn get_mode(&self) -> &str {
        self.0.mode.into()
    }

    #[setter]
    fn set_mode(&mut self, value: &str) -> Result<()> {
        self.0.mode = TransportMode::try_from(value)?;
        match self.0.mode {
            Backward => match self.0.compton_mode {
                Direct => {
                    self.0.compton_mode = Adjoint;
                },
                _ => (),
            },
            Forward => match self.0.compton_mode {
                Adjoint | Inverse => {
                    self.0.compton_mode = Direct;
                },
                _ => (),
            },
        }
        Ok(())
    }

    #[getter]
    fn get_absorption(&self) -> Option<&str> {
        to_optstr!(AbsorptionMode, self.0.absorption)
    }

    #[setter]
    fn set_absorption(&mut self, value: Option<&str>) -> Result<()> {
        from_optstr!(AbsorptionMode, self.0.absorption, value);
        Ok(())
    }

    #[getter]
    fn get_boundary(&self) -> Option<usize> {
        match self.0.boundary {
            TransportBoundary::None => None,
            TransportBoundary::Sector(index) => Some(index),
        }
    }

    #[setter]
    fn set_boundary(&mut self, value: Option<usize>) -> Result<()> {
        match value {
            None => self.0.boundary = TransportBoundary::None,
            Some(index) => self.0.boundary = TransportBoundary::Sector(index),
        };
        Ok(())
    }

    #[getter]
    fn get_compton_method(&self) -> &str {
        self.0.compton_method.into()
    }

    #[setter]
    fn set_compton_method(&mut self, value: &str) -> Result<()> {
        self.0.compton_method = ComptonMethod::try_from(value)?;
        Ok(())
    }

    #[getter]
    fn get_compton_mode(&self) -> Option<&str> {
        to_optstr!(ComptonMode, self.0.compton_mode)
    }

    #[setter]
    fn set_compton_mode(&mut self, value: Option<&str>) -> Result<()> {
        from_optstr!(ComptonMode, self.0.compton_mode, value);
        match self.0.compton_mode {
            Adjoint => {
                self.0.mode = Backward;
            },
            Direct => {
                self.0.mode = Forward;
            },
            Inverse => {
                self.0.mode = Backward;
                self.0.compton_method = ComptonMethod::InverseTransform;
            },
            ComptonMode::None => (),
        }
        Ok(())
    }

    #[getter]
    fn get_compton_model(&self) -> &str {
        self.0.compton_model.into()
    }

    #[setter]
    fn set_compton_model(&mut self, value: &str) -> Result<()> {
        self.0.compton_model = ComptonModel::try_from(value)?;
        Ok(())
    }

    #[getter]
    fn get_constraint(&self) -> Option<Float> {
        self.0.constraint
    }

    #[setter]
    fn set_constraint(&mut self, value: Option<Float>) -> Result<()> {
        self.0.constraint = value;
        Ok(())
    }

    #[getter]
    fn get_rayleigh(&self) -> Option<&str> {
        to_optstr!(RayleighMode, self.0.rayleigh)
    }

    #[setter]
    fn set_rayleigh(&mut self, value: Option<&str>) -> Result<()> {
        from_optstr!(RayleighMode, self.0.rayleigh, value);
        Ok(())
    }

    #[getter]
    fn get_energy_min(&self) -> Option<Float> {
        self.0.energy_min
    }

    #[setter]
    fn set_energy_min(&mut self, value: Option<Float>) -> Result<()> {
        self.0.energy_min = value;
        Ok(())
    }

    #[getter]
    fn get_energy_max(&self) -> Option<Float> {
        self.0.energy_max
    }

    #[setter]
    fn set_energy_max(&mut self, value: Option<Float>) -> Result<()> {
        self.0.energy_max = value;
        Ok(())
    }

    #[getter]
    fn get_length_max(&self) -> Option<Float> {
        self.0.length_max
    }

    #[setter]
    fn set_length_max(&mut self, value: Option<Float>) -> Result<()> {
        self.0.length_max = value;
        Ok(())
    }
}


// ===============================================================================================
// Main transport engine.
// ===============================================================================================

#[pyclass(name = "TransportEngine", module = "goupil")]
pub struct PyTransportEngine {
    #[pyo3(get)]
    geometry: Option<PyGeometryDefinition>,
    #[pyo3(get)]
    random: Py<PyRandomStream>,
    #[pyo3(get)]
    registry: Py<PyMaterialRegistry>,
    #[pyo3(get)]
    settings: Py<PyTransportSettings>,
}

#[derive(FromPyObject)]
enum GeometryArg {
    Object(PyGeometryDefinition),
    Path(String),
}

#[pymethods]
impl PyTransportEngine {
    #[new]
    fn new(
        py: Python,
        geometry: Option<GeometryArg>,
        random: Option<Py<PyRandomStream>>,
        registry: Option<Py<PyMaterialRegistry>>,
        settings: Option<Py<PyTransportSettings>>,
    ) -> Result<Self> {
        let geometry = match geometry {
            None => None,
            Some(geometry) => {
                let geometry = match geometry {
                    GeometryArg::Object(geometry) => geometry,
                    GeometryArg::Path(path) => {
                        let external = PyExternalGeometry::new(&path)?;
                        let external = Py::new(py, external)?;
                        PyGeometryDefinition::External(external)
                    },
                };
                Some(geometry)
            },
        };
        let random: Py<PyRandomStream> = match random {
            None => Py::new(py, PyRandomStream::new(None)?)?,
            Some(random) => random.into(),
        };
        let registry: Py<PyMaterialRegistry> = match registry {
            None => Py::new(py, PyMaterialRegistry::new(vec![])?)?,
            Some(registry) => registry.into(),
        };
        let settings: Py<PyTransportSettings> = match settings {
            None => Py::new(py, PyTransportSettings::new())?,
            Some(settings) => settings.into(),
        };
        Ok(Self { geometry, random, registry, settings })
    }

    fn __getattr__(&self, py: Python, name: &PyString) -> Result<PyObject> {
        Ok(self.settings.getattr(py, name)?)
    }

    fn __setattr__(&mut self, py: Python, name: &str, value: PyObject) -> Result<()> {
        match name {
            "geometry" => {
                if value.is_none(py) {
                    self.geometry = None;
                } else {
                    let geometry: PyGeometryDefinition = value.extract(py)?;
                    self.geometry = Some(geometry);
                }
            },
            "random" => self.random = value.extract(py)?,
            "registry" => self.registry = value.extract(py)?,
            "settings" => self.settings = value.extract(py)?,
            _ => self.settings.setattr(py, name, value)?,
        }
        Ok(())
    }

    // Implementation of pickling protocol.
    pub fn __setstate__(&mut self, py: Python, state: &PyBytes) -> Result<()> {
        let mut deserializer = Deserializer::new(state.as_bytes());

        let mut random = self.random.borrow_mut(py);
        *random = Deserialize::deserialize(&mut deserializer)?;

        let registry = &mut self.registry.borrow_mut(py).inner;
        *registry = Deserialize::deserialize(&mut deserializer)?;

        let settings = &mut self.settings.borrow_mut(py).0;
        *settings = Deserialize::deserialize(&mut deserializer)?;

        Ok(())
    }

    fn __getstate__<'py>(&self, py: Python<'py>) -> Result<&'py PyBytes> {
        let mut buffer = Vec::new();
        let mut serializer = Serializer::new(&mut buffer);

        let random = &self.random.borrow(py);
        random.serialize(&mut serializer)?;

        let registry = &self.registry.borrow(py).inner;
        registry.serialize(&mut serializer)?;

        let settings = &self.settings.borrow(py).0;
        settings.serialize(&mut serializer)?;

        Ok(PyBytes::new(py, &buffer))
    }

    #[pyo3(signature = (mode=None, atomic_data=None, **kwargs))]
    fn compile(
        &self,
        py: Python,
        mode: Option<&str>,
        atomic_data: Option<&str>,
        kwargs: Option<&PyDict>,
    ) -> Result<()> {
        enum CompileMode {
            All,
            Backward,
            Both,
            Forward,
        }

        let mode = match mode {
            None => match &self.settings.borrow(py).0.mode {
                TransportMode::Backward => CompileMode::Backward,
                TransportMode::Forward => CompileMode::Forward,
            },
            Some(mode) => match mode {
                "All" => CompileMode::All,
                "Backward" => CompileMode::Backward,
                "Both" => CompileMode::Both,
                "Forward" => CompileMode::Forward,
                _ => bail!(
                    "bad mode (expected 'All', 'Backward', 'Both' or 'Forward', found '{}')",
                    mode,
                ),
            }
        };

        {
            // Fetch material registry. Note that we scope this mutable borrow (see below).
            let registry = &mut self.registry.borrow_mut(py).inner;

            // Add current geometry materials to the registry.
            if let Some(geometry) = &self.geometry {
                match geometry {
                    PyGeometryDefinition::External(external) => {
                        self.update_with(&external.borrow(py).0, registry)?
                    },
                    PyGeometryDefinition::Simple(simple) => {
                        self.update_with(&simple.borrow(py).0, registry)?
                    },
                }
            }

            // Load atomic data.
            match atomic_data {
                None => if !registry.atomic_data_loaded() {
                    let mut path = prefix(py)?.clone();
                    path.push(PyMaterialRegistry::ELEMENTS_DATA);
                    registry.load_elements(&path)?;
                },
                Some(path) => registry.load_elements(&path)?,
            }
        }

        // Call the registry compute method through Python. This let us use keyword arguments,
        // thus avoiding to duplicate the registry.compute signature. However, we first need to
        // release the mutable borrow on the registry.
        match mode {
            CompileMode::All | CompileMode::Both | CompileMode::Forward => {
                let mut settings = self.settings.borrow(py).0.clone();
                settings.mode = Forward;
                match settings.compton_mode {
                    Adjoint | Inverse => settings.compton_mode = Direct,
                    _ =>(),
                }
                let args = (PyTransportSettings(settings),);
                self.registry.call_method(py, "compute", args, kwargs)?;
            },
            _ => (),
        }
        match mode {
            CompileMode::All | CompileMode::Both | CompileMode::Backward => {
                let mut settings = self.settings.borrow(py).0.clone();
                settings.mode = Backward;
                match settings.compton_mode {
                    Direct => settings.compton_mode = Adjoint,
                    _ =>(),
                }
                if let Inverse = settings.compton_mode {
                    settings.compton_method = ComptonMethod::InverseTransform;
                }
                let args = (PyTransportSettings(settings),);
                self.registry.call_method(py, "compute", args, kwargs)?;
            },
            _ => (),
        }
        match mode {
            CompileMode::All => {
                let mut settings = self.settings.borrow(py).0.clone();
                settings.mode = Backward;
                settings.compton_mode = Inverse;
                settings.compton_method = ComptonMethod::InverseTransform;
                let args = (PyTransportSettings(settings),);
                self.registry.call_method(py, "compute", args, kwargs)?;
            },
            _ => (),
        }
        Ok(())
    }

    fn transport(&self, py: Python, states: &PyArray<CState>) -> Result<PyObject> {
        match &self.geometry {
            None => type_error!(
                "bad geometry (expected an instance of 'ExternalGeometry' or 'SimpleGeometry' \
                 found 'none')"
            ),
            Some(geometry) => match geometry {
                PyGeometryDefinition::External(external) => {
                    self.transport_with::<_, ExternalTracer>(
                        py, &external.borrow(py).0, states
                    )
                },
                PyGeometryDefinition::Simple(simple) => {
                    self.transport_with::<_, SimpleTracer>(
                        py, &simple.borrow(py).0, states
                    )
                },
            },
        }
    }
}

impl PyTransportEngine {
    fn update_with<G>(&self, geometry: &G, registry: &mut MaterialRegistry) -> Result<()>
    where
        G: GeometryDefinition,
    {
        for material in geometry.materials().iter() {
            registry.add(material)?;
        }
        Ok(())
    }

    fn transport_with<'a, G, T>(
        &self,
        py: Python,
        geometry: &'a G,
        states: &PyArray<CState>,
    ) -> Result<PyObject>
    where
        G: GeometryDefinition,
        T: GeometryTracer<'a, G>,
    {
        // Create the status array.
        let status = PyArray::<i32>::empty(py, &states.shape())?;

        // Unpack registry and settings.
        let registry = &self.registry.borrow(py).inner;
        let settings = &self.settings.borrow(py).0;

        // Get a transport agent.
        let rng: &mut PyRandomStream = &mut self.random.borrow_mut(py);
        let mut agent = TransportAgent::<G, _, T>::new(geometry, registry, rng)?;

        // Do the Monte Carlo transport.
        let n = states.size();
        for i in 0..n {
            let mut state: PhotonState = states.get(i)?.into();
            let flag = agent.transport(settings, &mut state)?;
            states.set(i, state.into())?;
            status.set(i, flag.into())?;

            if i % 100 == 0 { // Check for a Ctrl+C interrupt, catched by Python.
                ctrlc_catched()?;
            }
        }

        let status: &PyAny = status;
        Ok(status.into())
    }
}


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
// Utility class for creating a numpy array of photon states.
// ===============================================================================================

#[pyclass(name = "PhotonState", module="goupil")]
pub struct PyPhotonState;

#[pymethods]
impl PyPhotonState {
    #[classattr]
    fn dtype(py: Python) -> Result<PyObject> {
        let dtype = CState::dtype(py)?;
        Ok(dtype)
    }

    #[staticmethod]
    fn empty(py: Python, shape: ShapeArg) -> Result<PyObject> {
        let shape: Vec<usize> = shape.into();
        let array: &PyAny = PyArray::<CState>::empty(py, &shape)?;
        Ok(array.into())
    }

    #[staticmethod]
    fn zeros(py: Python, shape: ShapeArg) -> Result<PyObject> {
        let shape: Vec<usize> = shape.into();
        let array: &PyAny = PyArray::<CState>::zeros(py, &shape)?;
        Ok(array.into())
    }
}

#[derive(FromPyObject)]
pub enum ShapeArg {
    Scalar(usize),
    Vector(Vec<usize>),
}

impl From<ShapeArg> for Vec<usize> {
    fn from(value: ShapeArg) -> Self {
        match value {
            ShapeArg::Scalar(value) => vec![value],
            ShapeArg::Vector(value) => value,
        }
    }
}


// ===============================================================================================
// Python class forwarding transport status codes.
// ===============================================================================================

#[pyclass(name = "TransportStatus", module="goupil")]
pub(crate) struct PyTransportStatus ();

#[allow(non_snake_case)]
#[pymethods]
impl PyTransportStatus {
    #[classattr]
    fn ABSORBED(py: Python<'_>) -> Result<PyObject> {
        Self::into_i32(py, TransportStatus::Absorbed)
    }

    #[classattr]
    fn BOUNDARY(py: Python<'_>) -> Result<PyObject> {
        Self::into_i32(py, TransportStatus::Boundary)
    }

    #[classattr]
    fn CONSTRAINT(py: Python<'_>) -> Result<PyObject> {
        Self::into_i32(py, TransportStatus::Constraint)
    }

    #[classattr]
    fn ENERGY_MAX(py: Python<'_>) -> Result<PyObject> {
        Self::into_i32(py, TransportStatus::EnergyMax)
    }

    #[classattr]
    fn ENERGY_MIN(py: Python<'_>) -> Result<PyObject> {
        Self::into_i32(py, TransportStatus::EnergyMin)
    }

    #[classattr]
    fn EXIT(py: Python<'_>) -> Result<PyObject> {
        Self::into_i32(py, TransportStatus::Exit)
    }

    #[classattr]
    fn LENGTH_MAX(py: Python<'_>) -> Result<PyObject> {
        Self::into_i32(py, TransportStatus::LengthMax)
    }

    /// Return the string representation of a `TransportStatus` integer code.
    #[staticmethod]
    fn str(code: i32) -> Result<String> {
        let status: TransportStatus = code.try_into()?;
        Ok(status.into())
    }
}

impl PyTransportStatus {
    fn into_i32(py: Python, status: TransportStatus) -> Result<PyObject> {
        let value: i32 = status.into();
        let scalar = PyScalar::new(py, value)?;
        Ok(scalar.into())
    }
}
