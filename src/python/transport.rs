use anyhow::Result;
use crate::numerics::float::Float;
use crate::physics::materials::MaterialRegistry;
use crate::physics::process::{
    absorption::AbsorptionMode,
    compton::{ComptonModel, ComptonMethod, ComptonMode::{self, Adjoint, Direct, Inverse}},
    rayleigh::RayleighMode,
};
use crate::transport::{
    agent::{TransportAgent, TransportStatus},
    boundary::TransportBoundary,
    geometry::{ExternalTracer, GeometryDefinition, GeometryTracer, SimpleTracer, StratifiedTracer},
    TransportMode::{self, Backward, Forward},
    TransportSettings, TransportVertex, VertexKind,
};
use derive_more::{AsMut, AsRef, From};
use pyo3::{
    prelude::*,
    gc::PyVisit,
    PyTraverseError,
    types::{PyBytes, PyDict, PyString},
};
use rmp_serde::{Deserializer, Serializer};
use serde::{Deserialize, Serialize};
use super::{
    boundary::PyTransportBoundary,
    ctrlc_catched,
    export::Export,
    geometry::{PyExternalGeometry, PyGeometryDefinition},
    macros::{type_error, value_error},
    materials::PyMaterialRegistry,
    namespace::Namespace,
    numpy::{ArrayOrFloat, PyArray, PyArrayMethods, PyScalar},
    rand::PyRandomStream,
    states::States,
    prefix,
};


// ===============================================================================================
// Python wrapper for a Goupil Monte Carlo engine.
// ===============================================================================================

#[pyclass(name = "TransportSettings", module = "goupil")]
pub(crate) struct PyTransportSettings {
    pub inner: TransportSettings,
    pub volume_sources: bool,
}

// Convert from raw type.
impl Into<PyTransportSettings> for TransportSettings {
    fn into(self) -> PyTransportSettings {
        let volume_sources = match self.constraint {
            None => false,
            Some(_) => true,
        };
        PyTransportSettings {
            inner: self,
            volume_sources
        }
    }
}

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
        let mut inner = TransportSettings::default();
        inner.constraint = Some(1.0);
        Self {
            inner,
            volume_sources: true,
        }
    }

    #[getter]
    fn get_mode(&self) -> &str {
        self.inner.mode.into()
    }

    #[setter]
    fn set_mode(&mut self, value: &str) -> Result<()> {
        self.inner.mode = TransportMode::try_from(value)?;
        match self.inner.mode {
            Backward => match self.inner.compton_mode {
                Direct => {
                    self.inner.compton_mode = Adjoint;
                },
                _ => (),
            },
            Forward => match self.inner.compton_mode {
                Adjoint | Inverse => {
                    self.inner.compton_mode = Direct;
                },
                _ => (),
            },
        }
        Ok(())
    }

    #[getter]
    fn get_absorption(&self) -> Option<&str> {
        to_optstr!(AbsorptionMode, self.inner.absorption)
    }

    #[setter]
    fn set_absorption(&mut self, value: Option<&str>) -> Result<()> {
        from_optstr!(AbsorptionMode, self.inner.absorption, value);
        Ok(())
    }

    #[getter]
    fn get_boundary(&self) -> TransportBoundary {
        self.inner.boundary
    }

    #[setter]
    fn set_boundary(&mut self, value: Option<PyTransportBoundary>) -> Result<()> {
        match value {
            None => self.inner.boundary = TransportBoundary::None,
            Some(boundary) => self.inner.boundary = boundary.into(),
        };
        Ok(())
    }

    #[getter]
    fn get_compton_method(&self) -> &str {
        self.inner.compton_method.into()
    }

    #[setter]
    fn set_compton_method(&mut self, value: &str) -> Result<()> {
        self.inner.compton_method = ComptonMethod::try_from(value)?;
        Ok(())
    }

    #[getter]
    fn get_compton_mode(&self) -> Option<&str> {
        to_optstr!(ComptonMode, self.inner.compton_mode)
    }

    #[setter]
    fn set_compton_mode(&mut self, value: Option<&str>) -> Result<()> {
        from_optstr!(ComptonMode, self.inner.compton_mode, value);
        match self.inner.compton_mode {
            Adjoint => {
                self.inner.mode = Backward;
            },
            Direct => {
                self.inner.mode = Forward;
            },
            Inverse => {
                self.inner.mode = Backward;
                self.inner.compton_method = ComptonMethod::InverseTransform;
            },
            ComptonMode::None => (),
        }
        Ok(())
    }

    #[getter]
    fn get_compton_model(&self) -> &str {
        self.inner.compton_model.into()
    }

    #[setter]
    fn set_compton_model(&mut self, value: &str) -> Result<()> {
        self.inner.compton_model = ComptonModel::try_from(value)?;
        Ok(())
    }

    #[getter]
    fn get_volume_sources(&self) -> bool {
        self.volume_sources
    }

    #[setter]
    fn set_volume_sources(&mut self, value: Option<bool>) -> Result<()> {
        let value = value.unwrap_or(false);
        self.volume_sources = value;
        if value {
            self.inner.constraint = Some(1.0);
        } else {
            self.inner.constraint = None;
        }
        Ok(())
    }

    #[getter]
    fn get_rayleigh(&self) -> bool {
        match self.inner.rayleigh {
            RayleighMode::FormFactor => true,
            RayleighMode::None => false,
        }
    }

    #[setter]
    fn set_rayleigh(&mut self, value: Option<bool>) -> Result<()> {
        let value = value.unwrap_or(false);
        if value {
            self.inner.rayleigh = RayleighMode::FormFactor;
        } else {
            self.inner.rayleigh = RayleighMode::None;
        }
        Ok(())
    }

    #[getter]
    fn get_energy_min(&self) -> Option<Float> {
        self.inner.energy_min
    }

    #[setter]
    fn set_energy_min(&mut self, value: Option<Float>) -> Result<()> {
        self.inner.energy_min = value;
        Ok(())
    }

    #[getter]
    fn get_energy_max(&self) -> Option<Float> {
        self.inner.energy_max
    }

    #[setter]
    fn set_energy_max(&mut self, value: Option<Float>) -> Result<()> {
        self.inner.energy_max = value;
        Ok(())
    }

    #[getter]
    fn get_length_max(&self) -> Option<Float> {
        self.inner.length_max
    }

    #[setter]
    fn set_length_max(&mut self, value: Option<Float>) -> Result<()> {
        self.inner.length_max = value;
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
    pub(crate) random: Py<PyRandomStream>,
    #[pyo3(get)]
    registry: Py<PyMaterialRegistry>,
    #[pyo3(get)]
    pub(crate) settings: Py<PyTransportSettings>,

    compiled: bool,
}

#[derive(FromPyObject)]
enum GeometryArg {
    Object(PyGeometryDefinition),
    Path(String),
}

#[pymethods]
impl PyTransportEngine {
    #[new]
    #[pyo3(signature = (geometry=None, *, random=None, registry=None, settings=None))]
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
                        let external = PyExternalGeometry::new(py, &path)?;
                        let external = Py::new(py, external)?;
                        PyGeometryDefinition::External(external)
                    },
                };
                Some(geometry)
            },
        };
        let random: Py<PyRandomStream> = match random {
            None => Py::new(py, PyRandomStream::new(None, None)?)?,
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
        Ok(Self { geometry, random, registry, settings, compiled: false })
    }

    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.random)?;
        visit.call(&self.registry)?;
        visit.call(&self.settings)?;
        Ok(())
    }

    fn __getattr__(&self, py: Python, name: &Bound<PyString>) -> Result<PyObject> {
        Ok(self.settings.getattr(py, name)?)
    }

    fn __setattr__(&mut self, py: Python, name: &str, value: PyObject) -> Result<()> {
        match name {
            "boundary" => {
                let boundary: TransportBoundary = if value.is_none(py) {
                    TransportBoundary::None
                } else {
                    let value: BoundaryArg = value.extract(py)?;
                    match value {
                        BoundaryArg::Description(description) => {
                            let index = match &self.geometry {
                                None => value_error!(
                                    "could not find any sector for '{}'",
                                    description
                                ),
                                Some(geometry) => geometry.sector_index(py, description.as_str())?,
                            };
                            TransportBoundary::Sector(index)
                        },
                        BoundaryArg::Explicit(boundary) => boundary.into(),
                    }
                };
                let settings = &mut self.settings.borrow_mut(py).inner;
                settings.boundary = boundary;
            },
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
    pub fn __setstate__(&mut self, state: &Bound<PyBytes>) -> Result<()> {
        let py = state.py();
        let mut deserializer = Deserializer::new(state.as_bytes());

        let mut random = self.random.borrow_mut(py);
        *random = Deserialize::deserialize(&mut deserializer)?;

        let registry = &mut self.registry.borrow_mut(py).inner;
        *registry = Deserialize::deserialize(&mut deserializer)?;

        let settings = &mut self.settings.borrow_mut(py);
        settings.inner = Deserialize::deserialize(&mut deserializer)?;
        match settings.inner.constraint {
            None => settings.volume_sources = false,
            Some(_) => settings.volume_sources = true,
        }

        self.compiled = Deserialize::deserialize(&mut deserializer)?;

        Ok(())
    }

    fn __getstate__<'py>(&self, py: Python<'py>) -> Result<Bound<'py, PyBytes>> {
        let mut buffer = Vec::new();
        let mut serializer = Serializer::new(&mut buffer);

        let random = &self.random.borrow(py);
        random.serialize(&mut serializer)?;

        let registry = &self.registry.borrow(py).inner;
        registry.serialize(&mut serializer)?;

        let settings = &self.settings.borrow(py).inner;
        settings.serialize(&mut serializer)?;

        self.compiled.serialize(&mut serializer)?;

        Ok(PyBytes::new_bound(py, &buffer))
    }

    #[pyo3(signature = (mode=None, *, atomic_data=None, **kwargs))]
    fn compile(
        &mut self,
        py: Python,
        mode: Option<&str>,
        atomic_data: Option<&str>,
        kwargs: Option<&Bound<PyDict>>,
    ) -> Result<()> {
        enum CompileMode {
            All,
            Backward,
            Both,
            Forward,
        }

        let mode = match mode {
            None => match &self.settings.borrow(py).inner.mode {
                TransportMode::Backward => CompileMode::Backward,
                TransportMode::Forward => CompileMode::Forward,
            },
            Some(mode) => match mode {
                "All" => CompileMode::All,
                "Backward" => CompileMode::Backward,
                "Both" => CompileMode::Both,
                "Forward" => CompileMode::Forward,
                _ => value_error!(
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
                        self.update_with(&external.borrow(py).inner, registry)?
                    },
                    PyGeometryDefinition::Simple(simple) => {
                        self.update_with(&simple.borrow(py).0, registry)?
                    },
                    PyGeometryDefinition::Stratified(stratified) => {
                        self.update_with(&stratified.borrow(py).inner, registry)?
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
                let mut settings = self.settings.borrow(py).inner.clone();
                settings.mode = Forward;
                match settings.compton_mode {
                    Adjoint | Inverse => settings.compton_mode = Direct,
                    _ =>(),
                }
                let args = (Into::<PyTransportSettings>::into(settings),);
                self.registry.call_method_bound(py, "compute", args, kwargs)?;
            },
            _ => (),
        }
        match mode {
            CompileMode::All | CompileMode::Both | CompileMode::Backward => {
                let mut settings = self.settings.borrow(py).inner.clone();
                settings.mode = Backward;
                match settings.compton_mode {
                    Direct => settings.compton_mode = Adjoint,
                    _ =>(),
                }
                if let Inverse = settings.compton_mode {
                    settings.compton_method = ComptonMethod::InverseTransform;
                }
                let args = (Into::<PyTransportSettings>::into(settings),);
                self.registry.call_method_bound(py, "compute", args, kwargs)?;
            },
            _ => (),
        }
        match mode {
            CompileMode::All => {
                let mut settings = self.settings.borrow(py).inner.clone();
                settings.mode = Backward;
                settings.compton_mode = Inverse;
                settings.compton_method = ComptonMethod::InverseTransform;
                let args = (Into::<PyTransportSettings>::into(settings),);
                self.registry.call_method_bound(py, "compute", args, kwargs)?;
            },
            _ => (),
        }

        // Record compilation step.
        self.compiled = true;

        Ok(())
    }

    #[pyo3(signature = (states, /, *, source_energies=None, random_index=None, vertices=None))]
    fn transport(
        &mut self,
        states: &Bound<PyAny>,
        source_energies: Option<ArrayOrFloat>,
        random_index: Option<RandomIndexArg>,
        vertices: Option<bool>,
    ) -> Result<PyObject> {
        // Extract Monte Carlo states.
        let py = states.py();
        let states = States::new(states)?;

        if self.settings.borrow(py).inner.mode == Backward {
            if !states.has_weights() {
                type_error!("bad states (expected 'weight' field, found None)")
            }
        }

        // Check constraints and states consistency.
        if let Some(constraints) = source_energies.as_ref() {
            if let ArrayOrFloat::Array(constraints) = constraints {
                if constraints.size() != states.size() {
                    value_error!(
                        "bad constraints (expected a scalar or a size {} array, \
                         found a size {} array)",
                        states.size(),
                        constraints.size(),
                    )
                }
            }
        }

        // Compile, if not already done.
        if !self.compiled {
            self.compile(py, Some("Both"), None, None)?;
        }

        // Run the Monte Carlo simulation.
        match &self.geometry {
            None => type_error!(
                "bad geometry (expected an instance of 'ExternalGeometry' or 'SimpleGeometry' \
                 found 'none')"
            ),
            Some(geometry) => match geometry {
                PyGeometryDefinition::External(external) => {
                    self.transport_with::<_, ExternalTracer>(
                        py, &external.borrow(py).inner, states, source_energies, random_index,
                        vertices,
                    )
                },
                PyGeometryDefinition::Simple(simple) => {
                    self.transport_with::<_, SimpleTracer>(
                        py, &simple.borrow(py).0, states, source_energies, random_index,
                        vertices,
                    )
                },
                PyGeometryDefinition::Stratified(stratified) => {
                    self.transport_with::<_, StratifiedTracer>(
                        py, &stratified.borrow(py).inner, states, source_energies, random_index,
                        vertices,
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
        states: States,
        constraints: Option<ArrayOrFloat>,
        random_index: Option<RandomIndexArg>,
        vertices: Option<bool>,
    ) -> Result<PyObject>
    where
        G: GeometryDefinition,
        T: GeometryTracer<'a, G>,
    {
        // Create the status array.
        let status = PyArray::<i32>::empty(py, &states.shape())?;

        // Unpack registry and settings.
        let registry = &self.registry.borrow(py).inner;
        let mut settings = self.settings.borrow(py).inner.clone();
        if constraints.is_none() {
            settings.constraint = None;
        }

        // Check consistency of settings with explicit constraints.
        if !constraints.is_none() {
            if settings.mode == TransportMode::Forward {
                value_error!("bad constraints (unused in 'Forward' mode)")
            } else {
                if settings.constraint.is_none() {
                    value_error!("bad constraints (disabled by transport settings)")
                }
            }
        }

        // Use registry energy limits if no explicit bound was specified.
        if settings.energy_min.is_none() {
            settings.energy_min = Some(registry.energy_min);
            settings.energy_max = Some(registry.energy_max);
        }

        // Get a transport agent.
        let rng: &mut PyRandomStream = &mut self.random.borrow_mut(py);
        let mut agent = TransportAgent::<G, _, T>::new(geometry, registry, rng)?;

        // Set random indices container.
        let random_index = random_index.unwrap_or(RandomIndexArg::Bool(false));
        let (indices_in, indices_out) = match random_index {
            RandomIndexArg::Bool(random_index) => match random_index {
                false => (None, None),
                true => {
                    #[cfg(not(feature = "f32"))]
                    let shape = {
                        let mut shape = states.shape();
                        shape.push(2);
                        shape
                    };

                    #[cfg(feature = "f32")]
                    let shape = states.shape();

                    (None, Some(PyArray::<u64>::empty(py, &shape)?))
                },
            },
            RandomIndexArg::Array(indices_in) => (Some(indices_in), None),
        };

        // Set vertices container.
        let mut vertices: Option<Vec<CVertex>> = match vertices.unwrap_or(false) {
            false => None,
            true => Some(Vec::new()),
        };

        // Do the Monte Carlo transport.
        let n = states.size();
        for i in 0..n {
            let mut state = states.get(i)?;
            if let Some(constraints) = constraints.as_ref() {
                let constraint = match constraints {
                    ArrayOrFloat::Array(constraints) => constraints.get(i)?,
                    ArrayOrFloat::Float(constraint) => *constraint,
                };
                settings.constraint = Some(constraint);
            }

            #[cfg(not(feature = "f32"))]
            if let Some(ref indices_in) = indices_in {
                let index = super::rand::arg::Index::Array([
                    indices_in.get(2 * i)?,
                    indices_in.get(2 * i + 1)?,
                ]);
                agent.rng_mut().set_index(Some(index))?;
            }

            #[cfg(feature = "f32")]
            if let Some(ref indices_in) = indices_in {
                let index = indices_in.get(i)?;
                agent.rng_mut().set_index(Some(index))?;
            }

            #[cfg(not(feature = "f32"))]
            if let Some(ref indices_out) = indices_out {
                let index = agent.rng().index_2u64();
                indices_out.set(2 * i, index[0])?;
                indices_out.set(2 * i + 1, index[1])?;
            }

            #[cfg(feature = "f32")]
            if let Some(ref indices_out) = indices_out {
                let index = agent.rng().index();
                indices_out.set(i, index)?;
            }

            let mut verts = vertices.as_ref().map(|_| Vec::new());
            let flag = agent.transport(&settings, &mut state, verts.as_mut())?;
            states.set(i, &state)?;
            status.set(i, flag.into())?;

            if let Some(mut verts) = verts {
                let vertices = vertices.as_mut().unwrap();
                vertices.reserve(verts.len());
                for vertex in verts.drain(..) {
                    vertices.push((i, flag, vertex).into());
                }
            }

            if i % 100 == 0 { // Check for a Ctrl+C interrupt, catched by Python.
                ctrlc_catched()?;
            }
        }

        let status = status.into_any().unbind();

        let result: PyObject = match indices_out {
            None => match vertices {
                None => status,
                Some(vertices) => {
                    let vertices = Export::export::<PyTransportVertices>(py, vertices)?;
                    Namespace::new(py, &[
                        ("status", status),
                        ("vertices", vertices),
                    ])?.unbind()
                }
            },
            Some(indices_out) => {
                let indices_out: PyObject = indices_out.into_py(py);
                match vertices {
                    None => Namespace::new(py, &[
                        ("status", status),
                        ("random_index", indices_out),
                    ])?.unbind(),
                    Some(vertices) => {
                        let vertices = Export::export::<PyTransportVertices>(py, vertices)?;
                        Namespace::new(py, &[
                            ("status", status),
                            ("random_index", indices_out),
                            ("vertices", vertices),
                        ])?.unbind()
                    }
                }
            },
        };

        Ok(result)
    }
}

#[derive(FromPyObject)]
enum BoundaryArg<'py> {
    Description(String),
    Explicit(PyTransportBoundary<'py>),
}

#[derive(FromPyObject)]
enum RandomIndexArg<'a> {
    Bool(bool),
    Array(&'a PyArray<u64>),
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
    fn ENERGY_CONSTRAINT(py: Python<'_>) -> Result<PyObject> {
        Self::into_i32(py, TransportStatus::EnergyConstraint)
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
        let status: &'static str = status.into();
        Ok(status.to_string())
    }
}

impl PyTransportStatus {
    fn into_i32(py: Python, status: TransportStatus) -> Result<PyObject> {
        let value: i32 = status.into();
        let scalar = PyScalar::new(py, value)?;
        Ok(scalar.into_py(py))
    }
}


// ===============================================================================================
// C-compliant transport vertices.
// ===============================================================================================
//
#[derive(AsMut, AsRef, From)]
#[pyclass(name = "TransportVertices", module="goupil")]
struct PyTransportVertices (Export<CVertex>);

#[repr(C)]
#[derive(Clone, Copy)]
pub(crate) struct CVertex {
    pub event: usize,
    pub sector: usize,
    pub energy: Float,
    pub position: [Float; 3],
    pub kind: [u8; 16],
}

impl From<(usize, TransportStatus, TransportVertex)> for CVertex {
    fn from(value: (usize, TransportStatus, TransportVertex)) -> Self {
        let (event, status, vertex) = value;
        let kind: [u8; 16] = {
            let kind = match vertex.kind {
                VertexKind::Compton => "Compton",
                VertexKind::Interface => "Interface",
                VertexKind::Rayleigh => "Rayleigh",
                VertexKind::Start => "Start",
                VertexKind::Stop => status.into(),
            };
            let kind = kind.as_bytes();
            let mut buffer = [0_u8; 16];
            let n = kind.len();
            for i in 0..n {
                if i == 15 {
                    break
                } else if (i == 14) && (n > 15) {
                    buffer[i] = u32::from('.') as u8;
                } else {
                    buffer[i] = kind[i];
                }
            }
            buffer
        };
        Self {
            event,
            sector: vertex.sector,
            energy: vertex.state.energy,
            position: vertex.state.position.into(),
            kind,
        }
    }
}
