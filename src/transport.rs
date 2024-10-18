use anyhow::{anyhow, Result};
use crate::numerics::float::{Float, Float3};
use crate::physics::{AbsorptionMode, ComptonMethod, ComptonMode, ComptonModel, RayleighMode};
use crate::pretty_enumerate;
use enum_iterator::{all, Sequence};
use serde_derive::{Deserialize, Serialize};
use std::fmt;

pub(crate) mod agent;
pub(crate) mod boundary;
pub(crate) mod density;
pub(crate) mod geometry;


// ===============================================================================================
// Public API.
// ===============================================================================================

pub use self::agent::{
    TransportAgent,
    TransportStatus
};
pub use self::boundary::{
    BoxShape,
    GeometryShape,
    SphereShape,
    TransportBoundary,
};
pub use self::density::DensityModel;
pub use self::geometry::{
    ExternalGeometry,
    ExternalTracer,
    GeometryDefinition,
    GeometrySector,
    GeometryTracer,
    SimpleGeometry,
    SimpleTracer,
    StratifiedGeometry,
    StratifiedTracer,
    TopographyMap,
    TopographySurface,
};


// ===============================================================================================
// Transport settings.
// ===============================================================================================

#[derive(Clone, Copy, Deserialize, Serialize)]
pub struct TransportSettings {
    // Global mode flag (forward or backward flow).
    pub mode: TransportMode,

    // Physics and sampling settings.
    pub absorption: AbsorptionMode,
    pub compton_method: ComptonMethod,
    pub compton_mode: ComptonMode,
    pub compton_model: ComptonModel,
    pub rayleigh: RayleighMode,

    // Transport boundary.
    pub boundary: TransportBoundary,

    // External limits.
    pub constraint: Option<Float>,
    pub energy_min: Option<Float>,
    pub energy_max: Option<Float>,
    pub length_max: Option<Float>,
}

impl Default for TransportSettings {
    fn default() -> Self {
        Self {
            mode: TransportMode::Forward,
            absorption: AbsorptionMode::Discrete,
            compton_method: ComptonMethod::RejectionSampling,
            compton_mode: ComptonMode::Direct,
            compton_model: ComptonModel::ScatteringFunction,
            rayleigh: RayleighMode::FormFactor,
            boundary: TransportBoundary::None,
            constraint: None,
            energy_min: None,
            energy_max: None,
            length_max: None,
        }
    }
}

// Private interface.
impl TransportSettings {
    #[inline]
    pub(crate) fn is_forward(&self) -> bool {
        match self.mode {
            TransportMode::Backward => false,
            TransportMode::Forward => true,
        }
    }
}


// ===============================================================================================
// Transport mode.
// ===============================================================================================

#[derive(Clone, Copy, Default, Deserialize, PartialEq, Sequence, Serialize)]
pub enum TransportMode {
    Backward,
    #[default]
    Forward,
}

impl TransportMode {
    const BACKWARD: &str = "Backward";
    const FORWARD: &str = "Forward";

    fn pretty_variants() -> String {
        let variants: Vec<_> = all::<Self>()
            .map(|e| format!("'{}'", e))
            .collect();
        pretty_enumerate(&variants)
    }
}

impl fmt::Display for TransportMode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s: &str = (*self).into();
        write!(f, "{}", s)
    }
}

impl TryFrom<&str> for TransportMode {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            Self::BACKWARD => Ok(Self::Backward),
            Self::FORWARD => Ok(Self::Forward),
            _ => Err(anyhow!(
                "bad transport mode (expected {}, found '{}')",
                Self::pretty_variants(),
                value,
            )),
        }
    }
}

impl From<TransportMode> for &str {
    fn from(value: TransportMode) -> Self {
        match value {
            TransportMode::Backward => TransportMode::BACKWARD,
            TransportMode::Forward => TransportMode::FORWARD,
        }
    }
}


// ===============================================================================================
// Monte Carlo state of an unpolarized photon.
// ===============================================================================================

#[derive(Clone, Default)]
pub struct PhotonState {
    pub energy: Float,    // MeV
    pub position: Float3, // cm
    pub direction: Float3,
    pub length: Float,    // cm
    pub weight: Float,
}

impl PhotonState {
    pub fn new(energy: Float, position: Float3, direction: Float3, weight: Float) -> Self {
        Self {energy, position, direction, weight, length: 0.0}
    }
}


// ===============================================================================================
// Monte Carlo transport vertex.
// ===============================================================================================

pub struct TransportVertex {
    pub sector: usize,
    pub kind: VertexKind,
    pub state: PhotonState,
}

pub enum VertexKind {
    Compton,
    Interface,
    Rayleigh,
    Start,
    Stop,
}
