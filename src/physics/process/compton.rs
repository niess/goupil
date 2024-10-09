use anyhow::{anyhow, bail, Result};
use crate::pretty_enumerate;
use enum_iterator::{all, Sequence};
use serde_derive::{Deserialize, Serialize};
use std::fmt;

pub(crate) mod compute;
pub(crate) mod sample;
pub(crate) mod table;


// ===============================================================================================
// Public API.
// ===============================================================================================

pub use self::sample::ComptonSampler;
pub use self::table::{
    ComptonCrossSection,
    ComptonCDF,
    ComptonInverseCDF,
};


// ===============================================================================================
// Compton models.
// ===============================================================================================

#[derive(Clone, Copy, Default, Deserialize, PartialEq, Sequence, Serialize)]
pub enum ComptonModel {
    /// Klein-Nishina model.
    ///
    /// It is assumed that the electrons being targeted are both free and at rest. This model
    /// yields the Klein-Nishina cross-section, which is an approximation that neglects atomic
    /// binding effects and Doppler broadening.
    KleinNishina,

    /// Penelope model for incoherent (Compton) scattering.
    ///
    /// The target is modelled as an incoherent superposition of electronic shells, following the
    /// Impulse Approximation. The resulting electronic structure can be reduced to a single
    /// longitudinal parameter, J_i, also known as the Compton profile, from which the DDCS is
    /// obtained in terms of outgoing energy and scattering angle. This model takes into account
    /// both atomic binding effects and Doppler broadening, striking a good balance between speed
    /// and accuracy when compared to a full computation. It is, however, only implemented in the
    /// forward mode.
    Penelope,

    /// Effective model, using Klein-Nishina DCS with a Scattering Function corrective factor.
    ///
    /// This model is an intermediate between [`Penelope`](Self::Penelope) and
    /// [`Klein-Nishina`](Self::KleinNishina). The differential cross section with respect to
    /// energy is calculated using Penelope's scattering function, which depends on the electronic
    /// structure, to recover Penelope's total cross section. However, the scattering angle is
    /// deterministic, assuming a collision with a free electron. As a result, the Doppler
    /// broadening caused by the motion of the target electrons is not included in the DCS, but
    /// only in the total cross section. Atomic binding effects are included in both the DCS and
    /// the total cross section.
    #[default]
    ScatteringFunction,
}

impl ComptonModel {
    const KLEIN_NISHINA: &str = "Klein-Nishina";
    const PENELOPE: &str = "Penelope";
    const SCATTERING_FUNCTION: &str = "Scattering Function";

    fn pretty_variants() -> String {
        let variants: Vec<_> = all::<Self>()
            .map(|e| format!("'{}'", e))
            .collect();
        pretty_enumerate(&variants)
    }
}

impl fmt::Display for ComptonModel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s: &str = (*self).into();
        write!(f, "{}", s)
    }
}

impl TryFrom<&str> for ComptonModel {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            Self::KLEIN_NISHINA => Ok(Self::KleinNishina),
            Self::PENELOPE => Ok(Self::Penelope),
            Self::SCATTERING_FUNCTION => Ok(Self::ScatteringFunction),
            _ => Err(anyhow!(
                "bad Compton model (expected {}, found '{}')",
                Self::pretty_variants(),
                value,
            )),
        }
    }
}

impl From<ComptonModel> for &str {
    fn from(value: ComptonModel) -> Self {
        match value {
            ComptonModel::KleinNishina => ComptonModel::KLEIN_NISHINA,
            ComptonModel::Penelope => ComptonModel::PENELOPE,
            ComptonModel::ScatteringFunction => ComptonModel::SCATTERING_FUNCTION,
        }
    }
}


// ===============================================================================================
// Sampling modes.
// ===============================================================================================

#[derive(Clone, Copy, Default, Deserialize, PartialEq, Sequence, Serialize)]
pub enum ComptonMode {
    Adjoint,
    #[default]
    Direct,
    Inverse,
    None,
}

impl ComptonMode {
    const ADJOINT: &str = "Adjoint";
    const DIRECT: &str = "Direct";
    const INVERSE: &str = "Inverse";
    const NONE: &str = "None";

    fn pretty_variants() -> String {
        let variants: Vec<_> = all::<Self>()
            .map(|e| format!("'{}'", e))
            .collect();
        pretty_enumerate(&variants)
    }
}

impl fmt::Display for ComptonMode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s: &str = (*self).into();
        write!(f, "{}", s)
    }
}

impl TryFrom<&str> for ComptonMode {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            Self::ADJOINT => Ok(Self::Adjoint),
            Self::DIRECT => Ok(Self::Direct),
            Self::INVERSE => Ok(Self::Inverse),
            Self::NONE => Ok(Self::None),
            _ => Err(anyhow!(
                "bad sampling mode (expected {}, found '{}')",
                Self::pretty_variants(),
                value,
            )),
        }
    }
}

impl From<ComptonMode> for &str {
    fn from(value: ComptonMode) -> Self {
        match value {
            ComptonMode::Adjoint => ComptonMode::ADJOINT,
            ComptonMode::Direct => ComptonMode::DIRECT,
            ComptonMode::Inverse => ComptonMode::INVERSE,
            ComptonMode::None => ComptonMode::NONE,
        }
    }
}


// ===============================================================================================
// Sampling methods.
// ===============================================================================================

#[derive(Default, Clone, Copy, Deserialize, PartialEq, Sequence, Serialize)]
pub enum ComptonMethod {
    InverseTransform,
    #[default]
    RejectionSampling,
}

impl ComptonMethod {
    const INVERSE_TRANSFORM: &str = "Inverse Transform";
    const REJECTION_SAMPLING: &str = "Rejection Sampling";

    fn pretty_variants() -> String {
        let variants: Vec<_> = all::<Self>()
            .map(|e| format!("'{}'", e))
            .collect();
        pretty_enumerate(&variants)
    }
}

impl fmt::Display for ComptonMethod {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s: &str = (*self).into();
        write!(f, "{}", s)
    }
}

impl TryFrom<&str> for ComptonMethod {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            Self::INVERSE_TRANSFORM => Ok(Self::InverseTransform),
            Self::REJECTION_SAMPLING => Ok(Self::RejectionSampling),
            _ => Err(anyhow!(
                "bad sampling method (expected {}, found '{}')",
                Self::pretty_variants(),
                value,
            )),
        }
    }
}

impl From<ComptonMethod> for &str {
    fn from(value: ComptonMethod) -> Self {
        match value {
            ComptonMethod::InverseTransform => ComptonMethod::INVERSE_TRANSFORM,
            ComptonMethod::RejectionSampling => ComptonMethod::REJECTION_SAMPLING,
        }
    }
}


// ===============================================================================================
// Validity matrix.
// ===============================================================================================

pub(crate) fn validate(
    model: ComptonModel,
    mode: ComptonMode,
    method: ComptonMethod,
) -> Result<()> {
    if let ComptonMode::None = mode {
        return Ok(())
    }
    match model {
        ComptonModel::Penelope => {
            match mode {
                ComptonMode::Direct => (),
                _ => bail!(
                    "bad sampling mode for '{}' Compton model (expected '{}', found '{}')",
                    model,
                    ComptonMode::Direct,
                    mode,
                )
            };
            match method {
                ComptonMethod::InverseTransform => bail!(
                    "bad sampling method for '{}' Compton model (expected '{}', found '{}')",
                    model,
                    ComptonMethod::RejectionSampling,
                    mode,
                ),
                ComptonMethod::RejectionSampling => (),
            };
        },
        ComptonModel::ScatteringFunction | ComptonModel::KleinNishina => {
            match method {
                ComptonMethod::InverseTransform => (),
                ComptonMethod::RejectionSampling => match mode {
                    ComptonMode::Inverse => bail!(
                        "bad sampling mode for '{}:{}' Compton process \
                            (expected '{}' or '{}', found '{}')",
                        model,
                        method,
                        ComptonMode::Adjoint,
                        ComptonMode::Direct,
                        mode,
                    ),
                    _ => (),
                },
            };
        },
    };
    Ok(())
}
