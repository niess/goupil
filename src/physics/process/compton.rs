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
    ComptonTable,
    ComptonCrossSection,
    ComptonCDF,
    ComptonInverseCDF,
};


// ===============================================================================================
// Compton models.
// ===============================================================================================

#[derive(Clone, Copy, Default, Deserialize, PartialEq, Sequence, Serialize)]
pub enum ComptonModel {
    /// The target is modeled as an incoherent superposition of electronic shells, following
    /// the Impulse Approximation.
    ///
    /// This is the most accurate model, but also the most CPU-intensive. It is implemented only in
    /// forward mode. The Impulse Approximation (IA) is used for representing the electronic
    /// structure of the target, with Penelope's analytical parametrisation of momentum
    /// distributions. While the sampling of collisions is rather straightforward, computing the
    /// total cross-section is not. Thus, the latter is estimated by Monte Carlo integration.
    ImpulseApproximation,

    /// Klein-Nishina model.
    ///
    /// Target electrons are assumed to be free and at rest. This is the most approximate model. It
    /// yields the klein-Nishina cross-section, neglecting atomic binding effects, and Doppler
    /// broadening.
    KleinNishina,

    /// Penelope model for incoherent (Compton) scattering.
    ///
    /// As an approximation, the target electronic structure can be reduced to a single
    /// longitudinal parameter, J_i, the so called Compton profile, from which the DDCS w.r.t.
    /// outgoing energy and scattering angle is obtained. This model renders both atomic binding
    /// effects and Doppler broadening, with a good balance between speed and accuracy, as compared
    /// to the full computation (i.e. [`ImpulseApproximation`](Self::ImpulseApproximation)).
    /// However, it is implemented only in forward mode.
    Penelope,

    /// Effective model, using Klein-Nishina DCS with a Scattering Function corrective factor.
    ///
    /// This model is intermediate between [`Penelope`](Self::Penelope) and
    /// [`KleinNishina`](Self::KleinNishina). The DCS w.r.t. energy is computed following
    /// Penelope's Scattering Function, depending on the electronic structure, such that Penelope's
    /// total cross-section is recovered. However, the scattering angle is deterministic, assuming
    /// a collision with a free electron. Consequently, Doppler broadening (due to target electrons
    /// motion) is not rendered in the DCS, but only in the total cross-section. Atomic binding
    /// effects are however, taken into account in both the DCS and total cross-section.
    #[default]
    ScatteringFunction,
}

impl ComptonModel {
    const IMPULSE_APPROXIMATION: &str = "Impulse Approximation"; // XXX Remove this.
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
            Self::IMPULSE_APPROXIMATION => Ok(Self::ImpulseApproximation),
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
            ComptonModel::ImpulseApproximation => ComptonModel::IMPULSE_APPROXIMATION,
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
    InverseCDF,
    #[default]
    RejectionSampling,
}

impl ComptonMethod {
    const INVERSE_CDF: &str = "Inverse CDF";
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
            Self::INVERSE_CDF => Ok(Self::InverseCDF),
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
            ComptonMethod::InverseCDF => ComptonMethod::INVERSE_CDF,
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
        ComptonModel::ImpulseApproximation | ComptonModel::Penelope => {
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
                ComptonMethod::InverseCDF => bail!(
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
                ComptonMethod::InverseCDF => (),
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
