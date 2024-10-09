use anyhow::{anyhow, Result};
use crate::pretty_enumerate;
use enum_iterator::{all, Sequence};
use serde_derive::{Deserialize, Serialize};
use std::fmt;

pub(crate) mod sample;
pub(crate) mod table;


// ===============================================================================================
// Public API.
// ===============================================================================================

pub use self::table::RayleighCrossSection;
pub use self::table::RayleighFormFactor;


// ===============================================================================================
// Rayleigh scattering mode.
// ===============================================================================================

#[derive(Clone, Copy, Default, Deserialize, PartialEq, Sequence, Serialize)]
pub enum RayleighMode {
    #[default]
    FormFactor,
    None,
}

impl RayleighMode {
    const FORM_FACTOR: &str = "Form Factor";
    const NONE: &str = "None";

    fn pretty_variants() -> String {
        let variants: Vec<_> = all::<Self>()
            .map(|e| format!("'{}'", e))
            .collect();
        pretty_enumerate(&variants)
    }
}

impl fmt::Display for RayleighMode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s: &str = (*self).into();
        write!(f, "{}", s)
    }
}

impl TryFrom<&str> for RayleighMode {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            Self::FORM_FACTOR => Ok(Self::FormFactor),
            Self::NONE => Ok(Self::None),
            _ => Err(anyhow!(
                "bad Rayleigh scattering mode (expected {}, found '{}')",
                Self::pretty_variants(),
                value,
            )),
        }
    }
}

impl From<RayleighMode> for &str {
    fn from(value: RayleighMode) -> Self {
        match value {
            RayleighMode::FormFactor => RayleighMode::FORM_FACTOR,
            RayleighMode::None => RayleighMode::NONE,
        }
    }
}
