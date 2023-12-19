use anyhow::{anyhow, Error, Result};
use crate::numerics::float::Float;
use crate::physics::{
    consts::AVOGADRO_NUMBER,
    materials::MaterialRecord,
};
use crate::pretty_enumerate;
use enum_iterator::{all, Sequence};
use serde_derive::{Deserialize, Serialize};
use std::fmt;

pub(crate) mod table;


// ===============================================================================================
// Public API.
// ===============================================================================================

pub use self::table::AbsorptionCrossSection;


// ===============================================================================================
// Photo-absorption mode.
// ===============================================================================================

#[derive(Clone, Copy, Default, Deserialize, PartialEq, Sequence, Serialize)]
pub enum AbsorptionMode {
    /// Treat photo-absorption as a continuous process.
    ///
    /// Monte Carlo events are weighted by the photo-absorption probability over the path length.
    /// Note that this might sometimes result in very unlikely tracks, with tiny weights.
    #[default]
    Continuous,

    /// Treat photo-absorption as a discrete process.
    ///
    /// This is the most detailled method for simulating photo-absorption processes. Note however
    /// that secondary photo-electrons, or pair production, are not generated.
    Discrete,

    /// Photo absorptive processes are not simulated.
    None,
}

impl AbsorptionMode {
    const CONTINUOUS: &str = "Continuous";
    const DISCRETE: &str = "Discrete";
    const NONE: &str = "None";

    pub(crate) fn transport_cross_section(
        &self,
        energy: Float,
        material: &MaterialRecord
    ) -> Result<Float> {
        match self {
            Self::Discrete => match &material.table.absorption {
                None => Err(Self::no_table(material)),
                Some(table) => Ok(table.interpolate(energy)),
            },
            _ => Ok(0.0),
        }
    }

    pub(crate) fn transport_weight(
        energy: Float,
        column_depth: Float,
        material: &MaterialRecord
    ) -> Result<Float> {
        match &material.table.absorption {
            None => Err(Self::no_table(material)),
            Some(table) => {
                let cs = table.interpolate(energy);
                let weight = if cs <= 0.0 {
                    1.0
                } else {
                    (-column_depth * cs * AVOGADRO_NUMBER / material.definition.mass()).exp()
                };
                Ok(weight)
            },
        }
    }

    fn no_table(material: &MaterialRecord) -> Error {
        anyhow!(
            "{}: no table for absorption cross-section",
            material.definition.name()
        )
    }

    fn pretty_variants() -> String {
        let variants: Vec<_> = all::<Self>()
            .map(|e| format!("'{}'", e))
            .collect();
        pretty_enumerate(&variants)
    }
}

impl fmt::Display for AbsorptionMode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s: &str = (*self).into();
        write!(f, "{}", s)
    }
}

impl TryFrom<&str> for AbsorptionMode {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            Self::CONTINUOUS => Ok(Self::Continuous),
            Self::DISCRETE => Ok(Self::Discrete),
            Self::NONE => Ok(Self::None),
            _ => Err(anyhow!(
                "bad absorption mode (expected {}, found '{}')",
                Self::pretty_variants(),
                value,
            )),
        }
    }
}

impl From<AbsorptionMode> for &str {
    fn from(value: AbsorptionMode) -> Self {
        match value {
            AbsorptionMode::Continuous => AbsorptionMode::CONTINUOUS,
            AbsorptionMode::Discrete => AbsorptionMode::DISCRETE,
            AbsorptionMode::None => AbsorptionMode::NONE,
        }
    }
}
