use anyhow::{anyhow, bail, Result};
use crate::numerics::consts::PI;
use crate::numerics::float::{Float, Float3};
use crate::numerics::rand::FloatRng;
use crate::physics::materials::MaterialRecord;
use crate::physics::process::rayleigh::table::RayleighFormFactor;
use super::RayleighMode::{self, FormFactor};

#[cfg(feature = "python")]
use crate::physics::consts::ELECTRON_RADIUS;


// ===============================================================================================
// Sampler utility for Rayleigh process.
// ===============================================================================================

#[derive(Clone, Copy)]
pub struct RayleighSampler {
    mode: RayleighMode,
}

// Public interface.
impl RayleighSampler {

    #[cfg(feature = "python")]
    pub(crate) fn dcs(
        &self,
        energy: Float,
        cos_theta: Float,
        form_factor: &RayleighFormFactor,
    ) -> Result<Float> {
        let dcs = match self.mode {
            FormFactor => {
                let q = energy * (2.0 * (1.0 - cos_theta)).sqrt();
                let r = ELECTRON_RADIUS * form_factor.interpolate(q);
                PI * r * r * (1.0 + cos_theta * cos_theta)
            },
            RayleighMode::None => 0.0,
        };
        Ok(dcs)
    }

    pub fn new(mode: RayleighMode) -> Self {
        Self { mode }
    }

    pub fn sample<R: FloatRng>(
        &self,
        rng: &mut R,
        energy_in: Float,
        direction_in: Float3,
        material: &MaterialRecord,
    ) -> Result<Float3> {
        match self.mode {
            FormFactor => {
                let form_factor = Self::form_factor(material)?;
                let cos_theta = self.sample_angle(rng, energy_in, &form_factor)?;
                let direction_out = if cos_theta == 0.0 {
                    direction_in
                } else {
                    let phi = 2.0 * PI * rng.uniform01();
                    let mut direction_out = direction_in.clone();
                    direction_out.rotate(cos_theta, phi);
                    direction_out
                };
                Ok(direction_out)
            },
            RayleighMode::None => Ok(direction_in.clone()),
        }
    }

    pub fn transport_cross_section(
        &self,
        energy: Float,
        material: &MaterialRecord
    ) -> Result<Float> {
        match self.mode {
            FormFactor => match material.rayleigh_cross_section() {
                None => Err(anyhow!(
                    "{}: no cross-section table for Rayleigh process",
                    material.definition.name(),
                )),
                Some(cross_section) => Ok(cross_section.interpolate(energy)),
            },
            RayleighMode::None => Ok(0.0),
        }
    }
}

// Private interface.
impl RayleighSampler {
    fn form_factor(material: &MaterialRecord) -> Result<&RayleighFormFactor> {
        match material.rayleigh_form_factor() {
            None => bail!(
                "{}: no form-factor table for Rayleigh process",
                material.definition.name(),
            ),
            Some(ff) => Ok(ff),
        }
    }

    pub(crate) fn sample_angle<R: FloatRng>(
        &self,
        rng: &mut R,
        energy: Float,
        form_factor: &RayleighFormFactor,
    ) -> Result<Float> {
        if let RayleighMode::None = self.mode {
            return Ok(0.0);
        }
        let q2max = 4.0 * energy * energy;
        if q2max <= 0.0 {
            return Ok(0.0);
        }
        loop {
            let u = rng.uniform01();
            let q2 = form_factor.inverse_cdf(u, q2max);
            let enveloppe = form_factor.enveloppe(q2);
            if enveloppe <= 0.0 {
                bail!(
                    "bad form factor enveloppe (expected a positive value, found {})",
                    enveloppe,
                )
            }
            let ff = form_factor.interpolate(q2.sqrt());
            let rf = ff / enveloppe;
            let cos_theta = (1.0 - 2.0 * q2 / q2max).clamp(-1.0, 1.0);
            let r = 0.5 * rf * rf * (1.0 + cos_theta * cos_theta);
            if rng.uniform01() < r {
                break Ok(cos_theta);
            }
        }
    }
}
