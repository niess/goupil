use anyhow::Result;
use crate::numerics::{
    consts::PI,
    float::Float,
    integrate::{GQIntegrator, LogSubstitution}
};
use crate::physics::{
    consts::{ELECTRON_MASS, ELECTRON_RADIUS},
    materials::electronic::ElectronicStructure
};
use serde_derive::{Deserialize, Serialize};
use super::{
    ComptonModel::{self, KleinNishina, Penelope, ScatteringFunction},
    ComptonMode::{self, Adjoint, Direct, Inverse}
};


// ===============================================================================================
// Computer utility for Compton process.
// ===============================================================================================

#[derive(Clone, Copy, Deserialize, Serialize)]
pub struct ComptonComputer {
    pub(crate) model: ComptonModel,
    pub(crate) mode: ComptonMode,
    pub(crate) precision: Float,
}

impl Default for ComptonComputer {
    fn default() -> Self {
        Self { model: ComptonModel::default(), mode: ComptonMode::default(), precision: 1.0 }
    }
}

impl ComptonComputer {

    pub fn new(model: ComptonModel, mode: ComptonMode) -> Self {
        Self { model, mode, precision: 1.0 }
    }

    // Compute the (restricted) cross-section, in cm^2.
    pub fn cross_section(&self, energy: Float, energy_min: Option<Float>,
        energy_max: Option<Float>, electrons: &ElectronicStructure) -> Result<Float> {

        if let ComptonMode::None = self.mode {
            return Ok(0.0)
        }

        #[cfg(feature = "python")]
        crate::python::ctrlc_catched()?;

        let result = match self.model {
            ScatteringFunction | Penelope => {
                let n = (1000.0 / self.precision) as usize;
                self.effective_cross_section(energy, energy_min, energy_max, electrons, n)
            },
            KleinNishina => electrons.charge() * self.free_cross_section(
                self.mode, energy, energy_min, energy_max),
        };
        Ok(result)
    }

    // Compute the differential cross-section w.r.t. outgoing energy, in cm^2 / MeV.
    pub fn dcs(&self, energy_in: Float, energy_out: Float, electrons: &ElectronicStructure)
        -> Float {

        if let ComptonMode::None = self.mode {
            return 0.0
        }

        match self.model {
            Penelope => 0.0,
            ScatteringFunction => self.effective_dcs(energy_in, energy_out, electrons),
            KleinNishina => electrons.charge() * self.free_dcs(energy_in, energy_out),
        }
    }

    // DCS support (for free and effective models, only).
    pub(crate) fn dcs_support(&self, energy: Float) -> (Float, Float) {
        match self.model {
            ScatteringFunction | KleinNishina => self.free_dcs_support(energy),
            _ => (energy, energy),
        }
    }

    // Pre-factor for cross-sections.
    const CS_FACTOR: Float = PI * ELECTRON_RADIUS * ELECTRON_RADIUS;

    // Adjoint sampling weight for KleinNishina model.
    //
    // Note that this is simply the ratio of cross-sections times the regularisation term, since
    // DCSs simplify out.
    pub(super) fn free_adjoint_weight(&self, energy_in: Float, energy_out: Float) -> Float {
        self.free_cross_section(Adjoint, energy_in, None, None) *  energy_out /
            (self.free_cross_section(Direct, energy_out, None, None) * energy_in)
    }

    // Total cross-section, in the rest frame of a free target electron.
    pub(crate) fn free_cross_section(&self, mode: ComptonMode,  energy: Float,
        energy_min: Option<Float>, energy_max: Option<Float>) -> Float {

        let x = energy / ELECTRON_MASS;
        if energy_min.is_none() && energy_max.is_none() {
            match mode {
                Adjoint => {
                    if x <= 0.5 {
                        let x2 = x * x;
                        8.0 * Self::CS_FACTOR / 3.0 * (x2 - x + 1.0)
                    } else {
                        let x4 = 4.0 * x;
                        Self::CS_FACTOR * (x4 * (x4 - 1.0) + 1.0) / (3.0 * x4 * x * x)
                    }
                },
                Inverse | Direct => {
                    let tmp0 = 1.0 + 2.0 * x;
                    let tmp1 = 1.0 / x;
                    Self::CS_FACTOR / x * ((1.0 - 2.0 * tmp1 - 2.0 * tmp1 * tmp1) * tmp0.ln() +
                        0.5 + 4.0 * tmp1 - 0.5 / (tmp0 * tmp0))
                },
                ComptonMode::None => unreachable!(),
            }
        } else {
            let (emin, emax) = self.free_dcs_support(energy);
            let xmin = {
                let emin = match energy_min {
                    None => emin,
                    Some(e) => e.max(emin),
                };
                emin / ELECTRON_MASS
            };
            let xmax = {
                let emax = match energy_max {
                    None => emax,
                    Some(e) => e.min(emax),
                };
                emax / ELECTRON_MASS
            };

            let f = match mode {
                Adjoint => |x: Float, y: Float| -> Float {
                    -1.0 / y + (2.0 * x - 1.0) / (2.0 * y.powi(2) * x) -
                        (x.powi(2) + 2.0 * x - 2.0) / (3.0 * y.powi(3)) - x / (4.0 * y.powi(4))
                },
                Inverse | Direct => |x: Float, y: Float| -> Float {
                    (y * y / (2.0 * x) + (1.0 + 2.0 * x) * y / (x * x) +
                        (x - 2.0 - 2.0 / x) * y.ln() - 1.0 / y) / (x * x)
                },
                ComptonMode::None => unreachable!(),
            };
            Self::CS_FACTOR * (f(x, xmax) - f(x, xmin))
        }
    }

    // Differential Cross-Section (DCS), in the rest frame of a free target electron.
    fn free_dcs(&self, energy_in: Float, energy_out: Float) -> Float {

        // Lorentz invariant DCS for collisions with free electrons.
        let invariant_dcs = |x0: Float, x1: Float| -> Float {
            let tmp0 = 1.0 / x1 - 1.0 / x0 - 1.0;
            let tmp1 = x0 / x1 + x1 / x0 + tmp0 * tmp0 - 1.0;
            Self::CS_FACTOR * tmp1 / (x0 * x0)
        };

        let (emin, emax) = self.free_dcs_support(energy_in);
        if (energy_out < emin) || (energy_out > emax) { return 0.0 }

        match self.mode {
            Adjoint => {
                invariant_dcs(energy_out / ELECTRON_MASS, energy_in / ELECTRON_MASS) *
                    energy_in / (energy_out * ELECTRON_MASS)
                },
            Inverse | Direct => {
                invariant_dcs(energy_in / ELECTRON_MASS, energy_out / ELECTRON_MASS) /
                    ELECTRON_MASS
            },
            ComptonMode::None => unreachable!(),
        }
    }

    fn free_dcs_support(&self, energy: Float) -> (Float, Float) {
        match self.mode {
            Adjoint | Inverse => if energy < 0.5 * ELECTRON_MASS {
                (energy, energy * ELECTRON_MASS / (ELECTRON_MASS - 2.0 * energy))
            } else {
                (energy, (energy * 1e+3).max(1e+3))
            },
            Direct => {
                (energy * ELECTRON_MASS / (ELECTRON_MASS + 2.0 * energy), energy)
            },
            ComptonMode::None => unreachable!(),
        }
    }

    // Restricted cross-section for effective model.
    //
    // Note that the unresricted total cross-section matches Penelope result.
    fn effective_cross_section(
        &self, energy_in: Float, energy_min: Option<Float>, energy_max: Option<Float>,
        electrons: &ElectronicStructure, n: usize) -> Float  {

        let (mut xmin, mut xmax) = match self.mode {
            Adjoint | Direct => self.free_dcs_support(energy_in),
            Inverse => {
                let computer = Self::new(ScatteringFunction, Direct);
                computer.free_dcs_support(energy_in)
            },
            ComptonMode::None => unreachable!(),
        };

        if let Some(v) = energy_min {
            if xmin < v {
                if v >= xmax { return 0.0 }
                xmin = v;
            }
        }

        if let Some(v) = energy_max {
            if xmax > v {
                if v <= xmin { return 0.0 }
                xmax = v;
            }
        }

        let integrator = GQIntegrator::<12>::new();
        let integrand = |x| { self.effective_dcs(energy_in, x, electrons) };
        integrator.integrate(integrand, xmin, xmax, n, LogSubstitution)
    }

    pub(super) fn effective_charge(
        &self, energy_in: Float, energy_out: Float, electrons: &ElectronicStructure) -> Float {

        let (w, v) = match self.mode {
            Adjoint => (energy_out, energy_in),
            Inverse | Direct => (energy_in, energy_out),
            ComptonMode::None => unreachable!(),
        };

        let effective_fraction = |shell_momentum: Float, pz: Float| -> Float {
            let z = pz.abs();
            let x = 1.0 + 2.0 * z / shell_momentum;
            let tmp = if x < 11.0 { 0.5 * (0.5 * (1.0 - x * x)).exp() } else { 0.0 };
            if pz >= 0.0 { 1.0 - tmp } else { tmp }
        };

        let mut z = 0.0;
        for shell in electrons.iter() {
            let u = shell.energy;
            if w <= u { continue }
            let tmp = w * (w - u) * (ELECTRON_MASS / v - ELECTRON_MASS / w);
            let pmax = (tmp - ELECTRON_MASS * u) / (2.0 * tmp + u * u).sqrt();
            z += effective_fraction(shell.momentum, pmax) * shell.occupancy as Float;
        }
        z
    }

    fn effective_dcs(&self, energy_in: Float, energy_out: Float, electrons: &ElectronicStructure)
        -> Float {

        self.free_dcs(energy_in, energy_out) *
            self.effective_charge(energy_in, energy_out, electrons)
    }
}
