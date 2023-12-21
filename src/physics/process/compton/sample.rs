use anyhow::{anyhow, bail, Error, Result};
use crate::numerics::{
    consts::{FRAC_1_SQRT_2, PI, SQRT_2},
    float::{Float, Float3},
    rand::FloatRng
};
use crate::physics::{
    consts::ELECTRON_MASS,
    materials::{electronic::ElectronicStructure, MaterialRecord}
};
use super::{
    ComptonModel::{self, KleinNishina, Penelope, ScatteringFunction},
    compute::ComptonComputer,
    ComptonMode::{self, Adjoint, Direct, Inverse},
    ComptonMethod::{self, InverseTransform, RejectionSampling}
};


// ===============================================================================================
// Sampler utility for Compton process.
// ===============================================================================================

#[derive(Clone, Copy)]
pub struct ComptonSampler {
    pub model: ComptonModel,
    pub mode: ComptonMode,
    pub method: ComptonMethod,
}

#[derive(Clone, Copy, Default)]
pub struct ComptonSample {
    pub momentum_out: Float3,
    pub weight: Float,
    pub constrained: bool,
}

impl ComptonSample {
    pub fn new(momentum_out: Float3, weight: Float) -> Self {
        Self { momentum_out, weight, constrained: false }
    }
}

impl Default for ComptonSampler {
    fn default() -> Self {
        Self { model: ScatteringFunction, mode: Direct, method: RejectionSampling }
    }
}

impl ComptonSampler {
    pub fn new(model: ComptonModel, mode: ComptonMode, method: ComptonMethod) -> Self {
        Self { model, mode, method }
    }

    pub fn sample<R: FloatRng>(
        &self,
        rng: &mut R,
        momentum_in: Float3,
        target: &MaterialRecord,
        energy_constrain: Option<Float>
    ) -> Result<ComptonSample> {
        if let ComptonMode::None = self.mode {
            return Ok(self.none_sample(momentum_in));
        }

        if !energy_constrain.is_none() {
            match self.model {
                Penelope => bail!(
                    "bad Compton model for energy constraint (expected '{}' or '{}', found '{}')",
                    ScatteringFunction,
                    KleinNishina,
                    self.model,
                ),
                _ => (),
            }
        }

        match self.method {
            InverseTransform => match self.model {
                Penelope =>
                    Err(self.bad_sampling_method(ComptonMethod::RejectionSampling)),
                ScatteringFunction | KleinNishina =>
                    self.inverse1(rng, momentum_in, target, energy_constrain),
            },
            RejectionSampling => match self.model {
                ScatteringFunction =>
                    self.effective_rejection(rng, momentum_in, target, energy_constrain),
                KleinNishina => self.free_rejection(rng, momentum_in, target, energy_constrain),
                Penelope => self.penelope_sample(rng, momentum_in, target),
            },
        }
    }

    pub fn transport_cross_section(
        &self,
        energy: Float,
        material: &MaterialRecord
    ) -> Result<Float> {
        if let ComptonMode::None = self.mode {
            return Ok(0.0);
        }
        match self.model {
            KleinNishina => {
                let computer = ComptonComputer::new(self.model, Direct);
                let electrons = self.get_electrons(material)?;
                let charge = electrons.charge();
                Ok(charge * computer.free_cross_section(Direct, energy, None, None))
            },
            _ => match &material.table.compton
                .get(self.model)
                .direct
                .cross_section {

                None => Err(anyhow!(
                    "{}: no cross-section table for {}:{} Compton process",
                    material.definition.name(),
                    self.model,
                    Direct,
                )),
                Some(table) => Ok(table.interpolate(energy)),
            }
        }
    }

    fn none_sample(&self, momentum_in: Float3) -> ComptonSample {
        ComptonSample::new(momentum_in, 1.0)
    }

    // Sample momentum direction rotatiuon in a free (or effective) collision.
    fn free_sample_direction<R: FloatRng>(&self, rng: &mut R, momentum_in: Float3,
        energy_in: Float, energy_out: Float) -> Float3 {

        let tmp = match self.mode {
            Adjoint | Inverse => ELECTRON_MASS / energy_in - ELECTRON_MASS / energy_out,
            Direct => ELECTRON_MASS / energy_out - ELECTRON_MASS / energy_in,
            ComptonMode::None => unreachable!(),
        };
        let cos_theta = (1.0 - tmp).clamp(-1.0, 1.0);
        let phi = 2.0 * PI * rng.uniform01();
        let mut momentum_out = momentum_in * (energy_out / energy_in);
        momentum_out.rotate(cos_theta, phi);
        momentum_out
    }

    // Sample a collision with a free electron (in the rest frame of the electron).
    fn free_rejection<R: FloatRng>(
        &self,
        rng: &mut R,
        momentum_in: Float3,
        target: &MaterialRecord,
        energy_constrain: Option<Float>
    ) -> Result<ComptonSample> {

        match self.mode {
            Adjoint => self.free_rejection_adjoint(rng, momentum_in, target, energy_constrain),
            Inverse => return Err(self.bad_sampling_mode(
                ComptonMode::Adjoint,
                Some(ComptonMode::Direct)
            )),
            Direct => Ok(self.free_rejection_forward(rng, momentum_in)),
            ComptonMode::None => unreachable!(),
        }
    }

    fn free_rejection_adjoint<R: FloatRng>(
        &self,
        rng: &mut R,
        momentum_in: Float3,
        target: &MaterialRecord,
        energy_constrain: Option<Float>
    ) -> Result<ComptonSample> {

        let mut constrain = EnergyConstrainChecker::from_option(self, energy_constrain);
        loop {
            let mut momentum_out = self.free_rejection_adjoint_raw(rng, momentum_in);

            let energy_in = momentum_in.norm();
            let energy_out = momentum_out.norm();

            let (constrained, weight) = match constrain.check(energy_in, energy_out, target)? {
                EnergyConstrain::None => {
                    let computer = ComptonComputer::default();
                    let weight = computer.free_adjoint_weight(energy_in, energy_out);
                    (false, weight)
                },
                EnergyConstrain::Some((energy_constrain, forced_weight)) => {
                    momentum_out = self.free_sample_direction(
                        rng, momentum_in, energy_in, energy_constrain);
                    (true, forced_weight)
                },
                EnergyConstrain::TryAgain => { continue },
            };
            return Ok(ComptonSample { momentum_out, weight, constrained })
        }
    }

    fn free_rejection_adjoint_raw<R: FloatRng>(&self, rng: &mut R, momentum_in: Float3)
        -> Float3  {

        let energy_in = momentum_in.norm();
        let x1 = energy_in / ELECTRON_MASS;
        let (a1, q) = if x1 < 0.5 {
            let x12 = x1 * x1;
            let a = (8.0 * x12 - 12.0 * x1 + 6.0) / x12;
            let b = x12 * a;
            (a, b / (b + 6.0))
        } else {
            (0.0, 0.0)
        };

        let generator = |x1: Float, rng: &mut R| -> Float {
            if x1 < 0.5 {
                let u = rng.uniform01();
                if u < q {
                    let uu = u / q;
                    (1.0 / (x1 * x1 * x1) - uu * a1).powf(-1.0 / 3.0)
                } else {
                    let uu = (u - q) / (1.0 - q);
                    x1 / (1.0 - 2.0 * uu * x1)
                }
            } else {
                let u = rng.uniform01();
                if u <= 0.25 {
                    let uu = u / 0.25;
                    x1 * uu.powf(-1.0 / 3.0)
                } else {
                    let uu = (u - 0.25) / 0.75;
                    x1 / uu
                }
            }
        };

        self.reject_sample(generator, x1, momentum_in, rng)
    }

    fn free_rejection_forward<R: FloatRng>(&self, rng: &mut R, momentum_in: Float3)
        -> ComptonSample {

        let momentum_out = self.free_rejection_forward_raw(rng, momentum_in);
        ComptonSample::new(momentum_out, 1.0)
    }

    fn free_rejection_forward_raw<R: FloatRng>(&self, rng: &mut R, momentum_in: Float3)
        -> Float3 {

        let x0 = momentum_in.norm() / ELECTRON_MASS;
        let tmp0 = 1.0 + 2.0 * x0;
        let pa = x0 * tmp0.ln();
        let pb = 2.0 * x0 * x0 * (1.0 + x0) / (tmp0 * tmp0);

        let generator = |x0: Float, rng: &mut R| -> Float {
            let u = rng.uniform01();
            let r = pa / (pa + pb);
            if u < r {
                    let v = u / r;
                    x0 * tmp0.powf(-v)
            } else {
                let tmp1 = 1.0 / (tmp0 * tmp0);
                let v = (u - r) / (1.0 - r);
                x0 * (tmp1 + v * (1.0 - tmp1)).sqrt()
            }
        };

        self.reject_sample(generator, x0, momentum_in, rng)
    }

    fn reject_sample<F, R>(&self, generator: F, x: Float, momentum_in: Float3, rng: &mut R)
        -> Float3
    where
        F: Fn(Float, &mut R) -> Float,
        R: FloatRng,
    {
        loop {
            let (x0, x1, xalt) = match self.mode {
                Adjoint => {
                    let tmp = (generator)(x, rng);
                    (tmp, x, tmp)
                },
                Inverse => { unreachable!() },
                Direct => {
                    let tmp = (generator)(x, rng);
                    (x, tmp, tmp)
                },
                ComptonMode::None => unreachable!(),
            };
            let den = x1 / x0 + x0 / x1;
            let tmp = 1.0 / x1 - 1.0 / x0;
            let num = den - tmp * (2.0 - tmp);
            if rng.uniform01() * den <= num {
                let cos_theta = 1.0 - tmp;
                let phi = 2.0 * PI * rng.uniform01();
                let mut momentum_out = momentum_in;
                momentum_out *= xalt / x;
                momentum_out.rotate(cos_theta, phi);
                return momentum_out;
            }
        }
    }

    fn effective_rejection<R: FloatRng>(
        &self,
        rng: &mut R,
        momentum_in: Float3,
        target: &MaterialRecord,
        energy_constrain: Option<Float>
    ) -> Result<ComptonSample> {

        let electrons = self.get_electrons(target)?;
        let computer = ComptonComputer::new(self.model, self.mode);
        let energy_in = momentum_in.norm();
        let zmax = match self.mode {
            Adjoint => if energy_in < 0.5 * ELECTRON_MASS {
                let max_out = energy_in * ELECTRON_MASS / (ELECTRON_MASS - 2.0 * energy_in);
                computer.effective_charge(energy_in, max_out, electrons)
            } else {
                electrons.charge()
            },
            Inverse => return Err(self.bad_sampling_mode(
                ComptonMode::Adjoint,
                Some(ComptonMode::Direct),
            )),
            Direct => {
                let min_out = energy_in * ELECTRON_MASS / (ELECTRON_MASS + 2.0 * energy_in);
                computer.effective_charge(energy_in, min_out, electrons)
            },
            ComptonMode::None => unreachable!(),
        };

        let mut constrain = EnergyConstrainChecker::from_option(self, energy_constrain);
        loop {
            let mut momentum_out = match self.mode {
                Adjoint => self.free_rejection_adjoint_raw(rng, momentum_in),
                Direct => self.free_rejection_forward_raw(rng, momentum_in),
                _ => unreachable!(),
            };
            let energy_out = momentum_out.norm();
            let r = computer.effective_charge(energy_in, energy_out, electrons) / zmax;
            if rng.uniform01() <= r {
                let (constrained, weight) = match constrain.check(energy_in, energy_out, target)? {
                    EnergyConstrain::None => {
                        let weight = match self.mode {
                            Adjoint => target.table.compton.effective.adjoint_weight(
                                ScatteringFunction,
                                energy_in,
                                energy_out
                            )?,
                            Direct => 1.0,
                            _ => unreachable!(),
                        };
                        (false, weight)
                    },
                    EnergyConstrain::Some((energy_constrain, forced_weight)) => {
                        momentum_out = self.free_sample_direction(
                            rng, momentum_in, energy_in, energy_constrain);
                        (true, forced_weight)
                    },
                    EnergyConstrain::TryAgain => { continue },
                };
                return Ok(ComptonSample { momentum_out, weight, constrained })
            }
        }
    }

    // Inverse sampling for a 1d distribution, depending only on the energy.
    //
    // Note: this method is valid only for the ScatteringFunction and KleinNishina models.
    fn inverse1<R: FloatRng>(
        &self,
        rng: &mut R,
        momentum_in: Float3,
        material: &MaterialRecord,
        energy_constrain: Option<Float>
    ) -> Result<ComptonSample> {

        let subtable = material.table.compton.get(self.model);
        let table = {
            let tmp = subtable.get(self.mode).inverse_cdf.as_ref();
            match tmp {
                None => bail!(
                    "{}: no inverse CDF table for {}:{} Compton process",
                    material.definition.name(),
                    self.model,
                    self.mode,
                ),
                Some(v) => v,
            }
        };
        let energy_in = momentum_in.norm();
        let mut constrain = EnergyConstrainChecker::from_option(self, energy_constrain);
        loop {
            let cdf = rng.uniform01();
            let mut constrained = false;
            let (energy_out, weight) = {
                let (mut energy_out, mut weight) = table.interpolate(energy_in, cdf);
                match constrain.check(energy_in, energy_out, material)? {
                    EnergyConstrain::None => if self.mode == Adjoint {
                        weight = subtable.adjoint_weight(
                            self.model,
                            energy_in,
                            energy_out
                        )?;
                    },
                    EnergyConstrain::Some((energy_constrain, forced_weight)) => {
                        energy_out = energy_constrain;
                        weight = forced_weight;
                        constrained = true;
                    },
                    EnergyConstrain::TryAgain => { continue },
                }
                (energy_out, weight)
            };
            let momentum_out = self.free_sample_direction(rng, momentum_in, energy_in, energy_out);

            return Ok(ComptonSample{ momentum_out, weight, constrained })
        }
    }

    // Penelope like sampler for Compton collisions with an atom.
    //
    // This is an implementation of Compton collisions according to the Penelope-2014 manual
    // (F.Salvat (2015), Penelope 2014 Workshop, Barcelona, section 2.3.2).
    #[allow(non_snake_case)]
    fn penelope_sample<R: FloatRng>(
        &self,
        rng: &mut R,
        momentum_in: Float3,
        material: &MaterialRecord
    ) -> Result<ComptonSample> {

        let electrons = self.get_electrons(material)?;

        if self.mode != Direct {
            return Err(self.bad_sampling_mode(
                ComptonMode::Direct,
                None
            ))
        }
        if self.method != RejectionSampling {
            return Err(self.bad_sampling_method(ComptonMethod::RejectionSampling))
        }

        // Sample cos(theta).
        let E = momentum_in.norm();
        let kappa = E / ELECTRON_MASS;
        let tmp0 = 1.0 + 2.0 * kappa;
        let tau_min = 1.0 / tmp0;
        let tau_min2 = tau_min * tau_min;
        let a1 = tmp0.ln();
        let a2 = 2.0 * kappa * (1.0 + kappa) * tau_min2;
        let p = a1 / (a1 + a2);
        let (Smax, _) = self.penelope_scattering_terms(E, -1.0, electrons);
        loop {
            let u = rng.uniform01();
            let tau = if u < p {
                let v = u / p;
                tau_min.powf(v)
            } else {
                let v = (u - p) / (1.0 - p);
                (tau_min2 + v * (1.0 - tau_min2)).sqrt()
            };
            let cos_theta = 1.0 - (1.0 - tau) / (kappa * tau);
            let (S, terms) = self.penelope_scattering_terms(E, cos_theta, electrons);
            let r = (1.0 - (1.0 - tau) * (tmp0 * tau - 1.0) /
                (kappa * kappa * tau * (1.0 + tau * tau))) * S / Smax;
            if rng.uniform01() <= r {
                // Sample target pz.
                let pz = self.penelope_sample_pz_raw(rng, E, cos_theta, electrons, S, &terms);

                // Compute final energy from pz.
                let Ep = self.penelope_final_energy(E, cos_theta, pz);

                // Build outgoing momentum.
                let mut momentum_out = momentum_in;
                momentum_out *= Ep / E;
                let phi = 2.0 * PI * rng.uniform01();
                momentum_out.rotate(cos_theta, phi);
                return Ok(ComptonSample::new(momentum_out, 1.0));
            }
        }
    }

    #[allow(non_snake_case)]
    fn penelope_final_energy(&self, E: Float, cos_theta: Float, pz: Float) -> Float {
        let t = (pz * pz) / (ELECTRON_MASS * ELECTRON_MASS);
        let kappa = E / ELECTRON_MASS;
        let tau = 1.0 / (1.0 + kappa * (1.0 - cos_theta));
        let sgn = if pz >= 0.0 { 1.0 } else { -1.0 };
        let tmp1 = 1.0 - t * tau * tau;
        let tmp2 = 1.0 - t * tau * cos_theta;
        E * tau / tmp1 * (tmp2 + sgn * (tmp2 * tmp2 - tmp1 * (1.0 - t)).sqrt())
    }

    #[allow(non_snake_case)]
    fn penelope_max_momentum(&self, E: Float, cos_theta: Float, U: Float) -> Float {
        let tmp = E * (E - U) * (1.0 - cos_theta);
        (tmp - ELECTRON_MASS * U) / (2.0 * tmp + U * U).sqrt()
    }

    fn penelope_neff(&self, shell_momentum: Float, pz: Float) -> Float {
        let z = pz.abs();
        let x = 1.0 + 2.0 * z / shell_momentum;
        let tmp = 0.5 * (0.5 * (1.0 - x * x)).exp();
        if pz >= 0.0 { 1.0 - tmp } else { tmp }
    }

    #[allow(non_snake_case)]
    fn penelope_sample_pz_raw<R: FloatRng>(&self, rng: &mut R, E: Float, cos_theta: Float,
        electrons: &ElectronicStructure, S: Float, terms: &[PenelopeScatteringTerm]) -> Float {

        let try_sample_pz = |rng: &mut R| -> Float {
            // Sample targeted shell.
            let (pmax, shell) = 'outer: loop {
                let v = S * rng.uniform01();
                let mut sum = 0.0;
                for (i, term) in terms.iter().enumerate() {
                    sum += term.S;
                    if v <= sum {
                        break 'outer (term.pmax, &electrons[i]);
                    }
                }
            };

            // Sample target momentum (z-component, actually).
            let A = self.penelope_neff(shell.momentum, pmax) * rng.uniform01();
            const D1: Float = FRAC_1_SQRT_2;
            const D2: Float = SQRT_2;
            if A < 0.5 {
                shell.momentum / D2 * (D1 - (D1 * D1 - (2.0 * A).ln()).sqrt())
            } else {
                shell.momentum / D2 * ((D1 * D1 - (2.0 * (1.0 - A)).ln()).sqrt() - D1)
            }
        };

        // Rejection sampling of pz.
        let kappa = E / ELECTRON_MASS;
        let Ec = E / (1.0 + kappa * (1.0 - cos_theta));
        let qc = (E * E + Ec * Ec - 2.0 * E * Ec * cos_theta).sqrt();
        let aF = qc / E * (1.0 + Ec * (Ec - E * cos_theta) / (qc * qc));
        loop {
            let pz = try_sample_pz(rng);
            if pz < -ELECTRON_MASS { continue }

            // Rejection sampling step, as F / Fmax.
            let F = 1.0 + aF * pz.clamp(-0.2 * ELECTRON_MASS, 0.2 * ELECTRON_MASS) / ELECTRON_MASS;
            let Fmax = 1.0 + aF.abs() * 0.2;
            if rng.uniform01() * Fmax <= F { break pz }
        }
    }

    #[allow(non_snake_case)]
    fn penelope_scattering_terms(&self, E: Float, cos_theta: Float,
        electrons: &ElectronicStructure) -> (Float, Vec<PenelopeScatteringTerm>) {

        let mut terms = vec![PenelopeScatteringTerm::default(); electrons.len()];
        let mut S = 0.0;
        for (shell, term) in std::iter::zip(electrons.iter(), &mut terms) {
            let U = shell.energy;
            if E <= U { continue }
            let pmax = self.penelope_max_momentum(E, cos_theta, U);
            let Si = self.penelope_neff(shell.momentum, pmax) * shell.occupancy as Float;
            S += Si;
            term.pmax = pmax;
            term.S = Si;
        }
        (S, terms)
    }

    // Formated errors.
    fn bad_sampling_method(&self, expected: ComptonMethod) -> Error {
        anyhow!(
            "bad sampling method for {}:{} Compton process (expected '{}', found '{}')",
            self.model,
            self.mode,
            expected,
            self.method,
        )
    }

    fn bad_sampling_mode(&self, expected0: ComptonMode, expected1: Option<ComptonMode>)
        -> Error {

        let expected = match expected1 {
            None => format!(
                "'{}'",
                expected0
            ),
            Some(expected1) => format!(
                "'{}' or '{}'",
                expected0,
                expected1
            ),
        };

        anyhow!("bad sampling mode for {}:{} Compton process (expected '{}', found '{}')",
            self.model,
            self.method,
            expected,
            self.mode,
        )
    }

    fn get_electrons<'a>(
        &self, material: &'a MaterialRecord
    ) -> Result<&'a ElectronicStructure> {
        material
            .electrons()
            .ok_or_else(|| anyhow!(
                "{}: no electronic structure for {} Compton process",
                material.definition.name(),
                self.model,
            ))
    }
}

#[allow(non_snake_case)]
#[derive(Clone, Copy, Default)]
struct PenelopeScatteringTerm {
    pmax: Float,
    S: Float,
}


// ===============================================================================================
// LOw level utility for checking energy constrain.
// ===============================================================================================

enum EnergyConstrainChecker {
    None,
    Some(ComptonSampler, Float, usize),
}

enum EnergyConstrain {
    None,
    Some((Float, Float)),
    TryAgain,
}

impl EnergyConstrainChecker {

    const MAX_TRIALS: usize = 20;

    fn check(
        &mut self,
        energy_in: Float,
        energy_out: Float,
        material: &MaterialRecord,
    ) -> Result<EnergyConstrain> {

        match self {
            Self::None => Ok(EnergyConstrain::None),
            Self::Some(sampler, energy_constrain, trials) => match sampler.mode {
                Direct => bail!(
                    "bad sampling mode for energy constraint (expected '{}' or '{}', found '{}')",
                     Adjoint,
                     Inverse,
                     sampler.mode,
                ),
                _ => if energy_out >= *energy_constrain {
                    // Compute weight of forced collision.
                    let (energy_in, energy_out, mode) = match sampler.mode {
                        Adjoint => (energy_in, *energy_constrain, Adjoint),
                        Inverse => (*energy_constrain, energy_in, Direct),
                        _ => unreachable!(),
                    };
                    let cdf = match sampler.model {
                        KleinNishina => {
                            // Compute CDF from analytical expression.
                            let computer = ComptonComputer::new(sampler.model, sampler.mode);
                            let cs0 = computer.free_cross_section(
                                sampler.mode, energy_in, None, None);
                            let cs1 = computer.free_cross_section(
                                sampler.mode, energy_in, None, Some(energy_out));
                            (cs1 / cs0).clamp(0.0, 1.0)
                        },
                        _ => match &material.table.compton
                            .get(sampler.model)
                            .get(mode)
                            .cdf {

                            None => bail!(
                                "{}: no forward CDF table for {}:{} Compton process",
                                material.definition.name(),
                                sampler.model,
                                mode,
                            ),
                            Some(table) => table.interpolate(energy_in, energy_out),
                        },
                    };

                    if cdf >= 1.0 || cdf <= 0.0 {
                        if *trials < Self::MAX_TRIALS {
                            *trials += 1;
                            Ok(EnergyConstrain::TryAgain)
                        } else {
                            Err(anyhow!(
                                "bad Compton CDF (expected a value in (0, 1), found {})",
                                cdf
                            ))
                        }
                    } else {
                        let (energy_in, energy_out, prob) = match sampler.mode {
                            Adjoint => (energy_out, energy_in, 1.0 - cdf),
                            Inverse => (energy_in, energy_out, cdf),
                            _ => unreachable!(),
                        };
                        let cross_section = sampler.transport_cross_section(energy_in, material)?;
                        let pdf = if cross_section <= 0.0 {
                            0.0
                        } else {
                            let computer = ComptonComputer::new(sampler.model, Direct);
                            let electrons = sampler.get_electrons(material)?;
                            let dcs = computer.dcs(energy_in, energy_out, electrons);
                            dcs / cross_section
                        };
                        let weight = pdf / prob;
                        Ok(EnergyConstrain::Some((*energy_constrain, weight)))
                    }
                } else {
                    Ok(EnergyConstrain::None)
                },
            },
        }
    }

    fn from_option(
        sampler: &ComptonSampler,
        energy_constrain: Option<Float>,
    ) -> Self {
        match energy_constrain {
            None => Self::None,
            Some(energy_constrain) => Self::Some(*sampler, energy_constrain, 0),
        }
    }
}
