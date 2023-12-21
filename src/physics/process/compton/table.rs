use anyhow::{anyhow, bail, Result};
use crate::numerics::{
    float::Float,
    grids::{Grid, LinearGrid, LogGrid, UnstructuredGrid},
    interpolate::{BilinearInterpolator, CubicInterpolator},
};
use crate::physics::materials::electronic::ElectronicStructure;
use crate::physics::process::compton::{
    ComptonModel::{self, KleinNishina, Penelope, ScatteringFunction},
    ComptonMode::{self, Adjoint, Direct, Inverse},
};
use serde_derive::{Deserialize, Serialize};
use super::compute::ComptonComputer;


// ===============================================================================================
// Material specific data table for Compton process.
// ===============================================================================================

#[derive(Default, Deserialize, Serialize)]
pub struct ComptonTable {
    pub(crate) effective: ComptonSubTable,
    pub(crate) free: ComptonSubTable,
}

#[derive(Default, Deserialize, Serialize)]
pub(crate) struct ComptonSubTable {
    pub adjoint:  ComptonSubSubTable,
    pub direct: ComptonSubSubTable,
    pub inverse:  ComptonSubSubTable,
    pub none:  ComptonSubSubTable,
}

#[derive(Default, Deserialize, Serialize)]
pub(crate) struct ComptonSubSubTable {
    pub cdf: Option<ComptonCDF>,
    pub cross_section: Option<ComptonCrossSection>,
    pub inverse_cdf: Option<ComptonInverseCDF>,
}

impl ComptonTable {

    pub(crate) fn get(&self, model: ComptonModel) -> &ComptonSubTable {
        match model {
            ScatteringFunction | Penelope => &self.effective,
            KleinNishina => &self.free,
        }
    }

    pub(crate) fn get_mut(&mut self, model: ComptonModel) -> &mut ComptonSubTable {
        match model {
            ScatteringFunction | Penelope => &mut self.effective,
            KleinNishina => &mut self.free,
        }
    }

    #[allow(dead_code)] // XXX remove if not used.
    pub fn weighted_sum(tables: &[(Float, &Self)]) -> Self {

        let effective: Vec<_> = tables
            .iter()
            .map(|(fraction, table)| (*fraction, &table.effective))
            .collect();

        let free: Vec<_> = tables
            .iter()
            .map(|(fraction, table)| (*fraction, &table.free))
            .collect();

        Self {
            effective: ComptonSubTable::weighted_sum(&effective),
            free: ComptonSubTable::weighted_sum(&free),
        }
    }
}

impl ComptonSubTable {

    pub fn get(&self, mode: ComptonMode) -> &ComptonSubSubTable {
        match mode {
            Adjoint => &self.adjoint,
            Direct => &self.direct,
            Inverse => &self.inverse,
            ComptonMode::None => &self.none,
        }
    }

    pub fn get_mut(&mut self, mode: ComptonMode) -> &mut ComptonSubSubTable {
        match mode {
            Adjoint => &mut self.adjoint,
            Direct => &mut self.direct,
            Inverse => &mut self.inverse,
            ComptonMode::None => &mut self.none,
        }
    }

    #[allow(dead_code)] // XXX remove if not used.
    pub fn weighted_sum(tables: &[(Float, &Self)]) -> Self {

        let adjoint: Vec<_> = tables
            .iter()
            .map(|(fraction, table)| (*fraction, &table.adjoint))
            .collect();

        let direct: Vec<_> = tables
            .iter()
            .map(|(fraction, table)| (*fraction, &table.direct))
            .collect();

        let inverse: Vec<_> = tables
            .iter()
            .map(|(fraction, table)| (*fraction, &table.inverse))
            .collect();

        Self {
            adjoint: ComptonSubSubTable::weighted_sum(&adjoint),
            direct: ComptonSubSubTable::weighted_sum(&direct),
            inverse: ComptonSubSubTable::weighted_sum(&inverse),
            none: ComptonSubSubTable::default(),
        }
    }

    // Adjoint sampling weight (valid for effective & free models only).
    //
    // Note that this is simply the ratio of cross-sections times the regularisation term, since
    // DCSs simplify out, for effective and free models.
    pub(crate) fn adjoint_weight(&self, model: ComptonModel, energy_in: Float, energy_out: Float)
        -> Result<Float> {

        let no_table = |mode| {
            anyhow!("no cross-section table for {}:{} Compton process",
                model,
                mode,
            )
        };
        match &self.direct.cross_section {
            None => Err(no_table(ComptonMode::Direct)),
            Some(csf) => match &self.adjoint.cross_section {
                None => Err(no_table(ComptonMode::Adjoint)),
                Some(csa) => {
                    Ok(csa.interpolate(energy_in) * energy_out /
                        (csf.interpolate(energy_out) * energy_in))
                },
            },
        }
    }
}

impl ComptonSubSubTable {
    #[allow(dead_code)] // XXX remove if not used.
    pub fn weighted_sum(tables: &[(Float, &Self)]) -> Self {
        let cs: Vec<_> = tables
            .iter()
            .filter(|(fraction, table)| { *fraction > 0.0 && table.cross_section.is_some() })
            .map(|(fraction, table)| (*fraction, table.cross_section.as_ref().unwrap()))
            .collect();

        let cdf: Vec<_> = tables
            .iter()
            .filter(|(fraction, table)| {
                *fraction > 0.0 && table.cdf.is_some() && table.cross_section.is_some() })
            .map(|(fraction, table)| (
                *fraction,
                table.cdf.as_ref().unwrap(),
                table.cross_section.as_ref().unwrap()
            ))
            .collect();

        Self {
            cross_section: ComptonCrossSection::weighted_sum(&cs),
            cdf: ComptonCDF::weighted_sum(&cdf),
            inverse_cdf: None,
        }
    }
}


// ===============================================================================================
// Cross-section tabulation, and related data.
// ===============================================================================================

#[derive(Default, Deserialize, Serialize)]
pub struct ComptonCrossSection {
    pub(crate) energies: LogGrid,
    pub(crate) values: Vec<Float>,
    pub(crate) interpolator: CubicInterpolator,
    pub(super) computer: ComptonComputer,
}

// Public API.
impl ComptonCrossSection {
    pub fn interpolate(&self, energy: Float) -> Float {
        let (i, d) = self
            .energies
            .transform(energy)
            .clamp();
        let dx = self.energies.width(i);
        self.interpolator.interpolate_raw(i, d, dx, &self.values)
    }

    pub fn energies(&self) -> &[Float] {
        self.energies.as_ref()
    }

    pub fn len(&self) -> usize { self.values.len() }

    pub fn mode(&self) -> ComptonMode {
        self.computer.mode
    }

    pub fn model(&self) -> ComptonModel {
        self.computer.model
    }

    pub fn values(&self) -> &[Float] {
        &self.values
    }
}


// Private API.
impl ComptonCrossSection {
    pub(crate) fn compute(
        &mut self,
        computer: &ComptonComputer,
        electrons: &ElectronicStructure
    ) -> Result<()> {
        for (i, ei) in self.energies.as_ref().iter().enumerate() {
            self.values[i] = computer.cross_section(*ei, None, None, electrons)?;
        }
        self.initialise();
        self.computer = *computer;
        Ok(())
    }

    pub(crate) fn initialise(&mut self) {
        self.interpolator.initialise(&self.energies, &self.values, false);
    }

    pub(crate) fn new(energy_min: Float, energy_max: Float, n: usize) -> Self {
        let energies = LogGrid::new(energy_min, energy_max, n);
        let values = vec![0.0; n];
        let interpolator = CubicInterpolator::new(n);
        let computer = ComptonComputer::default();
        Self { energies, values, interpolator, computer }
    }

    // Merge a weighted sum of cross-section data into a single table.
    #[allow(dead_code)] // XXX remove if not used.
    pub(crate) fn weighted_sum(data: &[(Float, &Self)]) -> Option<Self> {

        if data.is_empty() {
            return None
        }

        // Get the intersection of supports.
        let emin = data
            .iter()
            .map(|(_, table)| table.energies[0])
            .reduce(Float::max)
            .unwrap();

        let emax = data
            .iter()
            .map(|(_, table)| table.energies[table.len() - 1])
            .reduce(Float::min)
            .unwrap();

        // Use the finest resolution for the merged grid.
        let resolution = data
            .iter()
            .map(|(_, table)| {
                (table.energies[table.len() - 1] - table.energies[0]) / (table.len() as Float)
            })
            .reduce(Float::min)
            .unwrap();

        // Merge cross-section data over a new grid.
        let n = ((emax - emin) / resolution).round() as usize;
        let mut cs = ComptonCrossSection::new(emin, emax, n);
        for (i, ei) in cs.energies.as_ref().iter().enumerate() {
            cs.values[i] = data
                .iter()
                .map(|(fj, csj)| fj * csj.interpolate(*ei))
                .sum();
        }
        cs.initialise();

        Some(cs)
    }
}


// ===============================================================================================
// Compton CDF tabulation, and related data.
// ===============================================================================================

#[derive(Default, Deserialize, Serialize)]
pub struct ComptonCDF {
    pub(crate) energies_in: LogGrid,
    pub(crate) x: LinearGrid,
    pub(crate) values: BilinearInterpolator,
    pub(super) computer: ComptonComputer,
}

// Public API.
impl ComptonCDF {
    pub fn energies_in(&self) -> &[Float] {
        self.energies_in.as_ref()
    }

    pub fn energy_out(&self, i: usize, j: usize) -> Float {
        let energy_in = self.energies_in[i];
        let m = self.x.len();
        let x = (j as Float) / ((m - 1) as Float);
        self.from_x(energy_in, x)
    }

    pub fn interpolate(&self, energy_in: Float, energy_out: Float) -> Float {
        let x = self.to_x(energy_in, energy_out);
        self.values.interpolate(&self.x, x, &self.energies_in, energy_in)
    }

    pub fn mode(&self) -> ComptonMode {
        self.computer.mode
    }

    pub fn model(&self) -> ComptonModel {
        self.computer.model
    }

    pub fn shape(&self) -> (usize, usize) {
        self.values.shape()
    }

    pub fn values(&self, i: usize) -> &[Float] {
        let m = self.x.len();
        &self.values.as_ref()[(i * m)..((i + 1) * m)]
    }
}

// Private API.
impl ComptonCDF {
    pub(crate) fn compute(
        &mut self,
        computer: &ComptonComputer,
        electrons: &ElectronicStructure
    ) -> Result<ComptonCrossSection> {

        let computer = match computer.model {
            ScatteringFunction | KleinNishina => {
                let mut tmp: ComptonComputer = computer.clone();
                tmp.precision *= self.energies_in.len() as Float;
                tmp
            },
            _ => bail!(
                "bad Compton model for CDF computation (expected '{}' or '{}', found '{}')",
                ScatteringFunction,
                KleinNishina,
                computer.model,
            ),
        };
        let energies: &Vec<Float> = self.energies_in.as_ref();
        let (n, m) = self.values.shape();
        let mut cross_section = ComptonCrossSection::new(energies[0], energies[n - 1], n);
        for (i, ei) in energies.iter().enumerate() {
            let (energy_min, energy_max) = computer.dcs_support(*ei);
            let energy_out = LogGrid::new(energy_min, energy_max, m);
            let mut sum = 0.0;
            for j in 0..(m - 1) {
                self.values[(i, j)] = sum;
                let emin = Some(energy_out[j]);
                let emax = Some(energy_out[j + 1]);
                sum += computer.cross_section(*ei, emin, emax, electrons)?;
            }
            for j in 0..(m - 1) {
                self.values[(i, j)] /= sum;
            }
            self.values[(i, m - 1)] = 1.0;
            cross_section.values[i] = sum;
        }
        cross_section.computer = computer.clone();
        self.computer = computer;
        cross_section.initialise();
        Ok(cross_section)
    }

    pub(crate) fn from_x(&self, energy_in: Float, x: Float) -> Float {
        let (energy_min, energy_max) = self.computer.dcs_support(energy_in);
        let x = x.clamp(0.0, 1.0);
        let lne = (energy_max / energy_min).ln();
        energy_min * (x * lne).exp()
    }

    pub(crate) fn new(energy_min: Float, energy_max: Float, n: usize, m: usize) -> Self {
        let energies_in = LogGrid::new(energy_min, energy_max, n);
        let values = BilinearInterpolator::new(n, m);
        let x = LinearGrid::new(0.0, 1.0, m);
        let computer = ComptonComputer::default();
        Self { energies_in, x, values, computer }
    }

    fn is_similar(&self, other: &Self) -> bool {
        self.shape() == other.shape() &&
            self.energies_in == other.energies_in &&
            self.computer.model == other.computer.model &&
            self.computer.mode == other.computer.mode
    }

    pub(crate) fn to_x(&self, energy_in: Float, energy_out: Float) -> Float {
        let (energy_min, energy_max) = self.computer.dcs_support(energy_in);
        let x = (energy_out / energy_min).ln() / (energy_max / energy_min).ln();
        x.clamp(0.0, 1.0)
    }

    // Merge a weighted sum of CDFs data into a single table.
    //
    // Note that input data must all have the same shape and support. Otherwise, `None` is
    // returned.
    #[allow(dead_code)] // XXX remove if not used.
    pub(crate) fn weighted_sum(data: &[(Float, &Self, &ComptonCrossSection)]) -> Option<Self> {

        if data.is_empty() {
            return None
        }
        let first = &data[0].1;
        if data
            .iter()
            .any(|(_, datum, _)| !first.is_similar(datum)) { return None }

        let (n, m) = first.shape();
        let mut values = BilinearInterpolator::new(n, m);
        for i in 0..n {
            let ei = first.energies_in[i];
            let weights = {
                let mut sum = 0.0;
                let mut weights: Vec<_> = data
                    .iter()
                    .map(|(fraction, _, cs)| {
                        let w = cs.interpolate(ei) * fraction;
                        sum += w;
                        w
                    })
                    .collect();
                for w in weights.iter_mut() { *w /= sum }
                weights
            };
            for j in 0..m {
                let index = (i, j);
                values[index] = data
                    .iter()
                    .enumerate()
                    .map(|(k, (_, datum, _))| weights[k] * datum.values[index])
                    .sum();
            }
        }

        Some(Self {
            energies_in: first.energies_in.clone(),
            values,
            x: first.x.clone(),
            computer: first.computer.clone(),
        })
    }
}


// ===============================================================================================
// Inverse CDF table.
// ===============================================================================================

#[derive(Default, Deserialize, Serialize)]
pub struct ComptonInverseCDF {
    pub(crate) energies: LogGrid,
    pub(crate) cdf: LinearGrid,
    pub(crate) values: BilinearInterpolator,
    pub(crate) weights: Option<BilinearInterpolator>,
    pub(super) computer: ComptonComputer,
}

// Public API.
impl ComptonInverseCDF {
    pub fn cdf(&self, j: usize) -> Float {
        let m = self.cdf.len();
        (j as Float) / ((m - 1) as Float)
    }

    pub fn energies(&self) -> &[Float] {
        self.energies.as_ref()
    }

    pub fn interpolate(&self, energy: Float, cdf: Float) -> (Float, Float) {
        let (energy_out, weight) = match &self.weights {
            None => (
                self.values.interpolate(&self.cdf, cdf, &self.energies, energy),
                1.0
            ),
            Some(weight) => {
                let (i, hi) = self.energies.transform(energy).clamp();
                let (j, hj) = self.cdf.transform(cdf).clamp();
                (
                    self.values.interpolate_raw(i, hi, j, hj),
                    weight.interpolate_raw(i, hi, j, hj),
                )
            },
        };
        let (energy_min, energy_max) = self.computer.dcs_support(energy);
        (energy_out.clamp(energy_min, energy_max), weight)
    }

    pub fn mode(&self) -> ComptonMode {
        self.computer.mode
    }

    pub fn model(&self) -> ComptonModel {
        self.computer.model
    }

    pub fn shape(&self) -> (usize, usize) {
        self.values.shape()
    }

    pub fn values(&self, i: usize) -> &[Float] {
        let m = self.cdf.len();
        &self.values.as_ref()[(i * m)..((i + 1) * m)]
    }

    pub fn weights(&self, i: usize) -> Option<&[Float]> {
        match self.weights.as_ref() {
            None => None,
            Some(weights) => {
                let m = self.cdf.len();
                Some(&weights.as_ref()[(i * m)..((i + 1) * m)])
            },
        }
    }
}

// Private API.
impl ComptonInverseCDF {
    pub(crate) fn compute(
        &mut self,
        computer: &ComptonComputer,
        electrons: &ElectronicStructure
    ) -> Result<Option<ComptonCrossSection>> {
        match computer.model {
            ScatteringFunction | KleinNishina => {},
            _ => bail!(
                "bad Compton model (expected '{}' or '{}', found '{}')",
                ScatteringFunction,
                KleinNishina,
                computer.model,
            ),
        };
        let (n, m) = self.values.shape();
        let mut cross_section = ComptonCrossSection::new(
            self.energies[0], self.energies[n - 1], n);
        for (i, ei) in self.energies.as_ref().iter().enumerate() {
            let (energy_min, energy_max) = computer.dcs_support(*ei);
            let energies_out: Vec<Float> = LogGrid::new(energy_min, energy_max, m).into();
            let cdf_values: UnstructuredGrid = match computer.mode {
                Adjoint | Direct => {
                    let cs0 = computer.cross_section(*ei, None, None, electrons)?;
                    cross_section.values[i] = cs0;
                    let mut cdf_values = Vec::<Float>::with_capacity(m);
                    for ej in energies_out.iter() {
                        cdf_values.push(
                            (computer.cross_section(*ei, None, Some(*ej), electrons)? / cs0)
                            .clamp(0.0, 1.0)
                        )
                    }
                    cdf_values.into()
                },
                Inverse => {
                    let mut cdf = Vec::<Float>::with_capacity(m);
                    for ej in energies_out.iter() {
                        let cs0 = computer.cross_section(*ej, None, Some(*ei), electrons)?;
                        let cs1 = computer.cross_section(*ej, Some(*ei), None, electrons)?;
                        cdf.push(cs1 / (cs0 + cs1));
                    }
                    cdf.into()
                },
                ComptonMode::None => unreachable!(),
            };
            let mut interpolator = CubicInterpolator::new(m);
            interpolator.initialise(&cdf_values, &energies_out, false);
            for (j, cj) in self.cdf.iter().enumerate() {
                let (k, d) = cdf_values
                    .transform(cj)
                    .clamp();
                let dx = cdf_values.width(k);
                self.values[(i, j)] =
                    interpolator.interpolate_raw(k, d, dx, &energies_out);
            }
        }

        self.weights = match computer.mode {
            Adjoint | Direct => None,
            Inverse => {
                let mut weights = BilinearInterpolator::new(n, m);
                for j in 0..m {
                    let mut interpolator = CubicInterpolator::new(n);
                    let energies_out: Vec<Float> = (0..n)
                        .map(|i| self.values[(i, j)])
                        .collect();
                    interpolator.initialise(&self.energies, &energies_out, false);
                    for i in 0..n {
                        weights[(i, j)] = interpolator.as_ref()[i].abs();
                    }
                }
                Some(weights)
            },
            ComptonMode::None => unreachable!(),
        };
        self.computer = computer.clone();

        match computer.mode {
            Adjoint | Direct => {
                cross_section.initialise();
                cross_section.computer = computer.clone();
                Ok(Some(cross_section))
            },
            Inverse => Ok(None),
            ComptonMode::None => unreachable!(),
        }
    }

    // Invert a forward CDF, yielding an inverse CDF.
    pub(crate) fn invert(cdf: &ComptonCDF) -> Self {
        let (n, m) = cdf.shape();
        let mut cdf_values = UnstructuredGrid::new(m);
        let mut energies_out = vec![0.0; m];
        let mut interpolator = CubicInterpolator::new(m);
        let mut result = Self::new(cdf.energies_in[0], cdf.energies_in[n - 1], n, m);
        result.computer = cdf.computer;
        for i in 0..n {
            let (energy_min, energy_max) = cdf.computer.dcs_support(cdf.energies_in[i]);
            let grid = LogGrid::new(energy_min, energy_max, m);
            for j in 0..m {
                cdf_values[j] = cdf.values[(i, j)];
                energies_out[j] = grid[j];
            }
            interpolator.initialise(&cdf_values, &energies_out, false);
            for j in 0..m {
                let (k, dk) = cdf_values
                    .transform(result.cdf.get(j))
                    .clamp();
                let eij = interpolator
                    .interpolate_raw(k, dk, cdf_values.width(k), &energies_out);
                result.values[(i, j)] = eij.clamp(energy_min, energy_max);
            }
        }
        result
    }

    pub(crate) fn new(energy_min: Float, energy_max: Float, n: usize, m: usize) -> Self {
        let energies = LogGrid::new(energy_min, energy_max, n);
        let cdf = LinearGrid::new(0.0, 1.0, m);
        let values = BilinearInterpolator::new(n, m);
        Self {
            energies, cdf, values, weights: None, computer: ComptonComputer::default()
        }
    }
}
