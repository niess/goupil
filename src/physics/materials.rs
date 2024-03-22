use anyhow::{anyhow, bail, Result};
use crate::numerics::float::Float;
use crate::numerics::table::FromFile;
use crate::physics::elements::AtomicElement;
use crate::physics::process::absorption::{AbsorptionMode, table::AbsorptionCrossSection};
use crate::physics::process::compton::{
    self,
    ComptonModel::{self, KleinNishina, Penelope, ScatteringFunction},
    compute::ComptonComputer,
    ComptonMethod::{InverseTransform, RejectionSampling},
    ComptonMode::{self, Adjoint, Direct, Inverse},
    table::{ComptonCrossSection, ComptonCDF, ComptonInverseCDF},
};
use crate::physics::process::rayleigh::{
    RayleighMode,
    table::{
        RayleighCrossSection,
        RayleighFormFactor,
        RayleighTable,
    },
};
use crate::transport::TransportSettings;
use regex::Regex;
use serde_derive::{Deserialize, Serialize};
use std::{
    collections::{hash_map::Entry::{Occupied, Vacant}, HashMap},
    fmt::{self, Display},
    fs,
    path::Path,
};

pub(crate) mod electronic;
pub(crate) mod table;


// ===============================================================================================
// Public API.
// ===============================================================================================

pub use self::{
    electronic::{ElectronicShell, ElectronicStructure},
    table::MaterialTable,
};


// ===============================================================================================
// Material description.
// ===============================================================================================

#[derive(Clone, Default, Deserialize, PartialEq, Serialize)]
pub struct MaterialDefinition {
    name: String,
    mass: Float,
    mass_composition: Vec<WeightedElement>,
    mole_composition: Vec<WeightedElement>,
}

pub type WeightedElement = (Float, &'static AtomicElement);
pub type WeightedMaterial<'a> = (Float, &'a MaterialDefinition);

// Public interface.
impl MaterialDefinition {
    /// Defines a new material from a chemical formula indicating its atomic content, e.g `H2O`.
    pub fn from_formula(formula: &str) -> Result<Self> {
        let re = Regex::new(r"([A-Z][a-z]?)([0-9]*)")?;
        let mut composition = Vec::<(Float, &AtomicElement)>::default();
        for captures in re.captures_iter(formula) {
            let symbol = captures.get(1).unwrap().as_str();
            let element = AtomicElement::from_symbol(symbol)?;
            let weight = captures.get(2).unwrap().as_str();
            let weight: Float = if weight.len() == 0 {
                1.0
            } else {
                weight.parse::<Float>()?
            };
            composition.push((weight, element));
        }
        let definition = MaterialDefinition::from_mole(formula, &composition);
        Ok(definition)
    }

    /// Defines a new material as a mixture of atomic elements. Note that a *mass* composition is
    /// expected.
    pub fn from_mass(name: &str, composition: &[WeightedElement]) -> Self {
        let n = composition.len();
        let mut sum = 0.0;
        let mut mole_composition = Vec::<WeightedElement>::with_capacity(n);
        for (wi, element) in composition.iter() {
            let xi = wi / element.A;
            mole_composition.push((xi, element));
            sum += xi;
        }
        for i in 0..n {
            mole_composition[i].0 /= sum;
        }
        let mass = 1.0 / sum;

        Self {
            name: name.to_string(),
            mass,
            mass_composition: composition.into(),
            mole_composition,
        }
    }

    /// Defines a new material as a compound or mixture of atomic elements. Note that a *mole*
    /// composition is expected.
    pub fn from_mole(name: &str, composition: &[WeightedElement]) -> Self {
        let n = composition.len();
        let mut mass = 0.0;
        let mut mass_composition = Vec::<WeightedElement>::with_capacity(n);
        for (xi, element) in composition.iter() {
            let wi = xi * element.A;
            mass_composition.push((wi, element));
            mass += wi;
        }
        for i in 0..n {
            mass_composition[i].0 /= mass;
        }

        Self {
            name: name.to_string(),
            mass,
            mass_composition,
            mole_composition: composition.into(),
        }
    }

    /// Defines a new material as a mixture of other materials. Note that a *mass* composition is
    /// expected.
    pub fn from_others(name: &str, composition: &[WeightedMaterial]) -> Self {
        let mut weights = HashMap::<&'static AtomicElement, Float>::default();
        let mut sum = 0.0;
        for (wi, material) in composition.iter() {
            let xi = wi / material.mass();
            for (xj, element) in material.mole_composition().iter() {
                let xij = xi * xj;
                weights
                    .entry(element)
                    .and_modify(|x| *x += xij)
                    .or_insert(xij);
            }
            sum += xi;
        }
        let n = composition.len();
        let mut composition = Vec::<WeightedElement>::with_capacity(n);
        for (element, weight) in weights.iter() {
            composition.push((weight / sum, element))
        }
        composition.sort_by(|a, b| a.1.Z.partial_cmp(&b.1.Z).unwrap());
        Self::from_mole(name, &composition)
    }

    /// Returns the molar mass of this material, in g / mole.
    pub fn mass(&self) -> Float {
        self.mass
    }

    /// Returns the material's mass composition.
    pub fn mass_composition(&self) -> &[WeightedElement] {
        &self.mass_composition
    }

    /// Returns the material's mole composition.
    pub fn mole_composition(&self) -> &[WeightedElement] {
        &self.mole_composition
    }

    /// Returns the name of this material.
    pub fn name(&self) -> &str {
        &self.name
    }
}

// Private interface.
impl MaterialDefinition {
    pub(crate) fn compute_electrons(&self) -> Result<ElectronicStructure> {
        let composition = {
            let mut composition = Vec::<(Float, &ElectronicStructure)>::default();
            for (weight, element) in self
                .mole_composition
                .iter() {

                let electrons = element.electrons()?;
                composition.push((*weight, electrons))
            }
            composition
        };
        let electrons = ElectronicStructure::from_others(&composition);
        Ok(electrons)
    }
}

impl Display for MaterialDefinition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let composition: Vec<_> = self.mole_composition
            .iter()
            .map(|(x, element)| format!("({}, {})", x, element.symbol))
            .collect();
        let composition = composition.join(", ");
        write!(f, "({}, {}, [{}])", self.name, self.mass, composition)
    }
}


// ===============================================================================================
// Materials registry.
// ===============================================================================================

#[derive(Default, Deserialize, Serialize)]
pub struct MaterialRegistry {
    pub(crate) absorption: HashMap<&'static AtomicElement, AbsorptionCrossSection>,
    elements: Vec<&'static AtomicElement>,
    elements_refs: HashMap<&'static AtomicElement, usize>,
    materials: HashMap<String, MaterialRecord>,
    pub(crate) scattering_cs: HashMap<&'static AtomicElement, RayleighCrossSection>,
    pub(crate) scattering_ff: HashMap<&'static AtomicElement, RayleighFormFactor>,
    pub(crate) energy_min: Float,
    pub(crate) energy_max: Float,
}

// Public interface.
impl MaterialRegistry {
    /// Adds a material to the registry.
    pub fn add(&mut self, definition: &MaterialDefinition) -> Result<()> {
        let key = definition.name();
        match self.materials.get(key) {
            None => {
                let mut update_elements = false;
                for (_, element) in definition.mole_composition().iter() {
                    *self.elements_refs
                        .entry(element)
                        .or_insert_with(|| {
                            update_elements = true;
                            1
                        }) += 1;
                }
                if update_elements {
                    self.update_elements();
                }
                let entry = MaterialRecord::new(definition);
                self.materials.insert(key.to_string(), entry);
            },
            Some(entry) => if entry.definition() != definition {
                bail!(
                    "bad definition for '{}' (expected {}, found {})",
                    key,
                    definition,
                    entry.definition(),
                )
            },
        };
        Ok(())
    }

    /// Computes physics tables for all registered materials.
    pub fn compute(
        &mut self,
        settings: &TransportSettings,
        length: Option<usize>,
        width: Option<usize>,
        precision: Option<Float>,
    ) -> Result<()> {
        // Defaults.
        const DEFAULT_ENERGY_MIN: Float = 1e-3;
        const DEFAULT_ENERGY_MAX: Float = 1e+1;
        const DEFAULT_LENGTH: usize = 321;
        const DEFAULT_WIDTH: usize = 201;

        let energy_min = settings.energy_min.unwrap_or(DEFAULT_ENERGY_MIN);
        let energy_max = settings.energy_max.unwrap_or(DEFAULT_ENERGY_MAX);
        if energy_min >= energy_max {
            bail!(
                "bad energy range (expected 'min' < 'max', found {} >= {})",
                energy_min,
                energy_max,
            )
        }
        let n = length.unwrap_or(DEFAULT_LENGTH);
        let m = width.unwrap_or(DEFAULT_WIDTH);

        // Check consistency of Compton options.
        compton::validate(settings.compton_model, settings.compton_mode, settings.compton_method)?;

        // Check validity of any constraint.
        match settings.constraint {
            None => (),
            Some(value) => match settings.compton_model {
                ScatteringFunction | KleinNishina  => (),
                _ => bail!(
                    "bad constraint for '{}' Compton model (expected 'None', found {})",
                    settings.compton_model,
                    value,
                ),
            },
        };

        // Collect computation items for Compton process.
        enum Entry {
            ComputeCrossSection,
            ComputeCDF,
            InvertCDF,
            ComputeComptonInverseCDF,
        }

        struct Item {
            entry: Entry,
            model: ComptonModel,
            mode: ComptonMode,
        }

        let model = settings.compton_model;
        let constraint = settings.constraint;
        let mode = settings.compton_mode;
        let method = settings.compton_method;
        let mut items = Vec::<Item>::default();

        match method {
            InverseTransform => match mode {
                Adjoint | Direct => match constraint {
                    None => {
                        items.push(Item {
                            entry: Entry::ComputeComptonInverseCDF, model, mode
                        });
                    },
                    Some(_) => {
                        items.push(Item { entry: Entry::ComputeCDF, model, mode });
                        items.push(Item { entry: Entry::InvertCDF, model, mode });
                    },
                },
                Inverse => match constraint {
                    None => {
                        items.push(Item {
                            entry: Entry::ComputeCrossSection, model, mode: Direct
                        });
                        items.push(Item {
                            entry: Entry::ComputeComptonInverseCDF, model, mode
                        });
                    },
                    Some(_) => {
                        items.push(Item { entry: Entry::ComputeCDF, model, mode: Direct });
                        items.push(Item {
                            entry: Entry::ComputeComptonInverseCDF, model, mode
                        });
                    },
                },
                ComptonMode::None => (),
            },
            RejectionSampling => match mode {
                ComptonMode::None => (),
                _ => match model {
                    KleinNishina => (),
                    _ => match constraint {
                        None => match model {
                            Penelope => items.push(Item {
                                entry: Entry::ComputeCrossSection,
                                model: ScatteringFunction,
                                mode: Direct
                            }),
                            _ => items.push(Item {
                                entry: Entry::ComputeCrossSection, model, mode
                            }),
                        },
                        Some(_) => items.push(Item { entry: Entry::ComputeCDF, model, mode }),
                    },
                },
            },
        };

        if let Adjoint = mode { // Always compute transport cross-section.
            items.push(Item {
                entry: Entry::ComputeCrossSection, model, mode: Direct
            });
        }

        // Initialise Compton computer.
        let mut computer = ComptonComputer::default();
        if let Some(precision) = precision {
            computer.precision = precision
        }

        for (_, material_record) in self.materials.iter_mut() {
            let table = &mut material_record.table;
            if settings.absorption != AbsorptionMode::None {
                let composition = {
                    let mut composition = Vec::<(Float, &AbsorptionCrossSection)>::default();
                    for (weight, element) in material_record
                        .definition
                        .mole_composition
                        .iter() {

                        let absorption = self
                            .absorption
                            .get(element)
                            .ok_or_else(|| anyhow!(
                                "missing absorption data for '{}'",
                                element.symbol,
                            ))?;
                        composition.push((*weight, absorption))
                    }
                    composition
                };
                table.absorption = AbsorptionCrossSection::from_others(&composition);
            }

            if let RayleighMode::FormFactor = settings.rayleigh {
                let compositions = {
                    let mut composition_cs = Vec::<(Float, &RayleighCrossSection)>::default();
                    let mut composition_ff = Vec::<(Float, &RayleighFormFactor)>::default();
                    for (weight, element) in material_record
                        .definition
                        .mole_composition
                        .iter() {

                        let cross_section = self
                            .scattering_cs
                            .get(element)
                            .ok_or_else(|| anyhow!(
                                "missing scattering cross-section for '{}'",
                                element.symbol,
                            ))?;
                        composition_cs.push((*weight, cross_section));

                        let form_factor = self
                            .scattering_ff
                            .get(element)
                            .ok_or_else(|| anyhow!(
                                "missing scattering form-factor for '{}'",
                                element.symbol,
                            ))?;
                        composition_ff.push((*weight, form_factor));
                    }
                    (composition_cs, composition_ff)
                };
                table.rayleigh = RayleighTable {
                    cross_section: RayleighCrossSection::from_others(&compositions.0),
                    form_factor: RayleighFormFactor::from_others(&compositions.1),
                };
            }

            let electrons = material_record.definition.compute_electrons()?;

            for item in items.iter() {
                computer.model = item.model;
                computer.mode = item.mode;
                let table = table.compton.get_mut(item.model).get_mut(item.mode);
                match item.entry {
                    Entry::ComputeCrossSection => if table.cross_section.is_none() {
                        let mut cs = ComptonCrossSection::new(
                            energy_min, energy_max, n);
                        cs.compute(&computer, &electrons)?;
                        table.cross_section = Some(cs);
                    },
                    Entry::ComputeCDF => if table.cdf.is_none() {
                        let mut cdf = ComptonCDF::new(
                            energy_min, energy_max, n, m);
                        let cs = cdf.compute(&computer, &electrons)?;
                        table.cdf = Some(cdf);
                        if table.cross_section.is_none() {
                            table.cross_section = Some(cs);
                        }
                    },
                    Entry::InvertCDF => if table.inverse_cdf.is_none() {
                        let cdf = table.cdf.as_ref().unwrap();
                        table.inverse_cdf = Some(ComptonInverseCDF::invert(&cdf));
                    },
                    Entry::ComputeComptonInverseCDF => if table.inverse_cdf.is_none() {
                        let mut icdf = ComptonInverseCDF::new(energy_min, energy_max, n, m);
                        let cs = icdf.compute(&computer, &electrons)?;
                        table.inverse_cdf = Some(icdf);
                        if table.cross_section.is_none() {
                            if let Some(cs) = cs {
                                table.cross_section = Some(cs);
                            }
                        }
                    },
                }
            }
            material_record.electrons = Some(electrons);
        }

        self.energy_min = energy_min;
        self.energy_max = energy_max;

        Ok(())
    }

    /// Returns recorded data for the given `key`.
    pub fn get(&self, key: &str) -> Result<&MaterialRecord> {
        match self.materials.get(key) {
            None => bail!(
                "bad key (no such entry '{}' in registry)",
                key,
            ),
            Some(entry) => Ok(entry),
        }
    }

    /// Returns the number of material records.
    pub fn len(&self) -> usize {
        self.materials.len()
    }

    /// Loads atomic elements data from a directory.
    pub fn load_elements<P>(&mut self, path: &P) -> Result<()>
    where
        P: AsRef<Path>
    {
        let path: &Path = path.as_ref();
        if path.is_dir() {
            let subdir = path.join("absorption");
            self.absorption = self.load_elements_into(subdir.as_path())?;
            let subdir = path.join("scattering").join("cross-section");
            self.scattering_cs = self.load_elements_into(subdir.as_path())?;
            let subdir = path.join("scattering").join("form-factor");
            self.scattering_ff = self.load_elements_into(subdir.as_path())?;
        } else {
            bail!(
                "bad path (expected a directory, found '{}')",
                path.display(),
            );
        }
        Ok(())
    }

    /// Returns a sorted copy of all material keys.
    pub fn keys(&self) -> Vec<String> {
        let mut keys: Vec<_> = self.materials
            .keys()
            .map(|key| key.clone())
            .collect();
        keys.sort();
        keys
    }

    /// Removes a material from the registry, returning its data record.
    pub fn remove(&mut self, key: &str) -> Result<MaterialRecord> {
        let record = match self.materials.remove(key) {
            None => bail!("bad key ('{}')", key),
            Some(entry) => {
                for (_, element) in entry.definition().mole_composition().iter() {
                    let mut update_elements = false;
                    match self.elements_refs.entry(element) {
                        Occupied(mut entry) => {
                            let count = entry.get_mut();
                            if *count == 1 {
                                update_elements = true;
                                entry.remove();
                            } else {
                                *count -= 1;
                            }
                        },
                        Vacant(_) => unreachable!(),
                    }
                    if update_elements {
                        self.update_elements();
                    }
                }
                entry
            },
        };
        Ok(record)
    }
}

// Private interface.
impl MaterialRegistry {
    fn load_elements_into<T, P>(
        &mut self,
        path: &P
    ) -> Result<HashMap<&'static AtomicElement, T>>
    where
        T: FromFile,
        P: AsRef<Path> + ?Sized,
    {
        let path: &Path = path.as_ref();
        let map = if path.is_dir() {
            let mut map = HashMap::<&'static AtomicElement, T>::default();
            for entry in fs::read_dir(path)? {
                let path = entry?.path();
                if let Some(extension) = path.extension() {
                    if extension == "txt" {
                        if let Some(symbol) = path
                            .file_stem()
                            .and_then(|stem| stem.to_str()) {

                            let element = AtomicElement::from_symbol(symbol)?;
                            map.insert(
                                element,
                                T::from_file(path)?,
                            );
                        }
                    }
                }
            }
            map
        } else {
            bail!(
                "bad path (expected a directory, found '{}')",
                path.display(),
            );
        };
        Ok(map)
    }

    fn update_elements(&mut self) {
        self.elements = self.elements_refs
            .keys()
            .map(|x| *x)
            .collect();
        self.elements.sort_by(|a, b| a.Z.cmp(&b.Z));
    }

    #[cfg(feature = "python")]
    pub(crate) fn atomic_data_loaded(&self) -> bool {
        self.absorption.len() > 0 ||
        self.scattering_cs.len() > 0 ||
        self.scattering_ff.len() > 0
    }
}


// ===============================================================================================
// Material record.
// ===============================================================================================

#[derive(Default, Deserialize, Serialize)]
pub struct MaterialRecord {
    pub(crate) definition: MaterialDefinition,
    pub(crate) electrons: Option<ElectronicStructure>,
    pub(crate) table: MaterialTable,
}

// Public interface.
impl MaterialRecord {
    pub fn absorption_cross_section(&self) -> Option<&AbsorptionCrossSection> {
        self.table.absorption.as_ref()
    }

    pub fn compton_cdf(
        &self,
        model: ComptonModel,
        mode: ComptonMode,
    ) -> Option<&ComptonCDF> {
        self.table.compton
            .get(model)
            .get(mode)
            .cdf
            .as_ref()
    }

    pub fn compton_cross_section(
        &self,
        model: ComptonModel,
        mode: ComptonMode,
    ) -> Option<&ComptonCrossSection> {
        self.table.compton
            .get(model)
            .get(mode)
            .cross_section
            .as_ref()
    }

    pub fn compton_inverse_cdf(
        &self,
        model: ComptonModel,
        mode: ComptonMode,
    ) -> Option<&ComptonInverseCDF> {
        self.table.compton
            .get(model)
            .get(mode)
            .inverse_cdf
            .as_ref()
    }

    /// Returns the sampling weight for Compton collisions. Depending on the specified `model` and
    /// `mode`, pre-computed material tables might be required (see e.g. the registry
    /// [`compute`](MaterialRegistry::compute) method).
    ///
    pub fn compton_weight(
        &self,
        model: ComptonModel,
        mode: ComptonMode,
        energy_in: Float,
        energy_out: Float,
    ) -> Result<Float> {
        let weight = match mode {
            Adjoint => {
                let subtable = self.table.compton.get(model);
                subtable.adjoint_weight(model, energy_in, energy_out)?
            },
            Inverse => match self.compton_inverse_cdf(model, mode) {
                None => bail!(
                    "no inverse CDF table for {}:{} Compton process",
                    model,
                    mode,
                ),
                Some(table) => match self.compton_cdf(model, Direct) {
                    None => bail!(
                        "no CDF table for {}:{} Compton process",
                        model,
                        Direct,
                    ),
                    Some(cdf) => {
                        let u = 1.0 - cdf.interpolate(energy_out, energy_in);
                        table.interpolate(energy_in, u).1
                    },
                },
            },
            _ => 1.0,
        };
        Ok(weight)
    }

    pub fn definition(&self) -> &MaterialDefinition {
        &self.definition
    }

    /// Returns the electronic structure of this material, computed from atomic elements using the
    /// [`compute`](MaterialRegistry::compute) method.
    pub fn electrons(&self) -> Option<&ElectronicStructure> {
        self.electrons.as_ref()
    }

    /// Returns the cross-section table for Rayleigh scattering.
    pub fn rayleigh_cross_section(&self) -> Option<&RayleighCrossSection> {
        self.table.rayleigh
            .cross_section
            .as_ref()
    }

    /// Returns the form-factor table for Rayleigh scattering.
    pub fn rayleigh_form_factor(&self) -> Option<&RayleighFormFactor> {
        self.table.rayleigh
            .form_factor
            .as_ref()
    }
}

// Private interface.
impl MaterialRecord {
    fn new(definition: &MaterialDefinition) -> Self {
        Self {
            definition: definition.clone(),
            electrons: None,
            table: MaterialTable::default(),
        }
    }
}
