use anyhow::{anyhow, bail, Result};
use crate::numerics::{
    float::{Float, Float3},
    rand::FloatRng,
};
use crate::physics::{
    consts::AVOGADRO_NUMBER,
    materials::{MaterialRegistry, MaterialRecord},
    process::{
        absorption::AbsorptionMode,
        compton::{
            self,
            sample::ComptonSampler,
            ComptonMode::{self, Adjoint, Direct, Inverse},
        },
        rayleigh::sample::RayleighSampler,
    },
};
use std::fmt;
use super::{
    boundary::TransportBoundary,
    density::DensityModel,
    geometry::{GeometryDefinition, GeometrySector, GeometryTracer},
    PhotonState,
    TransportMode::{Backward, Forward},
    TransportSettings, TransportVertex, VertexKind,
};


// ===============================================================================================
// Monte Carlo transport settings.
// ===============================================================================================

pub struct TransportAgent<'a, 'b, 'c, G, R, T>
where
    G: GeometryDefinition,
    R: FloatRng,
    T: GeometryTracer<'a, G>,
{
    geometry: &'a G,
    tracer: T,
    records: Vec<&'b MaterialRecord>,
    rng: &'c mut R,
}

// Public interface.
impl<'a, 'b, 'c, G, R, T> TransportAgent<'a, 'b, 'c, G, R, T>
where
    G: GeometryDefinition,
    R: FloatRng,
    T: GeometryTracer<'a, G>,
{
    pub fn new(
        geometry: &'a G,
        materials: &'b MaterialRegistry,
        rng: &'c mut R,
    ) -> Result<Self> {
        // Resolve materials.
        let mut records = Vec::<&MaterialRecord>::default();
        for material in geometry.materials().iter() {
            let material = materials.get(material.name())?;
            records.push(material);
        }

        // Get a ray tracer.
        let tracer = T::new(geometry)?;

        Ok(Self {
            geometry,
            tracer,
            records,
            rng,
        })
    }

    pub fn rng(&self) -> &R {
        self.rng
    }

    pub fn rng_mut(&mut self) -> &mut R {
        self.rng
    }

    pub fn transport(
        &mut self,
        settings: &TransportSettings,
        state: &mut PhotonState,
        vertices: Option<&mut Vec<TransportVertex>>,
    ) -> Result<TransportStatus> {
        // Regularise direction.
        let norm_data = {
            let norm2 = state.direction.norm2();
            if (norm2 - 1.0).abs() > 10.0 * Float::EPSILON {
                let norm = norm2.sqrt();
                if norm2.abs() == 0.0 {
                    bail!(
                        "bad direction norm (expected a strictly positive value, found {})",
                        norm,
                    );
                }
                let initial_direction = state.direction;
                state.direction /= norm;
                Some((norm, initial_direction))
            } else {
                None
            }
        };

        // Do the transport.
        let initial_direction = state.direction;
        let status = self.regularised_transport(settings, state, vertices)?;

        // Unregularise direction.
        if let Some(norm_data) = norm_data {
            if state.direction == initial_direction {
                state.direction = norm_data.1;
            } else {
                state.direction *= norm_data.0;
            }
        }

        Ok(status)
    }
}

// Private interface.
impl<'a, 'b, 'c, G, R, T> TransportAgent<'a, 'b, 'c, G, R, T>
where
    G: GeometryDefinition,
    R: FloatRng,
    T: GeometryTracer<'a, G>,
{
    fn regularised_transport(
        &mut self,
        settings: &TransportSettings,
        state: &mut PhotonState,
        mut vertices: Option<&mut Vec<TransportVertex>>,
    ) -> Result<TransportStatus> {
        // Check configuration.
        compton::validate(settings.compton_model, settings.compton_mode, settings.compton_method)?;
        match settings.mode {
            Backward => match settings.compton_mode {
                Direct => bail!(
                    "bad compton mode for backward transport (expected '{}', '{}' or '{}', \
                        found '{}')",
                    Adjoint, Inverse, ComptonMode::None,
                    settings.compton_mode,
                ),
                _ => (),
            }
            Forward => match settings.compton_mode {
                Adjoint | Inverse => bail!(
                    "bad compton mode for forward transport (expected '{}', or '{}', \
                        found '{}')",
                    Forward, ComptonMode::None,
                    settings.compton_mode,
                ),
                _ => (),
            },
        }

        // Get samplers.
        let compton_sampler = ComptonSampler::new(
            settings.compton_model, settings.compton_mode, settings.compton_method);
        let rayleigh_sampler = RayleighSampler::new(settings.rayleigh);

        // Check for a valid energy constraint, in reverse mode.
        let energy_constraint = match settings.compton_mode {
            Adjoint | Inverse => match settings.constraint {
                None => None,
                Some(value) => {
                    if value <= 0.0 {
                        bail!(
                            "bad constraint (expected a strictly positive value, found {:})",
                            value,
                        )
                    }
                    Some(value)
                },
            },
            _ => None,
        };

        // Initial energy should be less or equal constraint, if any.
        let mut energy_in = state.energy;
        if let Some(energy_constraint) = energy_constraint {
            if energy_in > energy_constraint {
                bail!(
                    "bad initial energy (expected a value below constraint ({}), found {})",
                    energy_constraint,
                    energy_in,
                )
            }
        }

        // Initialise stepping.
        let mut direction = self.get_direction(&state, &settings);
        self.tracer.reset(state.position, direction)?;
        let mut properties = match self.tracer.sector() {
            None => return Ok(TransportStatus::Exit),
            Some(index) => {
                if settings.boundary.inside(state.position, index) {
                    return Ok(TransportStatus::Boundary)
                }
                LocalProperties::new(self.geometry, index, &self.records)?
            },
        };

        let mut status = SteppingStatus::new(settings, state.length, energy_in);
        if let SteppingStatus::Stop(status) = status {
            return Ok(status);
        }

        if let Some(ref mut vertices) = vertices {
            let vertex = TransportVertex {
                sector: properties.index,
                kind: VertexKind::Start,
                state: state.clone(),
            };
            vertices.push(vertex);
        }

        // Transport loop.
        let mut interaction_length: Option<Float>;
        loop {
            // Randomise the distance to the next collision. First, let us compute the
            // corresponding column depth.
            let absorption_cross_section = settings.absorption.transport_cross_section(
                energy_in, properties.material)?;
            let compton_cross_section = compton_sampler.transport_cross_section(
                energy_in,
                properties.material
            )?;
            interaction_length = if compton_cross_section <= 0.0 {
                None
            } else {
                let lambda = properties.material.definition.mass() /
                    (compton_cross_section * AVOGADRO_NUMBER);
                Some(lambda)
            };
            let rayleigh_cross_section = rayleigh_sampler.transport_cross_section(
                energy_in,
                properties.material
            )?;
            let cross_section =
                absorption_cross_section +
                compton_cross_section +
                rayleigh_cross_section;
            let column_depth = if cross_section <= 0.0 {
                Float::INFINITY
            } else {
                let lambda = properties.material.definition.mass() /
                    (cross_section * AVOGADRO_NUMBER);
                -lambda * self.rng.uniform01().ln()
            };

            // Check for initial energy constraint.
            if status.is_first() {
                status = SteppingStatus::Next;
                if !settings.is_forward() {
                    if let Some(energy_constraint) = energy_constraint {
                        if energy_in == energy_constraint {
                            status = SteppingStatus::Last;
                        }
                    }
                }
            }

            // Convert the column depth to an euclidian distance.
            let physical_length = if column_depth == Float::INFINITY { Float::INFINITY } else {
                    properties.density.range(state.position, direction, column_depth)
            };
            if physical_length <= 0.0 {
                bail!(
                    "bad physical step length (expected a positive value, found {})",
                    physical_length
                )
            }

            // Compute the geometry step length.
            let geometry_length = self.tracer.trace(physical_length)?;
            if geometry_length <= 0.0 {
                bail!(
                    "bad geometry step length (expected a positive value, found {})",
                    geometry_length
                )
            }

            // Compute the boundary step length.
            let boundary_length = settings.boundary.distance(state.position, direction);

            // Compute the effective step length.
            let mut kind: Option<VertexKind> = None;
            let length = {
                let length = physical_length
                    .min(geometry_length)
                    .min(boundary_length);
                match settings.length_max {
                    None => length,
                    Some(length_max) => {
                        let step_max = length_max - state.length;
                        if length > step_max {
                            status = SteppingStatus::Stop(TransportStatus::LengthMax);
                            step_max
                        } else {
                            length
                        }
                    },
                }
            };

            // Apply any continuous processes.
            if settings.absorption == AbsorptionMode::Continuous {
                let depth = if physical_length == length {
                    column_depth
                } else {
                    let tmp = properties.density.column_depth(state.position, direction, length);
                    if tmp > column_depth {
                        bail!(
                            "bad column depth (expected {} or less, found {})",
                            column_depth,
                            tmp
                        )
                    }
                    tmp
                };
                state.weight *= AbsorptionMode::transport_weight(
                    energy_in,
                    depth,
                    properties.material
                )?;
            }

            // Simulate any collision.
            if physical_length == length {
                let sigma = if cross_section > compton_cross_section {
                    self.rng.uniform01() * cross_section
                } else {
                    cross_section
                };
                if absorption_cross_section > 0. && sigma <= absorption_cross_section {
                    state.energy = 0.0;
                    status = SteppingStatus::Stop(TransportStatus::Absorbed);
                } else if rayleigh_cross_section > 0. &&
                    sigma <= absorption_cross_section + rayleigh_cross_section {
                    state.direction = rayleigh_sampler.sample(
                        self.rng, state.energy, state.direction, properties.material)?;
                    kind = Some(VertexKind::Rayleigh);
                } else if status.is_last() {
                    status = SteppingStatus::Stop(TransportStatus::EnergyConstraint);
                } else {
                    let momentum_in = state.direction * state.energy;
                    let mut sample = compton_sampler.sample(
                        self.rng,
                        momentum_in,
                        properties.material,
                        energy_constraint
                    )?;

                    state.energy = sample.momentum_out.norm();
                    state.direction = sample.momentum_out / state.energy;
                    if let Some(energy_constraint) = energy_constraint {
                        if state.energy >= energy_constraint {
                            sample.constrained = true;
                        }
                    }
                    state.weight *= sample.weight;
                    if sample.constrained {
                        state.energy = energy_constraint.unwrap();
                        status = SteppingStatus::Last;
                    }

                    if !settings.is_forward() {
                        // Apply the transport weight and correct for the process sampling bias.
                        // I.e., the final energy was used instead of the initial one, when
                        // selecting the active process.
                        let energy_out = state.energy;
                        let csi_com = compton_sampler.transport_cross_section(
                            energy_out,
                            properties.material
                        )?;
                        state.weight *= csi_com / compton_cross_section;
                    }

                    direction = self.get_direction(&state, &settings);
                    energy_in = state.energy;
                    status.update(settings, energy_in);
                    kind = Some(VertexKind::Compton);
                }
            }

            // Notify back to the geometry tracer, and update the Monte Carlo state accordingly.
            // Note that following this call the propagation medium might change.
            self.tracer.update(length, direction)?;
            state.position = self.tracer.position();
            state.length += length;

            // Check for any termination condition. Else, update the physical properties of the
            // local geometry.
            if status.is_stop() {
                break
            } else if length == boundary_length {
                status = SteppingStatus::Stop(TransportStatus::Boundary);
                break
            } else if length == geometry_length {
                match self.tracer.sector() {
                    None => {
                        status = SteppingStatus::Stop(TransportStatus::Exit);
                        break
                    },
                    Some(index) => {
                        if let TransportBoundary::Sector(sector) = settings.boundary {
                            if index == sector {
                                status = SteppingStatus::Stop(TransportStatus::Boundary);
                                break
                            }
                        }

                        if index != properties.index {
                            // Update physical properties of the local geometry.
                            properties.update(index, &self.records)?;
                            kind = Some(VertexKind::Interface);
                        }
                    },
                }
            }

            // Record the Monte Carlo step.
            if let Some(ref mut vertices) = vertices {
                if let Some(kind) = kind {
                    let vertex = TransportVertex {
                        sector: properties.index,
                        kind,
                        state: state.clone(),
                    };
                    vertices.push(vertex);
                }
            }
        }

        if let Some(ref mut vertices) = vertices {
            let vertex = TransportVertex {
                sector: properties.index,
                kind: VertexKind::Stop,
                state: state.clone(),
            };
            vertices.push(vertex);
        }

        // Unpack stop condition and process Monte Carlo weight accordingly.
        match status {
            SteppingStatus::Stop(status) => {
                match status {
                    TransportStatus::EnergyMin |
                    TransportStatus::EnergyMax |
                    TransportStatus::Absorbed => {
                        state.weight = 0.0;
                    },
                    TransportStatus::EnergyConstraint => {
                        // Apply volume weight in reverse mode.
                        if let Some(lambda) = interaction_length {
                            let density_value = Self::get_density(
                                properties.density,
                                state.position
                            )?;
                            state.weight *= lambda / density_value;
                        }
                    },
                    _ => (),
                }
                Ok(status)
            },
            _ => unreachable!(),
        }
    }

    #[inline]
    fn get_density(density: &DensityModel, position: Float3) -> Result<Float> {
        let density = density.value(position);
        if density > 0.0 {
            Ok(density)
        } else {
            Err(anyhow!(
                "bad density (expected a strictly positive value, found {})",
                density,
            ))
        }
    }

    #[inline]
    fn get_direction(&self, state: &PhotonState, settings: &TransportSettings) -> Float3 {
        if settings.is_forward() {
            state.direction
        } else {
            -state.direction
        }
    }
}

// ===============================================================================================
// Temporary container for managing the local properties of a geometry sector.
// ===============================================================================================

struct LocalProperties<'a, 'b, G: GeometryDefinition> {
    pub material: &'b MaterialRecord,
    pub density: &'a DensityModel,

    geometry: &'a G,
    index: usize,
}

impl<'a, 'b, G: GeometryDefinition> LocalProperties<'a, 'b, G> {
    fn new(geometry: &'a G, index: usize, records: &[&'b MaterialRecord]) -> Result<Self> {
        let sector = Self::get_sector(geometry, index)?;
        let material = Self::get_material(records, sector.material)?;
        let properties = Self { material, density: &sector.density, geometry, index };
        Ok(properties)
    }

    #[inline]
    fn update(&mut self, index: usize, records: &[&'b MaterialRecord]) -> Result<()> {
        if index != self.index {
            let sector = Self::get_sector(self.geometry, index)?;
            self.material = Self::get_material(records, sector.material)?;
            self.density = &sector.density;
            self.index = index;
        }
        Ok(())
    }

    #[inline]
    fn get_material(records: &[&'b MaterialRecord], index: usize) -> Result<&'b MaterialRecord> {
        records
            .get(index)
            .ok_or_else(|| anyhow!(
                "bad material index (expected a value in [0, {}), found {})",
                records.len(),
                index,
            ))
            .copied()
    }

    #[inline]
    fn get_sector(geometry: &'a G, index: usize) -> Result<&'a GeometrySector> {
        geometry
            .sectors()
            .get(index)
            .ok_or_else(|| anyhow!(
                "bad sector index (expected a value in [0, {}), found {})",
                geometry.sectors().len(),
                index,
            ))
    }
}


// ===============================================================================================
// Flag indicating the transport termination condition.
// ===============================================================================================

#[derive(Clone, Copy)]
pub enum TransportStatus {
    Absorbed,
    Boundary,
    EnergyConstraint,
    EnergyMax,
    EnergyMin,
    Exit,
    LengthMax,
}

impl TransportStatus {
    const ABSORBED: &str = "Absorbed";
    const BOUNDARY: &str = "Boundary";
    const ENERGY_CONSTRAINT: &str = "Energy Constraint";
    const ENERGY_MAX: &str = "Energy Max";
    const ENERGY_MIN: &str = "Energy Min";
    const EXIT: &str = "Exit";
    const LENGTH_MAX: &str = "Length Max";
}

impl From<TransportStatus> for i32 {
    fn from(status: TransportStatus) -> Self {
        match status {
            TransportStatus::Absorbed => 0,
            TransportStatus::Boundary => 1,
            TransportStatus::EnergyConstraint => 2,
            TransportStatus::EnergyMax => 3,
            TransportStatus::EnergyMin => 4,
            TransportStatus::Exit => 5,
            TransportStatus::LengthMax => 6,
        }
    }
}

impl TryFrom<i32> for TransportStatus {
    type Error = anyhow::Error;

    fn try_from(value: i32) -> Result<Self> {
        match value {
            0 => Ok(Self::Absorbed),
            1 => Ok(Self::Boundary),
            2 => Ok(Self::EnergyConstraint),
            3 => Ok(Self::EnergyMax),
            4 => Ok(Self::EnergyMin),
            5 => Ok(Self::Exit),
            6 => Ok(Self::LengthMax),
            _ => Err(anyhow!(
                "bad transport status (expected a value in [0, 6], found {})",
                value,
            )),
        }
    }
}

impl fmt::Display for TransportStatus {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let status: &'static str = (*self).into();
        write!(f, "{}", status)
    }
}

impl TryFrom<&str> for TransportStatus {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            Self::ABSORBED => Ok(Self::Absorbed),
            Self::BOUNDARY => Ok(Self::Boundary),
            Self::ENERGY_CONSTRAINT => Ok(Self::EnergyConstraint),
            Self::ENERGY_MAX => Ok(Self::EnergyMax),
            Self::ENERGY_MIN => Ok(Self::EnergyMin),
            Self::EXIT => Ok(Self::Exit),
            Self::LENGTH_MAX => Ok(Self::LengthMax),
            _ => Err(anyhow!("bad transport status ({value})")),
        }
    }
}

impl From<TransportStatus> for &'static str {
    fn from(value: TransportStatus) -> Self {
        match value {
            TransportStatus::Absorbed => TransportStatus::ABSORBED,
            TransportStatus::Boundary => TransportStatus::BOUNDARY,
            TransportStatus::EnergyConstraint => TransportStatus::ENERGY_CONSTRAINT,
            TransportStatus::EnergyMax => TransportStatus::ENERGY_MAX,
            TransportStatus::EnergyMin => TransportStatus::ENERGY_MIN,
            TransportStatus::Exit => TransportStatus::EXIT,
            TransportStatus::LengthMax => TransportStatus::LENGTH_MAX,
        }
    }
}


// ===============================================================================================
// Internal flag for the transport stepping status.
// ===============================================================================================

enum SteppingStatus {
    First,
    Next,
    Last,
    Stop(TransportStatus),
}

impl SteppingStatus {
    #[inline]
    fn is_first(&self) -> bool {
        match self {
            Self::First => true,
            _ => false,
        }
    }

    #[inline]
    fn is_last(&self) -> bool {
        match self {
            Self::Last => true,
            _ => false,
        }
    }

    #[inline]
    fn is_stop(&self) -> bool {
        match self {
            Self::Stop(_) => true,
            _ => false,
        }
    }

    fn new(
        settings: &TransportSettings,
        length: Float,
        energy: Float,
    ) -> Self {
        if let Some(length_max) = settings.length_max {
            if length >= length_max {
                Self::Stop(TransportStatus::LengthMax)
            } else {
                Self::First
            }
        } else {
            let mut status = Self::First;
            status.update(settings, energy);
            status
        }
    }

    fn update(&mut self, settings: &TransportSettings, energy: Float) {
        match self {
            Self::Last => {},
            _ => {
                if let Some(energy_min) = settings.energy_min {
                    if energy <= energy_min {
                        *self = Self::Stop(TransportStatus::EnergyMin);
                        return
                    }
                }
                if let Some(energy_max) = settings.energy_max {
                    if energy >= energy_max {
                        *self = Self::Stop(TransportStatus::EnergyMax);
                        return
                    }
                }
            }
        }
    }
}
