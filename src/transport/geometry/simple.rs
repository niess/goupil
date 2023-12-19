use anyhow::Result;
use crate::numerics::float::{Float, Float3};
use crate::physics::materials::MaterialDefinition;
use crate::transport::density::DensityModel;
use super::{GeometryDefinition, GeometrySector, GeometryTracer};


// ===============================================================================================
// Simple geometry containing a single medium of infinite extension.
// ===============================================================================================

pub struct SimpleGeometry {
    pub(crate) materials: [MaterialDefinition; 1],
    pub(crate) sectors: [GeometrySector; 1],
}

impl SimpleGeometry {
    /// Creates a simple geometry with the given `material` and bulk `density` model.
    pub fn new(material: &MaterialDefinition, density: DensityModel) -> Self {
        let materials = [material.clone()];
        let sectors = [GeometrySector { density, material: 0, description: None }];
        Self { materials, sectors }
    }

    /// Creates a simple geometry with the given `material` and a uniform bulk `density`, in
    /// g/cm<sup>3</sup>.
    pub fn uniform(material: &MaterialDefinition, density: Float) -> Result<Self> {
        let density = DensityModel::uniform(density)?;
        Ok(Self::new(material, density))
    }
}

impl GeometryDefinition for SimpleGeometry {
    #[inline]
    fn materials(&self)-> &[MaterialDefinition] {
        &self.materials
    }

    #[inline]
    fn sectors(&self)-> &[GeometrySector] {
        &self.sectors
    }
}


// ===============================================================================================
// Ray tracer for a simple geometry.
// ===============================================================================================

pub struct SimpleTracer<'a> {
    definition: &'a SimpleGeometry,
    position: Float3,
    direction: Float3,
}

impl<'a> GeometryTracer<'a, SimpleGeometry> for SimpleTracer<'a> {
    #[inline]
    fn definition(&self) -> &'a SimpleGeometry {
        self.definition
    }

    fn new(definition: &'a SimpleGeometry) -> Result<Self> {
        Ok(Self {
            definition,
            position: Float3::default(),
            direction: Float3::default(),
        })
    }

    #[inline]
    fn position(&self) -> Float3 {
        self.position
    }

    fn reset(&mut self, position: Float3, direction: Float3) -> Result<()> {
        self.position = position;
        self.direction = direction;
        Ok(())
    }

    #[inline]
    fn sector(&self) -> Option<usize> {
        Some(0)
    }

    #[allow(unused_variables)]
    fn trace(&mut self, physical_length: Float) -> Result<Float> {
        Ok(Float::INFINITY)
    }

    fn update(&mut self, length: Float, direction: Float3) -> Result<()> {
        self.position += length * self.direction;
        self.direction = direction;
        Ok(())
    }
}
