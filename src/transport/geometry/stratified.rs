use anyhow::{bail, Result};
use crate::numerics::Float;
use crate::numerics::grids::{Grid, GridCoordinate, LinearGrid};
use crate::numerics::interpolate::BilinearInterpolator;
use crate::physics::materials::MaterialDefinition;
use crate::transport::density::DensityModel;
use std::rc::Rc;
use super::{GeometryDefinition, GeometrySector};


// ===============================================================================================
// Topography data using a digital elevation model (DEM).
// ===============================================================================================

pub struct TopographyMap {
    x: LinearGrid,
    y: LinearGrid,
    pub(crate) z: BilinearInterpolator,
}

impl TopographyMap {
    pub fn new(xmin: Float, xmax: Float, nx: usize, ymin: Float, ymax: Float, ny: usize) -> Self {
        let x = LinearGrid::new(xmin, xmax, nx);
        let y = LinearGrid::new(ymin, ymax, ny);
        let z = BilinearInterpolator::new(ny, nx);
        Self { x, y, z }
    }

    pub fn z(&self, x: Float, y: Float) -> Option<Float> {
        let (i, hi) = match self.y.transform(y) {
            GridCoordinate::Inside(i, hi) => (i, hi),
            _ => return None,
        };
        let (j, hj) = match self.x.transform(x) {
            GridCoordinate::Inside(j, hj) => (j, hj),
            _ => return None,
        };
        let zij = self.z.interpolate_raw(i, hi, j, hj);
        Some(zij)
    }
}


// ===============================================================================================
// Representations of topography data.
// ===============================================================================================

#[derive(Clone)]
pub enum TopographyData<'a> {
    Constant(Float),
    Map(&'a Rc<TopographyMap>),
    Offset(Float, &'a Rc<TopographyMap>),
}

impl<'a> TopographyData<'a> {
    fn resolve(&self, maps: &mut Vec<Rc<TopographyMap>>) -> ResolvedData {
        let mut get_index = |map: &Rc<TopographyMap>| -> usize  {
            for (i, mi) in maps.iter().enumerate() {
                if Rc::ptr_eq(map, mi) {
                    return i
                }
            }
            maps.push(Rc::clone(map));
            maps.len() - 1
        };

        match self {
            Self::Constant(z) => ResolvedData::Constant(*z),
            Self::Map(m) => ResolvedData::Map(get_index(m)),
            Self::Offset(o, m) => ResolvedData::Offset(*o, get_index(m)),
        }
    }

    fn resolve_all(
        interface: &[Self],
        maps: &mut Vec<Rc<TopographyMap>>
    ) -> TopographyInterface {
        let interface: Vec<_> = interface
            .iter()
            .map(|data| data.resolve(maps))
            .collect();
        interface
    }
}

type TopographyInterface = Vec<ResolvedData>;

enum ResolvedData {
    Constant(Float),
    Map(usize),
    Offset(Float, usize),
}


// ===============================================================================================
// Stratified geometry containing a stack of geological layers.
// ===============================================================================================

pub struct StratifiedGeometry {
    interfaces: Vec<TopographyInterface>,
    maps: Vec<Rc<TopographyMap>>,
    materials: Vec<MaterialDefinition>,
    sectors: Vec<GeometrySector>,
}

// Public interface.
impl StratifiedGeometry {
    /// Creates a new stratified geometry initialised with the given `material` and bulk `density`
    /// model.
    pub fn new(
        material: &MaterialDefinition,
        density: DensityModel,
        description: Option<&str>
    ) -> Self {
        let interfaces = Vec::new();
        let maps = Vec::new();
        let materials = vec![material.clone()];
        let sectors = vec![Self::new_sector(0, density, description)];
        Self { interfaces, maps, materials, sectors }
    }

    pub fn push_layer(
        &mut self,
        interface: &[TopographyData],
        material: &MaterialDefinition,
        density: DensityModel,
        description: Option<&str>,
    ) -> Result<()> {
        let material = match self.find_material(material)? {
            None => {
                self.materials.push(material.clone());
                self.materials.len() - 1
            },
            Some(material) => material,
        };
        let sector = Self::new_sector(material, density, description);
        self.sectors.push(sector);
        let interface = TopographyData::resolve_all(interface, &mut self.maps);
        self.interfaces.push(interface);
        Ok(())
    }
}

// Private interface.
impl StratifiedGeometry {
    fn find_material(&self, material: &MaterialDefinition) -> Result<Option<usize>> {
        for (i, mi) in self.materials.iter().enumerate() {
            if material.name() == mi.name() {
                if material != mi {
                    return Ok(Some(i))
                } else {
                    bail!(
                        "material '{}' already exists with a different definition",
                        material.name()
                    )
                }
            }
        }
        Ok(None)
    }

    fn new_sector(
        material: usize,
        density: DensityModel,
        description: Option<&str>
    ) -> GeometrySector {
        let description = description
            .unwrap_or("Layer 0")
            .to_string();
        let description = Some(description);
        GeometrySector { density, material, description }
    }
}

impl GeometryDefinition for StratifiedGeometry {
    #[inline]
    fn materials(&self)-> &[MaterialDefinition] {
        &self.materials
    }

    #[inline]
    fn sectors(&self)-> &[GeometrySector] {
        &self.sectors
    }
}
