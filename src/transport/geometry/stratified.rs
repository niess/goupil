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
pub enum TopographyData {
    Constant(Float),
    Map(Rc<TopographyMap>),
    Offset(Rc<TopographyMap>, Float),
}

impl TopographyData {
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
            Self::Offset(m, o) => ResolvedData::Offset(get_index(m), *o),
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
    Offset(usize, Float),
}

impl ResolvedData {
    fn get_z(&self, z: &[Option<Float>]) -> Option<Float> {
        match self {
            Self::Constant(value) => Some(*value),
            Self::Map(index) => z[*index],
            Self::Offset(index, value) => z[*index].map(|v| v + *value),
        }
    }
}


// ===============================================================================================
// Stratified geometry containing a stack of geological layers.
// ===============================================================================================

pub struct StratifiedGeometry {
    interfaces: Vec<TopographyInterface>,
    maps: Vec<Rc<TopographyMap>>,
    pub(crate) materials: Vec<MaterialDefinition>,
    pub(crate) sectors: Vec<GeometrySector>,
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
        let interfaces = vec![Vec::new(), Vec::new()];
        let maps = Vec::new();
        let materials = vec![material.clone()];
        let sectors = vec![Self::new_sector(0, density, description)];
        Self { interfaces, maps, materials, sectors }
    }

    /// Adds a new layer on top of the geometry, separated by the provided interface.
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
        let last = self.interfaces.len() - 1;
        self.interfaces.insert(last, interface);
        Ok(())
    }

    /// Sets the geometry bottom interface. By default, the geometry is not bounded from below.
    pub fn set_bottom(&mut self, interface: &[TopographyData]) {
        let interface = TopographyData::resolve_all(interface, &mut self.maps);
        self.interfaces[0] = interface;
    }

    /// Sets the geometry to interface. By default, the geometry is not bounded from above.
    pub fn set_top(&mut self, interface: &[TopographyData]) {
        let interface = TopographyData::resolve_all(interface, &mut self.maps);
        let last = self.interfaces.len() - 1;
        self.interfaces[last] = interface;
    }

    /// Returns the interfaces' elevation values `z` at (`x`, `y`) coordinates.
    pub fn z(&self, x: Float, y: Float) -> Vec<Option<Float>> {
        let z_map: Vec<_> = self.maps
            .iter()
            .map(|m| m.z(x, y))
            .collect();
        let get_z = |data: &[ResolvedData]| -> Option<Float> {
            for d in data.iter() {
                let value = d.get_z(&z_map);
                if value.is_some() {
                    return value
                }
            }
            None
        };
        let mut result = Vec::<Option<Float>>::with_capacity(self.interfaces.len());
        for interface in self.interfaces.iter() {
            let zi = get_z(interface);
            result.push(zi);
        }
        result
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
