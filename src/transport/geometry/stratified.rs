use anyhow::{bail, Result};
use crate::numerics::{Float, Float3};
use crate::numerics::grids::{Grid, GridCoordinate, LinearGrid};
use crate::numerics::interpolate::BilinearInterpolator;
use crate::physics::materials::MaterialDefinition;
use crate::transport::density::DensityModel;
use std::rc::Rc;
use super::{GeometryDefinition, GeometrySector, GeometryTracer};


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


// ===============================================================================================
// Stratified geometry tracer.
// ===============================================================================================

pub struct StratifiedTracer<'a> {
    definition: &'a StratifiedGeometry,

    position: Float3,
    direction: Float3,
    current_sector: Option<usize>,
    next_sector: Option<usize>,
    length: Float,
    cache: Vec<CachedValue<'a>>,
    delta_min: Option<Float>,
}

struct CachedValue<'a> {
    x: Float,
    y: Float,
    z: Option<Float>,
    map: &'a TopographyMap,
}

impl<'a> CachedValue<'a> {
    fn new(map: &'a TopographyMap) -> Self {
        let x = 0.0;
        let y = 0.0;
        let z = map.z(x, y);
        Self { x, y, z, map }
    }

    #[inline]
    fn update(&mut self, x: Float, y: Float) -> Option<Float> {
        if (x != self.x) || (y != self.y) {
            self.z = self.map.z(x, y);
        }
        self.z
    }
}

impl ResolvedData {
    fn compute_z(
        &self,
        cache: &mut [CachedValue],
        x: Float,
        y: Float
    ) -> Option<Float> {
        match self {
            Self::Constant(value) => Some(*value),
            Self::Map(index) => cache[*index].update(x, y),
            Self::Offset(index, value) => {
                cache[*index]
                    .update(x, y)
                    .map(|v| v + *value)
            },
        }
    }

    fn interface_z(
        interface: &[Self],
        cache: &mut [CachedValue],
        x: Float,
        y: Float
    ) -> Option<Float> {
        for data in interface.iter() {
            let value = data.compute_z(cache, x, y);
            if value.is_some() {
                return value;
            }
        }
        None
    }
}

impl<'a> StratifiedTracer<'a> {
    fn locate(&mut self, r: Float3) -> (Option<usize>, Float) {
        let interfaces = &self.definition.interfaces;
        let n = interfaces.len();
        let mut delta = Float::INFINITY;

        let bound = |x: Float| -> Float {
            match self.delta_min {
                None => x,
                Some(delta_min) => x.max(delta_min),
            }
        };

        // Check bottom layer.
        let zb = ResolvedData::interface_z(&interfaces[0], &mut self.cache, r.0, r.1);
        if let Some(zb) = zb {
            if r.2 < zb {
                return (None, bound(zb - r.2))
            } else {
                delta = r.2 - zb;
            }
        }

        for i in 1..(n-1) {
            let zi = ResolvedData::interface_z(&interfaces[i], &mut self.cache, r.0, r.1);
            if let Some(zi) = zi {
                let d = (r.2 - zi).abs();
                if d < delta { delta = d }
                if r.2 < zi {
                    return (Some(i - 1), bound(delta))
                }
            }
        }
        (None, bound(delta))
    }
}

impl<'a> GeometryTracer<'a, StratifiedGeometry> for StratifiedTracer<'a> {
    #[inline]
    fn definition(&self) -> &'a StratifiedGeometry {
        self.definition
    }

    fn new(definition: &'a StratifiedGeometry) -> Result<Self> {
        // Initialise local state.
        let position = Float3::default();
        let direction = Float3::default();
        let current_sector = None;
        let next_sector = None;
        let length = 0.0;
        let cache: Vec<_> = definition.maps.iter()
            .map(|map| CachedValue::new(map))
            .collect();

        let delta_min = {
            let mut delta: Option<Float> = None;
            for map in definition.maps.iter() {
                let d = map.x.width(0)
                    .min(map.y.width(0));
                if d > 0.0 {
                    match delta {
                        None => { delta = Some(d) },
                        Some(value) => if value > d {
                            delta = Some(d);
                        },
                    }
                }
            }
            delta
        };

        Ok(Self {
            definition,
            position,
            direction,
            current_sector,
            next_sector,
            length,
            cache,
            delta_min
        })
    }

    fn position(&self) -> Float3 {
        self.position
    }

    fn reset(&mut self, position: Float3, direction: Float3) -> Result<()> {
        self.position = position;
        self.direction = direction;
        let (sector, _) = self.locate(position); // XXX Cache delta?
        self.current_sector = sector;
        Ok(())
    }

    fn sector(&self) -> Option<usize> {
        self.current_sector
    }

    fn trace(&mut self, physical_length: Float) -> Result<Float> {
        // XXX HERE I AM. Implement this.
        Ok(self.length)
    }

    fn update(&mut self, length: Float, direction: Float3) -> Result<()> {
        self.position += self.direction * length;
        self.direction = direction;
        if length == self.length {
            self.current_sector = self.next_sector;
        }
        Ok(())
    }
}
