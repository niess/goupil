use anyhow::{bail, Result};
use crate::numerics::{Float, Float3};
use crate::numerics::grids::{Grid, GridCoordinate, LinearGrid};
use crate::numerics::interpolate::BilinearInterpolator;
use crate::physics::materials::MaterialDefinition;
use crate::transport::density::DensityModel;
use std::fmt;
use std::ops::Deref;
use std::rc::Rc;
use super::{GeometryDefinition, GeometrySector, GeometryTracer};


// ===============================================================================================
// Topography data using a digital elevation model (DEM).
// ===============================================================================================

pub struct TopographyMap {
    x: LinearGrid,
    y: LinearGrid,
    pub(crate) z: MapData,
    pub(crate) zmin: Float,
    pub(crate) zmax: Float,
}

pub(crate) enum MapData {
    Interpolator(BilinearInterpolator),
    Scalar(Float),
}

// Public interface.
impl TopographyMap {
    pub fn new(xrange: &[Float; 2], yrange: &[Float; 2], shape: Option<&[usize; 2]>) -> Self {
        let zmin = 0.0;
        let zmax = 0.0;
        match shape {
            None => {
                let x = LinearGrid::new(xrange[0], xrange[1], 2);
                let y = LinearGrid::new(yrange[0], yrange[1], 2);
                let z = MapData::Scalar(0.0);
                Self { x, y, z, zmin, zmax }
            },
            Some(shape) => {
                let [ny, nx] = shape;
                let x = LinearGrid::new(xrange[0], xrange[1], *nx);
                let y = LinearGrid::new(yrange[0], yrange[1], *ny);
                let z = MapData::Interpolator(BilinearInterpolator::new(*ny, *nx));
                Self { x, y, z, zmin, zmax }
            }
        }
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
        let zij = match &self.z {
            MapData::Interpolator(z) => z.interpolate_raw(i, hi, j, hj),
            MapData::Scalar(z) => *z,
        };
        if zij.is_nan() {
            None
        } else {
            Some(zij)
        }
    }
}

// Private interface.
impl TopographyMap {
    pub(crate) fn get_box(&self) -> MapBox {
        let get_bounds = |grid: &LinearGrid| -> (Float, Float) {
            let x0 = grid.get(0);
            let x1 = grid.get(grid.len() - 1);
            if x0 < x1 {
                (x0, x1)
            } else {
                (x1, x0)
            }
        };
        let (xmin, xmax) = get_bounds(&self.x);
        let (ymin, ymax) = get_bounds(&self.y);
        MapBox { xmin, xmax, ymin, ymax }
    }
}

#[derive(Clone, Copy)]
pub(crate) struct MapBox {
    pub(crate) xmin: Float,
    pub(crate) xmax: Float,
    pub(crate) ymin: Float,
    pub(crate) ymax: Float,
}

impl MapBox {
    fn new() -> Self {
        let xmin = -Float::INFINITY;
        let xmax = Float::INFINITY;
        let ymin = -Float::INFINITY;
        let ymax = Float::INFINITY;
        Self { xmin, xmax, ymin, ymax }
    }

    fn clip(&mut self, other: &Self) {
        if other.xmin > self.xmin {
            self.xmin = other.xmin;
        }
        if other.xmax < self.xmax {
            self.xmax = other.xmax;
        }
        if other.ymin > self.ymin {
            self.ymin = other.ymin;
        }
        if other.ymax < self.ymax {
            self.ymax = other.ymax;
        }
    }

    fn inside(&self, x: Float, y: Float) -> bool {
        (x >= self.xmin) &&
        (x <= self.xmax) &&
        (y >= self.ymin) &&
        (y <= self.ymax)
    }
}

impl fmt::Display for MapBox {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{({}, {}), ({}, {})}}", self.xmin, self.xmax, self.ymin, self.ymax)
    }
}


// ===============================================================================================
// Representation of a topography surface.
// ===============================================================================================

#[derive(Clone)]
pub struct TopographySurface {
    maps: Vec<Rc<TopographyMap>>,
    pub offset: Float,
}

impl TopographySurface {
    pub fn new(maps: &[&Rc<TopographyMap>]) -> Result<Self> {
        let offset = 0.0;
        let result = if maps.len() == 0 {
            let maps = Vec::new();
            Self { maps, offset }
        } else {
            if maps.len() > 1 {
                let mut box0: MapBox = maps[0].get_box();
                for i in 1..maps.len() {
                    let box1 = maps[i].get_box();
                    if (box0.xmin < box1.xmin) ||
                       (box0.xmax > box1.xmax) ||
                       (box0.ymin < box1.ymin) ||
                       (box0.ymax > box1.ymax) {
                        bail!(
                            "bad maps size (expected {} to be included in {})",
                            box1,
                            box0,
                        )
                    }
                    box0 = box1;
                }
            }
            let maps: Vec<_> = maps
                .iter()
                .map(|m| Rc::clone(m))
                .collect();
            Self { maps, offset }
        };
        Ok(result)
    }

    pub fn z(&self, x: Float, y: Float) -> Option<Float> {
        for map in self.maps.iter() {
            let value = map.z(x, y);
            if let Some(value) = value {
                return Some(value + self.offset)
            }
        }
        None
    }
}

impl Deref for TopographySurface {
    type Target = [Rc<TopographyMap>];

    fn deref(&self) -> &Self::Target {
        &self.maps
    }
}

impl From<&Rc<TopographyMap>> for TopographySurface {
    fn from(value: &Rc<TopographyMap>) -> Self {
        let maps = vec![Rc::clone(value)];
        let offset = 0.0;
        Self { maps, offset }
    }
}


// ===============================================================================================
// Internal type storing a resolved topography surface.
// ===============================================================================================

#[derive(Default)]
struct ResolvedSurface {
    maps: Vec<usize>,
    offset: Float,
}

impl ResolvedSurface {
    fn new(interface: &TopographySurface, maps: &mut Vec<Rc<TopographyMap>>) -> Self {
        let mut get_index = |map: &Rc<TopographyMap>| -> usize  {
            for (i, mi) in maps.iter().enumerate() {
                if Rc::ptr_eq(map, mi) {
                    return i
                }
            }
            maps.push(Rc::clone(map));
            maps.len() - 1
        };
        let maps: Vec<_> = interface.maps
            .iter()
            .map(|m| get_index(m))
            .collect();
        let offset = interface.offset;
        Self { maps, offset }
    }
}


// ===============================================================================================
// Stratified geometry containing a stack of geological layers.
// ===============================================================================================

pub struct StratifiedGeometry {
    interfaces: Vec<ResolvedSurface>,
    maps: Vec<Rc<TopographyMap>>,
    pub(crate) materials: Vec<MaterialDefinition>,
    pub(crate) sectors: Vec<GeometrySector>,
    size: MapBox,
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
        let interfaces = vec![ResolvedSurface::default(), ResolvedSurface::default()];
        let maps = Vec::new();
        let materials = vec![material.clone()];
        let sectors = vec![Self::new_sector(0, density, 0, description)];
        let size = MapBox::new();
        Self { interfaces, maps, materials, sectors, size }
    }

    /// Adds a new layer on top of the geometry, separated by the provided interface.
    pub fn push_layer(
        &mut self,
        interface: &TopographySurface,
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
        let sector = Self::new_sector(material, density, self.sectors.len(), description);
        self.sectors.push(sector);
        if let Some(map) = interface.maps.last() {
            self.size.clip(&map.get_box());
        }
        let interface = ResolvedSurface::new(interface, &mut self.maps);
        let last = self.interfaces.len() - 1;
        self.interfaces.insert(last, interface);
        Ok(())
    }

    /// Sets the geometry bottom interface. By default, the geometry is not bounded from below.
    pub fn set_bottom(&mut self, interface: &TopographySurface) {
        if let Some(map) = interface.maps.last() {
            self.size.clip(&map.get_box());
        }
        let interface = ResolvedSurface::new(interface, &mut self.maps);
        self.interfaces[0] = interface;
    }

    /// Sets the geometry to interface. By default, the geometry is not bounded from above.
    pub fn set_top(&mut self, interface: &TopographySurface) {
        if let Some(map) = interface.maps.last() {
            self.size.clip(&map.get_box());
        }
        let interface = ResolvedSurface::new(interface, &mut self.maps);
        let last = self.interfaces.len() - 1;
        self.interfaces[last] = interface;
    }

    /// Returns the interfaces' elevation values `z` at (`x`, `y`) coordinates.
    pub fn z(&self, x: Float, y: Float) -> Vec<Option<Float>> {
        let z_map: Vec<_> = self.maps
            .iter()
            .map(|m| m.z(x, y))
            .collect();
        let get_z = |interface: &ResolvedSurface, i: usize, n: usize| -> Option<Float> {
            if interface.maps.len() == 0 {
                return if i == 0 {
                    if self.size.inside(x, y) {
                        Some(-Float::INFINITY)
                    } else {
                        None
                    }
                } else if i == n - 1 {
                    if self.size.inside(x, y) {
                        Some(Float::INFINITY)
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
            for index in interface.maps.iter() {
                let value = z_map[*index];
                if let Some(value) = value {
                    return Some(value + interface.offset)
                }
            }
            None
        };
        let mut result = Vec::<Option<Float>>::with_capacity(self.interfaces.len());
        let n = self.interfaces.len();
        for (i, interface) in self.interfaces.iter().enumerate() {
            let zi = get_z(interface, i, n);
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
                if material == mi {
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
        sectors: usize,
        description: Option<&str>
    ) -> GeometrySector {
        let description = match description {
            None => format!("Layer {}", sectors),
            Some(description) => description.to_string(),
        };
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
    inner_length: Float,
    outer_length: Float,
    next_delta: Float,
    cache: Vec<CachedValue<'a>>,
    delta_min: Float,
    zmin: Float,
    zmax: Float,
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

impl ResolvedSurface {
    fn z(
        &self,
        x: Float,
        y: Float,
        cache: &mut [CachedValue],
    ) -> Option<Float> {
        for index in self.maps.iter() {
            let value = cache[*index].update(x, y);
            if let Some(value) = value {
                return Some(value + self.offset);
            }
        }
        None
    }
}

impl<'a> StratifiedTracer<'a> {
    fn locate(&mut self, r: Float3) -> (Option<usize>, Float) {
        let size = &self.definition.size;
        let eps = 10.0 * Float::EPSILON;
        if (r.0 <= size.xmin + eps) || (r.0 >= size.xmax - eps) ||
           (r.1 <= size.ymin + eps) || (r.1 >= size.ymax - eps) {
               return (None, Float::INFINITY)
        }

        let interfaces = &self.definition.interfaces;
        let n = interfaces.len();
        let mut delta = Float::INFINITY;

        let bound = |x: Float| -> Float {
            x.max(self.delta_min)
        };

        // Check bottom layer.
        let zb = interfaces[0].z(r.0, r.1, &mut self.cache);
        if let Some(zb) = zb {
            if r.2 < zb {
                return (None, bound(zb - r.2))
            } else {
                delta = r.2 - zb;
            }
        }

        for i in 1..n {
            let zi = interfaces[i].z(r.0, r.1, &mut self.cache);
            match zi {
                None => {
                    if i == n - 1 {
                        return (Some(i - 1), bound(delta))
                    } else {
                        unreachable!("zi = {:?}, i = {}, r = {:}", zi, i, r);
                    }
                },
                Some(zi) => {
                    let d = (r.2 - zi).abs();
                    if d < delta { delta = d }
                    if r.2 < zi {
                        return (Some(i - 1), bound(delta))
                    }
                },
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
        let inner_length = 0.0;
        let outer_length = 0.0;
        let next_delta = 0.0;
        let cache: Vec<_> = definition.maps.iter()
            .map(|map| CachedValue::new(map))
            .collect();

        let delta_min = {
            let mut delta: Option<Float> = None;
            for map in definition.maps.iter() {
                if let MapData::Interpolator(_) = &map.z {
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
            }
            delta
        };
        let delta_min = delta_min.unwrap_or(1E+02);

        let (zmin, zmax) = {
            let mut zmin = Float::INFINITY;
            let mut zmax = -Float::INFINITY;
            for interface in definition.interfaces.iter() {
                for i in interface.maps.iter() {
                    let z = definition.maps[*i].zmin + interface.offset;
                    if z < zmin { zmin = z }
                    let z = definition.maps[*i].zmax + interface.offset;
                    if z > zmax { zmax = z }
                }
            }
            (zmin, zmax)
        };

        Ok(Self {
            definition,
            position,
            direction,
            current_sector,
            next_sector,
            inner_length,
            outer_length,
            next_delta,
            cache,
            delta_min,
            zmin,
            zmax,
        })
    }

    fn position(&self) -> Float3 {
        self.position
    }

    fn reset(&mut self, position: Float3, direction: Float3) -> Result<()> {
        self.position = position;
        self.direction = direction;
        let (sector, delta) = self.locate(position);
        self.inner_length = 0.0;
        self.outer_length = 0.0;
        self.next_delta = delta;
        self.current_sector = sector;
        self.next_sector = None;
        Ok(())
    }

    fn sector(&self) -> Option<usize> {
        self.current_sector
    }

    fn trace(&mut self, physical_length: Float) -> Result<Float> {
        // Compute length to xy sides.
        let side_length = {
            let size = self.definition.size;
            let compute_length = |x, ux, xmin, xmax| {
                if x < xmin {
                    if ux > 0.0 { (xmin - x) / ux } else { Float::INFINITY }
                } else if x > xmax {
                    if ux < 0.0 { (xmax - x) / ux } else { Float::INFINITY }
                } else {
                    if ux > 0.0 {
                        (xmax - x) / ux
                    } else if ux < 0.0 {
                        (xmin - x) / ux
                    } else {
                        Float::INFINITY
                    }
                }
            };
            let dx = compute_length(self.position.0, self.direction.0, size.xmin, size.xmax);
            let dy = compute_length(self.position.1, self.direction.1, size.ymin, size.ymax);
            dx.min(dy)
        };

        // Check unbounded cases.
        if self.current_sector.is_some() {
            let mut unbounded = false;
            if self.direction.2 < 0.0 {
                if self.position.2 < self.zmin {
                    unbounded = true;
                }
            } else if self.direction.2 > 0.0 {
                if self.position.2 > self.zmax {
                    unbounded = true;
                }
            }

            if unbounded {
                self.inner_length = side_length;
                self.outer_length = side_length;
                self.next_sector = None;
                return Ok(side_length)
            }
        }

        // Compute tentative step length.
        let factor = self.direction.2.abs().max(0.4);
        let mut inner_length = (self.next_delta * factor)
            .min(physical_length)
            .min(side_length);

        // Check corresponding location.
        let r1 = self.position + inner_length * self.direction;
        let (mut sector, mut delta) = self.locate(r1);
        if sector == self.current_sector {
            self.outer_length = inner_length;
        } else {
            let mut s0 = 0.0;
            let mut s1 = inner_length;
            while s1 - s0 > 1E-04 {
                let si = 0.5 * (s0 + s1);
                let ri = self.position + si * self.direction;
                let (sector_i, delta_i) = self.locate(ri);
                if sector_i == self.current_sector {
                    s0 = si;
                } else {
                    s1 = si;
                    sector = sector_i;
                    delta = delta_i;
                }
            }
            if s0 == 0.0 { s0 = s1 }
            inner_length = s0;
            self.outer_length = s1;
        }

        // Side check.
        if self.outer_length == side_length {
            sector = None;
        }

        // Update and return.
        self.inner_length = inner_length;
        self.next_sector = sector;
        self.next_delta = delta;
        Ok(inner_length)
    }

    fn update(&mut self, length: Float, direction: Float3) -> Result<()> {
        if length == self.inner_length {
            self.position += self.direction * self.outer_length;
            self.current_sector = self.next_sector;
        } else {
            self.position += self.direction * length;
        }
        self.direction = direction;
        Ok(())
    }
}
