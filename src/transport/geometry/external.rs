use anyhow::{anyhow, bail, Context, Result};
use anyhow::Error as AnyhowError;
use crate::numerics::{Float, Float3};
use crate::physics::elements::AtomicElement;
use crate::physics::materials::{MaterialDefinition, WeightedElement};
use crate::transport::density::DensityModel;
use libloading::Library;
use std::ffi::{c_int, CStr, OsStr};
use std::fmt::Display;
use super::{GeometryDefinition, GeometrySector, GeometryTracer};


// ===============================================================================================
// External geometry, dynamicaly loaded.
// ===============================================================================================

pub struct ExternalGeometry {
    #[allow(dead_code)]
    lib: Library, // for keeping the library alive.
    interface: CInterface,
    ptr: *mut CGeometry,
    pub(crate) materials: Vec<MaterialDefinition>,
    pub(crate) sectors: Vec<GeometrySector>,
}

#[repr(C)]
struct CInterface {
    new_definition: Option<extern "C" fn() -> *mut CGeometry>,
    new_tracer: Option<extern "C" fn(*const CGeometry) -> *mut CTracer>,
}

#[repr(C)]
struct CGeometry {
    destroy: Option<extern "C" fn(*mut CGeometry)>,
    get_material: Option<extern "C" fn(*const CGeometry, index: usize) -> *const CMaterial>,
    get_sector: Option<extern "C" fn(*const CGeometry, index: usize) -> CSector>,
    materials_len: Option<extern "C" fn(*const CGeometry) -> usize>,
    sectors_len: Option<extern "C" fn(*const CGeometry) -> usize>,
}

// Public interface.
impl ExternalGeometry {
    pub unsafe fn new<P>(path: P) -> Result<Self>
    where
        P: AsRef<OsStr> + Display
    {
        // Fetch geometry description from entry point.
        type Initialise = unsafe fn() -> CInterface;
        const INITIALISE: &[u8] = b"goupil_initialise\0";

        let library = Library::new(&path)?;
        let initialise = library.get::<Initialise>(INITIALISE)
            .with_context(|| format!(
                "{}: could not load geometry",
                path,
            ))?;
        let interface = unsafe { initialise() };
        let new_definition = as_fun_ok(&interface.new_definition)?;
        let geometry_ptr = new_definition();
        let geometry = as_ref_ok(geometry_ptr)?;

        // Build material definitions.
        let materials_len = as_fun_ok(&geometry.materials_len)?;
        let size = materials_len(geometry);
        let mut materials = Vec::<MaterialDefinition>::with_capacity(size);
        let get_material = as_fun_ok(&geometry.get_material)?;
        for i in 0..size {
            let material: &CMaterial = get_material(geometry, i)
                .as_ref()
                .ok_or_else(|| anyhow!(
                    "bad pointer for material {} (expected Some(address), found None)",
                    i,
                ))?;
            let material: MaterialDefinition = material
                .try_into()?;
            materials.push(material);
        }

        // Build geometry sectors.
        let sectors_len = as_fun_ok(&geometry.sectors_len)?;
        let size = sectors_len(geometry);
        let mut sectors = Vec::<GeometrySector>::with_capacity(size);
        let get_sector = as_fun_ok(&geometry.get_sector)?;
        for i in 0..size {
            let sector: GeometrySector = get_sector(geometry, i)
                .try_into()?;
            sectors.push(sector);
        }

        // Bundle the geometry definition.
        let geometry = Self {
            lib: library,
            interface,
            ptr: geometry_ptr,
            materials,
            sectors,
        };
        Ok(geometry)
    }

    pub fn update_material(
        &mut self,
        index: usize,
        material: &MaterialDefinition,
    ) -> Result<()> {
        self.check_material_index(index)?;
        self.materials[index] = material.clone();
        Ok(())
    }

    pub fn update_sector(
        &mut self,
        index: usize,
        material: Option<usize>,
        density: Option<&DensityModel>,
    ) -> Result<()> {
        self.check_sector_index(index)?;
        if let Some(material) = material {
            self.check_material_index(material)?;
            self.sectors[index].material = material;
        }
        if let Some(density) = density {
            self.sectors[index].density = density.clone();
        }
        Ok(())
    }
}

// Private interface.
impl ExternalGeometry {
    fn check_material_index(&self, index: usize) -> Result<()> {
        if index >= self.materials.len() {
            bail!(
                "bad material index (expected a value in [0, {}), found {}",
                self.materials.len(),
                index
            )
        }
        Ok(())
    }

    fn check_sector_index(&self, index: usize) -> Result<()> {
        if index >= self.sectors.len() {
            bail!(
                "bad sector index (expected a value in [0, {}), found {}",
                self.sectors.len(),
                index
            )
        }
        Ok(())
    }
}

impl Drop for ExternalGeometry {
    fn drop(&mut self) {
        if let Some(geometry) = unsafe { self.ptr.as_mut() } {
            if let Some(destroy) = geometry.destroy { destroy(self.ptr) }
        }
    }
}

impl GeometryDefinition for ExternalGeometry {
    #[inline]
    fn materials(&self) -> &[MaterialDefinition] {
        &self.materials
    }

    #[inline]
    fn sectors(&self) -> &[GeometrySector] {
        &self.sectors
    }
}

unsafe impl Send for ExternalGeometry {}

unsafe impl Sync for ExternalGeometry {}


// ===============================================================================================
// External geometry tracer.
// ===============================================================================================

pub struct ExternalTracer<'a> {
    definition: &'a ExternalGeometry,
    ptr: *mut CTracer,
    size: usize,
}

#[repr(C)]
struct CTracer {
    geometry: *const CGeometry,

    destroy: Option<extern "C" fn(*mut CTracer)>,
    position: Option<extern "C" fn(*const CTracer) -> CFloat3>,
    reset: Option<extern "C" fn(*mut CTracer, CFloat3, CFloat3)>,
    sector: Option<extern "C" fn(*const CTracer) -> usize>,
    trace: Option<extern "C" fn(*mut CTracer, Float) -> Float>,
    update: Option<extern "C" fn(*mut CTracer, Float, CFloat3)>,
}

impl<'a> Drop for ExternalTracer<'a> {
    fn drop(&mut self) {
        if let Some(tracer) = unsafe { self.ptr.as_mut() } {
            if let Some(destroy) = tracer.destroy { destroy(self.ptr) }
        }
    }
}

unsafe impl<'a> Send for ExternalTracer<'a> {}

impl<'a> GeometryTracer<'a, ExternalGeometry> for ExternalTracer<'a> {
    #[inline]
    fn definition(&self) -> &'a ExternalGeometry {
        self.definition
    }

    fn new(definition: &'a ExternalGeometry) -> Result<Self> {
        let tracer = {
            let new_tracer = as_fun_ok(&definition.interface.new_tracer)?;
            as_mut_ok(new_tracer(definition.ptr))?
        };

        if tracer.geometry != definition.ptr {
            Err(anyhow!(
                "bad geometry definition (expected {:?}, found {:?})",
                definition.ptr,
                tracer.geometry,
            ))
        } else {
            let tracer = ExternalTracer {
                definition,
                ptr: tracer,
                size: definition.sectors.len(),
            };
            Ok(tracer)
        }
    }

    fn position(&self) -> Float3 {
        let tracer = as_ref_ex(self.ptr);
        let position = as_fun_ex(&tracer.position);
        let position: Float3 = position(self.ptr).into();
        position
    }

    fn reset(&mut self, position: Float3, direction: Float3) -> Result<()> {
        let tracer = as_mut_ok(self.ptr)?;
        let reset = as_fun_ok(&tracer.reset)?;
        let position: CFloat3 = position.into();
        let direction: CFloat3 = direction.into();
        reset(self.ptr, position, direction);
        Ok(())
    }

    fn sector(&self) -> Option<usize> {
        let tracer = as_ref_ex(self.ptr);
        let sector = as_fun_ex(&tracer.sector);
        let sector = sector(self.ptr);
        if sector >= self.size { None } else { Some(sector) }
    }

    fn trace(&mut self, physical_length: Float) -> Result<Float> {
        let tracer = as_mut_ok(self.ptr)?;
        let trace = as_fun_ok(&tracer.trace)?;
        let length = trace(self.ptr, physical_length);
        Ok(length)
    }

    fn update(&mut self, length: Float, direction: Float3) -> Result<()> {
        let tracer = as_mut_ok(self.ptr)?;
        let update = as_fun_ok(&tracer.update)?;
        let direction: CFloat3 = direction.into();
        update(self.ptr, length, direction);
        Ok(())
    }
}

// ===============================================================================================
// External Float3 (C representation).
// ===============================================================================================

#[repr(C)]
struct CFloat3 {
    x: Float,
    y: Float,
    z: Float,
}

impl From<CFloat3> for Float3 {
    fn from(value: CFloat3) -> Self {
        Self::new(
            value.x,
            value.y,
            value.z,
        )
    }
}

impl From<Float3> for CFloat3 {
    fn from(value: Float3) -> Self {
        Self {
            x: value.0,
            y: value.1,
            z: value.2
        }
    }
}

// ===============================================================================================
// External geometry sector (C representation).
// ===============================================================================================

#[repr(C)]
struct CSector {
    material: usize,
    density: Float,
    description: *const i8,
}

impl TryInto<GeometrySector> for CSector {
    type Error = AnyhowError;

    fn try_into(self) -> Result<GeometrySector> {
        let material = self.material;
        let density = DensityModel::uniform(self.density)?;
        let description = if self.description.is_null() {
            None
        } else {
            let description = unsafe { CStr::from_ptr(self.description) }
                .to_str()?
                .to_string();
            Some(description)
        };
        let sector = GeometrySector {
            material,
            density,
            description,
        };
        Ok(sector)
    }
}


// ===============================================================================================
// External material definition (C representation).
// ===============================================================================================

#[repr(C)]
struct CMaterial {
    composition_len: Option<extern "C" fn(*const CMaterial) -> usize>,
    get_composition: Option<extern "C" fn(*const CMaterial, usize) -> CElement>,
    name: Option<extern "C" fn(*const CMaterial) -> *const i8>,
}

impl TryInto<MaterialDefinition> for &CMaterial {
    type Error = AnyhowError;

    fn try_into(self) -> Result<MaterialDefinition> {
        // Parse composition.
        let ptr = self as *const CMaterial;
        let composition_len = as_fun_ex(&self.composition_len);
        let composition_len = composition_len(ptr);
        let mut composition = Vec::<WeightedElement>::with_capacity(composition_len);
        let get_composition = as_fun_ex(&self.get_composition);
        for i in 0..composition_len {
            let weighted: WeightedElement = get_composition(ptr, i)
                .try_into()?;
            composition.push(weighted);
        }

        // Parse name.
        let name = as_fun_ex(&self.name);
        let name = name(ptr);
        let name = unsafe { CStr::from_ptr(name) }
                .to_str()?
                .to_string();

        Ok(MaterialDefinition::from_mole(&name, &composition))
    }
}


// ===============================================================================================
// External weighted element (C representation).
// ===============================================================================================

#[allow(non_snake_case)]
#[repr(C)]
struct CElement {
    Z: c_int,
    weight: Float,
}

impl TryInto<WeightedElement> for CElement {
    type Error = AnyhowError;

    fn try_into(self) -> Result<WeightedElement> {
        let z = i32::try_from(self.Z)?;
        let element = AtomicElement::from_Z(z)?;
        Ok((self.weight, element))
    }
}


// ===============================================================================================
// Shortcuts for converting pointers to references.
// ===============================================================================================

#[inline]
fn as_fun_ex<T>(ptr: &Option<T>) -> &T {
    ptr.as_ref()
        .expect(BAD_POINTER)
}

#[inline]
fn as_fun_ok<T>(ptr: &Option<T>) -> Result<&T> {
    match ptr.as_ref() {
        None => Err(anyhow!(BAD_POINTER)),
        Some(reference) => Ok(reference),
    }
}

#[inline]
fn as_mut<'a, T>(ptr: *mut T) -> Option<&'a mut T> {
    unsafe { ptr.as_mut() }
}

#[inline]
fn as_mut_ok<'a, T>(ptr: *mut T) -> Result<&'a mut T> {
    match as_mut(ptr) {
        None => Err(anyhow!(BAD_POINTER)),
        Some(reference) => Ok(reference),
    }
}

#[inline]
fn as_ref<'a, T>(ptr: *const T) -> Option<&'a T> {
    unsafe { ptr.as_ref() }
}

#[inline]
fn as_ref_ex<'a, T>(ptr: *const T) -> &'a T {
    as_ref(ptr)
        .expect(BAD_POINTER)
}

#[inline]
fn as_ref_ok<'a, T>(ptr: *const T) -> Result<&'a T> {
    match as_ref(ptr) {
        None => Err(anyhow!(BAD_POINTER)),
        Some(reference) => Ok(reference),
    }
}

static BAD_POINTER: &'static str = "bad pointer (null)";
