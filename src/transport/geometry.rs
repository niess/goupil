use anyhow::Result;
use crate::numerics::float::{Float, Float3};
use crate::physics::materials::MaterialDefinition;
use super::density::DensityModel;

pub(crate) mod external;
pub(crate) mod simple;
pub(crate) mod stratified;

pub use external::{ExternalGeometry, ExternalTracer};
pub use simple::{SimpleGeometry, SimpleTracer};
pub use stratified::{
    StratifiedGeometry,
    StratifiedTracer,
    TopographyMap,
    TopographySurface,
};


// ===============================================================================================
// Sectorised geometry interface.
// ===============================================================================================

pub trait GeometryDefinition
{
    /// Returns the materials composing this geometry.
    fn materials(&self) -> &[MaterialDefinition];

    /// Returns the sectors of this geometry.
    fn sectors(&self) -> &[GeometrySector];
}


// ===============================================================================================
// Ray tracing interface.
// ===============================================================================================

pub trait GeometryTracer<'a, D>
where
    Self: Sized,
    D: GeometryDefinition
{
    /// Returns the geometry definition associated to this ray tracer.
    fn definition(&self) -> &'a D;

    /// Constructs a new ray tracer for the given geometry `definition`.
    fn new(definition: &'a D) -> Result<Self>;

    /// Returns the current `position`.
    fn position(&self) -> Float3;

    /// Resets the geometry tracer for a new traversal. The `position` and pointing `direction` are
    /// initialised with the provided values.
    fn reset(&mut self, position: Float3, direction: Float3) -> Result<()>;

    /// Returns the current sector index, or `None` if the tracer left the geometry.
    fn sector(&self) -> Option<usize>;

    /// Computes the length to the next geometric boundary, from the current `position` and along
    /// the current `direction`.
    ///
    /// At input, the physical step length is provided. At output, this function must return a step
    /// length smaller than or equal to the geometric distance to the next medium.
    fn trace(&mut self, physical_length: Float) -> Result<Float>;

    /// Moves the geometry tracer by `length` along the current `direction`, and then updates the
    /// pointing `direction`.
    fn update(&mut self, length: Float, direction: Float3) -> Result<()>;
}


// ===============================================================================================
// Definition of a geometry sector.
// ===============================================================================================

pub struct GeometrySector {
    pub density: DensityModel,
    pub material: usize,
    pub description: Option<String>,
}
