use anyhow::Result;
use crate::transport::{
    density::DensityModel,
    geometry::{ExternalGeometry, SimpleGeometry},
};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use super::materials::PyMaterialDefinition;


// ===============================================================================================
// Python wrapper for a simple geometry object.
// ===============================================================================================

#[pyclass(name = "SimpleGeometry", module = "goupil")]
pub struct PySimpleGeometry (pub SimpleGeometry);

#[pymethods]
impl PySimpleGeometry {
    #[new]
    fn new(
        material: PyRef<PyMaterialDefinition>,
        density: DensityModel,
    ) -> Result<Self> {
        let geometry = SimpleGeometry::new(&material.0, density);
        Ok(Self(geometry))
    }

    #[getter]
    fn get_density(&self, py: Python) -> PyObject {
        self.0.sectors[0].density.into_py(py)
    }

    #[setter]
    fn set_density(&mut self, value: DensityModel) -> Result<()> {
        self.0.sectors[0].density = value;
        Ok(())
    }

    #[getter]
    fn get_material(&self) -> PyMaterialDefinition {
        PyMaterialDefinition(self.0.materials[0].clone())
    }
}


// ===============================================================================================
// Python wrapper for an external geometry object.
// ===============================================================================================

#[pyclass(name = "ExternalGeometry", module = "goupil")]
pub struct PyExternalGeometry (pub ExternalGeometry);

#[pymethods]
impl PyExternalGeometry {
    #[new]
    pub fn new(path: &str) -> Result<Self> {
        let geometry = unsafe { ExternalGeometry::new(path)? };
        Ok(Self(geometry))
    }

    #[getter]
    fn get_materials<'p>(&self, py: Python<'p>) -> &'p PyTuple {
        let mut materials = Vec::<PyObject>::with_capacity(self.0.materials.len());
        for material in self.0.materials.iter() {
            let material = PyMaterialDefinition(material.clone());
            materials.push(material.into_py(py));
        }
        PyTuple::new(py, materials)
    }

    #[getter]
    fn get_sectors<'p>(&self, py: Python<'p>) -> &'p PyTuple {
        let sectors: Vec<_> = self.0
            .sectors
            .iter()
            .map(|sector| (
                sector.material,
                sector.density,
                sector.description
                    .as_ref()
                    .map(|description| description.to_string()),
            ))
            .collect();
        PyTuple::new(py, sectors)
    }

    fn update_material(
        &mut self,
        index: usize,
        material: PyRef<PyMaterialDefinition>
    ) -> Result<()> {
        self.0.update_material(index, &material.0)
    }

    fn update_sector(
        &mut self,
        index: usize,
        material: Option<usize>,
        density: Option<DensityModel>,
    ) -> Result<()> {
        self.0.update_sector(index, material, density.as_ref())
    }
}


// ===============================================================================================
// Unresolved geometry definition.
// ===============================================================================================

#[derive(Clone, FromPyObject)]
pub enum PyGeometryDefinition {
    External(Py<PyExternalGeometry>),
    Simple(Py<PySimpleGeometry>),
}

impl IntoPy<PyObject> for PyGeometryDefinition {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            Self::External(external) => external.into_py(py),
            Self::Simple(simple) => simple.into_py(py),
        }
    }
}
