use crate::numerics::float::{Float, Float3, Float3x3};
use crate::physics::materials::electronic::ElectronicShell;
// PyO3 interface.
use pyo3::conversion::{FromPyObject, IntoPy, ToPyObject};
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::{ffi, FromPyPointer};
use pyo3::marker::Python;
use pyo3::once_cell::GILOnceCell;
use pyo3::{Py, PyErr, PyNativeType, PyObject, PyResult};
use pyo3::type_object::PyTypeInfo;
use pyo3::types::{PyAny, PyCapsule, PyModule};
// Standard library.
use std::cell::UnsafeCell;
use std::ffi::{c_char, c_int, c_uchar, c_void};
use std::marker::PhantomData;
use std::ops::Deref;
// Local Python interface.
use super::transport::CState;


//================================================================================================
// Numpy floating point type (defined at compilation).
//================================================================================================
#[cfg(not(feature = "f32"))]
pub(super) const FLOAT_FORMAT: &str = "f8";
#[cfg(feature = "f32")]
pub(super) const FLOAT_FORMAT: &str = "f4";


// ===============================================================================================
//
// Numpy array interface.
//
// ===============================================================================================

struct ArrayInterface {
    // Keep the capsule alive.
    #[allow(dead_code)]
    capsule: PyObject,
    // Type objects.
    dtype_bool: PyObject,
    dtype_float: PyObject,
    dtype_int32: PyObject,
    dtype_shell: PyObject,
    dtype_state: PyObject,
    dtype_usize: PyObject,
    type_ndarray: PyObject,
    // Functions.
    empty: *const PyArray_Empty,
    equiv_types: *const PyArray_EquivTypes,
    new_from_descriptor: *const PyArray_NewFromDescriptor,
    scalar: *const PyArray_Scalar,
    scalar_as_ctype: *const PyArray_ScalarAsCtype,
    set_base_object: *const PyArray_SetBaseObject,
    zeros: *const PyArray_Zeros,
}

#[allow(non_camel_case_types)]
pub type npy_intp = ffi::Py_intptr_t;

#[allow(non_camel_case_types)]
type PyArray_Empty = extern "C" fn(
    nd: c_int,
    dims: *const npy_intp,
    dtype: *mut ffi::PyObject,
    fortran: c_int,
) -> *mut ffi::PyObject;

#[allow(non_camel_case_types)]
type PyArray_EquivTypes = extern "C" fn(
    type1: *mut ffi::PyObject,
    type2: *mut ffi::PyObject,
) -> c_uchar;

#[allow(non_camel_case_types)]
type PyArray_NewFromDescriptor = extern "C" fn(
    subtype: *mut ffi::PyObject,
    descr: *mut ffi::PyObject,
    nd: c_int,
    dims: *const npy_intp,
    strides: *const npy_intp,
    data: *mut c_void,
    flags: c_int,
    obj: *mut ffi::PyObject,
) -> *mut ffi::PyObject;

#[allow(non_camel_case_types)]
type PyArray_Scalar = extern "C" fn(
    data: *mut c_void,
    dtype: *mut ffi::PyObject,
    base: *mut ffi::PyObject,
) -> *mut ffi::PyObject;

#[allow(non_camel_case_types)]
type PyArray_ScalarAsCtype = extern "C" fn(
    scalar: *mut ffi::PyObject,
    ctypeptr: *mut c_void,
);

#[allow(non_camel_case_types)]
type PyArray_SetBaseObject = extern "C" fn(
    arr: *mut ffi::PyObject,
    obj: *mut ffi::PyObject
) -> c_int;

#[allow(non_camel_case_types)]
type PyArray_Zeros = extern "C" fn(
    nd: c_int,
    dims: *const npy_intp,
    dtype: *mut ffi::PyObject,
    fortran: c_int,
) -> *mut ffi::PyObject;

unsafe impl Send for ArrayInterface {}
unsafe impl Sync for ArrayInterface {}

static ARRAY_INTERFACE: GILOnceCell<ArrayInterface> = GILOnceCell::new();

fn api(py: Python) -> &ArrayInterface {
    ARRAY_INTERFACE
        .get(py)
        .expect("Numpy Array API not initialised")
}

pub fn initialise(py: Python) -> PyResult<()> {
    if let Some(_) = ARRAY_INTERFACE.get(py) {
        return Err(PyValueError::new_err("Numpy Array API already initialised"))
    }

    // Import interfaces.
    let numpy = PyModule::import(py, "numpy")?;
    let capsule = PyModule::import(py, "numpy.core.multiarray")?
        .getattr("_ARRAY_API")?;

    // Cache used dtypes, generated from numpy Python interface.
    let dtype = numpy.getattr("dtype")?;

    let dtype_bool: PyObject = dtype
        .call1(("bool",))?
        .into_py(py);

    let dtype_float: PyObject = dtype
        .call1((FLOAT_FORMAT,))?
        .into_py(py);

    let dtype_int32: PyObject = dtype
        .call1(("i4",))?
        .into_py(py);

    let dtype_shell: PyObject = dtype
        .call1(([
            ("occupancy", FLOAT_FORMAT),
            ("energy", FLOAT_FORMAT),
            ("momentum", FLOAT_FORMAT),
        ],))?
        .into_py(py);

    let dtype_state: PyObject = {
        let arg: [PyObject; 5] = [
            ("energy", FLOAT_FORMAT).into_py(py),
            ("position", FLOAT_FORMAT, 3).into_py(py),
            ("direction", FLOAT_FORMAT, 3).into_py(py),
            ("length", FLOAT_FORMAT).into_py(py),
            ("weight", FLOAT_FORMAT).into_py(py),
        ];
        dtype
            .call1((arg,))?
            .into_py(py)
    };

    let dtype_usize: PyObject = dtype
        .call1((format!("u{}", std::mem::size_of::<usize>()),))?
        .into_py(py);

    // Parse C interface.
    // See e.g. numpy/_core/code_generators/numpy_api.py for API mapping.
    let ptr = capsule
        .downcast::<PyCapsule>()?
        .pointer() as *const *const c_void;

    let object = |offset: isize| -> PyObject {
        unsafe {
            PyAny::from_borrowed_ptr(py, *ptr.offset(offset) as *mut ffi::PyObject)
                .into_py(py)
        }
    };

    let function = |offset: isize| unsafe {
        ptr.offset(offset)
    };

    let api = ArrayInterface {
        capsule: capsule.into(),
        // Type objects.
        dtype_bool,
        dtype_float,
        dtype_int32,
        dtype_shell,
        dtype_state,
        dtype_usize,
        type_ndarray: object(2),
        // Functions.
        empty:               function(184) as *const PyArray_Empty,
        equiv_types:         function(182) as *const PyArray_EquivTypes,
        new_from_descriptor: function( 94) as *const PyArray_NewFromDescriptor,
        scalar:              function(60)  as *const PyArray_Scalar,
        scalar_as_ctype:     function(62)  as *const PyArray_ScalarAsCtype,
        set_base_object:     function(282) as *const PyArray_SetBaseObject,
        zeros:               function(183) as *const PyArray_Zeros,
    };

    // Initialise static data and return.
    match ARRAY_INTERFACE.set(py, api) {
        Err(_) => unreachable!(),
        Ok(_) => (),
    }
    Ok(())
}


// ===============================================================================================
//
// Generic (untyped) array.
//
// ===============================================================================================

#[repr(transparent)]
pub struct PyUntypedArray(UnsafeCell<PyArrayObject>);

#[repr(C)]
pub struct PyArrayObject {
    pub object: ffi::PyObject,
    pub data: *mut c_char,
    pub nd: c_int,
    pub dimensions: *mut npy_intp,
    pub strides: *mut npy_intp,
    pub base: *mut ffi::PyObject,
    pub descr: *mut ffi::PyObject,
    pub flags: c_int,
}

// Public interface.
impl PyUntypedArray {
    #[inline]
    pub fn dtype(&self) -> &PyAny {
        unsafe { PyAny::from_borrowed_ptr(self.py(), self.as_ptr()) }
    }

    #[inline]
    pub fn ndim(&self) -> usize {
        let obj: &PyArrayObject = self.as_ref();
        obj.nd as usize
    }

    pub fn readonly(&self) {
        let obj = unsafe { &mut *(self.as_ptr() as *mut PyArrayObject) };
        obj.flags &= !PyArrayFlags::WRITEABLE;
    }

    #[inline]
    pub fn shape(&self) -> Vec<usize> {
        self.shape_slice()
            .iter()
            .map(|v| *v as usize)
            .collect()
    }

    #[inline]
    pub fn size(&self) -> usize {
        self.shape_slice()
            .iter()
            .product::<npy_intp>() as usize
    }
}

// Private interface.
impl PyUntypedArray {
    fn data(&self, index: usize) -> PyResult<*mut c_char> {
        let size = self.size();
        if index >= size {
            Err(PyIndexError::new_err(format!(
                "ndarray index out of range (expected an index in [0, {}), found {})",
                size,
                index
            )))
        } else {
            let offset = self.offset_of(index);
            let obj: &PyArrayObject = self.as_ref();
            let data = unsafe { obj.data.offset(offset as isize) };
            Ok(data)
        }
    }

    fn offset_of(&self, index: usize) -> isize {
        let shape = self.shape_slice();
        let strides = self.strides_slice();
        let n = shape.len();
        if n == 0 {
            0
        } else {
            let mut remainder = index;
            let mut offset = 0_isize;
            for i in (0..n).rev() {
                let m = shape[i] as usize;
                let j = remainder % m;
                remainder = (remainder - j) / m;
                offset += (j as isize) * strides[i];
            }
            offset
        }
    }

    #[inline]
    fn shape_slice(&self) -> &[npy_intp] {
        let obj: &PyArrayObject = self.as_ref();
        unsafe { std::slice::from_raw_parts(obj.dimensions, obj.nd as usize) }
    }

    #[inline]
    fn strides_slice(&self) -> &[npy_intp] {
        let obj: &PyArrayObject = self.as_ref();
        unsafe { std::slice::from_raw_parts(obj.strides, obj.nd as usize) }
    }

    #[inline]
    fn transmute(&self) -> &PyAny {
        unsafe { &*(self as * const Self as *const PyAny) }
    }
}

// Trait implementations.
impl AsRef<PyArrayObject> for PyUntypedArray {
    #[inline]
    fn as_ref(&self) -> &PyArrayObject {
        unsafe { &*self.0.get() }
    }
}

unsafe impl PyNativeType for PyUntypedArray {}

unsafe impl PyTypeInfo for PyUntypedArray {
    type AsRefTarget = PyUntypedArray;

    const NAME: &'static str = "numpy.ndarray";
    const MODULE: Option<&'static str> = Some("numpy");

    fn type_object_raw(py: Python<'_>) -> *mut ffi::PyTypeObject {
        api(py)
            .type_ndarray
            .as_ptr() as *mut ffi::PyTypeObject
    }
}

impl AsRef<PyAny> for PyUntypedArray {
    #[inline]
    fn as_ref(&self) -> &PyAny { self.transmute() }
}

impl Deref for PyUntypedArray {
    type Target = PyAny;

    #[inline]
    fn deref(&self) -> &Self::Target { self.transmute() }
}

impl<'py> FromPyObject<'py> for &'py PyUntypedArray {
    fn extract(obj: &'py PyAny) -> PyResult<Self> {
        obj.downcast().map_err(std::convert::Into::into)
    }
}

impl IntoPy<Py<PyUntypedArray>> for &'_ PyUntypedArray {
    #[inline]
    fn into_py(self, py: Python) -> Py<PyUntypedArray> {
        unsafe { Py::from_borrowed_ptr(py, self.as_ptr()) }
    }
}

impl ToPyObject for PyUntypedArray {
    #[inline]
    fn to_object(&self, py: Python) -> PyObject {
        unsafe { PyObject::from_borrowed_ptr(py, self.as_ptr()) }
    }
}


// ===============================================================================================
//
// Typed array.
//
// ===============================================================================================

#[repr(transparent)]
pub struct PyArray<T>(PyUntypedArray, PhantomData<T>);

// Public interface.
impl<T> PyArray<T>
where
    T: Copy + Dtype,
{
    pub fn empty<'py>(py: Python<'py>, shape: &[usize]) -> PyResult<&'py Self> {
        let api = api(py);
        let empty = unsafe { *api.empty };
        let dtype = T::dtype(py)?;
        let (ndim, shape) = Self::try_shape(shape)?;
        let array = empty(
            ndim,
            shape.as_ptr() as *const npy_intp,
            dtype.as_ptr(),
            0,
        );
        if PyErr::occurred(py) {
            match PyErr::take(py) {
                None => unreachable!(),
                Some(err) => return Err(err),
            }
        }
        let array = unsafe { &*(array as *const Self) };
        Ok(array)
    }

    pub fn from_data<'py>(
        py: Python<'py>,
        data: &[T],
        base: &PyAny,
        flags: PyArrayFlags,
        shape: Option<&[usize]>,
    ) -> PyResult<&'py Self> {
        let api = api(py);
        let new_from_descriptor = unsafe { *api.new_from_descriptor };
        let dtype = T::dtype(py)?;
        let (ndim, shape) = match shape {
            None => {
                let size = Self::try_size(data.len())?;
                (1, vec![size as npy_intp])
            },
            Some(shape) => {
                let size = shape.iter().product::<usize>();
                if size != data.len() {
                    return Err(PyValueError::new_err(format!(
                        "bad ndarray size (expected {}, found {})",
                        data.len(),
                        size,
                    )))
                }
                Self::try_shape(shape)?
            },
        };
        let array = new_from_descriptor(
            api.type_ndarray.as_ptr(),
            dtype.as_ptr(),
            ndim,
            shape.as_ptr() as *const npy_intp,
            std::ptr::null_mut(),
            data.as_ptr() as *mut c_void,
            flags.into(),
            std::ptr::null_mut(),
        );
        if PyErr::occurred(py) {
            match PyErr::take(py) {
                None => unreachable!(),
                Some(err) => return Err(err),
            }
        }
        let set_base_object = unsafe { *api.set_base_object };
        let ptr = base.as_ptr();
        set_base_object(array, ptr);
        unsafe { pyo3::ffi::Py_INCREF(ptr); }
        let array = unsafe { &*(array as *const Self) };
        Ok(array)
    }

    pub fn from_iter<'py, I>(py: Python<'py>, shape: &[usize], iter: I) -> PyResult<&'py Self>
    where
        I: Iterator<Item=T>,
    {
        let array = Self::empty(py, shape)?;
        let data = unsafe { array.slice_mut()? };
        for (xi, val) in std::iter::zip(data.iter_mut(), iter) {
            *xi = val;
        }
        Ok(array)
    }

    pub fn get(&self, index: usize) -> PyResult<T> {
        let data = self.data(index)?;
        let value = unsafe { *(data as *const T) };
        Ok(value)
    }

    pub fn set(&self, index: usize, value: T) -> PyResult<()> {
        self.is_writeable()?;
        let data = self.data(index)?;
        let element = unsafe { &mut *(data as *mut T) };
        *element = value;
        Ok(())
    }

    pub unsafe fn slice(&self) -> PyResult<&[T]> {
        self.is_contiguous()?;
        let obj: &PyArrayObject = self.as_ref();
        let ptr = obj.data as *const T;
        let size = self.size();
        let slice = unsafe { std::slice::from_raw_parts(ptr, size) };
        Ok(slice)
    }

    pub unsafe fn slice_mut(&self) -> PyResult<&mut [T]> {
        self.is_contiguous()?;
        self.is_writeable()?;
        let obj: &PyArrayObject = self.as_ref();
        let ptr = obj.data as *mut T;
        let size = self.size();
        let slice = unsafe { std::slice::from_raw_parts_mut(ptr, size) };
        Ok(slice)
    }

    pub fn zeros<'py>(py: Python<'py>, shape: &[usize]) -> PyResult<&'py Self> {
        let api = api(py);
        let zeros = unsafe { *api.zeros };
        let dtype = T::dtype(py)?;
        let (ndim, shape) = Self::try_shape(shape)?;
        let array = zeros(
            ndim,
            shape.as_ptr() as *const npy_intp,
            dtype.as_ptr(),
            0,
        );
        if PyErr::occurred(py) {
            match PyErr::take(py) {
                None => unreachable!(),
                Some(err) => return Err(err),
            }
        }
        let array = unsafe { &*(array as *const Self) };
        Ok(array)
    }
}

// Private interface.
impl<T> PyArray<T> {
    fn is_contiguous(&self) -> PyResult<()> {
        let obj: &PyArrayObject = self.as_ref();
        if obj.flags & PyArrayFlags::C_CONTIGUOUS == 0 {
            Err(PyValueError::new_err("memory is not C-contiguous"))
        } else {
            Ok(())
        }
    }

    fn is_writeable(&self) -> PyResult<()> {
        let obj: &PyArrayObject = self.as_ref();
        if obj.flags & PyArrayFlags::WRITEABLE == 0 {
            Err(PyValueError::new_err("assignment destination is read-only"))
        } else {
            Ok(())
        }
    }

    fn try_shape(shape: &[usize]) -> PyResult<(i32, Vec<npy_intp>)> {
        let ndim = match i32::try_from(shape.len()) {
            Err(_) => return Err(PyValueError::new_err(format!(
                "bad i32 value ({})",
                shape.len(),
            ))),
            Ok(ndim) => ndim,
        };
        let mut raw_shape = Vec::<npy_intp>::with_capacity(shape.len());
        for v in shape.iter() {
            let v = Self::try_size(*v)?;
            raw_shape.push(v);
        }
        Ok((ndim, raw_shape))
    }

    fn try_size(size: usize) -> PyResult<npy_intp> {
        match npy_intp::try_from(size) {
            Err(_) => Err(PyValueError::new_err(format!(
                "bad npy_intp value ({})",
                size,
            ))),
            Ok(size) => Ok(size),
        }
    }
}

// Traits implementations.
impl<T> AsRef<PyArrayObject> for PyArray<T> {
    #[inline]
    fn as_ref(&self) -> &PyArrayObject {
        unsafe { &*self.0.0.get() }
    }
}

impl<T> Deref for PyArray<T> {
    type Target = PyUntypedArray;

    #[inline]
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl<'a, T> From<&'a PyArray<T>> for &'a PyUntypedArray {
    #[inline]
    fn from(ob: &'a PyArray<T>) -> &'a PyUntypedArray {
        unsafe { &*(ob as *const PyArray<T> as *const PyUntypedArray) }
    }
}

impl<'a, T> TryFrom<&'a PyUntypedArray> for &'a PyArray<T>
where
    T: Dtype,
{
    type Error = PyErr;

    #[inline]
    fn try_from(ob: &'a PyUntypedArray) -> Result<&'a PyArray<T>, Self::Error> {
        let dtype = T::dtype(ob.py())?;
        let descr = unsafe { &*ob.0.get() }.descr;
        let mut same = descr as * const ffi::PyObject == dtype.as_ptr();
        if !same {
            let api = api(ob.py());
            let equiv_types = unsafe { *api.equiv_types };
            same = equiv_types(descr as * mut ffi::PyObject, dtype.as_ptr()) != 0;
        }
        if same {
            Ok(unsafe { &*(ob as *const PyUntypedArray as *const PyArray<T>) })
        } else {
            Err(PyTypeError::new_err(format!(
                "bad dtype (expected {:?}, found {:?})",
                dtype.as_ref(ob.py()),
                unsafe { &*(descr as *mut PyAny) },
            )))
        }
    }
}

impl<'py, T> FromPyObject<'py> for &'py PyArray<T>
where
    T: Dtype,
{
    fn extract(obj: &'py PyAny) -> PyResult<Self> {
        let untyped: &PyUntypedArray = FromPyObject::extract(obj)?;
        let typed: &PyArray<T> = std::convert::TryFrom::try_from(untyped)?;
        Ok(typed)
    }
}


// ===============================================================================================
//
// D-types.
//
// ===============================================================================================

pub trait Dtype {
    fn dtype(py: Python) -> PyResult<PyObject>;
}

impl Dtype for bool {
    #[inline]
    fn dtype(py: Python) -> PyResult<PyObject> {
        Ok(api(py).dtype_bool.clone_ref(py))
    }
}

impl Dtype for Float {
    #[inline]
    fn dtype(py: Python) -> PyResult<PyObject> {
        Ok(api(py).dtype_float.clone_ref(py))
    }
}

impl Dtype for i32 {
    #[inline]
    fn dtype(py: Python) -> PyResult<PyObject> {
        Ok(api(py).dtype_int32.clone_ref(py))
    }
}

impl Dtype for ElectronicShell {
    #[inline]
    fn dtype(py: Python) -> PyResult<PyObject> {
        Ok(api(py).dtype_shell.clone_ref(py))
    }
}

impl Dtype for CState {
    #[inline]
    fn dtype(py: Python) -> PyResult<PyObject> {
        Ok(api(py).dtype_state.clone_ref(py))
    }
}

impl Dtype for usize {
    #[inline]
    fn dtype(py: Python) -> PyResult<PyObject> {
        Ok(api(py).dtype_usize.clone_ref(py))
    }
}


//================================================================================================
// Control flags for Numpy arrays.
//================================================================================================

pub enum PyArrayFlags {
    ReadOnly,
    ReadWrite,
}

impl PyArrayFlags {
    pub const C_CONTIGUOUS: c_int = 0x0001;
    pub const WRITEABLE:    c_int = 0x0400;
}

impl From<PyArrayFlags> for c_int {
    fn from(value: PyArrayFlags) -> Self {
        match value {
            PyArrayFlags::ReadOnly =>  PyArrayFlags::C_CONTIGUOUS,
            PyArrayFlags::ReadWrite => PyArrayFlags::C_CONTIGUOUS | PyArrayFlags::WRITEABLE,
        }
    }
}

// ===============================================================================================
//
// Typed scalar.
//
// ===============================================================================================

#[repr(transparent)]
pub struct PyScalar<T>(ffi::PyObject, PhantomData<T>);

// Public interface.
impl<T> PyScalar<T>
where
    T: Copy + Default + Dtype + Sized,
{
    pub fn get(&self) -> PyResult<T> {
        let py = self.py();
        let api = api(py);
        let scalar_as_ctype = unsafe { *api.scalar_as_ctype };
        let mut data = T::default();
        scalar_as_ctype(
            &self.0 as *const ffi::PyObject as *mut ffi::PyObject,
            &mut data as *mut T as *mut c_void,
        );
        if PyErr::occurred(py) {
            match PyErr::take(py) {
                None => unreachable!(),
                Some(err) => return Err(err),
            }
        }
        Ok(data)
    }

    pub fn new<'py>(py: Python<'py>, value: T) -> PyResult<&'py Self> {
        let api = api(py);
        let scalar = unsafe { *api.scalar };
        let dtype = T::dtype(py)?;
        let scalar = scalar(
            &value as *const T as *mut c_void,
            dtype.as_ptr(),
            std::ptr::null_mut(),
        );
        if PyErr::occurred(py) {
            match PyErr::take(py) {
                None => unreachable!(),
                Some(err) => return Err(err),
            }
        }
        let scalar = unsafe { &*(scalar as *const Self) };
        Ok(scalar)
    }
}

// Private interface.
impl<T> PyScalar<T> {
    #[inline]
    fn transmute(&self) -> &PyAny {
        unsafe { &*(self as * const Self as *const PyAny) }
    }
}

// Trait implementations.
unsafe impl<T> PyNativeType for PyScalar<T> {}

impl<T> AsRef<PyAny> for PyScalar<T> {
    #[inline]
    fn as_ref(&self) -> &PyAny { self.transmute() }
}

impl<T> Deref for PyScalar<T> {
    type Target = PyAny;

    #[inline]
    fn deref(&self) -> &Self::Target { self.transmute() }
}

impl<T> IntoPy<Py<PyScalar<T>>> for &'_ PyScalar<T> {
    #[inline]
    fn into_py(self, py: Python) -> Py<PyScalar<T>> {
        unsafe { Py::from_borrowed_ptr(py, self.as_ptr()) }
    }
}

impl<T> ToPyObject for PyScalar<T> {
    #[inline]
    fn to_object(&self, py: Python) -> PyObject {
        unsafe { PyObject::from_borrowed_ptr(py, self.as_ptr()) }
    }
}


// ===============================================================================================
//
// Arguments conversion.
//
// ===============================================================================================

#[derive(pyo3::FromPyObject)]
pub enum ArrayOrFloat<'a> {
    Array(&'a PyArray<Float>),
    Float(Float),
}

impl<'a> ArrayOrFloat<'a> {
    pub fn get(&self, index: usize) -> PyResult<Float> {
        match self {
            Self::Array(a) => a.get(index),
            Self::Float(s) => Ok(*s),
        }
    }

    pub fn is_float(&self) -> bool {
        match self {
            Self::Array(_) => false,
            Self::Float(_) => true,
        }
    }

    pub fn size(&self) -> usize {
        match self {
            Self::Array(a) => a.size(),
            Self::Float(_) => 1,
        }
    }
}

#[derive(pyo3::FromPyObject)]
pub enum ArrayOrFloat3<'a> {
    Array(&'a PyArray<Float>),
    Float3(Float3),
}

impl IntoPy<PyObject> for Float3 {
    fn into_py(self, py: Python) -> PyObject {
        let result = PyArray::<Float>::empty(py, &[3]).unwrap();
        result.set(0, self.0).unwrap();
        result.set(1, self.1).unwrap();
        result.set(2, self.2).unwrap();
        result.readonly();
        result.into_py(py)
    }
}

impl IntoPy<PyObject> for Float3x3 {
    fn into_py(self, py: Python) -> PyObject {
        let result = PyArray::<Float>::empty(py, &[3, 3]).unwrap();
        let data: &[Float] = self.as_ref();
        for i in 0..9 {
            result.set(i, data[i]).unwrap();
        }
        result.readonly();
        result.into_py(py)
    }
}

#[derive(pyo3::FromPyObject)]
pub enum FloatOrFloat3 {
    Float(Float),
    Float3(Float3),
}

#[derive(pyo3::FromPyObject)]
pub enum ShapeArg {
    Scalar(usize),
    Vector(Vec<usize>),
}

impl From<ShapeArg> for Vec<usize> {
    fn from(value: ShapeArg) -> Self {
        match value {
            ShapeArg::Scalar(value) => vec![value],
            ShapeArg::Vector(value) => value,
        }
    }
}
