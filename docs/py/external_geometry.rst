.. _ExternalGeometry:

`ExternalGeometry`_
===================

This class encapsulates an external geometry, usually implemented through the
:doc:`C interface <../c/index>`. The geometry is represented as a set of sectors
with a specific atomic composition, as defined by a :doc:`material_definition`,
and a density model.


Constructor
-----------

.. py:class:: ExternalGeometry(path)

   Loads a Monte Carlo geometry from a shared library located at the specified
   *path*. For example, on a Linux system,

   >>> geometry = goupil.ExternalGeometry("/path/to/libgeometry.so") # doctest: +SKIP

.. warning::

   One must specify an absolute or explicit relative path to the geometry
   library, unless it is known system wide (e.g. located under
   :c:`LD_LIBRARY_PATH` on a Linux system).


Attributes
----------

.. note::

   :doc:`external_geometry` attributes are read-only. However, physical
   properties can be edited using the :py:meth:`update_material
   <ExternalGeometry.update_material>` and :py:meth:`update_sector
   <ExternalGeometry.update_sector>` methods described below.

.. py:attribute:: ExternalGeometry.lib
   :type: ~ctypes.CDLL

   This attribute provides a :external:py:mod:`ctypes`  representation of the
   geometry library.

.. tip::

   The :py:attr:`lib <ExternalGeometry.lib>` attribute grants access to
   user-defined functions that are embedded in the geometry library. These
   functions can be used for various purposes, such as initializing Monte Carlo
   states. This use case is exemplified by `examples/geant4/
   <https://github.com/niess/goupil/tree/master/examples/geant4>`_\ .

.. py:attribute:: ExternalGeometry.materials
   :type: tuple[MaterialDefinition]

   This attribute lists all geometry materials as a tuple.

.. py:attribute:: ExternalGeometry.path
   :type: str

   This attribute contains the path to the geometry library.

.. py:attribute:: ExternalGeometry.sectors
   :type: tuple[GeometrySector]

   This attribute lists all geometry sectors as a tuple.


Methods
-------

.. py:method:: ExternalGeometry.locate(states) -> numpy.ndarray

   Locates the specified *states* within the geometry. The input *states* must
   be a structured :external:py:class:`numpy.ndarray` containing the
   :python:`"position"` field, e.g. as returned by the :py:func:`states
   <states>` function. Upon completion, the function returns a
   :external:py:class:`numpy.ndarray` of sector indices.

.. py:method:: ExternalGeometry.material_index(name) -> int

   Returns the index of a material in the list of geometry :py:attr:`materials
   <ExternalGeometry.materials>` based on its :py:attr:`name
   <MaterialDefinition.name>`. For instance

   >>> geometry.material_index("CaCO3") # doctest: +SKIP
   1

.. py:method:: ExternalGeometry.sector_index(description) -> int

   Returns the index of a sector in the list of geometry :py:attr:`sectors
   <ExternalGeometry.sectors>` based on its :py:attr:`description
   <GeometrySector.description>`. For instance

   >>> geometry.sector_index("Atmosphere") # doctest: +SKIP
   0

.. py:method:: ExternalGeometry.trace(states, lengths=None, density=None) -> numpy.ndarray

   Casts rays through the geometry, starting from the specified *states*. The
   *states* must be a structured :external:py:class:`numpy.ndarray` containing
   the :python:`"position"` and :python:`"direction"` fields, e.g. as returned
   by the :py:func:`states <states>` function. Upon completion, this function
   returns a :external:py:class:`numpy.ndarray` containing the path length of
   rays in each geometry sector. Optionally, you can provide a *lengths*
   :external:py:class:`numpy.ndarray` of floats, or a single float, indicating
   the lengths of rays. If no *lengths* are specified, rays are traced until the
   geometry outer boundary.

   If the *density* parameter is set to :python:`True`, this function will
   return the column depth (grammage) along rays, in each sector, rather than
   the path length.

.. py:method:: ExternalGeometry.update_material(material, definition)

   Replaces the *material* at the specified :external:py:class:`int` index in
   the list of geometry :py:attr:`materials <ExternalGeometry.materials>`.
   Alternatively, the *material* argument could be a :external:py:class:`str`
   indicating the :py:attr:`name <MaterialDefinition.name>` of the material to
   be replaced. The *definition* argument should be in line with a
   :doc:`material_definition`. For instance

   >>> geometry.update_material("CaCO3", "SiO2") # doctest: +SKIP


.. py:method:: ExternalGeometry.update_sector(sector, material=None, density=None)

   Alters the *material* or *density* model of a *sector*, identified by its
   :external:py:class:`int` index in the list of geometry :py:attr:`sectors
   <ExternalGeometry.sectors>`. Alternatively, the *sector* argument could be
   a :external:py:class:`str` describing the sector to be altered. For instance

   >>> geometry.update_sector("Atmosphere", density=1.205E-03) # doctest: +SKIP
