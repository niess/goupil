:class:`ExternalGeometry`
=========================

.. _ExternalGeometry:

----

This class encapsulates an external geometry, usually implemented through the
:doc:`C interface <../interface/index>`. The geometry is represented as a set of
sectors with a specific atomic composition, as defined by a
:doc:`material_definition`, and a density model.


Constructor
-----------

.. py:class:: ExternalGeometry(path)

   Loads a Monte Carlo geometry from a shared library located at the specified
   *path*. For example, on a Linux system,

   >>> geometry = goupil.ExternalGeometry("/path/to/libgeometry.so")


Attributes
----------

.. note::

   Geometry attributes are read-only. However, the corresponding data can be
   edited using the methods_ described below.

.. py:attribute:: ExternalGeometry.materials
   :type: tuple[MaterialDefinition]

   This attribute lists all geometry materials as a tuple.

.. py:attribute:: ExternalGeometry.sectors
   :type: tuple[int, float | DensityGradient, str]

   This attribute lists all geometry sectors as a tuple. The description of each
   sector is presented as a sub-tuple that includes the material index (as an
   :external:py:class:`int`), the density model, and a brief explanation (as a
   :external:py:class:`str`).


Methods
-------

.. py:method:: ExternalGeometry.update_material(index, material)

   Replaces the material at the specified *index* in the list of geometry
   :py:attr:`materials <ExternalGeometry.materials>`. The new *material* must be
   a :doc:`material_definition` object.

.. py:method:: ExternalGeometry.update_sector(index, material=None, density=None)

   Alters the *material* or *density* model of a sector, identified by its
   *index* in the list of geometry :py:attr:`sectors
   <ExternalGeometry.sectors>`.
