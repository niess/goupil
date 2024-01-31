.. _StratifiedGeometry:

`StratifiedGeometry`_
=====================

This class defines a stratified geometry, such as geological layers, represented
as :doc:`geometry_sector`. The layers are separated by
:doc:`topography_surface`. Within each layer, a uniform atomic composition is
assumed, as defined by a :doc:`material_definition`. However, the density may
vary e.g., following a :doc:`density_gradient`.


Constructor
-----------

.. py:class:: StratifiedGeometry(*args: GeometrySector | TopographySurface)

   Creates a Monte Carlo geometry that is stratified by using a sequence of
   alternating :doc:`geometry_sector` and :doc:`topography_surface`. 
   The elements are organised in reading order, with the first element of the
   sequence located on top of the geometry. For instance, the following

   >>> geometry = goupil.StratifiedGeometry(
   ...     goupil.GeometrySector("N2", 1.205E-03),
   ...     water_surface,
   ...     goupil.GeometrySector("H2O", 1.0),
   ...     soil_surface,
   ...     goupil.GeometrySector("SiO2", 2.0)
   ... )

   describes a vertical section of water covered by a nitrogen atmosphere and
   bounded below by a sandy soil.


Attributes
----------

.. note::

   The attributes of a :doc:`stratified_geometry` are read-only. Any
   modifications require rebuilding the geometry object.

.. py:attribute:: StratifiedGeometry.materials
   :type: tuple[MaterialDefinition]

   This attribute lists all geometry materials as a tuple.

.. py:attribute:: StratifiedGeometry.sectors
   :type: tuple[int, float | DensityGradient, str]

   This attribute lists all geometry sectors as a tuple. The description of each
   sector is presented as a sub-tuple that includes the material index (as an
   :external:py:class:`int`), the density model, and a brief explanation (as a
   :external:py:class:`str`).


Methods
-------

.. py:method:: StratifiedGeometry.locate(states) -> numpy.ndarray

   Locates the specified *states* within the geometry. The input *states* must
   be a structured :external:py:class:`numpy.ndarray` as returned by the
   :py:func:`states <states>` function. The function returns a
   :external:py:class:`numpy.ndarray` of sector indices.

.. py:method:: StratifiedGeometry.trace(states, lengths=None, density=None) -> numpy.ndarray

   Casts rays through the geometry, starting from the specified *states*. The
   *states* must be a structured :external:py:class:`numpy.ndarray` as returned
   by the :py:func:`states <states>` function. This function returns a
   :external:py:class:`numpy.ndarray` containing the path length of rays in each
   geometry sector. Optionally, you can provide a *lengths*
   :external:py:class:`numpy.ndarray` of floats, or a single float, indicating
   the lengths of rays. If no *lengths* are specified, rays are traced until the
   geometry outer boundary.

   If the *density* parameter is set to :python:`True`, this function will
   return the column depth (grammage) along rays, in each sector, rather than
   the path length.
