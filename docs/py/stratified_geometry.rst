.. _StratifiedGeometry:

`StratifiedGeometry`_
=====================

This class defines a stratified geometry, such as geological layers, represented
by a collection of :doc:`geometry_sector`. The layers are separated by
:doc:`topography_surface` objects. Within each layer, a uniform atomic
composition is assumed, specified by a :doc:`material_definition`. The
density might vary as a :doc:`density_gradient`.


Constructor
-----------

.. py:class:: StratifiedGeometry(*args: GeometrySector | TopographySurface)

   Creates a stratified Monte Carlo geometry from a sequence of
   alternating :doc:`geometry_sector` and :doc:`topography_surface` objects.
   The geometry is specified in reading order, with the first element of the
   sequence located on top of the geometry. For instance, the following

   .. doctest::
      :hide:

      >>> water_surface = goupil.TopographyMap((-1e5, 1e5), (-1e5, 1e5))
      >>> soil_surface = goupil.TopographyMap((-1e5, 1e5), (-1e5, 1e5), z=-1e2)

   >>> geometry = goupil.StratifiedGeometry(
   ...     goupil.GeometrySector("N2", 1.205E-03, "Atmosphere"),
   ...     water_surface,
   ...     goupil.GeometrySector("H2O", 1.0, "Water"),
   ...     soil_surface,
   ...     goupil.GeometrySector("SiO2", 2.0, "Soil")
   ... )

   defines a vertical section of water covered by a nitrogen atmosphere and
   bounded below by a sandy soil.

.. note::

   In the previous example, the :python:`soil_surface` may extend above the
   :python:`water_surface`. This means that lower layers take precedence over
   higher ones in case of overlaps.


Attributes
----------

.. note::

   The attributes of a :doc:`stratified_geometry` are read-only. Any
   modifications require rebuilding the geometry object.

.. py:attribute:: StratifiedGeometry.materials
   :type: tuple[MaterialDefinition]

   This attribute lists all geometry materials as a tuple.

.. py:attribute:: StratifiedGeometry.sectors
   :type: tuple[GeometrySector]

   This attribute lists all geometry sectors as a tuple.

.. warning::

   Geometry sectors are stored in indexing order. That is, :python:`sectors[0]`
   corresponds to the bottom layer.


Methods
-------

.. py:method:: StratifiedGeometry.locate(states) -> numpy.ndarray

   Locates the specified *states* within the geometry. The input *states* must
   be a structured :external:py:class:`numpy.ndarray` as returned by the
   :py:func:`states <states>` function. The function returns a
   :external:py:class:`numpy.ndarray` of sector indices.

.. py:method:: StratifiedGeometry.material_index(name) -> int

   Returns the index of a material in the list of geometry :py:attr:`materials
   <StratifiedGeometry.materials>` based on its :py:attr:`name
   <MaterialDefinition.name>`. For instance

   >>> geometry.material_index("SiO2")
   0

.. py:method:: StratifiedGeometry.sector_index(description) -> int

   Returns the index of a sector in the list of geometry :py:attr:`sectors
   <StratifiedGeometry.sectors>` based on its :py:attr:`description
   <GeometrySector.description>`. For instance

   >>> geometry.sector_index("Atmosphere")
   2

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

.. note::

   The `Turtle <https://niess.github.io/turtle-pages/>`_ algorithm is used to
   perform ray tracing. For more information, refer to [Niess20]_.

.. py:method:: StratifiedGeometry.z(x, y, grid=None) -> numpy.ndarray

   Returns the elevation values of each :doc:`topography_surface` at coordinates
   :math:`(x, y)`. The `x` and `y` arguments can be :external:py:class:`float`
   or :external:py:class:`numpy.ndarray` with consistent sizes. If `grid` is set
   to :python:`True`, elevation values are computed over a grid that corresponds
   to the outer product of `x` and `y`, similar to the
   :py:meth:`TopographyMap.__call__` method.
