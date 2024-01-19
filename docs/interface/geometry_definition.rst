.. _goupil_geometry_definition:

`goupil_geometry_definition`_
=============================

----

This objects provides a high level definition of a Monte Carlo geometry,
according to Goupil. For tracing operations on the geometry, refer to the
related :c:struct:`goupil_geometry_tracer` object.

.. c:struct:: goupil_geometry_definition

   .. c:function:: void destroy(struct goupil_geometry_definition * self)

      Destroys the Monte Carlo geometry, e.g. freeing any dynamically allocated
      memory.

      .. warning::

         The *definition* object should not be accessed anymore after invocation
         of this function.

   .. c:function:: const struct goupil_material_definition * get_material(const struct goupil_geometry_definition * self, size_t index)

      Provides information about the material that corresponds to the specified
      *index* in the list of geometry materials.

   .. c:function:: const struct goupil_geometry_sector get_sector(const struct goupil_geometry_definition * self, size_t index)

      Provides information about the sector that corresponds to the specified
      *index* in the list of geometry sectors.

   .. c:function:: size_t materials_len(const struct goupil_geometry_definition * self)

      Returns the total number of materials defined by the Monte Carlo geometry.

   .. c:function:: size_t sectors_len(const struct goupil_geometry_definition * self)

      Returns the total number of sectors defined by the Monte Carlo geometry.

