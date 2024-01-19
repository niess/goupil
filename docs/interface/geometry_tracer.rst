.. _goupil_geometry_tracer:

`goupil_geometry_tracer`_
=========================

----

This objects enables tracing operations on a Monte Carlo geometry. For
information on the geometry, refer to the related
:c:struct:`goupil_geometry_definition` object.


.. c:struct:: goupil_geometry_tracer

   .. c:var:: const struct goupil_geometry_definition * definition

      The related geometry definition.

   .. c:function:: void destroy(struct goupil_geometry_definition * self)

      Destroys the geometry tracer, e.g. freeing any dynamically allocated
      memory.

      .. warning::

         The *tracer* object should not be accessed anymore after invocation
         of this function.

   .. c:function:: struct goupil_float3 position(const struct goupil_geometry_tracer * self)

      Returns the current position of the ray tracer.

   .. c:function:: void reset(struct goupil_geometry_tracer * self, struct goupil_float3 position, struct goupil_float3 direction)

      Reset the ray tracer in preparation for a new geometry traversal using
      the *position* and *direction* arguments to set the starting position and
      direction of the ray.

   .. c:function:: size_t sector(const struct goupil_geometry_tracer * self)

      Returns the index of the current geometry sector.

   .. c:function:: goupil_float_t trace(struct goupil_geometry_tracer * self, goupil_float_t max_length)

      Returns the length :math:`s` of the next tracing step, which must satisfy

      .. math::

         s \leq \min\left(s_\text{geo}, s_\text{max}\right),

      where :math:`s_\text{geo}` represents the distance to the closest geometry
      interface along the ray, and :math:`s_\text{max}` represents the maximum
      allowed step length given by *max_length*.

   .. c:function:: void update(struct goupil_geometry_tracer * self, goupil_float_t length, struct goupil_float3 direction)

      Applies a tracing step with the specified *length*. After completion, the
      *direction* of the ray is updated to the specified value.
