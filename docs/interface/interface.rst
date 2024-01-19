.. _goupil_interface:

`goupil_interface`_
===================

----

This is an initialised C interface to Goupil, containing constructors for
objects that interact with a Monte Carlo geometry.

.. c:struct:: goupil_interface

   .. c:function:: struct goupil_geometry_definition * new_geometry_definition(void)

      Returns a high level definition (Goupil-wise) of a Monte Carlo geometry.

      .. note::

         The Monte Carlo geometry associated with the returned *definition*
         object is expected to remain unchanged and valid throughout its
         lifetime.

   .. c:function:: struct goupil_geometry_tracer * new_geometry_tracer(struct goupil_geometry_definition * definition)

      Returns a tracer for the given geometry *definition*.

      .. note::

         The geometry *definition* associated with the returned *tracer*
         object is expected to remain unchanged and valid throughout its
         lifetime.
