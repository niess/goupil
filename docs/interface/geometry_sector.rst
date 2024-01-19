.. _goupil_geometry_sector:

`goupil_geometry_sector`_
=========================

----

This structure stores information related to a geometry sector.

.. c:struct:: goupil_geometry_sector

   .. c:var:: size_t material

      Index of the constitutive material.

   .. c:var:: goupil_float_t density

      Bulk density of this geometry sector, in g/cm\ :sup:`3`.

   .. c:var:: const char * definition

      A brief definition of this geometry sector.

      .. note::

         Providing a definition is optional. This field may be :c:`NULL`.
