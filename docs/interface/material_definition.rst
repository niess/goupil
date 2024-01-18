:c:struct:`goupil_material_definition`
======================================

.. _material_definition:

----

This object provides information on a Monte Carlo material, as per Goupil.

.. c:struct:: goupil_material_definition

   .. c:function:: size_t composition_len(const struct goupil_material_definition * self)

      Indicates the number of distinct atomic elements that make up the
      material.

   .. c:function:: const struct goupil_weighted_element get_composition(const struct goupil_material_definition * self, size_t index)

      Provides information about the atomic element corresponding to the
      specified *index* in the list of elements that make up the material.

   .. c:function:: const char * name(const struct goupil_material_definition * self)

      Returns the material's name.
