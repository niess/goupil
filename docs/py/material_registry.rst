.. _MaterialRegistry:

`MaterialRegistry`_
===================

This class acts as a container for simulation materials and handles the
precomputation of properties, such as cross-sections, that are important for
Monte Carlo transport. To add new materials to the registry, a
:doc:`material_definition` must be provided. The relevant physical properties of
the registered materials are then calculated based on the Monte Carlo settings
provided in the form of :doc:`transport_settings`. These properties can be
retrieved as a :doc:`material_record`.


Constructor
-----------

.. py:class:: MaterialRegistry(*args: MaterialDefinition)

   To initialise a new registry with some materials, specify their definitions
   as positional arguments. For example

   >>> goupil.MaterialRegistry(H2O, SiO2)
   {H20, SiO2}


Methods
-------

.. py:method:: MaterialRegistry.__len__()

   Returns the number of registered materials. For example

   >>> len(registry)
   2

.. py:method:: MaterialRegistry.__getitem__(name: str)

   Returns the :doc:`material_record` corresponding to the registered material
   name. For example

   >>> registry["H2O"]
   H2O

.. py:method:: MaterialRegistry.add(material: MaterialDefinition)

   Registers a material. The material must be an instance of a
   :doc:`material_definition`. The registry uses the material name as an index.

   .. note::

      If a material with the same name but with a different definition has
      already been registered, an error will be raised.

.. py:method:: MaterialRegistry.compute(settings=None, shape=None, precision=None, **kwargs)

   Calculates physical properties of all registered materials based on the
   provided `settings` which must be a :doc:`transport_settings` instance.
   Optionally, transport settings can also be detailed explicitly as named
   arguments (`kwargs`). Physical properties are tabulated over grids (1d or
   2d), the `shape` of which can be optionally specified. The `precision`
   parameter allows one to specify the numerical accuracy of computations,
   relative to 1 being the default.

.. py:method:: MaterialRegistry.load_elements(path: str=None)

   Loads atomic elements data from a specific `path`. If no path is given
   default data are loaded. These data are necessary for computing material
   properties.
