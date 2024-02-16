.. _MaterialDefinition:

`MaterialDefinition`_
=====================

This class is an abstract representation of the target material of a collision.
It is used, for example, in defining the geometry of Monte Carlo simulations.
Materials are treated as ideal gases in :mod:`goupil`, meaning they are mixtures
of non-interacting atoms. Therefore, a material is fully defined by its atomic
composition.


Constructor
-----------

.. py:class:: MaterialDefinition(name, mass_composition=None, mole_composition=None)

   Basic materials can be defined by their chemical formula, for example

   >>> goupil.MaterialDefinition("H2O")
   H2O

   More sophisticated constructions are possible by explicitly defining the mass
   or mole composition as a sequence of weighted :doc:`atomic_element` (or
   their corresponding :external:py:class:`str` symbols), e.g

   >>> goupil.MaterialDefinition(
   ...     name = "Water",
   ...     mole_composition = ((2, "H"), (1, "O"))
   ... )
   ...
   Water

   Composites, which are mixtures of different materials, may also be defined as

   >>> goupil.MaterialDefinition(
   ...     name = "WetSand",
   ...     mass_composition = ((0.7, "SiO2"), (0.3, "H2O"))
   ... )
   ...
   WetSand

   .. note::

      When explicitly defining the mass or mole composition of the material, it
      is also necessary to specify a name.

Attributes
----------

.. py:attribute:: MaterialDefinition.mass
   :type: float

   The molar mass of the substance, expressed in grams per mole. For compounds,
   it is the sum of the molar masses of its constituent parts, and for mixtures,
   it is the average molar mass of its constituent parts. For instance

   .. testsetup::
      >>> H2O = goupil.MaterialDefinition("H2O")

   >>> H2O.mass
   18.0167

.. py:attribute:: MaterialDefinition.mass_composition
   :type: tuple[float, AtomicElement]

   The composition of the material in terms of mass, represented as a tuple of
   :doc:`atomic_element` with corresponding mass fractions.

.. py:attribute:: MaterialDefinition.mole_composition
   :type: tuple[float, AtomicElement]

   The material's atomic composition, represented as a tuple of
   :doc:`atomic_element` with corresponding mole fractions.

.. py:attribute:: MaterialDefinition.name
   :type: str

   A name that identifies the material and is used to index it in a
   :doc:`material_registry`.


Methods
-------

.. py:method:: MaterialDefinition.electrons() -> ElectronicStructure

   Calculates the electronic structure based on the atomic composition of the
   material.
