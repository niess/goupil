.. _GeometrySector:

`GeometrySector`_
=================

This class represents a sector of a Monte Carlo geometry. A sector is
characterised by a uniform atomic composition, specified by a
:doc:`material_definition`. However, the material density might vary
continuously within a sector, for example following a :doc:`density_gradient`.

.. note::

   :doc:`geometry_sector` objects provide immutable representations of a sector
   of the geometry. That is, all attributes below are read-only.


Constructor
-----------

.. py:class:: GeometrySector(material, density, description=None)

   The *material* argument must be consistent with a :doc:`material_definition`.
   The `density` can be an instance of a :external:py:class:`float`, in which
   case a uniform density is assumed (in g/cm\ :sup:`3`), or a
   :doc:`density_gradient`. Optionaly, a `description` can be provided as a
   :external:py:class:`str`. For instance,

   >>> sector = goupil.GeometrySector("H2O", 1.0, "Water")


Attributes
----------

.. py:attribute:: GeometrySector.density
   :type: DensityGradient | float

   The density model of the geometry sector. A :external:py:class:`float` value
   indicates a uniform density, expressed in g/cm\ :sup:`3`.

.. py:attribute:: GeometrySector.description
   :type: None | str

   A textual description of the geometry sector.


.. py:attribute:: GeometrySector.material
   :type: MaterialDefinition

   The uniform atomic composition of the geometry sector.
