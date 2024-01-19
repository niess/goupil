.. _SimpleGeometry:

`SimpleGeometry`_
=================

----

This class represents a basic unbounded geometry, composed of a single material
and extending infinitely.  It is intended solely for testing purposes. For
practical geometry implementations, please refer instead to the
:doc:`external_geometry` class.


Constructor
-----------

.. py:class:: SimpleGeometry(material, density)

   Generates a Monte Carlo geometry consisting of a single sector with the
   specified *material* and *density*. The argument for the *material* must be a
   :doc:`material_definition`. The density can be specified as a
   :external:py:class:`float`, indicating a uniform density value (in g/cm\
   :sup:`3`), or as a :doc:`density_gradient` object.


Attributes
----------

.. py:attribute:: SimpleGeometry.density
   :type: DensityGradient | float

   The medium density model. This attribute is mutable. For instance,

   >>> geometry.density = 2.65 # g/cm^3

.. py:attribute:: SimpleGeometry.material
   :type: MaterialDefinition

   The material that makes up the propagation medium. This attribute is **not**
   mutable.
