.. _AbsorptionProcess:

`AbsorptionProcess`_
====================

This class provides access to the implementation of absorption processes, that
is photoelectric interactions and :math:`e^+ e^-` pair production.

.. warning::

   All methods described below are defined at the class level.


Constructor
-----------

.. py:class:: AbsorptionProcess

   .. note::

      This class serves as a namespace and is not intended for instantiation.
      Any attempt to do so will result in a :external:py:class:`TypeError`.


Methods
-------

.. py:method:: AbsorptionProcess.cross_section(energy, material)

   Computes the total cross-section for absorption of a photon with a specified
   *energy* (in MeV) by an atom of a given *material*. The *energy* can be a
   :external:py:class:`float` or a :external:py:class:`numpy.ndarray` of floats.
   The *material* must be consistent with a :doc:`material_definition`. For
   instance,

   >>> goupil.AbsorptionProcess.cross_section(0.1, "H2O")
   8.264...e-26
