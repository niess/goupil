.. _CrossSection:

`CrossSection`_
===============

This class provides read-only views of tables that relate to the cross-section
of a collision process. The cross-section values are tabulated over a grid of
:math:`m` nodes, which depend on the energy :math:`\nu_i` of the incoming
photon.


Constructor
-----------

.. py:class:: CrossSection(*args, **kwargs)

.. warning::

   Cross-section views should not be instantiated directly, but only from a
   :doc:`material_record`. Direct instantiation will result in a
   :external:py:class:`TypeError`.


Attributes
----------

.. py:attribute:: CrossSection.energies
   :type: numpy.ndarray

   The :math:`m` energy values of the incoming photon (:math:`\nu_i`), for which
   the cross-section was pre-computed.

.. py:attribute:: CrossSection.material
   :type: MaterialRecord

   The collision's target material.

.. py:attribute:: CrossSection.process
   :type: str

   A description of the interaction process.

.. py:attribute:: CrossSection.values
   :type: numpy.ndarray

   The :math:`m` pre-computed cross-section values.


Methods
-------

.. py:method:: CrossSection.__call__(energy)

   Returns interpolated cross-section value(s) for the specified *energy*
   (:math:`\nu_i`) which can be either a float or a
   :external:py:class:`numpy.ndarray`.
