.. _FormFactor:

`FormFactor`_
=============

This class provides read-only views of tables that relate to the form-factor of
a collision process. The form-factor values are tabulated over a grid of
:math:`m` nodes, which depend on the transfered momentum :math:`q`.


Constructor
-----------

.. py:class:: FormFactor(*args, **kwargs)

.. warning::

   Form-factor views should not be instantiated directly, but only from a
   :doc:`material_record`. Direct instantiation will result in a
   :external:py:class:`TypeError`.


Attributes
----------

.. py:attribute:: FormFactor.material
   :type: MaterialRecord

   The collision's target material.

.. py:attribute:: FormFactor.momenta
   :type: numpy.ndarray

   The :math:`m` momenta values (:math:`q`), for which the form-factor was
   computed.

.. py:attribute:: FormFactor.process
   :type: str

   A description of the interaction process.

.. py:attribute:: FormFactor.values
   :type: numpy.ndarray

   The :math:`m` computed form-factor values.


Methods
-------

.. py:method:: FormFactor.__call__(momentum)

   Returns interpolated form-factor value(s) for the specified *momentum*
   (:math:`q`) which can be either a float or a
   :external:py:class:`numpy.ndarray`.
