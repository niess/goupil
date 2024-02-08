.. _TransportStatus:

`TransportStatus`_
==================

This class serves as a namespace for return codes from the :py:meth:`transport
<TransportEngine.transport>` method. All attributes and methods described below
are defined at the class level. For instance,

>>> goupil.TransportStatus.EXIT
5


Constructor
-----------

.. py:class:: TransportStatus

   .. note::

      This class is not meant to be instantiated. Attempting to do so will
      result in a :external:py:class:`TypeError`.


Attributes
----------

.. py:attribute:: TransportStatus.ABSORBED

   Indicates that the photon was absorbed, e.g. in a photo-electric interaction.

.. py:attribute:: TransportStatus.BOUNDARY

   Indicates that the photon reached an inner transport :py:attr:`boundary
   <TransportSettings.boundary>`.

.. py:attribute:: TransportStatus.ENERGY_CONSTRAINT

   Indicates that the photon reached the energy constraint of a potential volume
   source, during backward transport.

.. py:attribute:: TransportStatus.ENERGY_MAX

   Indicates that the photon energy exceeded its maximum allowed value, e.g.
   according to material tables.

.. py:attribute:: TransportStatus.ENERGY_MIN

   Indicates that the photon energy dropped below its minimum allowed value,
   e.g. according to material tables.

.. py:attribute:: TransportStatus.EXIT

   Indicates that the photon exit the simulation :py:attr:`geometry
   <TransportEngine.geometry>`.

.. py:attribute:: TransportStatus.LENGTH_MAX

   Indicates that the photon path length reached a limit value, as specified by
   :py:attr:`TransportSettings.length_max`.


Methods
-------

.. py:method:: TransportStatus.str(code: int) -> str

   Returns a string representation of a transport status *code*. For instance,

   >>> goupil.TransportStatus.str(5)
   'Exit'
