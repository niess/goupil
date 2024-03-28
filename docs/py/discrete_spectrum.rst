.. _DiscreteSpectrum:

`DiscreteSpectrum`_
===================

This class represents a spectrum of discrete emission lines, for instance
corresponding to one or more radioactive isotopes. It can be used for sourcing
the energies of Monte Carlo states.


Constructor
-----------

.. py:class:: DiscreteSpectrum(energies, intensities=None, background=None, energy_min=None)

   This function creates a spectrum of distinct emission lines at specified
   *energies*, in MeV, which must be a sequence of :external:py:class:`float`.
   The corresponding *intensities* can also be specified as a sequence of
   floats. By default, all emission lines have equal probability.

   The optional *background* and *energy_min* arguments are only applicable when
   sampling the spectrum in backward mode (refer to the corresponding
   attributes_ below).

Attributes
----------

.. py:attribute:: DiscreteSpectrum.background
   :type: float

   The fraction (in [0,1]) of background events when sampling backwards.

.. py:attribute:: DiscreteSpectrum.energies
   :type: numpy.ndarray

   The energies, in MeV, of emission lines.

.. py:attribute:: DiscreteSpectrum.energy_min
   :type: numpy.ndarray

   The minimum energy of background events, in MeV, when sampling backwards.

.. py:attribute:: DiscreteSpectrum.intensities
   :type: numpy.ndarray | None

   The relative intensities of emissin lines.



Methods
-------

.. py:method:: DiscreteSpectrum.sample(states, engine=None, mode=None)

   Samples Monte Carlo state energies according to the spectrum emission lines.
   The *mode* argument indicates wether a conventional (forward) sampling should
   be performed or not. Possible values are :python:`"Backward"` or
   :python:`"Forward"`.

   If a :doc:`transport_engine` is provided (using the *engine* argument), then
   the sampling is configured based one the simulation :py:attr:`mode
   <TransportSettings.mode>`.

.. note::

   At return from the :py:meth:`sample <DiscreteSpectrum.sample>` method, the
   energies (and weights in backward mode) of Monte Carlo states are modified
   in-place.
