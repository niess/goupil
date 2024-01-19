.. _TransportSettings:

`TransportSettings`_
====================

----

This class collects settings related to Monte Carlo transport. These settings
can be accessed as mutable attributes, which are described below.


Constructor
-----------

.. py:class:: TransportSettings()

   Generates a fresh set of configuration parameters for Monte Carlo transport.
   The set will be initialised with the default values listed below.

   .. list-table::

      * - absorption
        - :python:`"Discrete"`

      * - mode
        - :python:`"Forward"`

      * - rayleigh
        - :python:`True`

      * - volume_sources
        - :python:`True`


Attributes
----------

.. py:attribute:: TransportSettings.absorption
   :type: str | None

   This attribute controls the simulation mode for absorption processes. Below
   are the possible values.

   .. list-table::

      * - :python:`"Continuous"`
        - Absorption is considered a continuous process by weighting Monte Carlo
          trajectories.

      * - :python:`"Discrete"`
        - Absorption is modelled as a discrete interaction process, resulting in
          the random termination of Monte Carlo trajectories.

      * - :python:`None`
        - Absorption is not considered.

.. py:attribute:: TransportSettings.boundary
   :type: int | None

   This attribute represents an inner geometry boundary for Monte Carlo
   transport, which is specified as a sector index, for example, in an
   :doc:`external_geometry`. The Monte Carlo transport stops whenever a
   trajectory enters the corresponding sector.

.. py:attribute:: TransportSettings.compton_method
   :type: str

   This attribute determines the sampling method for Compton collisions. Refer
   to the :py:attr:`ComptonProcess.method` section for a list of available
   values.

.. py:attribute:: TransportSettings.compton_mode
   :type: str | None

   This attribute determines the simulation mode for Compton collisions. Refer
   to the :py:attr:`ComptonProcess.mode` section for a list of available
   values.

.. py:attribute:: TransportSettings.compton_model
   :type: str

   This attribute determines the physics model for Compton collisions. Refer
   to the :py:attr:`ComptonProcess.model` section for a list of available
   values.

.. py:attribute:: TransportSettings.energy_max
   :type: float | None

   This attribute specifies the maximum energy of transported photons, if
   not :python:`None`.

.. py:attribute:: TransportSettings.energy_min
   :type: float | None

   This attribute specifies the minimum energy of transported photons, if
   not :python:`None`.

.. py:attribute:: TransportSettings.length_max
   :type: float | None

   This attribute specifies the maximum path length of transported photons, if
   not :python:`None`.

.. py:attribute:: TransportSettings.mode
   :type: str

   This flag controls the direction of flow for Monte Carlo transport. Switching
   this flag between :python:`"Backward"` and :python:`"Forward"` results in a
   default set of settings being selected. Use the :py:attr:`compton_mode
   <TransportSettings.compton_mode>` attribute instead if thinner control is
   needed.

   .. list-table::

      * - :python:`"Backward"`
        - Reverse Monte Carlo transport.

      * - :python:`"Forward"`
        - Standard (forward) Monte Carlo transport.

.. py:attribute:: TransportSettings.rayleigh
   :type: bool

   Enable (true) or disable (false) Rayleigh collisions during the Monte Carlo
   transport.

.. py:attribute:: TransportSettings.volume_sources
   :type: bool

   This flag controls whether the backward Monte Carlo transport considers
   volume sources with discrete energies or not.
