.. _TransportEngine:

`TransportEngine`_
==================

This class is the primary interface for the Monte Carlo simulation. It allows
one to convert a geometric description of materials distributions to an
operational Monte-Carlo simulation of photons transport.


Constructor
-----------

.. py:class:: TransportEngine(*args, **kwargs)

   Generates a new Monte Carlo transport engine. Optional arguments can be used
   to initialise the attributes described below (following lexicographic order
   for positional syntax).

   >>> engine = goupil.TransportEngine(geometry)


Attributes
----------

.. py:attribute:: TransportEngine.geometry
   :type: ExternalGeometry | SimpleGeometry

   The geometry seen by Monte Carlo trajectories.

.. py:attribute:: TransportEngine.random
   :type: RandomStream

   The pseudo-random stream utilised by the Monte Carlo engine.

.. py:attribute:: TransportEngine.registry
   :type: MaterialRegistry

   The registry stores pre-computed material tables for use in Monte Carlo
   transport. It is populated from the :py:attr:`geometry
   <TransportEngine.geometry>` description when the :py:meth:`compile
   <TransportEngine.compile>` method is called.

.. py:attribute:: TransportEngine.settings
   :type: TransportSettings

   Settings that control the Monte Carlo simulation. These settings are also
   accessible directly from the engine for convenience, such as

   >>> engine.compton_model
   "Scattering Function"


Methods
-------

.. py:method:: TransportEngine.compile(mode=None, atomic_data=None, **kwargs)

   Compiles material tables based on the current engine settings. This function
   fills the engine's material :py:attr:`registry <TransportEngine.registry>`
   based on the :py:attr:`geometry <TransportEngine.geometry>` description. The
   *mode* parameter determines the strategy used to construct the material
   tables. Available options include:

   .. list-table::

      * - :python:`"All"`
        - Compute all possible tables (not recommended for standard usage).

      * - :python:`"Backward"`
        - Compute only tables needed for backward transport.

      * - :python:`"Both"`
        - Compute tables needed for both forward and backward transport.

      * - :python:`"Forward"`
        - Compute only tables needed for forward transport.

   If no explicit *mode* is specified, then the engine will follow either the
   :python:`"Backward"` or :python:`"Forward"` strategy, depending on the
   configured :py:attr:`mode <TransportSettings.mode>`.

   The optional argument *atomic_data* can be used to specify non-default atomic
   data. Additional build options can be provided as keyword arguments
   (*kwargs*). For more details, refer to the :py:meth:`compute
   <MaterialRegistry.compute>` method in the :doc:`material_registry`.

.. py:method:: TransportEngine.transport(states, sources_energies=None) -> numpy.ndarray

   Performs a Monte Carlo transport of photon *states*, which were, e.g.,
   initially generated using the :doc:`states` function. The optional
   *sources_energies* argument should be used to specify the energies of volume
   sources, either as a :external:py:class:`float` or a
   :external:py:class:`numpy.ndarray` of floats, for backward mode. Upon
   completion, the function returns a :external:py:class:`numpy.ndarray` of
   integers that map the end-condition of each simulated trajectory. Please
   refer to :doc:`transport_status` for the precise definition of these numbers.

   .. warning::

      The *states* array is modified in-place. This means that the final states
      overwrite the initial ones. If the initial values need to be conserved,
      then a copy of the states array must be made before calling this method.
