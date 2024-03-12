Overview
========

.. _description:

This section provides a practical overview of Goupil illustrated with examples.
Refer to the :doc:`Python <py/index>` and :doc:`C <py/index>` interface sections
of this document for a detailed description of the user interfaces.

.. note::

   The following instructions assume technical familiarity with Python and prior
   experience in using :external:py:class:`numpy.ndarray`.


Installing goupil
-----------------

Goupil is distributed as a Python 3 module, available from PyPI. For instance,
it can be installed as:

.. code:: console

   $ pip3 install goupil


Specifying a geometry
---------------------

To begin with, we will define a Monte Carlo geometry. For this overview we will
use a :doc:`py/stratified_geometry` object with two layers: a lower layer of
limestone and an upper layer filled with air.

.. note::

   Goupil also allows an :doc:`py/external_geometry` to be plugged in via its
   :doc:`c/index`. In particular, Goupil includes a :doc:`geant4`.

We start by importing the :doc:`goupil <py/index>` module as

>>> import goupil

Then, we specify the air composition as

>>> air = goupil.MaterialDefinition(
...     "Air",
...     mass_composition = (
...         (0.76, "N"),
...         (0.23, "O"),
...         (0.01, "Ar"),
...     )
... )

Let us consider an exponential :doc:`py/density_gradient` for the air
density as

>>> density = goupil.DensityGradient(1.225E-03, 1.04E+06)

which would describe the lower part of the Earth's atmosphere (i.e. the
troposphere).

We also need to define the interface between the lower and upper layers. This is
done using a Digital Elevation Model (DEM). For this example, we will create a
flat interface at :math:`z = 0` m, covering :math:`[-1, 1] \times [-1,
1]` km\ :sup:`2` in :math:`(x, y)`. This is done as

>>> interface = goupil.TopographyMap((-1e5, 1e5), (-1e5, 1e5), z=0)

.. note::

   More complex interfaces can be specified. See the :doc:`py/topography_map`
   and :doc:`py/topography_surface` objects for additional information.

Finally, we define the Monte Carlo geometry as

>>> geometry = goupil.StratifiedGeometry(
...     goupil.GeometrySector(air, density, "Atmosphere"),
...     interface,
...     goupil.GeometrySector("CaCO3", 2.8, "Ground")
... )

where we specify the composition of the limestone by its chemical formula
(CaCO3) and assume a uniform density of 2.8 g/cm\ :sup:`3`.

.. note::

   The geometry defined previously lacks upper and lower bounds. Specifically,
   the Atmosphere sector extends to :math:`z \to +\infty` and the Ground sector
   extends to :math:`z \to -\infty`.


Running a simulation
--------------------

The Monte Carlo transport of photons is managed by a :doc:`py/transport_engine`
taking in charge a specific geometry. A :doc:`py/transport_engine` is created
as:

>>> engine = goupil.TransportEngine(geometry)

Each engine has its own :doc:`py/random_stream`, which can be accessed through
the :py:attr:`random <TransportEngine.random>` attribute. By default, this
stream is seeded from the system entropy. For example purposes, let us set a
specific seed value.

>>> engine.random.seed = 123456789

.. note::

   Setting a seed has the effect of reseting the pseudo-random stream.

The transport engine is set to perform a conventional (forward) Monte Carlo
simulation by default. Let us instead configure the engine for backwards
transport. This is done as:

>>> engine.mode = "Backwards"

.. tip::

   See :doc:`py/transport_settings` for a summary of configurable parameters.


Then, let us define a set of :python:`100` Monte Carlo states representing
photons with an energy of :python:`0.5` MeV, located at :math:`z =
100` cm, that is 1 m above the ground. This is done with the
:doc:`py/states` function as

>>> states = goupil.states(100, energy=0.5, position=(0,0,1e2))

The :doc:`py/states` function returns a `numpy structured array
<https://numpy.org/doc/stable/user/basics.rec.html>`_ of states, containing the
photons energies, their locations, etc. Since we perform a backwards simulation,
these states represent expected final states, e.g., at a particular observation
point. In practice, it is also necessary to specify the arrival directions of
these photons. However, for the purposes of this overview, default values will
be used. That is

>>> states["direction"]
array([[0., 0., 1.],
...

Then, let us backwards propagate the expected photons through the geometry. This
is done with the :py:meth:`transport <TransportEngine.transport>` method, as:

>>> status = engine.transport(states, source_energies=1.0)

.. warning::

   The :py:meth:`transport <TransportEngine.transport>` method modifies the
   *states* array in-place. Thus, after completion, the *states* array will
   contain the propagated photons instead of the original ones.

The second argument, *source_energies*, requires further explanation. When
running a backwards Monte Carlo simulation, information about sources is needed
to correctly terminate the transport. Goupil considers two types of sources:

- Surface sources with a distributed energy spectrum, such as an external flux
  of gamma-rays.
- Volume sources with a discrete energy spectrum, such as scattered
  radio-isotopes.

In the previous example, a constant value of :python:`1.0` MeV was assumed for
the energy of volume sources.

.. note::

   The *source_energies* argument should be omitted if there are no volume
   sources or in the case of a forward Monte Carlo.

.. tip::

   In a backwards transport, contained surface sources (i.e. not located on an
   outer boundary of the geometry) can be specified as a sector
   :py:attr:`boundary <TransportSettings.boundary>` at the level of the
   :doc:`py/transport_engine`.


Inspecting results
------------------

The :py:meth:`transport <TransportEngine.transport>` method returns an array of
integer codes (:doc:`py/transport_status`) which indicate the termination
condition for each propagated photon. For instance, backwards propagated photons
that are consistent with a volume source can be selected as follows:

>>> constrained = (status == goupil.TransportStatus.ENERGY_CONSTRAINT)

These photons should have an energy of :python:`1.0` MeV, as requested:

>>> states[constrained]["energy"]
array([1., 1., ...])

The corresponding geometry sectors can be located as:

>>> geometry.locate(states[constrained])
array([0, 0, ...], dtype=uint64)


Computing an estimate
---------------------

An important property that you will use is the transport weight (hereafter noted
:math:`\omega`) associated with each backwards propagated photon. These weights
are given as:

>>> weights = states["weight"]

A backwards Monte Carlo estimate of the gamma-ray flux for the expected
state :math:`\mathcal{S}_f` is given by

.. math::

   \phi(\mathcal{S}_f) \simeq \frac{1}{N} \sum_{i=1}^N {
        \omega\left(\mathcal{S}_f,\mathcal{S_i}\right)
        S(\mathcal{S}_i)
   },

where the :math:`\mathcal{S}_i` denote the :math:`N` backwards sampled
photon states, and where the source term :math:`S` depends on the
termination condition of each Monte Carlo event, as

.. math::

   S(\mathcal{S}_i) = \begin{cases}
        \mathcal{A}(\mathcal{S}_i) & \text{on }\scriptstyle{ENERGY\_CONSTRAINT} \\
        \phi_0(\mathcal{S}_i) & \text{on }{\scriptstyle{BOUNDARY}}\text{ or }\scriptstyle{EXIT} \\
        0 & \text{otherwise} \\
   \end{cases}.

In the previous equation, :math:`\mathcal{A}` is the activity per unit volume
and solid angle of volume sources, while :math:`\phi_0` is an external flux
associated with surface sources.

.. note::

   In case of an :python:`ENERGY_CONSTRAINT` termination, transport weights have
   units cm |nbsp| MeV\ :sup:`-1`, if :math:`\nu_f < \nu_i` or cm, if
   :math:`\nu_f = \nu_i`, where :math:`\nu_f` (:math:`\nu_i`) is the final
   (initial) photon energy. In other cases, transport weights are unitless.

.. note::

   In the case of a forward Monte Carlo simulation, Goupil's transport weights
   are all equal to one, i.e., Goupil's forward transport is *analogue*.

As an example, consider only volume sources with a uniform activity
:math:`\mathcal{A}_0` per unit volume and solid angle. Then the expected flux
can be written as

.. math::

   \phi(\mathcal{S}_f) = K({\mathcal{S}_f}) \mathcal{A}_0, \quad
   K({\mathcal{S}_f}) \simeq \frac{1}{N} \sum_i{
        \omega\left(\mathcal{S}_f,\mathcal{S_i}\right)
   },

where it should be understood that the sum only runs over events with an
:python:`ENERGY_CONSTRAINT` termination, but the normalisation :math:`N`
considers all simulated events. The quantity :math:`K` can be interpreted as a
sensitivity factor to volume sources. It is estimated as

>>> K = sum(weights[constrained]) / weights.size

This section concludes the current overview of Goupil. For further insight,
please refer to the `examples/
<https://github.com/niess/goupil/tree/master/examples>`_ folder that is
distributed with Goupil's source.
