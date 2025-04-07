#! /usr/bin/env python3
import argparse
import gzip
import numpy
from pathlib import Path
import pickle
import sys
import time
try:
    import tomli as toml
except ImportError:
    import tomllib as toml
from types import SimpleNamespace

import calzone
import goupil


ARGS = None

PREFIX = Path(__file__).parent


class Chrono:
    """Measure CPU-time consumption."""

    def __init__(self):
        self._t0 = None

    def __enter__(self):
        self._t0 = time.process_time_ns()
        return self

    def __exit__(self, type, value, traceback):
        self.duration = (time.process_time_ns() - self._t0) * 1E-09


def load_sources(path):
    """Load the radio-isotope(s) emission lines."""

    with open(path, "rb") as f:
        sources = toml.load(f)

    energies, intensities = [], []
    for source in sources.values():
        energies += source["energies"]
        intensities += source["intensities"]
    energies, intensities = map(numpy.array, (energies, intensities))

    return SimpleNamespace(
        energies = energies,
        intensities = intensities
    )


def log(msg):
    """Logging wrapper."""

    t = time.strftime("%a %d %b %Y %H:%M:%S")
    sys.stdout.write(f"[{t}] {msg}\n")
    sys.stdout.flush()


class Runner:
    """Monte Carlo simulation runner."""

    def __init__(self):
        """Configure the Monte Carlo simulation."""

        # Create the Geant4 geometry (using Calzone).
        if ARGS.radius:
            self.geometry = calzone.GeometryBuilder(ARGS.geometry)           \
                .modify("Environment.Source", shape={"sphere": ARGS.radius}) \
                .build()
        else:
            self.geometry = calzone.Geometry(ARGS.geometry)

        volume = self.geometry.find(ARGS.generation_volume)

        # Create the forward transport engine.
        self.forward_engine = calzone.Simulation(self.geometry)
        self.forward_engine.physics = ARGS.physics
        self.forward_engine.physics.default_cut = ARGS.default_cut

        # Create the backward transport engine.
        self.backward_engine = goupil.TransportEngine(
            geometry = self.geometry.export()
        )
        self.backward_engine.energy_min = ARGS.energy_min
        self.backward_engine.mode = "Backward"

        self.backward_engine.boundary = \
            self.backward_engine.geometry.sector_index(volume.path)

        # Compile materials tables for the backward simulation.
        self.backward_engine.compile()

        # Create the emission spectrum.
        self.source = load_sources(ARGS.sources)
        self.spectrum = goupil.DiscreteSpectrum(
            energies = self.source.energies,
            intensities = self.source.intensities,
            background = ARGS.background_fraction,
            energy_min = ARGS.energy_min
        )


    def run(self, *, data=None):
        """Run a mixed Monte Carlo simulation."""

        # Reset the pseudo-random streams.
        self.forward_engine.random.seed = None
        self.backward_engine.random.seed = None

        # Generate primary particles entering the detector volume.
        volume = self.geometry.find(ARGS.generation_volume)
        initial_states = self.forward_engine.particles(weight=True) \
            .pid(22)                                                \
            .on(volume, direction="ingoing")                        \
            .generate(ARGS.number_of_events)

        # Backward-generate the energies.
        source_energies = self.spectrum.sample(initial_states,
                                               engine=self.backward_engine)

        # Normalise weights to the radio-isotope activity.
        total_intensity = sum(self.source.intensities) / 100.0
        initial_states["weight"] *= total_intensity

        # Run the forward simulation.
        with Chrono() as chrono:
            forward_result = self.forward_engine.run(initial_states)

        try:
            volume = self.geometry.find("Scintillator")
            deposits = forward_result.deposits[volume.path]
        except KeyError:
            deposits = numpy.empty(0)
            weights = numpy.empty(0)
            dt0, dt1 = chrono.duration, 0.0
            log(f"collected no deposit in {dt0:.1f} s")
        else:
            sel = deposits["value"] > ARGS.energy_min
            deposits = deposits[sel]

            dt0 = chrono.duration
            log(f"collected {deposits.size} deposits in {dt0:.1f} s")

            # Run the backward simulation.
            index = deposits["event"]
            deposits = deposits["value"]
            initial_states = initial_states[index]
            source_energies = source_energies[index]

            source_states = initial_states.copy()
            EPSILON = 1E-04 # cm
            source_states["position"] -= EPSILON * source_states["direction"]

            with Chrono() as chrono:
                status = self.backward_engine.transport(
                    source_states,
                    source_energies = source_energies,
                )

            # Check sources locations.
            locations = self.backward_engine.geometry.locate(source_states)
            volume = self.geometry.find("Source")
            index = self.backward_engine.geometry.sector_index(volume.path)
            is_source = locations == index

            sel = (status == goupil.TransportStatus.ENERGY_CONSTRAINT) & \
                  is_source
            deposits = deposits[sel]
            source_states = source_states[sel]

            dt1 = chrono.duration
            log(f"backward sampled {source_states.size} sources in {(dt1):.1f} s")

            weights = source_states["weight"]

        # Save the simulation results.
        outdir = ARGS.output_directory
        outdir.mkdir(parents=True, exist_ok=True)
        seed0 = self.forward_engine.random.seed
        seed1 = self.backward_engine.random.seed
        filename = f"{seed0:0X}-{seed1:0X}.pkl.gz"
        data = {
            "deposits": deposits,
            "weights": weights,
            "cpu": SimpleNamespace(forward = dt0, backward = dt1),
            "seeds": SimpleNamespace(forward = seed0, backward = seed1),
            "versions": SimpleNamespace(
                calzone = calzone.VERSION,
                geant4 = calzone.GEANT4_VERSION,
                goupil = goupil.VERSION
            )
        }
        data.update(ARGS.__dict__)
        data = SimpleNamespace(**data)

        out_path = outdir / filename
        with gzip.open(out_path, "wb") as f:
            pickle.dump(data, f)

        log(f"dumped events to {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = Runner.run.__doc__)
    parser.add_argument("-b", "--background-fraction",
        help = "fraction of background events",
        type = float,
        default = 0.75
    )
    parser.add_argument("-d", "--default-cut",
        help = "Geant4 default physics cut, in cm",
        type = float,
        default = 0.01
    )
    parser.add_argument("-e", "--energy-min",
        help = "minimum energy to simulate / record, in MeV",
        type = float,
        default = 1E-03
    )
    parser.add_argument("-g", "--geometry",
        help = "path to the geometry file",
        type = Path,
        default = PREFIX / "geometry.toml"
    )
    parser.add_argument("-G", "--generation-volume",
        help = "name of the generation volume",
        default = "Detector"
    )
    parser.add_argument("-l", "--loop",
        help = "loop Monte Carlo runs",
        nargs = "?",
        type = int,
        default = 1,
        const = 0,
    )
    parser.add_argument("-n", "--number-of-events",
        help = "number of Monte Carlo events (per run)",
        type = int,
        default = 1000000
    )
    parser.add_argument("-o", "--output-directory",
        help = "Monte Carlo events output directory",
        type = Path,
        default = PREFIX / "data/mixed"
    )
    parser.add_argument("-p", "--physics",
        help = "Geant4 EM physics",
        default = "penelope"
    )
    parser.add_argument("-r", "--radius",
        help = "source radius, in cm",
        type = float
    )
    parser.add_argument("-s", "--sources",
        help = "path to the sources file",
        type = Path,
        default = PREFIX / "sources.toml"
    )

    ARGS = parser.parse_args()

    with Chrono() as chrono:
        runner = Runner()

    dt = chrono.duration
    log(f"runner configured in {dt:.1f} s")

    if ARGS.loop <= 0:
        while True:
            runner.run()
    else:
        for _ in range(ARGS.loop):
            runner.run()
