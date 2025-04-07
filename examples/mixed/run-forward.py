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

        # Create the transport engine.
        self.engine = calzone.Simulation(self.geometry)
        self.engine.physics = ARGS.physics
        self.engine.physics.default_cut = ARGS.default_cut

        # Create the emission spectrum.
        self.source = load_sources(ARGS.sources)
        self.spectrum = list(zip(
            self.source.energies,
            self.source.intensities
        ))


    def run(self, *, data=None):
        """Run a Geant4 Monte Carlo simulation (using Calzone)."""

        # Reset the pseudo-random stream.
        self.engine.random.seed = None

        # Generate primary particles over the source volume.
        volume = self.geometry.find("Source")
        primaries = self.engine.particles(weight=True) \
            .pid(22)                                   \
            .inside(volume)                            \
            .solid_angle()                             \
            .spectrum(self.spectrum, weight=False)     \
            .generate(ARGS.number_of_events)

        # Normalise weights to the radio-isotope(s) activity.
        total_intensity = sum(self.source.intensities) / 100.0
        primaries["weight"] *= total_intensity

        # Run the forward simulation.
        with Chrono() as chrono:
            result = self.engine.run(primaries)

        try:
            volume = self.geometry.find("Scintillator")
            deposits = result.deposits[volume.path]
        except KeyError:
            deposits = numpy.empty(0)
            weights = numpy.empty(0)
            dt = chrono.duration
            log(f"collected no deposit in {dt:.1f} s")
        else:
            sel = deposits["value"] > ARGS.energy_min
            deposits = deposits[sel]

            dt = chrono.duration
            log(f"collected {deposits.size} deposits in {dt:.1f} s")

            # Run the backward simulation.
            index = deposits["event"]
            deposits = deposits["value"]
            primaries = primaries[index]
            weights = primaries["weight"]

        # Save the simulation results.
        outdir = ARGS.output_directory
        outdir.mkdir(parents=True, exist_ok=True)
        seed = self.engine.random.seed
        filename = f"{seed:0X}.pkl.gz"
        data = {
            "deposits": deposits,
            "weights": weights,
            "cpu": SimpleNamespace(forward=dt),
            "seeds": SimpleNamespace(forward=seed),
            "versions": SimpleNamespace(
                calzone = calzone.VERSION,
                geant4 = calzone.GEANT4_VERSION
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
        default = PREFIX / "data/forward"
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
