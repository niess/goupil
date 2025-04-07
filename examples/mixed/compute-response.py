#! /usr/bin/env python3
import argparse
import glob
import gzip
import numpy
from pathlib import Path
import pickle
from types import SimpleNamespace

try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda x: x


ARGS = None

CS_LINE = 0.662 # MeV


def compute():
    """Compute the detector energy response."""

    x = numpy.linspace(
        ARGS.energy_min,
        ARGS.energy_max,
        ARGS.nodes
    )
    y = numpy.zeros((ARGS.nodes, 2))
    events, collected, rate, cpu = 0, 0, 0.0, 0.0

    # Source activity, per unit volume and solid angle.
    m3 = 1E+06 # m3 to cm3
    activity = ARGS.activity / (4.0 * numpy.pi * m3)

    for path in tqdm(glob.glob(str(ARGS.path / "*.pkl.gz"))):
        with gzip.open(path, "rb") as f:
            data = pickle.load(f)

        weights = data.weights * activity

        sigma = ARGS.resolution * numpy.sqrt(data.deposits * CS_LINE) / 2.355
        for i, xi in enumerate(x):
            u = (xi - data.deposits) / sigma
            w = weights * numpy.exp(-0.5 * u**2) / sigma
            y[i, 0] += sum(w) / numpy.sqrt(2.0 * numpy.pi)
            y[i, 1] += sum(w**2) / (2.0 * numpy.pi)

        events += data.number_of_events
        collected += data.deposits.size
        rate += sum(weights)
        cpu += data.cpu.forward
        try:
            cpu += data.cpu.backward
        except AttributeError:
            pass

    s1 = y[:, 0] / events
    s2 = y[:, 1] / events
    y[:, 0] = s1
    y[:, 1] = numpy.sqrt(numpy.maximum(s2 - s1**2, 0.0) / events)
    rate /= events

    response = SimpleNamespace(
        x = x,
        y = y,
        n = events,
        collected = collected,
        rate = rate,
        cpu = cpu
    )

    if ARGS.output_file is None:
        ARGS.output_file = ARGS.path / "response.pkl"

    ARGS.output_file.parent.mkdir(parents=True, exist_ok=True)
    with ARGS.output_file.open("wb") as f:
        pickle.dump(response, f)

    print(f"simulated {events} events")
    print(f"collected {collected} events")
    print(f"total cpu: {cpu:.3E} s")

    if ARGS.plot:
        import matplotlib.pyplot as plt
        plt.fill_between(
            response.x,
            response.y[:,0] - response.y[:,1],
            response.y[:,0] + response.y[:,1],
            color = "k",
            alpha = 0.25
        )
        plt.plot(response.x, response.y[:,0], "k-")
        plt.yscale("log")
        plt.xlabel("reconstructed energy (MeV)")
        plt.ylabel("differential rate (Hz/MeV)")
        plt.ylim(2E-06, 5E-02)
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = compute.__doc__)
    parser.add_argument("path",
        help = "Monte Carlo events directory",
        type = Path,
        nargs = "?",
        default = Path("data/mixed")
    )
    parser.add_argument("-a", "--activity",
        help = "source volume activity, in Bq / m^3",
        type = float,
        default = 1.0
    )
    parser.add_argument("--energy-max",
        help = "maximum energy, in MeV",
        type = float,
        default = 2.5
    )
    parser.add_argument("--energy-min",
        help = "minimum energy, in MeV",
        type = float,
        default = 0.0
    )
    parser.add_argument("-n", "--nodes",
        help = "number of energy nodes",
        type = int,
        default = 251
    )
    parser.add_argument("-o", "--output-file",
        help = "energy response output file",
        type = Path,
    )
    parser.add_argument("-p", "--plot",
        help = "plot the computed response",
        action = "store_true",
    )
    parser.add_argument("-r", "--resolution",
        help = f"energy resolution (at {CS_LINE} MeV)",
        type = float,
        default = 6.7E-02
    )

    ARGS = parser.parse_args()
    compute()
