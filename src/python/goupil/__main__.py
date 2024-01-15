import argparse
import goupil


def main():
    parser = argparse.ArgumentParser(
        prog = "python3 -m goupil",
        description = "Configuration utility for Goupil.",
        epilog = "Copyright (C) Université Clermont Auvergne, CNRS/IN2P3, LPC"
    )
    parser.add_argument("-p", "--prefix",
        help = "Goupil installation prefix.",
        action = "store_true"
    )
    parser.add_argument("-v", "--version",
        help = "Goupil version.",
        action = "store_true"
    )

    args = parser.parse_args()

    result = []
    if args.prefix:
        result.append(goupil.PREFIX)
    if args.version:
        result.append(goupil.VERSION)

    if result:
        print(" ".join(result))


if __name__ == "__main__":
    main()
