# Goupil <img src="https://github.com/niess/goupil/blob/master/docs/goupil.svg" width="30"> [![Documentation Status][RTD_BADGE]][RTD]
(_**G**amma transp**O**rt **U**tility, a**P**proximate but revers**I**b**L**e_)

Goupil is a Monte Carlo transport engine designed for low-energy gamma-rays
(0.1-3 MeV), typically emitted from radioactive isotopes. By using a backward
transport algorithm, and making a few approximations, Goupil can significantly
accelerate Monte Carlo simulations, resulting in computation time savings of up
to 5 orders of magnitude in the air.

## Documentation

Documentation for Goupil can be found online at [Read the Docs][RTD].

## License

The Goupil library is  under the **GNU LGPLv3** license (see the provided
[COPYING](COPYING) and [COPYING.LESSER](COPYING.LESSER) files), except for
[examples][EXAMPLES] and [interfaces][INTERFACES] which are MIT-licensed.


[EXAMPLES]: https://github.com/niess/goupil/tree/master/examples
[INTERFACES]: https://github.com/niess/goupil/tree/master/src/interfaces
[RTD]: https://goupil.readthedocs.io/en/latest/
[RTD_BADGE]: https://readthedocs.org/projects/goupil/badge/?version=latest
