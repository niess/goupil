use super::AtomicElement;


//================================================================================================
// Atomic element data from https://pdg.lbl.gov/2023/AtomicNuclearProperties.
//================================================================================================

pub static ELEMENTS: [AtomicElement; 118] = [
    AtomicElement {
        name: "Hydrogen",
        symbol: "H",
        Z: 1,
        A: 1.0087,
    },
    AtomicElement {
        name: "Helium",
        symbol: "He",
        Z: 2,
        A: 4.0026022,
    },
    AtomicElement {
        name: "Lithium",
        symbol: "Li",
        Z: 3,
        A: 6.942,
    },
    AtomicElement {
        name: "Beryllium",
        symbol: "Be",
        Z: 4,
        A: 9.01218315,
    },
    AtomicElement {
        name: "Boron",
        symbol: "B",
        Z: 5,
        A: 10.817,
    },
    AtomicElement {
        name: "Carbon",
        symbol: "C",
        Z: 6,
        A: 12.01078,
    },
    AtomicElement {
        name: "Nitrogen",
        symbol: "N",
        Z: 7,
        A: 14.0072,
    },
    AtomicElement {
        name: "Oxygen",
        symbol: "O",
        Z: 8,
        A: 15.9993,
    },
    AtomicElement {
        name: "Fluorine",
        symbol: "F",
        Z: 9,
        A: 18.9984031636,
    },
    AtomicElement {
        name: "Neon",
        symbol: "Ne",
        Z: 10,
        A: 20.17976,
    },
    AtomicElement {
        name: "Sodium",
        symbol: "Na",
        Z: 11,
        A: 22.989769282,
    },
    AtomicElement {
        name: "Magnesium",
        symbol: "Mg",
        Z: 12,
        A: 24.3056,
    },
    AtomicElement {
        name: "Aluminium",
        symbol: "Al",
        Z: 13,
        A: 26.98153857,
    },
    AtomicElement {
        name: "Silicon",
        symbol: "Si",
        Z: 14,
        A: 28.08553,
    },
    AtomicElement {
        name: "Phosphorus",
        symbol: "P",
        Z: 15,
        A: 30.9737619985,
    },
    AtomicElement {
        name: "Sulfur",
        symbol: "S",
        Z: 16,
        A: 32.0655,
    },
    AtomicElement {
        name: "Chlorine",
        symbol: "Cl",
        Z: 17,
        A: 35.4532,
    },
    AtomicElement {
        name: "Argon",
        symbol: "Ar",
        Z: 18,
        A: 39.9481,
    },
    AtomicElement {
        name: "Potassium",
        symbol: "K",
        Z: 19,
        A: 39.09831,
    },
    AtomicElement {
        name: "Calcium",
        symbol: "Ca",
        Z: 20,
        A: 40.0784,
    },
    AtomicElement {
        name: "Scandium",
        symbol: "Sc",
        Z: 21,
        A: 44.9559085,
    },
    AtomicElement {
        name: "Titanium",
        symbol: "Ti",
        Z: 22,
        A: 47.8671,
    },
    AtomicElement {
        name: "Vanadium",
        symbol: "V",
        Z: 23,
        A: 50.94151,
    },
    AtomicElement {
        name: "Chromium",
        symbol: "Cr",
        Z: 24,
        A: 51.99616,
    },
    AtomicElement {
        name: "Manganese",
        symbol: "Mn",
        Z: 25,
        A: 54.9380443,
    },
    AtomicElement {
        name: "Iron",
        symbol: "Fe",
        Z: 26,
        A: 55.8452,
    },
    AtomicElement {
        name: "Cobalt",
        symbol: "Co",
        Z: 27,
        A: 58.9331944,
    },
    AtomicElement {
        name: "Nickel",
        symbol: "Ni",
        Z: 28,
        A: 58.69344,
    },
    AtomicElement {
        name: "Copper",
        symbol: "Cu",
        Z: 29,
        A: 63.5463,
    },
    AtomicElement {
        name: "Zinc",
        symbol: "Zn",
        Z: 30,
        A: 65.382,
    },
    AtomicElement {
        name: "Gallium",
        symbol: "Ga",
        Z: 31,
        A: 69.7231,
    },
    AtomicElement {
        name: "Germanium",
        symbol: "Ge",
        Z: 32,
        A: 72.6301,
    },
    AtomicElement {
        name: "Arsenic",
        symbol: "As",
        Z: 33,
        A: 74.9215956,
    },
    AtomicElement {
        name: "Selenium",
        symbol: "Se",
        Z: 34,
        A: 78.9718,
    },
    AtomicElement {
        name: "Bromine",
        symbol: "Br",
        Z: 35,
        A: 79.9041,
    },
    AtomicElement {
        name: "Krypton",
        symbol: "Kr",
        Z: 36,
        A: 83.7982,
    },
    AtomicElement {
        name: "Rubidium",
        symbol: "Rb",
        Z: 37,
        A: 85.46783,
    },
    AtomicElement {
        name: "Strontium",
        symbol: "Sr",
        Z: 38,
        A: 87.621,
    },
    AtomicElement {
        name: "Yttrium",
        symbol: "Y",
        Z: 39,
        A: 88.905842,
    },
    AtomicElement {
        name: "Zirconium",
        symbol: "Zr",
        Z: 40,
        A: 91.2242,
    },
    AtomicElement {
        name: "Niobium",
        symbol: "Nb",
        Z: 41,
        A: 92.906372,
    },
    AtomicElement {
        name: "Molybdenum",
        symbol: "Mo",
        Z: 42,
        A: 95.951,
    },
    AtomicElement {
        name: "Technetium",
        symbol: "Tc",
        Z: 43,
        A: 97.907223,
    },
    AtomicElement {
        name: "Ruthenium",
        symbol: "Ru",
        Z: 44,
        A: 101.072,
    },
    AtomicElement {
        name: "Rhodium",
        symbol: "Rh",
        Z: 45,
        A: 102.905502,
    },
    AtomicElement {
        name: "Palladium",
        symbol: "Pd",
        Z: 46,
        A: 106.421,
    },
    AtomicElement {
        name: "Silver",
        symbol: "Ag",
        Z: 47,
        A: 107.86822,
    },
    AtomicElement {
        name: "Cadmium",
        symbol: "Cd",
        Z: 48,
        A: 112.4144,
    },
    AtomicElement {
        name: "Indium",
        symbol: "In",
        Z: 49,
        A: 114.8183,
    },
    AtomicElement {
        name: "Tin",
        symbol: "Sn",
        Z: 50,
        A: 118.7107,
    },
    AtomicElement {
        name: "Antimony",
        symbol: "Sb",
        Z: 51,
        A: 121.7601,
    },
    AtomicElement {
        name: "Tellurium",
        symbol: "Te",
        Z: 52,
        A: 127.603,
    },
    AtomicElement {
        name: "Iodine",
        symbol: "I",
        Z: 53,
        A: 126.904473,
    },
    AtomicElement {
        name: "Xenon",
        symbol: "Xe",
        Z: 54,
        A: 131.2936,
    },
    AtomicElement {
        name: "Caesium",
        symbol: "Cs",
        Z: 55,
        A: 132.905451966,
    },
    AtomicElement {
        name: "Barium",
        symbol: "Ba",
        Z: 56,
        A: 137.3277,
    },
    AtomicElement {
        name: "Lanthanum",
        symbol: "La",
        Z: 57,
        A: 138.905477,
    },
    AtomicElement {
        name: "Cerium",
        symbol: "Ce",
        Z: 58,
        A: 140.1161,
    },
    AtomicElement {
        name: "Praseodymium",
        symbol: "Pr",
        Z: 59,
        A: 140.907662,
    },
    AtomicElement {
        name: "Neodymium",
        symbol: "Nd",
        Z: 60,
        A: 144.2423,
    },
    AtomicElement {
        name: "Promethium",
        symbol: "Pm",
        Z: 61,
        A: 144.912753,
    },
    AtomicElement {
        name: "Samarium",
        symbol: "Sm",
        Z: 62,
        A: 150.362,
    },
    AtomicElement {
        name: "Europium",
        symbol: "Eu",
        Z: 63,
        A: 151.9641,
    },
    AtomicElement {
        name: "Gadolinium",
        symbol: "Gd",
        Z: 64,
        A: 157.253,
    },
    AtomicElement {
        name: "Terbium",
        symbol: "Tb",
        Z: 65,
        A: 158.925352,
    },
    AtomicElement {
        name: "Dysprosium",
        symbol: "Dy",
        Z: 66,
        A: 162.5001,
    },
    AtomicElement {
        name: "Holmium",
        symbol: "Ho",
        Z: 67,
        A: 164.930332,
    },
    AtomicElement {
        name: "Erbium",
        symbol: "Er",
        Z: 68,
        A: 167.2593,
    },
    AtomicElement {
        name: "Thulium",
        symbol: "Tm",
        Z: 69,
        A: 168.934222,
    },
    AtomicElement {
        name: "Ytterbium",
        symbol: "Yb",
        Z: 70,
        A: 173.0545,
    },
    AtomicElement {
        name: "Lutetium",
        symbol: "Lu",
        Z: 71,
        A: 174.96681,
    },
    AtomicElement {
        name: "Hafnium",
        symbol: "Hf",
        Z: 72,
        A: 178.492,
    },
    AtomicElement {
        name: "Tantalum",
        symbol: "Ta",
        Z: 73,
        A: 180.947882,
    },
    AtomicElement {
        name: "Tungsten",
        symbol: "W",
        Z: 74,
        A: 183.841,
    },
    AtomicElement {
        name: "Rhenium",
        symbol: "Re",
        Z: 75,
        A: 186.2071,
    },
    AtomicElement {
        name: "Osmium",
        symbol: "Os",
        Z: 76,
        A: 190.233,
    },
    AtomicElement {
        name: "Iridium",
        symbol: "Ir",
        Z: 77,
        A: 192.2173,
    },
    AtomicElement {
        name: "Platinum",
        symbol: "Pt",
        Z: 78,
        A: 195.0849,
    },
    AtomicElement {
        name: "Gold",
        symbol: "Au",
        Z: 79,
        A: 196.9665695,
    },
    AtomicElement {
        name: "Mercury",
        symbol: "Hg",
        Z: 80,
        A: 200.5922,
    },
    AtomicElement {
        name: "Thallium",
        symbol: "Tl",
        Z: 81,
        A: 204.382,
    },
    AtomicElement {
        name: "Lead",
        symbol: "Pb",
        Z: 82,
        A: 207.21,
    },
    AtomicElement {
        name: "Bismuth",
        symbol: "Bi",
        Z: 83,
        A: 208.980401,
    },
    AtomicElement {
        name: "Polonium",
        symbol: "Po",
        Z: 84,
        A: 208.982432,
    },
    AtomicElement {
        name: "Astatine",
        symbol: "At",
        Z: 85,
        A: 209.987156,
    },
    AtomicElement {
        name: "Radon",
        symbol: "Rn",
        Z: 86,
        A: 222.017582,
    },
    AtomicElement {
        name: "Francium",
        symbol: "Fr",
        Z: 87,
        A: 223.019742,
    },
    AtomicElement {
        name: "Radium",
        symbol: "Ra",
        Z: 88,
        A: 226.025412,
    },
    AtomicElement {
        name: "Actinium",
        symbol: "Ac",
        Z: 89,
        A: 227.027752,
    },
    AtomicElement {
        name: "Thorium",
        symbol: "Th",
        Z: 90,
        A: 232.03774,
    },
    AtomicElement {
        name: "Protactinium",
        symbol: "Pa",
        Z: 91,
        A: 231.035882,
    },
    AtomicElement {
        name: "Uranium",
        symbol: "U",
        Z: 92,
        A: 238.028913,
    },
    AtomicElement {
        name: "Neptunium",
        symbol: "Np",
        Z: 93,
        A: 237.048172,
    },
    AtomicElement {
        name: "Plutonium",
        symbol: "Pu",
        Z: 94,
        A: 244.064204,
    },
    AtomicElement {
        name: "Americium",
        symbol: "Am",
        Z: 95,
        A: 243.061382,
    },
    AtomicElement {
        name: "Curium",
        symbol: "Cm",
        Z: 96,
        A: 247.070353,
    },
    AtomicElement {
        name: "Berkelium",
        symbol: "Bk",
        Z: 97,
        A: 247.070314,
    },
    AtomicElement {
        name: "Californium",
        symbol: "Cf",
        Z: 98,
        A: 251.079593,
    },
    AtomicElement {
        name: "Einsteinium",
        symbol: "Es",
        Z: 99,
        A: 252.08304,
    },
    AtomicElement {
        name: "Fermium",
        symbol: "Fm",
        Z: 100,
        A: 257.095105,
    },
    AtomicElement {
        name: "Mendelevium",
        symbol: "Md",
        Z: 101,
        A: 258.098433,
    },
    AtomicElement {
        name: "Nobelium",
        symbol: "No",
        Z: 102,
        A: 259.10107,
    },
    AtomicElement {
        name: "Lawrencium",
        symbol: "Lr",
        Z: 103,
        A: 262.109612,
    },
    AtomicElement {
        name: "Rutherfordium",
        symbol: "Rf",
        Z: 104,
        A: 267.121794,
    },
    AtomicElement {
        name: "Dubnium",
        symbol: "Db",
        Z: 105,
        A: 268.125674,
    },
    AtomicElement {
        name: "Seaborgium",
        symbol: "Sg",
        Z: 106,
        A: 269.128635,
    },
    AtomicElement {
        name: "Bohrium",
        symbol: "Bh",
        Z: 107,
        A: 270.133364,
    },
    AtomicElement {
        name: "Hassium",
        symbol: "Hs",
        Z: 108,
        A: 269.1337513,
    },
    AtomicElement {
        name: "Meitnerium",
        symbol: "Mt",
        Z: 109,
        A: 278.15637,
    },
    AtomicElement {
        name: "Darmstadtium",
        symbol: "Ds",
        Z: 110,
        A: 281.16456,
    },
    AtomicElement {
        name: "Roentgenium",
        symbol: "Rg",
        Z: 111,
        A: 282.16917,
    },
    AtomicElement {
        name: "Copernicium",
        symbol: "Cn",
        Z: 112,
        A: 285.177125,
    },
    AtomicElement {
        name: "Nihonium",
        symbol: "Nh",
        Z: 113,
        A: 286.182216,
    },
    AtomicElement {
        name: "Flerovium",
        symbol: "Fl",
        Z: 114,
        A: 289.190425,
    },
    AtomicElement {
        name: "Moscovium",
        symbol: "Mc",
        Z: 115,
        A: 289.193636,
    },
    AtomicElement {
        name: "Livermorium",
        symbol: "Lv",
        Z: 116,
        A: 293.20456,
    },
    AtomicElement {
        name: "Tennessine",
        symbol: "Ts",
        Z: 117,
        A: 294.21057,
    },
    AtomicElement {
        name: "Oganesson",
        symbol: "Og",
        Z: 118,
        A: 294.213928,
    },
];
