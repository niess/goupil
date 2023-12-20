use anyhow::{anyhow, bail, Result};
use crate::numerics::float::Float;
use once_cell::sync::OnceCell;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use super::materials::electronic::{ElectronicShell, ElectronicStructure};

pub(crate) mod data;
pub(crate) mod serialise;

use self::data::{ELEMENTS, SHELLS};


//================================================================================================
// Atomic element.
//================================================================================================

#[allow(non_snake_case)]
pub struct AtomicElement {
    pub name: &'static str,
    pub symbol: &'static str,
    pub Z: i32,
    pub A: Float,
}

// Public interface.
#[allow(non_snake_case)]
impl AtomicElement {
    pub fn electrons(&'static self) -> Result<&'static ElectronicStructure> {
        static STRUCTURES: OnceCell<Vec<ElectronicStructure>> = OnceCell::new();
        let structures = STRUCTURES.get_or_init(|| {
            let mut structures = Vec::<ElectronicStructure>::default();
            let mut z: u8 = 0;
            let mut shells = Vec::<ElectronicShell>::default();
            for shell in SHELLS.iter() {
                if (shell.0 != z) && (z > 0) {
                    structures.push(ElectronicStructure(shells.clone()));
                    shells.clear();
                }
                z = shell.0;
                shells.push(ElectronicShell::from(shell))
            };
            structures.push(ElectronicStructure(shells.clone()));
            structures
        });
        let index = usize::try_from(self.Z - 1)?;
        structures.get(index)
           .as_ref()
           .copied()
           .ok_or_else(|| anyhow!(
                "no electronic data for atomic element '{}'",
                self.symbol
            ))
    }

    pub fn from_symbol(symbol: &str) -> Result<&'static Self> {
        static MAP: OnceCell<HashMap<&'static str, &'static AtomicElement>> = OnceCell::new();
        let map = MAP.get_or_init(|| {
            let mut map = HashMap::<&'static str, &'static AtomicElement>::default();
            for element in ELEMENTS.iter() {
                map.insert(element.symbol, element);
            }
            map
        });
        map.get(symbol)
           .copied()
           .ok_or_else(|| anyhow!(
                "no such atomic element '{}'",
                symbol
            ))
    }

    pub fn from_Z(Z: i32) -> Result<&'static Self> {
        match ELEMENTS.get((Z - 1) as usize) {
            None => bail!(
                "bad atomic number (expected a value in [1, {}], found {})",
                ELEMENTS.len(),
                Z,
            ),
            Some(element) => Ok(element),
        }
    }
}

// Private interface.
#[cfg(feature = "python")]
impl AtomicElement {
    pub(crate) fn none() -> &'static Self {
        &NONE_ELEMENT
    }
}

#[cfg(feature = "python")]
static NONE_ELEMENT: AtomicElement = AtomicElement {
    name: "None",
    symbol: "",
    Z: 0,
    A: 0.0,
};

// Beware: traits implementations below are correct provided that all instances are static.
impl PartialEq for AtomicElement {
    fn eq(&self, other: &Self) -> bool {
        (self as *const Self) == (other as *const Self)
    }
}

impl Eq for AtomicElement {}

impl Hash for AtomicElement {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let addr = self as *const Self;
        addr.hash(state);
    }
}
