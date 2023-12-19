use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde::de::{Error, Unexpected};
use super::{AtomicElement, ELEMENTS};


// ===============================================================================================
//
// Serialization and deserialization for static reference to an AtomicElement.
//
// ===============================================================================================

impl<'de> Deserialize<'de> for &'static AtomicElement {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let z: i32 = Deserialize::deserialize(deserializer)?;
        AtomicElement::from_Z(z)
            .map_err(|_| D::Error::invalid_value(
                Unexpected::Signed(z.into()),
                &format!("a value in [1, {}]", ELEMENTS.len()).as_str(),
            ))
    }
}

impl Serialize for &'static AtomicElement {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_i32(self.Z)
    }
}
