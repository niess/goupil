use std::fmt::Display;


// ===============================================================================================
// Public API.
// ===============================================================================================

pub mod numerics;
pub mod physics;
pub mod transport;


// ===============================================================================================
// Python interface.
// ===============================================================================================

#[cfg(feature = "python")]
mod python;


// ===============================================================================================
// Utilities, e.g. for error messages.
// ===============================================================================================

fn pretty_enumerate<T>(elements: &[T]) -> String
where
    T: Display,
{
    let n = elements.len();
    match n {
        0 => unreachable!(),
        1 => format!("{}", elements[0]),
        2 => format!("{} or {}", elements[0], elements[1]),
        _ => {
            let elements: Vec<_> = elements
                .iter()
                .map(|e| format!("{}", e))
                .collect();
            format!(
                "{} or {}",
                elements[0..(n - 1)].join(", "),
                elements[n - 1],
            )
        },
    }
}
