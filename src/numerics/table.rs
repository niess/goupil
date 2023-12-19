use anyhow::{bail, Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use super::float::Float;


//================================================================================================
// A utility trait for abstracting a 1D tabulation of a physical property.
//================================================================================================

pub(crate) trait Table1D {
    fn interpolate(&self, x: Float) -> Float;
    fn len(&self) -> usize;
    fn x(&self, i: usize) -> Float;
    fn y(&self, i: usize) -> Float;

    // Merge (weighted) tabulated data taking care of preserving all x-values.
    fn merge(tables: &[(Float, &Self)]) -> Option<(Vec<Float>, Vec<Float>)> {
        // Reduce tables.
        let mut data: Vec<_> = tables
            .iter()
            .filter(|(weight, _)| { *weight > 0.0 })
            .collect();
        if data.is_empty() {
            return None
        }

        // Loop over x values and build the new data set.
        let n = data.len();
        let mut index = vec![0_usize; n];
        let mut x = Vec::<Float>::default();
        let mut y = Vec::<Float>::default();
        while !data.is_empty() {
            // Find smallest energy.
            let xmin = data
                .iter()
                .enumerate()
                .map(|(i, (_, table))| table.x(index[i]))
                .reduce(Float::min)
                .unwrap();

            // Update merged data.
            x.push(xmin);
            let y_value = data
                .iter()
                .enumerate()
                .map(|(i, (weight, table))| {
                    let ii = index[i];
                    let yi = if table.x(ii) == xmin {
                        table.y(ii)
                    } else {
                        table.interpolate(xmin)
                    };
                    weight * yi
                })
                .sum();
            y.push(y_value);

            // Update indices.
            let mut i = 0;
            while i < index.len() {
                if Table1D::x(data[i].1, index[i]) <= xmin {
                    // Consume table entry.
                    let new_index = index[i] + 1;
                    if new_index >= data[i].1.len() {
                        // Table is fully consumed. Let us remove it.
                        index.remove(i);
                        data.remove(i);
                        continue;
                    } else {
                        index[i] = new_index;
                    }
                }
                i += 1;
            }
        }
        Some((x, y))
    }
}


//================================================================================================
// Data loader for 1D tables.
//================================================================================================

pub(crate) struct Data1D {
    pub x: Vec<Float>,
    pub y: Vec<Float>,
}

impl Data1D {
    pub fn from_file<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let mut x = Vec::<Float>::default();
        let mut y = Vec::<Float>::default();
        let path: &Path = path.as_ref();
        let file = File::open(path)
            .with_context(|| format!("could not open {}", path.display()))?;
        for (lineno, line) in BufReader::new(file).lines().enumerate() {
            if let Ok(line) = line {
                if line.starts_with("#") || (line.len() == 0) { continue }

                let words: Vec<&str> = line.split_whitespace().collect();
                let n_words = words.len();
                const EXPECTED: usize = 2;
                if n_words != EXPECTED {
                    bail!(
                        "{}:{}: bad format (expected {} items, found {})",
                        path.display(),
                        lineno,
                        EXPECTED,
                        n_words,
                    );
                }
                let xi = words[0]
                    .parse::<Float>()
                    .with_context(|| format!(
                        "{}:{}: could not parse first column",
                        path.display(),
                        lineno,
                    ))?;
                let yi = words[1]
                    .parse::<Float>()
                    .with_context(|| format!(
                        "{}:{}: could not parse second column",
                        path.display(),
                        lineno,
                    ))?;
                x.push(xi);
                y.push(yi);
            }
        }

        Ok(Self { x, y })
    }
}


//================================================================================================
// Generic unpacking for Data1D like types.
//================================================================================================

pub(crate) trait FromFile {
    fn from_file<P>(path: P) -> Result<Self>
    where
        Self: Sized,
        P: AsRef<Path>;
}

impl<T> FromFile for T
where
    T: From<Data1D>,
{
    fn from_file<P>(path: P) -> Result<T>
    where
        P: AsRef<Path>,
    {
        let data = Data1D::from_file(path)?;
        Ok(data.into())
    }
}
