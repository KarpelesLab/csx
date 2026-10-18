//! Typed identifiers for property settings stored as numbers in files.
//!
//! Every enum keeps unrecognised numeric ids in an `Other(i32)` variant so
//! files round-trip losslessly. [`from_id`](ExcitationType::from_id) always
//! maps known ids to their named variant; an `Other` holding a known id
//! compares unequal to the named variant, so prefer `from_id` to build
//! values from numbers.

use core::fmt;
use core::str::FromStr;

use crate::error::Error;

macro_rules! id_enum {
    (
        $(#[$m:meta])*
        $name:ident, default = $def:ident {
            $( $(#[$vm:meta])* $v:ident = $id:literal => $s:literal ),* $(,)?
        }
    ) => {
        $(#[$m])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        #[non_exhaustive]
        pub enum $name {
            $( $(#[$vm])* $v, )*
            /// An id without a named variant, preserved verbatim.
            Other(i32),
        }

        impl $name {
            /// Every named variant.
            pub const ALL: &'static [$name] = &[$($name::$v),*];

            /// Variant for a numeric id (unknown ids become `Other`).
            pub const fn from_id(id: i32) -> Self {
                match id {
                    $( $id => $name::$v, )*
                    o => $name::Other(o),
                }
            }

            /// Numeric id as used in files.
            pub const fn id(self) -> i32 {
                match self {
                    $( $name::$v => $id, )*
                    $name::Other(o) => o,
                }
            }

            /// Symbolic name, `None` for `Other`.
            pub const fn name(self) -> Option<&'static str> {
                match self {
                    $( $name::$v => Some($s), )*
                    $name::Other(_) => None,
                }
            }
        }

        impl Default for $name {
            fn default() -> Self {
                $name::$def
            }
        }

        impl From<i32> for $name {
            fn from(id: i32) -> Self {
                $name::from_id(id)
            }
        }

        impl From<$name> for i32 {
            fn from(v: $name) -> i32 {
                v.id()
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                match self.name() {
                    Some(n) => f.write_str(n),
                    None => write!(f, "{}", self.id()),
                }
            }
        }

        impl FromStr for $name {
            type Err = Error;

            /// Accepts a symbolic name (case-insensitive, `-`/`_` ignored)
            /// or a numeric id (also float-formatted, e.g. `1.000000e+00`).
            ///
            /// # Errors
            ///
            /// [`Error::Format`] if `s` is neither a known name nor an
            /// integer id.
            fn from_str(s: &str) -> Result<Self, Error> {
                let key = normalize(s);
                $( if normalize($s) == key { return Ok($name::$v); } )*
                parse_id(s).map($name::from_id).ok_or_else(|| {
                    Error::Format(format!(concat!("invalid ", stringify!($name), " {:?}"), s))
                })
            }
        }
    };
}

fn normalize(s: &str) -> String {
    s.chars()
        .filter(|c| *c != '_' && *c != '-' && !c.is_whitespace())
        .flat_map(char::to_lowercase)
        .collect()
}

/// Parses an integer id, accepting float formatting (`2.000000e+00`).
pub(crate) fn parse_id(s: &str) -> Option<i32> {
    let t = s.trim();
    if let Ok(v) = t.parse::<i32>() {
        return Some(v);
    }
    let v = t.parse::<f64>().ok()?;
    (v.is_finite() && v.fract() == 0.0 && v.abs() <= i32::MAX as f64).then_some(v as i32)
}

id_enum! {
    /// Kind of field excitation (`Type` attribute of `<Excitation>`).
    ExcitationType, default = EFieldSoft {
        /// Soft electric field source.
        EFieldSoft = 0 => "EFieldSoft",
        /// Hard electric field source.
        EFieldHard = 1 => "EFieldHard",
        /// Soft magnetic field source.
        HFieldSoft = 2 => "HFieldSoft",
        /// Hard magnetic field source.
        HFieldHard = 3 => "HFieldHard",
        /// Plane wave (total-field / scattered-field).
        PlaneWave = 10 => "PlaneWave",
    }
}

impl ExcitationType {
    /// True for electric field sources (types 0 and 1).
    pub const fn is_electric(self) -> bool {
        matches!(
            self,
            ExcitationType::EFieldSoft | ExcitationType::EFieldHard
        )
    }

    /// True for magnetic field sources (types 2 and 3).
    pub const fn is_magnetic(self) -> bool {
        matches!(
            self,
            ExcitationType::HFieldSoft | ExcitationType::HFieldHard
        )
    }
}

id_enum! {
    /// Kind of probe (`Type` attribute of `<ProbeBox>` / `<DumpBox>`).
    ProbeType, default = Voltage {
        /// Voltage, integrated along a line.
        Voltage = 0 => "Voltage",
        /// Current, integrated around a face.
        Current = 1 => "Current",
        /// Electric field at a point.
        EField = 2 => "EField",
        /// Magnetic field at a point.
        HField = 3 => "HField",
        /// Waveguide voltage (mode matching).
        WaveguideVoltage = 10 => "WaveguideVoltage",
        /// Waveguide current (mode matching).
        WaveguideCurrent = 11 => "WaveguideCurrent",
    }
}

id_enum! {
    /// Quantity written by a field dump (`DumpType`).
    DumpType, default = EField {
        /// Electric field (time domain).
        EField = 0 => "EField",
        /// Magnetic field (time domain).
        HField = 1 => "HField",
        /// Electric current density (time domain).
        Current = 2 => "Current",
        /// Total current density rot(H) (time domain).
        TotalCurrentDensity = 3 => "TotalCurrentDensity",
        /// Electric flux density D (time domain).
        DField = 4 => "DField",
        /// Magnetic flux density B (time domain).
        BField = 5 => "BField",
        /// Electric field (frequency domain).
        EFieldFd = 10 => "EFieldFD",
        /// Magnetic field (frequency domain).
        HFieldFd = 11 => "HFieldFD",
        /// Electric current density (frequency domain).
        CurrentFd = 12 => "CurrentFD",
        /// Total current density rot(H) (frequency domain).
        TotalCurrentDensityFd = 13 => "TotalCurrentDensityFD",
        /// Electric flux density D (frequency domain).
        DFieldFd = 14 => "DFieldFD",
        /// Magnetic flux density B (frequency domain).
        BFieldFd = 15 => "BFieldFD",
        /// Local specific absorption rate.
        SarLocal = 20 => "SARLocal",
        /// SAR averaged over 1 g.
        Sar1g = 21 => "SAR1g",
        /// SAR averaged over 10 g.
        Sar10g = 22 => "SAR10g",
        /// Raw SAR data.
        SarRaw = 29 => "SARRaw",
    }
}

impl DumpType {
    /// True for frequency-domain dumps (including SAR), which need
    /// frequencies.
    pub const fn is_frequency_domain(self) -> bool {
        matches!(self.id(), 10..=15 | 20..=29)
    }
}

id_enum! {
    /// Interpolation applied by a field dump (`DumpMode`).
    DumpMode, default = Node {
        /// Raw values, no interpolation.
        NoInterpolation = 0 => "NoInterpolation",
        /// Interpolated to the mesh nodes.
        Node = 1 => "Node",
        /// Interpolated to the cell centres.
        Cell = 2 => "Cell",
    }
}

id_enum! {
    /// File format written by a field dump (`FileType`).
    DumpFileType, default = Vtk {
        /// VTK files.
        Vtk = 0 => "VTK",
        /// HDF5 files.
        Hdf5 = 1 => "HDF5",
    }
}

id_enum! {
    /// Connection of the R, L and C parts of a lumped element (`LEtype`).
    LumpedElementType, default = Parallel {
        /// R, L and C in parallel.
        Parallel = 0 => "Parallel",
        /// R, L and C in series.
        Series = 1 => "Series",
    }
}

id_enum! {
    /// Absorbing boundary condition algorithm (`AbsorbingBoundaryType`).
    /// The numeric ids are not confirmed by documentation.
    AbcType, default = Mur1st {
        /// First order Mur with one phase velocity.
        Mur1st = 0 => "MUR_1ST",
        /// First order Mur with one phase velocity, super-absorption variant.
        Mur1stSa = 1 => "MUR_1ST_SA",
    }
}

/// A coordinate axis.
///
/// Used where a file stores an axis index (lumped element `Direction`,
/// probe `NormDir`); ids outside 0..=2 are rejected when reading.
///
/// ```
/// use csx::properties::Axis;
///
/// assert_eq!(Axis::try_from(2usize)?, Axis::Z);
/// assert_eq!("rho".parse::<Axis>()?, Axis::X);
/// assert_eq!(Axis::Y.index(), 1);
/// let v = [10.0, 20.0, 30.0];
/// assert_eq!(v[Axis::Z.index()], 30.0);
/// # Ok::<(), csx::Error>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
pub enum Axis {
    /// x (or ρ in cylindrical coordinates); id 0.
    #[default]
    X,
    /// y (or α); id 1.
    Y,
    /// z; id 2.
    Z,
}

impl Axis {
    /// All axes in order.
    pub const ALL: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];

    /// Axis for an index 0..=2.
    pub const fn from_index(i: usize) -> Option<Axis> {
        match i {
            0 => Some(Axis::X),
            1 => Some(Axis::Y),
            2 => Some(Axis::Z),
            _ => None,
        }
    }

    /// Index 0..=2, usable to index `[T; 3]` arrays.
    pub const fn index(self) -> usize {
        self as usize
    }
}

impl From<Axis> for usize {
    fn from(a: Axis) -> usize {
        a.index()
    }
}

impl TryFrom<usize> for Axis {
    type Error = Error;

    /// # Errors
    ///
    /// [`Error::Invalid`] for indices above 2.
    fn try_from(i: usize) -> Result<Axis, Error> {
        Axis::from_index(i).ok_or_else(|| Error::Invalid(format!("axis {i} is not 0, 1 or 2")))
    }
}

impl TryFrom<i32> for Axis {
    type Error = Error;

    /// # Errors
    ///
    /// [`Error::Invalid`] for ids outside 0..=2.
    fn try_from(i: i32) -> Result<Axis, Error> {
        usize::try_from(i)
            .ok()
            .and_then(Axis::from_index)
            .ok_or_else(|| Error::Invalid(format!("axis {i} is not 0, 1 or 2")))
    }
}

impl fmt::Display for Axis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Axis::X => "x",
            Axis::Y => "y",
            Axis::Z => "z",
        })
    }
}

impl FromStr for Axis {
    type Err = Error;

    /// Accepts `x`/`y`/`z` (also `rho`, `alpha`) or an index 0..=2.
    ///
    /// # Errors
    ///
    /// [`Error::Format`] for anything else.
    fn from_str(s: &str) -> Result<Axis, Error> {
        let bad = || Error::Format(format!("invalid axis {s:?}"));
        match s.trim().to_ascii_lowercase().as_str() {
            "x" | "rho" | "r" => Ok(Axis::X),
            "y" | "a" | "alpha" => Ok(Axis::Y),
            "z" => Ok(Axis::Z),
            t => parse_id(t)
                .and_then(|i| Axis::try_from(i).ok())
                .ok_or_else(bad),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids_round_trip() {
        for v in DumpType::ALL {
            assert_eq!(DumpType::from_id(v.id()), *v);
            assert_eq!(v.to_string().parse::<DumpType>().unwrap(), *v);
        }
        assert_eq!(DumpType::from_id(7), DumpType::Other(7));
        assert_eq!(DumpType::Other(7).id(), 7);
        assert_eq!("7".parse::<DumpType>().unwrap(), DumpType::Other(7));
        assert_eq!(
            "0.000000e+00".parse::<LumpedElementType>().unwrap(),
            LumpedElementType::Parallel
        );
        assert_eq!(
            "plane_wave".parse::<ExcitationType>().unwrap(),
            ExcitationType::PlaneWave
        );
        assert_eq!("mur-1st-sa".parse::<AbcType>().unwrap(), AbcType::Mur1stSa);
        assert!("bogus".parse::<ProbeType>().is_err());
        assert!("1.5".parse::<ProbeType>().is_err());
        assert!(DumpType::Sar10g.is_frequency_domain());
        assert!(!DumpType::BField.is_frequency_domain());
        assert!(ExcitationType::HFieldHard.is_magnetic());
        assert_eq!(DumpMode::default(), DumpMode::Node);
    }

    #[test]
    fn axes() {
        assert_eq!(Axis::try_from(2i32).unwrap(), Axis::Z);
        assert!(Axis::try_from(-1i32).is_err());
        assert_eq!("2.000000e+00".parse::<Axis>().unwrap(), Axis::Z);
        assert_eq!("y".parse::<Axis>().unwrap(), Axis::Y);
        assert!("3".parse::<Axis>().is_err());
        assert_eq!(Axis::Z.index(), 2);
    }
}
