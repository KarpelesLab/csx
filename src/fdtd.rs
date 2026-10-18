//! openEMS simulation settings: time stepping, excitation signal and
//! boundary conditions (the `<FDTD>` element of an openEMS file).
//!
//! Every setting is optional so that files only carrying a subset of the
//! attributes are written back unchanged; [`FdtdSettings::gaussian`]
//! gives a complete set of commonly used values. With the `xml` feature,
//! [`FdtdSettings::from_xml_element`] and [`FdtdSettings::to_xml_element`]
//! convert from / to the XML representation.

use core::f64::consts::PI;

use crate::error::{Error, Result};
use crate::expr::Expr;

/// The time signal used by all excitations.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum ExcitationSignal {
    /// Gaussian-modulated sine (type `0`): centre frequency `f0` and 20 dB
    /// corner frequency `fc`, i.e. the spectrum covers `f0 ± fc`.
    Gaussian {
        /// Centre frequency (Hz).
        f0: f64,
        /// 20 dB corner frequency (Hz).
        fc: f64,
    },
    /// Continuous sine of frequency `f0` (type `1`).
    Sinus {
        /// Frequency (Hz).
        f0: f64,
    },
    /// Dirac pulse (type `2`); `f_max` is the highest frequency of interest.
    Dirac {
        /// Highest frequency of interest (Hz).
        f_max: f64,
    },
    /// Unit step (type `3`); `f_max` is the highest frequency of interest.
    Step {
        /// Highest frequency of interest (Hz).
        f_max: f64,
    },
    /// User defined signal (type `10`): an expression of the time `t` in
    /// seconds (the names `f0` and `f_max` may also be used).
    Custom {
        /// Expression of `t`.
        function: String,
        /// Characteristic frequency (Hz).
        f0: f64,
        /// Highest frequency of interest (Hz).
        f_max: f64,
    },
}

/// `ln(10)`: a Gaussian spectrum `exp(-(π Δf τ)²)` drops by 20 dB in
/// amplitude when `(π Δf τ)² = ln 10`.
const LN_10: f64 = core::f64::consts::LN_10;

impl ExcitationSignal {
    /// Numeric type id used in files.
    pub fn id(&self) -> i32 {
        match self {
            ExcitationSignal::Gaussian { .. } => 0,
            ExcitationSignal::Sinus { .. } => 1,
            ExcitationSignal::Dirac { .. } => 2,
            ExcitationSignal::Step { .. } => 3,
            ExcitationSignal::Custom { .. } => 10,
        }
    }

    /// Highest frequency contained in the signal (Hz): `f0 + fc` for a
    /// Gaussian pulse, `f0` for a sine, `f_max` otherwise.
    pub fn max_frequency(&self) -> f64 {
        match *self {
            ExcitationSignal::Gaussian { f0, fc } => f0 + fc,
            ExcitationSignal::Sinus { f0 } => f0,
            ExcitationSignal::Dirac { f_max } | ExcitationSignal::Step { f_max } => f_max,
            ExcitationSignal::Custom { f_max, .. } => f_max,
        }
    }

    /// Envelope width `τ` of the Gaussian pulse (s).
    ///
    /// The pulse is `cos(2π f0 (t − T0)) · exp(−((t − T0)/τ)²)`. Its amplitude
    /// spectrum is proportional to `exp(−(π (f − f0) τ)²)`, which is 20 dB
    /// below the peak at `f − f0 = fc` when `τ = √(ln 10) / (π fc)`.
    /// Returns `None` for other signal types.
    pub fn gaussian_width(&self) -> Option<f64> {
        match *self {
            ExcitationSignal::Gaussian { fc, .. } if fc > 0.0 => Some(LN_10.sqrt() / (PI * fc)),
            _ => None,
        }
    }

    /// Delay `T0` of the Gaussian pulse centre (s): chosen so the envelope
    /// is 10⁻⁶ of its peak at `t = 0`, i.e. `T0 = τ √(6 ln 10)`.
    pub fn gaussian_delay(&self) -> Option<f64> {
        self.gaussian_width().map(|tau| tau * (6.0 * LN_10).sqrt())
    }

    /// Duration of the signal (s): `2 T0` for the Gaussian pulse (the
    /// envelope is again below 10⁻⁶ afterwards), `None` for signals without
    /// a natural end.
    pub fn duration(&self) -> Option<f64> {
        self.gaussian_delay().map(|t0| 2.0 * t0)
    }

    /// Builds an evaluator of the signal amplitude at time `t` (s).
    /// Fails only for custom signals whose expression is invalid.
    pub fn evaluator(&self) -> Result<SignalFn> {
        Ok(match self {
            ExcitationSignal::Custom {
                function,
                f0,
                f_max,
            } => {
                let (f0, f_max) = (*f0, *f_max);
                let e = Expr::parse(function)?.compile(&["t"], |n| match n {
                    "f0" => Some(f0),
                    "f_max" | "fmax" => Some(f_max),
                    _ => None,
                })?;
                SignalFn(Inner::Custom(e))
            }
            other => SignalFn(Inner::Builtin(other.clone())),
        })
    }

    /// Signal amplitude at time `t` (s). Custom signals that fail to parse
    /// evaluate to NaN; use [`ExcitationSignal::evaluator`] for repeated
    /// evaluation.
    pub fn eval(&self, t: f64) -> f64 {
        self.evaluator().map_or(f64::NAN, |f| f.eval(t))
    }
}

#[derive(Debug, Clone)]
enum Inner {
    Builtin(ExcitationSignal),
    Custom(crate::expr::CompiledExpr),
}

/// A compiled excitation signal; see [`ExcitationSignal::evaluator`].
#[derive(Debug, Clone)]
pub struct SignalFn(Inner);

impl SignalFn {
    /// Amplitude at time `t` (s).
    pub fn eval(&self, t: f64) -> f64 {
        match &self.0 {
            Inner::Custom(e) => e.eval(&[t]),
            Inner::Builtin(s) => match *s {
                ExcitationSignal::Gaussian { f0, .. } => {
                    let (Some(tau), Some(t0)) = (s.gaussian_width(), s.gaussian_delay()) else {
                        return 0.0;
                    };
                    let x = (t - t0) / tau;
                    (2.0 * PI * f0 * (t - t0)).cos() * (-x * x).exp()
                }
                ExcitationSignal::Sinus { f0 } => (2.0 * PI * f0 * t).sin(),
                ExcitationSignal::Dirac { .. } => {
                    if t == 0.0 {
                        1.0
                    } else {
                        0.0
                    }
                }
                ExcitationSignal::Step { .. } => {
                    if t >= 0.0 {
                        1.0
                    } else {
                        0.0
                    }
                }
                ExcitationSignal::Custom { .. } => unreachable!("compiled separately"),
            },
        }
    }
}

/// Boundary condition of one face of the simulation domain.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Boundary {
    /// Perfect electric conductor (`PEC`, numeric `0`).
    #[default]
    Pec,
    /// Perfect magnetic conductor (`PMC`, numeric `1`).
    Pmc,
    /// Mur first order absorbing boundary (`MUR`, numeric `2`).
    Mur,
    /// Perfectly matched layer with the given number of cells
    /// (`PML_n`, numeric `3` = 8 cells).
    Pml(u32),
}

impl Boundary {
    fn parse_opt(s: &str) -> Option<Boundary> {
        let t = s.trim().to_ascii_uppercase();
        match t.as_str() {
            "PEC" | "0" => return Some(Boundary::Pec),
            "PMC" | "1" => return Some(Boundary::Pmc),
            "MUR" | "2" => return Some(Boundary::Mur),
            "PML" | "3" => return Some(Boundary::Pml(8)),
            _ => {}
        }
        let rest = t.strip_prefix("PML")?.trim_start_matches(['_', ' ', '-']);
        rest.parse().ok().map(Boundary::Pml)
    }
}

impl core::str::FromStr for Boundary {
    type Err = Error;

    /// Parses `PEC`, `PMC`, `MUR`, `PML_n` (also `PML n`, `PML-n`, `PML`)
    /// or the numeric forms `0`–`3` (case-insensitive).
    fn from_str(s: &str) -> Result<Boundary> {
        Boundary::parse_opt(s)
            .ok_or_else(|| Error::Format(format!("unknown boundary condition {s:?}")))
    }
}

impl core::fmt::Display for Boundary {
    /// The textual form used in files (`PEC`, `PMC`, `MUR`, `PML_n`).
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Boundary::Pec => f.write_str("PEC"),
            Boundary::Pmc => f.write_str("PMC"),
            Boundary::Mur => f.write_str("MUR"),
            Boundary::Pml(n) => write!(f, "PML_{n}"),
        }
    }
}

/// Boundary conditions of the six faces.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BoundaryConditions {
    /// Faces in the order `xmin, xmax, ymin, ymax, zmin, zmax`.
    pub faces: [Boundary; 6],
    /// Phase velocity used by Mur boundaries (m/s); file attribute
    /// `MUR_PhaseVelocity` (name not confirmed by documentation).
    pub mur_phase_velocity: Option<f64>,
    /// PML grading function of the variables `D, dl, W, N, Z`; file
    /// attribute `PML_Grading` (name not confirmed by documentation).
    pub pml_grading: Option<String>,
    /// Unrecognised attributes preserved from files.
    pub extra_attributes: Vec<(String, String)>,
}

/// Attribute names of the six faces.
pub const FACE_NAMES: [&str; 6] = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"];

impl BoundaryConditions {
    /// The same condition on all faces.
    pub fn all(b: Boundary) -> Self {
        BoundaryConditions {
            faces: [b; 6],
            ..Default::default()
        }
    }
}

/// The `<FDTD>` settings of an openEMS simulation.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FdtdSettings {
    /// Maximum number of time steps (`NumberOfTimesteps`).
    pub number_of_timesteps: Option<u64>,
    /// Energy decay end criterion, e.g. `1e-5` = −50 dB (`endCriteria`).
    pub end_criteria: Option<f64>,
    /// Highest frequency of interest (Hz) (`f_max`). Dirac, step and custom
    /// signals take their `f_max` from this attribute when read from a file.
    pub f_max: Option<f64>,
    /// Over-sampling factor of probes (`OverSampling`).
    pub over_sampling: Option<u32>,
    /// Time step method (`TimeStepMethod`).
    pub time_step_method: Option<TimeStepMethod>,
    /// Maximum simulated time (s) (`MaxTime`).
    pub max_time: Option<f64>,
    /// Fixed time step (s) (`TimeStep`).
    pub time_step: Option<f64>,
    /// Factor applied to the computed time step (`TimeStepFactor`).
    pub time_step_factor: Option<f64>,
    /// Use cell-constant materials (`CellConstantMaterial`).
    pub cell_constant_material: Option<bool>,
    /// Coordinate system id of the simulation (`CoordSystem`).
    pub coord_system: Option<i32>,
    /// Multi-grid radii for cylindrical meshes (`MultiGrid`, verbatim).
    pub multi_grid: Option<String>,
    /// Excitation signal (`<Excitation>` child).
    pub excitation: Option<ExcitationSignal>,
    /// Unrecognised attributes of `<Excitation>`.
    pub excitation_extra_attributes: Vec<(String, String)>,
    /// Boundary conditions (`<BoundaryCond>` child).
    pub boundary: Option<BoundaryConditions>,
    /// Unrecognised attributes of `<FDTD>`.
    pub extra_attributes: Vec<(String, String)>,
}

impl FdtdSettings {
    /// Settings with no attribute set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Typical settings for a Gaussian pulse excitation (centre `f0`, 20 dB
    /// half bandwidth `fc`): 10⁹ time steps, −50 dB end criterion and 8-cell
    /// PML on all faces.
    pub fn gaussian(f0: f64, fc: f64) -> Self {
        FdtdSettings {
            number_of_timesteps: Some(1_000_000_000),
            end_criteria: Some(1e-5),
            f_max: Some(f0 + fc),
            excitation: Some(ExcitationSignal::Gaussian { f0, fc }),
            boundary: Some(BoundaryConditions::all(Boundary::Pml(8))),
            ..Default::default()
        }
    }

    /// Checks the settings for obviously invalid values.
    pub fn validate(&self) -> Result<()> {
        let bad = |m: &str| Err(Error::Invalid(m.into()));
        if self.end_criteria.is_some_and(|e| e.is_nan() || e <= 0.0) {
            return bad("end criteria must be positive");
        }
        if let Some(ExcitationSignal::Gaussian { f0, fc }) = self.excitation
            && (fc.is_nan() || fc <= 0.0 || f0 < 0.0)
        {
            return bad("Gaussian excitation needs fc > 0 and f0 >= 0");
        }
        if let Some(e @ ExcitationSignal::Custom { .. }) = &self.excitation {
            e.evaluator()?;
        }
        Ok(())
    }
}

/// How the FDTD time step is chosen (`TimeStepMethod`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TimeStepMethod {
    /// Courant–Friedrichs–Lewy criterion; id 1.
    Cfl,
    /// Rennings' criterion (openEMS default); id 3.
    Rennings,
    /// An id without a named variant, preserved verbatim.
    Other(i32),
}

impl TimeStepMethod {
    /// Variant for a numeric id (unknown ids become `Other`).
    pub const fn from_id(id: i32) -> Self {
        match id {
            1 => TimeStepMethod::Cfl,
            3 => TimeStepMethod::Rennings,
            o => TimeStepMethod::Other(o),
        }
    }

    /// Numeric id as used in files.
    pub const fn id(self) -> i32 {
        match self {
            TimeStepMethod::Cfl => 1,
            TimeStepMethod::Rennings => 3,
            TimeStepMethod::Other(o) => o,
        }
    }
}

impl From<i32> for TimeStepMethod {
    fn from(id: i32) -> Self {
        TimeStepMethod::from_id(id)
    }
}

impl From<TimeStepMethod> for i32 {
    fn from(m: TimeStepMethod) -> i32 {
        m.id()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boundaries() {
        assert_eq!("PML_8".parse::<Boundary>().ok(), Some(Boundary::Pml(8)));
        assert_eq!("PML 12".parse::<Boundary>().ok(), Some(Boundary::Pml(12)));
        assert_eq!("pml".parse::<Boundary>().ok(), Some(Boundary::Pml(8)));
        assert_eq!("1".parse::<Boundary>().ok(), Some(Boundary::Pmc));
        assert_eq!("3".parse::<Boundary>().ok(), Some(Boundary::Pml(8)));
        assert_eq!(" MUR ".parse::<Boundary>().ok(), Some(Boundary::Mur));
        assert_eq!("PEC".parse::<Boundary>().ok(), Some(Boundary::Pec));
        assert_eq!("PML_x".parse::<Boundary>().ok(), None);
        assert_eq!("foo".parse::<Boundary>().ok(), None);
        assert_eq!(Boundary::Pml(16).to_string(), "PML_16");
    }

    #[test]
    fn gaussian() {
        let g = ExcitationSignal::Gaussian { f0: 1e9, fc: 5e8 };
        let tau = g.gaussian_width().unwrap();
        let t0 = g.gaussian_delay().unwrap();
        assert!((tau - 10f64.ln().sqrt() / (PI * 5e8)).abs() < 1e-20);
        // peak at T0, envelope ~1e-6 at t = 0 and t = 2 T0
        assert!((g.eval(t0) - 1.0).abs() < 1e-12);
        assert!(g.eval(0.0).abs() <= 1.0001e-6);
        assert!(g.eval(g.duration().unwrap()).abs() <= 1.0001e-6);
        assert_eq!(g.max_frequency(), 1.5e9);
        assert_eq!(g.id(), 0);
    }

    #[test]
    fn other_signals() {
        let s = ExcitationSignal::Sinus { f0: 1.0 };
        assert!((s.eval(0.25) - 1.0).abs() < 1e-12);
        assert_eq!(ExcitationSignal::Step { f_max: 1e9 }.eval(-1.0), 0.0);
        assert_eq!(ExcitationSignal::Step { f_max: 1e9 }.eval(1.0), 1.0);
        assert_eq!(ExcitationSignal::Dirac { f_max: 1e9 }.eval(0.0), 1.0);
        let c = ExcitationSignal::Custom {
            function: "sin(2*pi*f0*t)".into(),
            f0: 2.0,
            f_max: 4.0,
        };
        assert!((c.eval(0.125) - 1.0).abs() < 1e-12);
        let bad = ExcitationSignal::Custom {
            function: "q*t".into(),
            f0: 1.0,
            f_max: 1.0,
        };
        assert!(bad.eval(1.0).is_nan());
        assert!(bad.evaluator().is_err());
    }

    #[test]
    fn validation() {
        assert!(FdtdSettings::gaussian(1e9, 5e8).validate().is_ok());
        let mut s = FdtdSettings::new();
        s.excitation = Some(ExcitationSignal::Gaussian { f0: 1e9, fc: 0.0 });
        assert!(s.validate().is_err());
    }
}
