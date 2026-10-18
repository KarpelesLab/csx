//! Property kind payloads.

use crate::error::Result;
use crate::math::Vec3;
use crate::param::{ParameterSet, Scalar};

use super::weight::{Aniso, WeightFn, weight_vars};

/// Electric/magnetic material with (optionally anisotropic and position
/// dependent) constitutive parameters.
///
/// * `epsilon`: relative permittivity
/// * `mue`: relative permeability
/// * `kappa`: electric conductivity (S/m)
/// * `sigma`: magnetic conductivity (Ω/m)
/// * `density`: mass density (kg/m³)
#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    /// Relative permittivity.
    pub epsilon: Aniso,
    /// Relative permeability.
    pub mue: Aniso,
    /// Electric conductivity.
    pub kappa: Aniso,
    /// Magnetic conductivity.
    pub sigma: Aniso,
    /// Mass density.
    pub density: Scalar,
    /// Weighting function of the density.
    pub density_weight: WeightFn,
    /// If true only the first axis of each quantity is used.
    pub isotropic: bool,
}

impl Default for Material {
    fn default() -> Self {
        Material {
            epsilon: Aniso::iso(1.0),
            mue: Aniso::iso(1.0),
            kappa: Aniso::iso(0.0),
            sigma: Aniso::iso(0.0),
            density: Scalar::new(0.0),
            density_weight: WeightFn::default(),
            isotropic: true,
        }
    }
}

/// Selects one of the four main material quantities.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MaterialQuantity {
    /// Relative permittivity.
    Epsilon,
    /// Relative permeability.
    Mue,
    /// Electric conductivity.
    Kappa,
    /// Magnetic conductivity.
    Sigma,
}

impl MaterialQuantity {
    /// All quantities in file order.
    pub const ALL: [MaterialQuantity; 4] = [
        MaterialQuantity::Epsilon,
        MaterialQuantity::Mue,
        MaterialQuantity::Kappa,
        MaterialQuantity::Sigma,
    ];

    /// Attribute name used in files.
    pub fn name(self) -> &'static str {
        match self {
            MaterialQuantity::Epsilon => "Epsilon",
            MaterialQuantity::Mue => "Mue",
            MaterialQuantity::Kappa => "Kappa",
            MaterialQuantity::Sigma => "Sigma",
        }
    }
}

impl Material {
    /// Material with the given relative permittivity and conductivity.
    pub fn new(epsilon: f64, kappa: f64) -> Self {
        let mut m = Material::default();
        m.epsilon.set(epsilon);
        m.kappa.set(kappa);
        m
    }

    /// The quantity selected by `q`.
    pub fn quantity(&self, q: MaterialQuantity) -> &Aniso {
        match q {
            MaterialQuantity::Epsilon => &self.epsilon,
            MaterialQuantity::Mue => &self.mue,
            MaterialQuantity::Kappa => &self.kappa,
            MaterialQuantity::Sigma => &self.sigma,
        }
    }

    /// Mutable access to the quantity selected by `q`.
    pub fn quantity_mut(&mut self, q: MaterialQuantity) -> &mut Aniso {
        match q {
            MaterialQuantity::Epsilon => &mut self.epsilon,
            MaterialQuantity::Mue => &mut self.mue,
            MaterialQuantity::Kappa => &mut self.kappa,
            MaterialQuantity::Sigma => &mut self.sigma,
        }
    }

    #[inline]
    fn axis(&self, dir: usize) -> usize {
        if self.isotropic { 0 } else { dir.min(2) }
    }

    /// Weighted value of `q` along axis `dir` at Cartesian point `p`.
    pub fn value_at(&self, q: MaterialQuantity, dir: usize, p: Vec3) -> f64 {
        let a = self.quantity(q);
        let d = self.axis(dir);
        if a.weight[d].is_unity() {
            a.base(d)
        } else {
            a.at_vars(d, &weight_vars(p))
        }
    }

    /// Weighted values of `q` on all three axes at `p`.
    pub fn values_at(&self, q: MaterialQuantity, p: Vec3) -> [f64; 3] {
        let a = self.quantity(q);
        if !a.is_weighted() {
            return core::array::from_fn(|d| a.base(self.axis(d)));
        }
        let v = weight_vars(p);
        core::array::from_fn(|d| a.at_vars(self.axis(d), &v))
    }

    /// Unweighted value of `q` along `dir`.
    pub fn base_value(&self, q: MaterialQuantity, dir: usize) -> f64 {
        self.quantity(q).base(self.axis(dir))
    }

    /// Weighted density at `p`.
    pub fn density_at(&self, p: Vec3) -> f64 {
        self.density.value() * self.density_weight.eval(p)
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        for q in MaterialQuantity::ALL {
            self.quantity_mut(q).update(p)?;
        }
        self.density.update(p)?;
        self.density_weight.update(p)
    }
}

/// One pole of a Drude/Lorentz dispersion model.
///
/// For the electric part (similar for the magnetic part):
/// `ε(ω) = ε∞ − ω_p² / (ω² − ω_L² − iω/τ)`; a zero pole frequency gives a
/// Drude pole. Frequencies are in Hz (not rad/s).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LorentzPole {
    /// Electric plasma frequency.
    pub eps_plasma_frequency: Aniso,
    /// Electric Lorentz pole frequency (0 for Drude).
    pub eps_lorentz_pole_frequency: Aniso,
    /// Electric relaxation time.
    pub eps_relax_time: Aniso,
    /// Magnetic plasma frequency.
    pub mue_plasma_frequency: Aniso,
    /// Magnetic Lorentz pole frequency.
    pub mue_lorentz_pole_frequency: Aniso,
    /// Magnetic relaxation time.
    pub mue_relax_time: Aniso,
}

impl LorentzPole {
    fn quantities_mut(&mut self) -> [&mut Aniso; 6] {
        [
            &mut self.eps_plasma_frequency,
            &mut self.eps_lorentz_pole_frequency,
            &mut self.eps_relax_time,
            &mut self.mue_plasma_frequency,
            &mut self.mue_lorentz_pole_frequency,
            &mut self.mue_relax_time,
        ]
    }

    /// Creates an all-zero pole (zero relaxation time means lossless).
    pub fn zero() -> Self {
        LorentzPole {
            eps_plasma_frequency: Aniso::iso(0.0),
            eps_lorentz_pole_frequency: Aniso::iso(0.0),
            eps_relax_time: Aniso::iso(0.0),
            mue_plasma_frequency: Aniso::iso(0.0),
            mue_lorentz_pole_frequency: Aniso::iso(0.0),
            mue_relax_time: Aniso::iso(0.0),
        }
    }
}

/// One pole of a Debye dispersion model:
/// `ε(ω) = ε∞ + Δε / (1 + iωτ)`.
#[derive(Debug, Clone, PartialEq)]
pub struct DebyePole {
    /// Permittivity step Δε.
    pub eps_delta: Aniso,
    /// Relaxation time τ (s).
    pub eps_relax_time: Aniso,
}

impl Default for DebyePole {
    fn default() -> Self {
        DebyePole {
            eps_delta: Aniso::iso(0.0),
            eps_relax_time: Aniso::iso(0.0),
        }
    }
}

/// Drude/Lorentz dispersive material.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LorentzMaterial {
    /// Static (high frequency limit) material parameters.
    pub material: Material,
    /// Dispersion poles.
    pub poles: Vec<LorentzPole>,
}

/// Debye dispersive material.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DebyeMaterial {
    /// Static (high frequency limit) material parameters.
    pub material: Material,
    /// Dispersion poles.
    pub poles: Vec<DebyePole>,
}

impl LorentzMaterial {
    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.material.update(p)?;
        for pole in &mut self.poles {
            pole.quantities_mut()
                .into_iter()
                .try_for_each(|q| q.update(p))?;
        }
        Ok(())
    }
}

impl DebyeMaterial {
    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.material.update(p)?;
        for pole in &mut self.poles {
            pole.eps_delta.update(p)?;
            pole.eps_relax_time.update(p)?;
        }
        Ok(())
    }
}

/// A material defined on a discrete voxel grid stored in an external file.
///
/// Only the description is modelled; the voxel data itself (usually HDF5)
/// is not loaded by this crate.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DiscMaterial {
    /// Static parameters used outside the voxel data.
    pub material: Material,
    /// Data file name.
    pub file: String,
    /// Scale factor applied to the voxel grid coordinates.
    pub scale: Scalar,
    /// Remaining attributes, preserved verbatim.
    pub attributes: Vec<(String, String)>,
}

/// A thin conducting sheet with finite conductivity.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ConductingSheet {
    /// Conductivity (S/m).
    pub conductivity: Scalar,
    /// Sheet thickness (m).
    pub thickness: Scalar,
}

/// Field excitation.
#[derive(Debug, Clone, PartialEq)]
pub struct Excitation {
    /// Port / excitation number.
    pub number: u32,
    /// Whether the excitation is active.
    pub enabled: bool,
    /// Excitation type: 0 = soft E, 1 = hard E, 2 = soft H, 3 = hard H,
    /// 10 = plane wave.
    pub excitation_type: i32,
    /// Excitation vector (amplitude per axis).
    pub excite: [Scalar; 3],
    /// Spatial weighting per axis.
    pub weight: [WeightFn; 3],
    /// Propagation direction (plane waves).
    pub prop_dir: [Scalar; 3],
    /// Time delay (s).
    pub delay: Scalar,
    /// Frequency hint (Hz), 0 if unused.
    pub frequency: Scalar,
}

impl Default for Excitation {
    fn default() -> Self {
        Excitation {
            number: 0,
            enabled: true,
            excitation_type: 0,
            excite: Default::default(),
            weight: Default::default(),
            prop_dir: Default::default(),
            delay: Scalar::new(0.0),
            frequency: Scalar::new(0.0),
        }
    }
}

impl Excitation {
    /// Weighted excitation amplitude along `dir` at `p`.
    pub fn amplitude_at(&self, dir: usize, p: Vec3) -> f64 {
        self.excite[dir].value() * self.weight[dir].eval(p)
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.excite.iter_mut().try_for_each(|s| s.update(p))?;
        self.weight.iter_mut().try_for_each(|w| w.update(p))?;
        self.prop_dir.iter_mut().try_for_each(|s| s.update(p))?;
        self.delay.update(p)?;
        self.frequency.update(p)
    }
}

/// A probe recording a field integral or current. Waveguide mode functions
/// are stored as the user attributes `ModeFunctionX/Y/Z` of the property.
#[derive(Debug, Clone, PartialEq)]
pub struct ProbeBox {
    /// Probe number.
    pub number: u32,
    /// Probe type: 0 = voltage, 1 = current, 2 = E-field, 3 = H-field,
    /// 10 = waveguide voltage mode, 11 = waveguide current mode.
    pub probe_type: i32,
    /// Weight (sign / scaling) of the probe result.
    pub weight: Scalar,
    /// Normal direction (for current probes), if any.
    pub norm_dir: Option<i32>,
    /// Start time of the recording (s).
    pub start_time: Scalar,
    /// Stop time of the recording (s).
    pub stop_time: Scalar,
    /// Frequencies for frequency-domain probing (`FD_Samples`).
    pub frequencies: Vec<f64>,
}

impl Default for ProbeBox {
    fn default() -> Self {
        ProbeBox {
            number: 0,
            probe_type: 0,
            weight: Scalar::new(1.0),
            norm_dir: None,
            start_time: Scalar::new(0.0),
            stop_time: Scalar::new(0.0),
            frequencies: Vec::new(),
        }
    }
}

/// A field dump box.
#[derive(Debug, Clone, PartialEq)]
pub struct DumpBox {
    /// Common probe settings.
    pub probe: ProbeBox,
    /// What to dump: 0 = E, 1 = H, 2 = J, 3 = rot(H), 4 = D, 5 = B; 10–15 the
    /// frequency-domain variants; 20 = local SAR, 21 = 1 g, 22 = 10 g
    /// averaged SAR, 29 = raw SAR data.
    pub dump_type: i32,
    /// Interpolation mode: 0 = none, 1 = node (default), 2 = cell.
    pub dump_mode: i32,
    /// File type: 0 = VTK, 1 = HDF5.
    pub file_type: i32,
    /// Spatial sub-sampling per axis.
    pub sub_sampling: Option<[u32; 3]>,
    /// Optional target resolution per axis.
    pub opt_resolution: Option<[Scalar; 3]>,
    /// Multi-grid level.
    pub multi_grid_level: i32,
}

impl Default for DumpBox {
    fn default() -> Self {
        DumpBox {
            probe: ProbeBox::default(),
            dump_type: 0,
            dump_mode: 1,
            file_type: 0,
            sub_sampling: None,
            opt_resolution: None,
            multi_grid_level: 0,
        }
    }
}

/// Lumped R/L/C element. Unset `c` / `l` values are NaN.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedElement {
    /// Direction of the element (0, 1, 2).
    pub direction: i32,
    /// Whether end caps (PEC) are added (default true).
    pub caps: bool,
    /// Resistance (Ω).
    pub r: Scalar,
    /// Capacitance (F).
    pub c: Scalar,
    /// Inductance (H).
    pub l: Scalar,
    /// Element type: 0 = parallel, 1 = series (optional).
    pub le_type: Option<i32>,
}

impl Default for LumpedElement {
    fn default() -> Self {
        LumpedElement {
            direction: 0,
            caps: true,
            r: Scalar::new(f64::NAN),
            c: Scalar::new(f64::NAN),
            l: Scalar::new(f64::NAN),
            le_type: None,
        }
    }
}

impl LumpedElement {
    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.r.update(p)?;
        self.c.update(p)?;
        self.l.update(p)
    }
}

/// Absorbing boundary condition sheet.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct AbsorbingBC {
    /// Whether the outward normal points in positive direction.
    pub normal_sign_positive: bool,
    /// Phase velocity of the absorbed wave (m/s); 0 = speed of light.
    pub phase_velocity: Scalar,
    /// Boundary type id.
    pub abc_type: i32,
}
