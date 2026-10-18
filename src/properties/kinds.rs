//! Property kind payloads (re-exported from [`crate::properties`]).

use crate::error::Result;
use crate::math::Vec3;
use crate::param::{ParameterSet, Scalar};

use super::types::{
    AbcType, Axis, DumpFileType, DumpMode, DumpType, ExcitationType, LumpedElementType, ProbeType,
};
use super::weight::{Aniso, WeightFn, weight_vars};

/// Electric/magnetic material with (optionally anisotropic and position
/// dependent) constitutive parameters.
///
/// * `epsilon`: relative permittivity
/// * `mue`: relative permeability
/// * `kappa`: electric conductivity (S/m)
/// * `sigma`: magnetic conductivity (Ω/m)
/// * `density`: mass density (kg/m³)
///
/// The default is vacuum (εᵣ = μᵣ = 1, lossless), isotropic.
///
/// ```
/// use csx::properties::{Material, MaterialQuantity};
/// use csx::Vec3;
///
/// // isotropic FR4
/// let fr4 = Material::new(4.3, 0.02);
/// assert_eq!(fr4.base_value(MaterialQuantity::Epsilon, 2), 4.3);
///
/// // anisotropic (uniaxial) material: ε_z differs
/// let mut m = Material::new(9.4, 0.0);
/// m.isotropic = false;
/// m.epsilon.value[2].set(11.6);
/// assert_eq!(m.values_at(MaterialQuantity::Epsilon, Vec3::ZERO), [9.4, 9.4, 11.6]);
///
/// // graded permittivity along x
/// let mut g = Material::new(2.0, 0.0);
/// g.epsilon.set_weight("1 + x/10");
/// assert_eq!(g.value_at(MaterialQuantity::Epsilon, 0, Vec3::new(5.0, 0.0, 0.0)), 3.0);
/// ```
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

    /// Weighted value of `q` along axis `dir` (0..=2; higher values are
    /// clamped to 2) at Cartesian point `p`. Isotropic materials use axis 0
    /// for every direction.
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
/// Drude pole. Frequencies are in Hz (not rad/s): `ω_p = 2π·f_p`.
///
/// File attributes (1-based pole index `n`): `EpsilonPlasmaFrequency_n`,
/// `EpsilonLorPoleFrequency_n`, `EpsilonRelaxTime_n` and the `Mue…`
/// counterparts.
///
/// ```
/// use csx::properties::{Aniso, LorentzPole};
///
/// // Drude pole: plasma frequency 2 THz, relaxation time 10 fs
/// let pole = LorentzPole {
///     eps_plasma_frequency: Aniso::iso(2e12),
///     eps_relax_time: Aniso::iso(1e-14),
///     ..LorentzPole::zero()
/// };
/// assert_eq!(pole.eps_lorentz_pole_frequency.base(0), 0.0);
/// ```
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
///
/// File attributes (1-based pole index `n`): `EpsilonDelta_n`,
/// `EpsilonRelaxTime_n`.
///
/// ```
/// use csx::properties::DebyePole;
///
/// let p = DebyePole::new(74.0, 8.3e-12); // water-like relaxation
/// assert_eq!(p.eps_delta.base(0), 74.0);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct DebyePole {
    /// Permittivity step Δε.
    pub eps_delta: Aniso,
    /// Relaxation time τ (s).
    pub eps_relax_time: Aniso,
}

impl DebyePole {
    /// Isotropic pole with permittivity step `eps_delta` and relaxation time
    /// `relax_time` (s).
    pub fn new(eps_delta: f64, relax_time: f64) -> Self {
        DebyePole {
            eps_delta: Aniso::iso(eps_delta),
            eps_relax_time: Aniso::iso(relax_time),
        }
    }
}

impl Default for DebyePole {
    fn default() -> Self {
        DebyePole {
            eps_delta: Aniso::iso(0.0),
            eps_relax_time: Aniso::iso(0.0),
        }
    }
}

/// Drude/Lorentz dispersive material (`<LorentzMaterial>`): a static
/// [`Material`] (ε∞, μ∞, losses) plus any number of [`LorentzPole`]s.
///
/// ```
/// use csx::properties::{Aniso, LorentzMaterial, LorentzPole, Material};
/// use csx::Property;
///
/// let drude = LorentzMaterial {
///     material: Material::new(1.0, 0.0),
///     poles: vec![LorentzPole {
///         eps_plasma_frequency: Aniso::iso(2.1e15),
///         eps_relax_time: Aniso::iso(1e-14),
///         ..LorentzPole::zero()
///     }],
/// };
/// let p = Property::lorentz_material("metal_drude", drude);
/// assert_eq!(p.kind.type_name(), "LorentzMaterial");
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LorentzMaterial {
    /// Static (high frequency limit) material parameters.
    pub material: Material,
    /// Dispersion poles.
    pub poles: Vec<LorentzPole>,
}

/// Debye dispersive material (`<DebyeMaterial>`): a static [`Material`]
/// (ε∞) plus any number of [`DebyePole`]s.
///
/// ```
/// use csx::properties::{DebyeMaterial, DebyePole, Material};
/// use csx::Property;
///
/// let tissue = DebyeMaterial {
///     material: Material::new(4.0, 0.2),
///     poles: vec![DebyePole::new(50.0, 7.2e-12), DebyePole::new(1100.0, 3.2e-8)],
/// };
/// let p = Property::debye_material("muscle", tissue);
/// assert_eq!(p.kind.material().unwrap().epsilon.base(0), 4.0);
/// ```
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
/// is not loaded by this crate. In files the data source is a
/// `<DiscFile Type File UseDBBackground Scale/>` child element.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DiscMaterial {
    /// Static parameters used outside the voxel data.
    pub material: Material,
    /// Data file name.
    pub file: String,
    /// Scale factor applied to the voxel grid coordinates.
    pub scale: Scalar,
    /// Remaining attributes of the `<DiscFile>` element (`Type`,
    /// `UseDBBackground`, ...), preserved verbatim.
    pub attributes: Vec<(String, String)>,
}

/// A thin conducting sheet with finite conductivity (`<ConductingSheet>`),
/// applied to 2-D primitives. A thickness of 0 behaves like a perfect
/// conductor.
///
/// ```
/// use csx::properties::ConductingSheet;
///
/// let cu = ConductingSheet::new(5.8e7, 35e-6); // 35 µm copper
/// assert_eq!(cu.thickness.value(), 35e-6);
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ConductingSheet {
    /// Conductivity (S/m).
    pub conductivity: Scalar,
    /// Sheet thickness (m).
    pub thickness: Scalar,
}

impl ConductingSheet {
    /// Sheet with conductivity (S/m) and thickness (m).
    pub fn new(conductivity: impl Into<Scalar>, thickness: impl Into<Scalar>) -> Self {
        ConductingSheet {
            conductivity: conductivity.into(),
            thickness: thickness.into(),
        }
    }
}

/// Field excitation (`<Excitation>`): a source applied on the mesh edges
/// covered by the property's primitives.
///
/// The amplitude per axis is `excite[d] · weight[d](x, y, z)`; the time
/// signal itself is part of the simulation settings
/// ([`crate::fdtd::ExcitationSignal`]).
///
/// ```
/// use csx::properties::{Excitation, ExcitationType, WeightFn};
/// use csx::Vec3;
///
/// // soft E-field along z with a cosine taper across x
/// let mut e = Excitation::new([0.0, 0.0, 1.0]).with_number(1);
/// e.weight[2] = WeightFn::new("cos(pi*x/10)");
/// assert_eq!(e.excitation_type, ExcitationType::EFieldSoft);
/// assert!((e.amplitude_at(2, Vec3::new(0.0, 0.0, 0.0)) - 1.0).abs() < 1e-12);
/// assert!(e.amplitude_at(2, Vec3::new(5.0, 0.0, 0.0)).abs() < 1e-12);
/// assert_eq!(e.amplitude_at(0, Vec3::ZERO), 0.0);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Excitation {
    /// Port / excitation number.
    pub number: u32,
    /// Whether the excitation is active.
    pub enabled: bool,
    /// Excitation type.
    pub excitation_type: ExcitationType,
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
            excitation_type: ExcitationType::EFieldSoft,
            excite: Default::default(),
            weight: Default::default(),
            prop_dir: Default::default(),
            delay: Scalar::new(0.0),
            frequency: Scalar::new(0.0),
        }
    }
}

impl Excitation {
    /// Soft E-field excitation (type 0) with the given excitation vector.
    pub fn new(excite: [f64; 3]) -> Self {
        Excitation {
            excite: excite.map(Scalar::new),
            ..Default::default()
        }
    }

    /// Builder: excitation type.
    #[must_use]
    pub fn with_type(mut self, excitation_type: ExcitationType) -> Self {
        self.excitation_type = excitation_type;
        self
    }

    /// Builder: port / excitation number.
    #[must_use]
    pub fn with_number(mut self, number: u32) -> Self {
        self.number = number;
        self
    }

    /// Weighted excitation amplitude along `dir` (0..=2) at Cartesian point
    /// `p`.
    ///
    /// # Panics
    ///
    /// If `dir > 2`.
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

/// A probe recording a field integral or current (`<ProbeBox>`). Waveguide
/// mode functions are stored as the user attributes `ModeFunctionX/Y/Z` of
/// the property.
///
/// See [`ProbeType`] for what is measured and which primitive shape it
/// expects (line for voltage, face for current, point for fields).
///
/// ```
/// use csx::properties::{Axis, ProbeBox, ProbeType};
///
/// let u = ProbeBox::new(ProbeType::Voltage).with_weight(-1.0);
/// let i = ProbeBox::new(ProbeType::Current).with_norm_dir(Axis::Z).with_frequencies([1e9, 2e9]);
/// assert_eq!(u.weight.value(), -1.0);
/// assert_eq!(i.frequencies, vec![1e9, 2e9]);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct ProbeBox {
    /// Probe number.
    pub number: u32,
    /// Probe type.
    pub probe_type: ProbeType,
    /// Weight (sign / scaling) of the probe result.
    pub weight: Scalar,
    /// Normal direction (for current probes), if any (`NormDir`, `-1` in
    /// files when unset).
    pub norm_dir: Option<Axis>,
    /// Start time of the recording (s).
    pub start_time: Scalar,
    /// Stop time of the recording (s).
    pub stop_time: Scalar,
    /// Frequencies for frequency-domain probing (`FD_Samples`).
    pub frequencies: Vec<f64>,
}

impl ProbeBox {
    /// Probe of the given type.
    pub fn new(probe_type: ProbeType) -> Self {
        ProbeBox {
            probe_type,
            ..Default::default()
        }
    }

    /// Builder: result weight.
    #[must_use]
    pub fn with_weight(mut self, w: f64) -> Self {
        self.weight = Scalar::new(w);
        self
    }

    /// Builder: normal direction.
    #[must_use]
    pub fn with_norm_dir(mut self, d: Axis) -> Self {
        self.norm_dir = Some(d);
        self
    }

    /// Builder: frequency-domain samples.
    #[must_use]
    pub fn with_frequencies(mut self, f: impl IntoIterator<Item = f64>) -> Self {
        self.frequencies = f.into_iter().collect();
        self
    }
}

impl Default for ProbeBox {
    fn default() -> Self {
        ProbeBox {
            number: 0,
            probe_type: ProbeType::Voltage,
            weight: Scalar::new(1.0),
            norm_dir: None,
            start_time: Scalar::new(0.0),
            stop_time: Scalar::new(0.0),
            frequencies: Vec::new(),
        }
    }
}

/// A field dump box (`<DumpBox>`): writes field data over the bounding box
/// of its primitives.
///
/// Frequency-domain dump types ([`DumpType::is_frequency_domain`]) need
/// frequencies in [`ProbeBox::frequencies`] of `probe`.
///
/// ```
/// use csx::properties::{DumpBox, DumpFileType, DumpMode, DumpType};
///
/// let mut d = DumpBox::new(DumpType::EFieldFd).with_file_type(DumpFileType::Hdf5);
/// d.probe.frequencies = vec![2.4e9];
/// d.dump_mode = DumpMode::Cell;
/// assert!(d.dump_type.is_frequency_domain());
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct DumpBox {
    /// Common probe settings.
    pub probe: ProbeBox,
    /// What to dump.
    pub dump_type: DumpType,
    /// Interpolation mode (default: node).
    pub dump_mode: DumpMode,
    /// Output file format.
    pub file_type: DumpFileType,
    /// Spatial sub-sampling per axis.
    pub sub_sampling: Option<[u32; 3]>,
    /// Optional target resolution per axis.
    pub opt_resolution: Option<[Scalar; 3]>,
    /// Multi-grid level.
    pub multi_grid_level: i32,
}

impl DumpBox {
    /// Dump of the given type.
    pub fn new(dump_type: DumpType) -> Self {
        DumpBox {
            dump_type,
            ..Default::default()
        }
    }

    /// Builder: output file format.
    #[must_use]
    pub fn with_file_type(mut self, file_type: DumpFileType) -> Self {
        self.file_type = file_type;
        self
    }
}

impl Default for DumpBox {
    fn default() -> Self {
        DumpBox {
            probe: ProbeBox::default(),
            dump_type: DumpType::EField,
            dump_mode: DumpMode::Node,
            file_type: DumpFileType::Vtk,
            sub_sampling: None,
            opt_resolution: None,
            multi_grid_level: 0,
        }
    }
}

/// Lumped R/L/C element (`<LumpedElement>`) spread over the mesh edges along
/// `direction` inside its primitives. Unset `r` / `c` / `l` values are NaN
/// (and omitted in files).
///
/// ```
/// use csx::properties::{Axis, LumpedElement, LumpedElementType};
///
/// // 50 Ω in series with 1 nH
/// let mut le = LumpedElement::resistor(Axis::Z, 50.0).with_l(1e-9);
/// le.le_type = Some(LumpedElementType::Series);
/// assert!(le.c.value().is_nan());
/// assert!(le.caps);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedElement {
    /// Direction of the element.
    pub direction: Axis,
    /// Whether end caps (PEC) are added (default true).
    pub caps: bool,
    /// Resistance (Ω).
    pub r: Scalar,
    /// Capacitance (F).
    pub c: Scalar,
    /// Inductance (H).
    pub l: Scalar,
    /// How R, L and C are connected (`LEtype`, optional in files).
    pub le_type: Option<LumpedElementType>,
}

impl Default for LumpedElement {
    fn default() -> Self {
        LumpedElement {
            direction: Axis::X,
            caps: true,
            r: Scalar::new(f64::NAN),
            c: Scalar::new(f64::NAN),
            l: Scalar::new(f64::NAN),
            le_type: None,
        }
    }
}

impl LumpedElement {
    /// Element along axis `direction` with resistance `r` (Ω); C and L
    /// unset.
    pub fn resistor(direction: Axis, r: f64) -> Self {
        LumpedElement {
            direction,
            r: Scalar::new(r),
            ..Default::default()
        }
    }

    /// Builder: capacitance (F).
    #[must_use]
    pub fn with_c(mut self, c: f64) -> Self {
        self.c = Scalar::new(c);
        self
    }

    /// Builder: inductance (H).
    #[must_use]
    pub fn with_l(mut self, l: f64) -> Self {
        self.l = Scalar::new(l);
        self
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.r.update(p)?;
        self.c.update(p)?;
        self.l.update(p)
    }
}

/// Absorbing boundary condition sheet (`<AbsorbingBC>`), applied to 2-D
/// primitives.
///
/// ```
/// use csx::properties::{AbcType, AbsorbingBC};
/// use csx::Property;
///
/// let abc = AbsorbingBC { normal_sign_positive: true, abc_type: AbcType::Mur1st, ..Default::default() };
/// let p = Property::absorbing_bc("abc", abc);
/// assert_eq!(p.kind.type_name(), "AbsorbingBC");
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct AbsorbingBC {
    /// Whether the outward normal points in positive direction.
    pub normal_sign_positive: bool,
    /// Phase velocity of the absorbed wave (m/s); 0 = speed of light.
    pub phase_velocity: Scalar,
    /// Boundary algorithm.
    pub abc_type: AbcType,
}
