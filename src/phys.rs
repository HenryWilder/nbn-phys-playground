use std::ops::*;
use raylib::prelude as rl;
// use crate::VectorMath;

#[derive(Debug, Clone, Copy, Default)]
struct Vec3 {
    x: f32,
    y: f32,
    z: f32,
}

impl From<rl::Vector3> for Vec3 {
    fn from(rl::Vector3 { x, y, z }: rl::Vector3) -> Self {
        Self { x, y, z }
    }
}

impl From<Vec3> for rl::Vector3 {
    fn from(Vec3 { x, y, z }: Vec3) -> Self {
        Self { x, y, z }
    }
}

impl Vec3 {
    pub const fn zero() -> Self {
        Self { x: 0.0, y: 0.0, z: 0.0 }
    }

    pub const fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }

    pub fn dot(self, rhs: Self) -> f32 {
        self.x * rhs.x +
        self.y * rhs.y +
        self.z * rhs.z
    }

    pub fn length_sqr(self) -> f32 {
        self.dot(self)
    }

    pub fn length(self) -> f32 {
        self.length_sqr().sqrt()
    }

    pub fn distance_sqr(self, other: Self) -> f32 {
        (self - other).length_sqr()
    }

    pub fn distance(self, other: Self) -> f32 {
        (self - other).length()
    }

    pub fn normalized(self) -> Self {
        self / self.length()
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}
impl Add for Vec3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}
impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}
impl Mul<f32> for Vec3 {
    type Output = Self;
    fn mul(self, rhs: f32) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}
impl Div<f32> for Vec3 {
    type Output = Self;
    fn div(self, rhs: f32) -> Self::Output {
        let inv = rhs.recip();
        Self {
            x: self.x * inv,
            y: self.y * inv,
            z: self.z * inv,
        }
    }
}
impl Mul<Vec3> for f32 {
    type Output = Vec3;
    fn mul(self, rhs: Vec3) -> Self::Output {
        Vec3 {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}
impl Div<Vec3> for f32 {
    type Output = Vec3;
    fn div(self, rhs: Vec3) -> Self::Output {
        Vec3 {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}
impl Add<f32> for Vec3 {
    type Output = Self;
    /// Add magnitude
    fn add(self, rhs: f32) -> Self::Output {
        self * (1.0 + rhs / self.length())
    }
}
impl Sub<f32> for Vec3 {
    type Output = Self;
    /// Subtract magnitude
    fn sub(self, rhs: f32) -> Self::Output {
        self * (1.0 - rhs / self.length())
    }
}
impl Add<Vec3> for f32 {
    type Output = Vec3;
    /// Add magnitude
    fn add(self, rhs: Vec3) -> Self::Output {
        rhs * (1.0 + self / rhs.length())
    }
}
impl Sub<Vec3> for f32 {
    type Output = Vec3;
    /// Subtract magnitude
    fn sub(self, rhs: Vec3) -> Self::Output {
        -rhs * (1.0 + self / rhs.length())
    }
}
impl Mul for Vec3 {
    type Output = f32;
    /// Multiply magnitudes
    fn mul(self, rhs: Self) -> Self::Output {
        (self.length_sqr() * rhs.length_sqr()).sqrt()
    }
}
impl Div for Vec3 {
    type Output = f32;
    /// Divide magnitudes
    fn div(self, rhs: Self) -> Self::Output {
        (self.length_sqr() / rhs.length_sqr()).sqrt()
    }
}
impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.add(rhs);
    }
}
impl SubAssign for Vec3 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.sub(rhs);
    }
}
impl AddAssign<f32> for Vec3 {
    fn add_assign(&mut self, rhs: f32) {
        *self = self.add(rhs);
    }
}
impl SubAssign<f32> for Vec3 {
    fn sub_assign(&mut self, rhs: f32) {
        *self = self.sub(rhs);
    }
}
impl MulAssign<f32> for Vec3 {
    fn mul_assign(&mut self, rhs: f32) {
        *self = self.mul(rhs);
    }
}
impl DivAssign<f32> for Vec3 {
    fn div_assign(&mut self, rhs: f32) {
        *self = self.div(rhs);
    }
}

macro_rules! rule_op {
    (@impl $({$($impliedby:tt)+})? $trait:ident::$fn:ident($a:ident $op:tt $b:ident = $c:ident)) => {
        impl $trait<$b> for $a {
            type Output = $c;
            #[doc = stringify!($($($impliedby)+ => )? $c = $a $op $b)]
            fn $fn(self, rhs: $b) -> Self::Output {
                $c(self.0 $op rhs.0)
            }
        }
    };

    ($c:ident = $a:ident + $b:ident $(<== $($impliedby:tt)+)?) => { rule_op!{ @impl $({$($impliedby)+})? Add::add($a + $b = $c) } };
    ($c:ident = $a:ident - $b:ident $(<== $($impliedby:tt)+)?) => { rule_op!{ @impl $({$($impliedby)+})? Sub::sub($a - $b = $c) } };
    ($c:ident = $a:ident * $b:ident $(<== $($impliedby:tt)+)?) => { rule_op!{ @impl $({$($impliedby)+})? Mul::mul($a * $b = $c) } };
    ($c:ident = $a:ident / $b:ident $(<== $($impliedby:tt)+)?) => { rule_op!{ @impl $({$($impliedby)+})? Div::div($a / $b = $c) } };
}

macro_rules! calc_rule {
    ($self:ident += $rhs:ident) => {
        calc_rule!($self = $self + $rhs);
        impl AddAssign<$rhs> for $self { #[doc = stringify!($self = $self + $rhs)] fn add_assign(&mut self, rhs: $rhs) { self.0 += rhs.0; } }
        impl SubAssign<$rhs> for $self { #[doc = stringify!($self = $self + $rhs => $self = $self - $rhs)] fn sub_assign(&mut self, rhs: $rhs) { self.0 -= rhs.0; } }
    };

    ($self:ident *= $rhs:ident) => {
        calc_rule!($self = $self * $rhs);
        impl MulAssign<$rhs> for $self { #[doc = stringify!($self = $self * $rhs)] fn mul_assign(&mut self, rhs: $rhs) { self.0 *= rhs.0; } }
        impl DivAssign<$rhs> for $self { #[doc = stringify!($self = $self * $rhs => $self = $self / $rhs)] fn div_assign(&mut self, rhs: $rhs) { self.0 /= rhs.0; } }
    };

    ($c:ident = $a:ident + $b:ident, ..) => { rule_op!($c = $a + $b); rule_op!($a = $c - $b <== $c = $a + $b); rule_op!($b = $c - $a <== $c = $a + $b); rule_op!($c = $b + $a <== $c = $a + $b); };
    ($c:ident = $a:ident - $b:ident, ..) => { rule_op!($c = $a - $b); rule_op!($a = $c + $b <== $c = $a - $b); rule_op!($b = $a - $c <== $c = $a - $b); };
    ($c:ident = $a:ident * $b:ident, ..) => { rule_op!($c = $a * $b); rule_op!($a = $c / $b <== $c = $a * $b); rule_op!($b = $c / $a <== $c = $a * $b); rule_op!($c = $b * $a <== $c = $a * $b); };
    ($c:ident = $a:ident / $b:ident, ..) => { rule_op!($c = $a / $b); rule_op!($a = $c * $b <== $c = $a / $b); rule_op!($b = $a / $c <== $c = $a / $b); };
    ($c:ident = $a:ident + $b:ident    ) => { rule_op!($c = $a + $b); rule_op!($a = $c - $b <== $c = $a + $b); };
    ($c:ident = $a:ident - $b:ident    ) => { rule_op!($c = $a - $b); rule_op!($a = $c + $b <== $c = $a - $b); };
    ($c:ident = $a:ident * $b:ident    ) => { rule_op!($c = $a * $b); rule_op!($a = $c / $b <== $c = $a * $b); };
    ($c:ident = $a:ident / $b:ident    ) => { rule_op!($c = $a / $b); rule_op!($a = $c * $b <== $c = $a / $b); };
}

macro_rules! rules {
    () => {};

    // Prevent duplicate implementations
    ($c:ident = $a:ident $op:tt $b:ident, .. $($rest:tt)*) => {
        calc_rule!{ $c = $a $op $b, .. }
        rules!($($rest)*);
    };

    ($c:ident = $a:ident $op:tt $b:ident $($rest:tt)*) => {
        calc_rule!{ $c = $a $op $b }
        rules!($($rest)*);
    };
}

macro_rules! calc_wrapper {
    (
        $name:ident in $units:literal ~ $mag:ident
        $(+= $add_assign:ident,)*
        $(*= $mul_assign:ident,)*
        $(= $a:ident $op:tt $b:ident,)*
    ) => {
        #[doc = $units]
        #[derive(Debug, Clone, Copy, Default)]
        pub struct $name(pub Vec3);
        impl $name {
            pub const fn zero() -> Self {
                Self(Vec3::zero())
            }
            pub fn magnitude(self) -> $mag {
                $mag(self.0.length())
            }
            pub fn magnitude_sqr(self) -> $mag {
                $mag(self.0.length_sqr())
            }
            pub fn direction(self) -> Self {
                Self(self.0.normalized())
            }
            pub fn direction_and_magnitude(self) -> (Self, $mag) {
                let mag = self.0.length();
                (Self(self.0 * mag.recip()), $mag(mag))
            }
        }
        $(calc_rule!{ $name += $add_assign })*
        $(calc_rule!{ $name *= $mul_assign })*
        $(calc_rule!{ $name = $a $op $b })*
    };

    (
        $name:ident in $units:literal
        $(+= $add_assign:ident,)*
        $(*= $mul_assign:ident,)*
        $(= $a:ident $op:tt $b:ident,)*
    ) => {
        #[doc = $units]
        pub struct $name(pub f32);
        $(calc_rule!{ $name += $add_assign })*
        $(calc_rule!{ $name *= $mul_assign })*
        $(calc_rule!{ $name = $a $op $b })*
    };
}

macro_rules! calc {
    (
        $(
            $name:ident $(~ $mag:ident)? in $units:literal
            $(+= $add_assign:ident,)*
            $(*= $mul_assign:ident,)*
            $(= $a:ident $op:tt $b:ident,)*
        )+
    ) => {
        $(
            calc_wrapper!{
                $name in $units $(~ $mag)?
                $(+= $add_assign,)*
                $(*= $mul_assign,)*
                $(= $a $op $b,)*
            }
        )+
    };
}

calc!{
    Scale                             in "scalar"         += Self,                  *= Self,
    Time                              in "s"              += Self,                  *= Scale,
    Length                            in "m"              += Self,                  *= Scale,
    Mass                              in "kg"             += Self,                  *= Scale,
    Position      ~ Length            in "(m,m,m)"        += Self, += Length,       *= Scale,
    Area                              in "m²"             += Self,                  *= Scale,
    Volume                            in "m³"             += Self,                  *= Scale,
    Speed                             in "m/s"            += Self,                  *= Scale,
    Velocity      ~ Speed             in "(m,m,m)/s"      += Self, += Speed,        *= Scale,
    Acceleration                      in "m/s/s"          += Self,                  *= Scale,
    AccelerationVector ~ Acceleration in "(m,m,m)/s/s"    += Self, += Acceleration, *= Scale,
    Force                             in "kg*m/s/s"       += Self,                  *= Scale,
    ForceVector        ~ Force        in "kg*(m,m,m)/s/s" += Self, += Force,        *= Scale,
    Pressure                          in "kg/m²"          += Self,                  *= Scale,
}

rules!{
    Speed    = Length   / Time, ..
    Velocity = Position / Time, ..

    Acceleration       = Speed    / Time, ..
    AccelerationVector = Velocity / Time, ..

    Force       = Mass * Acceleration, ..
    ForceVector = Mass * AccelerationVector, ..

    Pressure = Force / Area, ..

    Area = Length * Length

    Volume = Area * Length, ..
}

pub struct MassPoint {
    pub position: Position,
    pub mass: Mass,
    pub velocity: Velocity,
    pub force: ForceVector,
}

impl AddAssign<Position> for MassPoint {
    fn add_assign(&mut self, rhs: Position) {
        self.position += rhs;
    }
}
impl AddAssign<Velocity> for MassPoint {
    fn add_assign(&mut self, rhs: Velocity) {
        self.velocity += rhs;
    }
}
impl AddAssign<ForceVector> for MassPoint {
    fn add_assign(&mut self, rhs: ForceVector) {
        self.force += rhs;
    }
}
impl AddAssign<AccelerationVector> for MassPoint {
    fn add_assign(&mut self, rhs: AccelerationVector) {
        self.force += rhs * self.mass;
    }
}

impl MassPoint {
    pub const fn init(position: Position, mass: Mass) -> MassPointBuilder {
        MassPointBuilder::new(position, mass)
    }

    // pub fn apply_point_force(&mut self, source: &PointForce) {
    //     let delta = self.position - source.position;
    //     let distance_sqr = delta.length_sqr();
    //     if distance_sqr != 0.0 {
    //         let force_direction = delta / distance_sqr.sqrt(); // normalize
    //         self.add_force(force_direction * source.force);
    //     }
    // }

    pub fn resolve(mut self, dt: Time) {
        self.velocity += self.force / self.mass * dt;
        self.force = ForceVector::zero();
        self.position += self.velocity * dt;
    }
}

pub struct MassPointBuilder {
    position: Position,
    mass: Mass,
    velocity: Option<Velocity>,
}

impl MassPointBuilder {
    pub const fn new(position: Position, mass: Mass) -> Self {
        Self {
            position,
            mass,
            velocity: None,
        }
    }

    pub const fn initial_velocity(&mut self, velocity: Velocity) -> &mut Self {
        self.velocity = Some(velocity);
        self
    }

    pub const fn build(self) -> MassPoint {
        MassPoint {
            position: self.position,
            mass: self.mass,
            velocity: self.velocity.unwrap_or(Velocity::zero()),
            force: ForceVector::zero(),
        }
    }
}

pub struct Fluid {
    /// density = mass / volume
    density: f32,
    // /// pressure = force / area
    // pressure: f32,
    // /// force = viscosity * area * rate of sheer deformation
    // viscosity: f32,
    // /// specific volume = volume / mass
    // specific_volume: f32,
}
const AIR: Fluid = Fluid { density: 1.225 };
const WATER: Fluid = Fluid { density: 1.0 };

/// Positionless shape
pub enum Form {
    Sphere { radius: f32 },
}

impl Form {
    pub const fn drag_coefficient(&self) -> f32 {
        match self {
            Form::Sphere { .. } => 0.47,
        }
    }

    pub fn volume(&self) -> f32 {
        match self {
            Form::Sphere { radius } => 4.0 * std::f32::consts::FRAC_PI_3 * radius.powi(3),
        }
    }

    pub fn surface_area(&self) -> f32 {
        match self {
            Form::Sphere { radius } => 4.0 * std::f32::consts::PI * radius.powi(2),
        }
    }
}

pub struct RigidBody {
    pub mp: MassPoint,
    pub form: Form,
}

impl RigidBody {
    // pub fn apply_drag(&mut self, coef: f32, exposed_surface_area: f32) {
    //     const AIR_DENSITY: f32 = 1.21;
    //     self.mp.force += (self.mp.force.powi(2)) * -(coef * exposed_surface_area * AIR_DENSITY * 0.5);
    // }
}

pub struct PointForce {
    pub position: Vec3,
    pub force: f32,
}
