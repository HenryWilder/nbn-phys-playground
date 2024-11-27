use std::ops::*;
use super::vec3::*;

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
    ($c:ident =        1 / $b:ident $(<== $($impliedby:tt)+)?) => {
        impl $b {
            #[doc = stringify!($($($impliedby)+ => )? $c = 1 / $b)]
            pub fn recip(self) -> $c {
                $c(1.0 / self.0)
            }
        }
    };
}

macro_rules! calc_rule {
    ($self:ident += $rhs:ident) => {
        calc_rule!($self = $self + $rhs);
        impl AddAssign<$rhs> for $self {
            #[doc = stringify!($self = $self + $rhs)]
            fn add_assign(&mut self, rhs: $rhs) { self.0 += rhs.0; }
        }
        impl SubAssign<$rhs> for $self {
            #[doc = stringify!($self = $self + $rhs => $self = $self - $rhs)]
            fn sub_assign(&mut self, rhs: $rhs) { self.0 -= rhs.0; }
        }
    };

    ($self:ident *= $rhs:ident) => {
        calc_rule!($self = $self * $rhs);
        impl MulAssign<$rhs> for $self {
            #[doc = stringify!($self = $self * $rhs)]
            fn mul_assign(&mut self, rhs: $rhs) { self.0 *= rhs.0; }
        }
        impl DivAssign<$rhs> for $self {
            #[doc = stringify!($self = $self * $rhs => $self = $self / $rhs)]
            fn div_assign(&mut self, rhs: $rhs) { self.0 /= rhs.0; }
        }
    };

    ($c:ident = $a:ident + $b:ident, ..) => { rule_op!($c = $a + $b); rule_op!($a = $c - $b <== $c = $a + $b); rule_op!($b = $c - $a <== $c = $a + $b); rule_op!($c = $b + $a <== $c = $a + $b); };
    ($c:ident = $a:ident - $b:ident, ..) => { rule_op!($c = $a - $b); rule_op!($a = $c + $b <== $c = $a - $b); rule_op!($b = $a - $c <== $c = $a - $b); };
    ($c:ident = $a:ident * $b:ident, ..) => { rule_op!($c = $a * $b); rule_op!($a = $c / $b <== $c = $a * $b); rule_op!($b = $c / $a <== $c = $a * $b); rule_op!($c = $b * $a <== $c = $a * $b); };
    ($c:ident = $a:ident / $b:ident, ..) => { rule_op!($c = $a / $b); rule_op!($a = $c * $b <== $c = $a / $b); rule_op!($b = $a / $c <== $c = $a / $b); };
    ($c:ident = $a:ident + $b:ident    ) => { rule_op!($c = $a + $b); rule_op!($a = $c - $b <== $c = $a + $b); };
    ($c:ident = $a:ident - $b:ident    ) => { rule_op!($c = $a - $b); rule_op!($a = $c + $b <== $c = $a - $b); };
    ($c:ident = $a:ident * $b:ident    ) => { rule_op!($c = $a * $b); rule_op!($a = $c / $b <== $c = $a * $b); };
    ($c:ident = $a:ident / $b:ident    ) => { rule_op!($c = $a / $b); rule_op!($a = $c * $b <== $c = $a / $b); };
    ($c:ident =        1 / $b:ident    ) => { rule_op!($c =  1 / $b); rule_op!($b =  1 / $c <== $c =  1 / $b); };
}

macro_rules! rules {
    () => {};

    ($c:ident = $a:ident $op:tt $b:ident, .. $($rest:tt)*) => {
        calc_rule!{ $c = $a $op $b, .. }
        rules!($($rest)*);
    };

    // Prevent duplicate implementations
    ($c:ident = $a:ident $op:tt $b:ident $($rest:tt)*) => {
        calc_rule!{ $c = $a $op $b }
        rules!($($rest)*);
    };

    ($c:ident = 1 / $b:ident $($rest:tt)*) => {
        calc_rule!{ $c = 1 / $b }
        rules!($($rest)*);
    };

    ($self:ident += $rhs:ident $($rest:tt)*) => {
        calc_rule!{ $self += $rhs }
        rules!($($rest)*);
    };

    ($self:ident *= $rhs:ident $($rest:tt)*) => {
        calc_rule!{ $self *= $rhs }
        rules!($($rest)*);
    };
}

pub trait CalcVector: Sized {
    type Magnitude;
    fn clear(&mut self);
    fn magnitude(self) -> Self::Magnitude;
    /// Get the square of the magnitude, saving a sqrt.
    fn magnitude_sqr(self) -> Self::Magnitude;
    fn direction(self) -> Direction;
    /// Decompose a vector into its direction (without magnitude) and its magnitude (without direction)
    fn direction_and_magnitude(self) -> (Direction, Self::Magnitude);
    fn direction_and_distance(self, other: Self) -> (Direction, Self::Magnitude) where Self: Sub<Output = Self>;
}

macro_rules! calc_wrapper {
    (
        $name:ident $(in $units:literal)? ~ $mag:ident
    ) => {
        $(#[doc = $units])?
        #[derive(Debug, Clone, Copy, Default)]
        pub struct $name(pub Vec3);
        impl $name {
            pub const fn zero() -> Self {
                Self(Vec3::zero())
            }
        }
        impl From<Vec3> for $name {
            fn from(value: Vec3) -> Self {
                Self(value)
            }
        }
        impl Neg for $name {
            type Output = Self;
            fn neg(self) -> Self::Output {
                Self(-self.0)
            }
        }
        impl CalcVector for $name {
            type Magnitude = $mag;
            fn clear(&mut self) {
                *self = Self::zero();
            }
            fn magnitude(self) -> $mag {
                $mag(self.0.length())
            }
            fn magnitude_sqr(self) -> $mag {
                $mag(self.0.length_sqr())
            }
            fn direction(self) -> Direction {
                Direction(self.0.normalized())
            }
            fn direction_and_magnitude(self) -> (Direction, $mag) {
                let len = self.0.length();
                (if len != 0.0 {
                    Direction(self.0 * len.recip())
                } else {
                    Direction::zero()
                }, $mag(len))
            }
            fn direction_and_distance(self, other: Self) -> (Direction, $mag)
            where Self: Sub<Output = Self> {
                (other - self).direction_and_magnitude()
            }
        }
    };

    (
        $name:ident $(in $units:literal)?
    ) => {
        $(#[doc = $units])?
        #[derive(Debug, Clone, Copy, Default, PartialEq, PartialOrd)]
        pub struct $name(pub f32);
        impl $name {
            pub const fn zero() -> Self {
                Self(0.0)
            }
        }
        impl Neg for $name {
            type Output = Self;
            fn neg(self) -> Self::Output {
                Self(-self.0)
            }
        }
    };
}

macro_rules! calc {
    // unitless scalar
    (@auto_rules $mag:ident) => {
        rules!{
            $mag += $mag
            $mag *= Scale
        }
    };
    // unitful scalar
    (@auto_rules $mag:ident in $mag_units:literal) => {
        rules!{
            $mag += $mag
            $mag *= Scale
        }
    };

    // unitless vector
    (@auto_rules $name:ident [$mag:ident $(in $mag_units:literal)?]) => {
        rules!{
            $name += $name
            $name *= Scale
        }
        calc!(@auto_rules $mag $(in $mag_units)?);
    };
    // unitful vector
    (@auto_rules $name:ident in $units:literal [$mag:ident $(in $mag_units:literal)?]) => {
        rules!{
            $name += $name
            $name *= Scale
            $name += $mag
            $name = Direction * $mag
            $name = $mag * Direction
        }
        calc!(@auto_rules $mag $(in $mag_units)?);
    };

    (
        $($name:ident $(in $units:literal)? $([ $mag:ident $(in $mag_units:literal)? ])?;)*
    ) => {
        $(
            $(
                calc_wrapper!{
                    $mag $(in $mag_units)?
                }
            )?
            calc_wrapper!{
                $name $(in $units)? $(~ $mag)?
            }
            calc!(@auto_rules $name $(in $units)? $([$mag $(in $mag_units)?])?);
        )*
    };
}

calc!{
    Direction [ Scale ];

    Time               in "s";
    Mass               in "kg";
    Area               in "m²";
    Volume             in "m³";
    Frequency          in "1/s";
    Displacement       in "(m,m,m)"       [ Length             in "m"       ];
    Velocity           in "(m,m,m)/s"     [ Speed              in "m/s"     ];
    AccelerationVector in "(m,m,m)/s²"    [ AccelerationScalar in "m/s²"    ];
    ForceVector        in "kg*(m,m,m)/s²" [ ForceScalar        in "kg*m/s²" ];
    MomentumVector     in "kg*(m,m,m)/s"  [ MomentumScalar     in "kg*m/s"  ];
    Pressure           in "kg*m/s²/m²";
    Density            in "kg/m³";
    SpecificVolume     in "m³/kg";
    Energy             in "kg*m²/s²";
    Power              in "kg*m²/s³";
}

rules!{
    Speed    = Length       / Time, ..
    Velocity = Displacement / Time, ..

    Frequency = 1 / Time

    AccelerationScalar = Speed    / Time, ..
    AccelerationVector = Velocity / Time, ..

    ForceScalar = Mass * AccelerationScalar, ..
    ForceVector = Mass * AccelerationVector, ..

    ForceScalar = MomentumScalar / Time, ..
    ForceVector = MomentumVector / Time, ..

    MomentumScalar = Mass * Speed, ..
    MomentumVector = Mass * Velocity, ..

    Pressure = ForceScalar / Area, ..
    Pressure = Energy / Volume, ..

    Area = Length * Length

    Volume = Area * Length, ..
    SpecificVolume = Volume / Mass, ..
    SpecificVolume = 1 / Density

    Density = Mass / Volume, ..

    Energy = ForceScalar * Length, ..
    Power = Energy / Time, ..
}
