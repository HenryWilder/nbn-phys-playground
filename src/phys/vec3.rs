use std::ops::*;
use raylib::prelude as rl;

#[derive(Debug, Clone, Copy, Default)]
pub struct Vec3 {
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
