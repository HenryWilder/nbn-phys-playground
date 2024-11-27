use std::ops::*;

pub mod vec3;
pub mod units;
use units::*;

pub struct MassPoint {
    pub position: Displacement,
    pub mass: Mass,
    pub velocity: Velocity,
    pub force: ForceVector,
}

impl AddAssign<Displacement> for MassPoint {
    fn add_assign(&mut self, rhs: Displacement) {
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
        let mass = self.mass;
        self.force += rhs * mass;
    }
}

impl MassPoint {
    pub const fn init(position: Displacement, mass: Mass) -> MassPointBuilder {
        MassPointBuilder::new(position, mass)
    }

    /// Apply a force from a point-source with inverse square falloff
    pub fn apply_point_force(&mut self, source: &PointForce) {
        let (direction, distance) = self.position.direction_and_distance(source.center);
        let area = distance * distance;
        let pressure = source.force / area;
        let force = pressure * Area(1.0) * direction;
        *self += force;
    }

    pub fn resolve(mut self, dt: Time) {
        self += self.force / self.mass * dt;
        self.force.clear();
        self += self.velocity * dt;
    }
}

pub struct MassPointBuilder {
    position: Displacement,
    mass: Mass,
    velocity: Velocity,
}

impl MassPointBuilder {
    pub const fn new(position: Displacement, mass: Mass) -> Self {
        Self {
            position,
            mass,
            velocity: Velocity::zero(),
        }
    }

    pub const fn initial_velocity(mut self, velocity: Velocity) -> Self {
        self.velocity = velocity;
        self
    }

    pub const fn build(self) -> MassPoint {
        MassPoint {
            position: self.position,
            mass: self.mass,
            velocity: self.velocity,
            force: ForceVector::zero(),
        }
    }
}

pub struct Fluid {
    density: Density,
    pressure: Pressure,
    /// force = viscosity * area * rate of sheer deformation
    viscosity: f32,
    /// specific volume = volume / mass
    specific_volume: SpecificVolume,
}
// const AIR: Fluid = Fluid { density: 1.225 };
// const WATER: Fluid = Fluid { density: 1.0 };

/// Positionless shape
pub enum Form {
    Sphere { radius: f32 },
}

impl Form {
    pub const fn drag_coefficient(&self) -> Scale {
        match self {
            Form::Sphere { .. } => Scale(0.47),
        }
    }

    pub fn volume(&self) -> Volume {
        match self {
            Form::Sphere { radius } => Volume(4.0 * std::f32::consts::FRAC_PI_3 * radius.powi(3)),
        }
    }

    pub fn surface_area(&self) -> Area {
        match self {
            Form::Sphere { radius } => Area(4.0 * std::f32::consts::PI * radius.powi(2)),
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
    pub center: Displacement,
    pub force: ForceScalar,
}
