// use raylib::prelude::*;
// use crate::{sim::*, phys::*};
// pub use crate::sim::Simulation;

// pub struct PropulsionSimulation {
//     radius: f32,
//     source: Vector3,
//     impulse: f32,
//     p_init: Vector3,
//     p_final: Vector3,
//     v_init: Vector3,
//     gravity: Vector3,
//     trail: Trail,
//     expected: Path3D,
//     p: MassPoint,
// }

// impl PropulsionSimulation {
//     pub fn new(source: Vector3, impulse: f32, radius: f32, mass: f32, p_init: Vector3, p_final: Vector3, gravity: Vector3) -> Self {
//         let v_init = (p_init - source).normalized() * impulse / mass;
//         Self {
//             radius,
//             source,
//             impulse,
//             p_init,
//             p_final,
//             v_init,
//             gravity,
//             trail: Trail::new(0.5, Some(256)),
//             expected: Path3D::from_fn(|t|
//                 Vector3 {
//                     //  Acceleration                Velocity         Position
//                     x:  0.5 * gravity.x * t * t  +  v_init.x * t  +  p_init.x,
//                     y:  0.5 * gravity.y * t * t  +  v_init.y * t  +  p_init.y,
//                     z:  0.5 * gravity.z * t * t  +  v_init.z * t  +  p_init.z,
//                 },
//                 0.0..=100.0,
//                 100,
//             ),
//             p: MassPoint::init(p_init, mass)
//                 .initial_velocity(v_init)
//                 .build(),
//         }
//     }
// }

// impl Simulation for PropulsionSimulation {
//     fn tick(&mut self, rl: &mut RaylibHandle, scale: f32) {
//         self.p.add_accel(self.gravity);
//         self.p.apply_drag(0.47, 2.0 * std::f32::consts::PI * self.radius * self.radius);
//         self.p.resolve(rl.get_frame_time() * scale);
//         // if self.p.position.y <= 0.0 {
//         //     self.p.position.y = 0.0;
//         //     self.p.velocity = Vector3::zero();
//         // }
//         self.trail.try_push_point(self.p.position);
//     }

//     fn draw(&self, d: &mut impl RaylibDraw3D) {
//         self.expected.draw(d, Color::MAGENTA);

//         // d.draw_debug_sphere(self.p.position, self.radius, Color::BLUEVIOLET);
//         // d.draw_line_3D(self.p.position, self.p.position + self.p.velocity, Color::BLUE);
//         // d.draw_line_3D(self.p.position, self.p.position + self.gravity, Color::ORANGERED);
//         // self.trail.draw(d, Color::BLUEVIOLET);

//         // d.draw_debug_sphere(self.source, 0.5, Color::YELLOW);
//         // d.draw_line_3D(self.source, self.source + self.v_init, Color::ORANGE);
//         // d.draw_debug_sphere(self.p_final, 0.5, Color::DODGERBLUE);
//     }
// }
