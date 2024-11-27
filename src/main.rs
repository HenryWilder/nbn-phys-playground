use raylib::prelude::*;

pub mod extras;
pub mod phys;
pub mod sim;
use extras::*;

pub mod propulsion_sim;
use propulsion_sim::*;

fn main() {
    let window_width = 1280;
    let window_height = 720;
    let (mut rl, thread) = init()
        .size(window_width, window_height)
        .title("Nothin but Nade Physics Playground")
        .build();

    rl.set_target_fps(60);
    rl.disable_cursor();

    let mut sim = PropulsionSimulation::new(
        Vector3::zero(),
        5.0,
        0.1,
        1.0,
        Vector3::new(0.1, 0.1, 0.0),
        Vector3::new(8.0, 0.0, 0.0),
        Vector3::new(0.0, -9.8, 0.0),
    );

    let mut camera = Camera3D::perspective(Vector3::new(10.0, 10.0, 10.0), Vector3::zero(), Vector3::up(), 45.0);

    while !rl.window_should_close() {
        rl.update_camera(&mut camera, CameraMode::CAMERA_FREE);
        // sim.tick(&mut rl, 0.5);
        let mut d = rl.begin_drawing(&thread);
        d.clear_background(Color::BLACK);
        let mut d = d.begin_mode3D(camera);
        d.draw_debug_world_grid(1.0, 5, Color::RED, Color::GREEN, Color::BLUE, Some(Color::WHITE));
        sim.draw(&mut d);
    }
}
