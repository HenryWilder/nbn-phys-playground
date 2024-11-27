use raylib::prelude::*;

pub trait MathWithVector<V> {
    fn add(self, v: V) -> V;
    fn sub(self, v: V) -> V;
    fn mul(self, v: V) -> V;
    fn div(self, v: V) -> V;
}

impl MathWithVector<Vector3> for f32 {
    fn add(self, v: Vector3) -> Vector3 {
        Vector3 {
            x: self + v.x,
            y: self + v.y,
            z: self + v.z,
        }
    }

    fn sub(self, v: Vector3) -> Vector3 {
        Vector3 {
            x: self - v.x,
            y: self - v.y,
            z: self - v.z,
        }
    }

    fn mul(self, v: Vector3) -> Vector3 {
        Vector3 {
            x: self * v.x,
            y: self * v.y,
            z: self * v.z,
        }
    }

    fn div(self, v: Vector3) -> Vector3 {
        Vector3 {
            x: self / v.x,
            y: self / v.y,
            z: self / v.z,
        }
    }
}

pub trait VectorMath {
    fn length_sqr(self) -> f32;
    fn distance_sqr(self, other: Self) -> f32;
    fn powi(self, n: i32) -> Self;
    fn powf(self, n: f32) -> Self;
}

impl VectorMath for Vector3 {
    fn length_sqr(self) -> f32 {
        self.x * self.x +
        self.y * self.y +
        self.z * self.z
    }

    fn distance_sqr(self, other: Self) -> f32 {
        (other - self).length_sqr()
    }

    fn powi(self, n: i32) -> Self {
        Self {
            x: self.x.powi(n),
            y: self.y.powi(n),
            z: self.z.powi(n),
        }
    }

    fn powf(self, n: f32) -> Self {
        Self {
            x: self.x.powf(n),
            y: self.y.powf(n),
            z: self.z.powf(n),
        }
    }
}

pub trait ColorMath {
    fn add(self, other: Self) -> Self;
    fn scale(self, scale: f32) -> Self;
}
impl ColorMath for Color {
    fn add(self, other: Self) -> Self {
        Self {
            r: self.r.saturating_add(other.r),
            g: self.g.saturating_add(other.g),
            b: self.b.saturating_add(other.b),
            a: self.a.saturating_add(other.a),
        }
    }

    fn scale(self, scale: f32) -> Self {
        Self {
            r: (self.r as f32 * scale) as u8,
            g: (self.g as f32 * scale) as u8,
            b: (self.b as f32 * scale) as u8,
            a: self.a,
        }
    }
}

pub trait DebugDraw3D: RaylibDraw3D {
    fn draw_debug_sphere(&mut self, center: Vector3, radius: f32, color: Color) {
        self.draw_sphere_wires(center, radius, 3, 6, color);
    }

    fn draw_debug_world_grid(
        &mut self,
        scale: f32,
        minor: u32,
        x_color: Color,
        y_color: Color,
        z_color: Color,
        o_color: Option<Color>,
    ) {
        let slices: i32 = 100;
        let extent = slices as f32 * scale;
        let minor_mod = minor + 1;

        let x_axis = (rvec3(1.0, 0.0, 0.0), x_color);
        let y_axis = (rvec3(0.0, 1.0, 0.0), y_color);
        let z_axis = (rvec3(0.0, 0.0, 1.0), z_color);

        for ((along, along_color), against) in [
            (x_axis, [y_axis, z_axis]),
            (y_axis, [x_axis, z_axis]),
            (z_axis, [x_axis, y_axis]),
        ] {
            for (against, against_color) in against {
                let major_color = along_color.add(against_color.scale(0.5));
                let minor_color = major_color.alpha(0.25);
                let mut start = (against + along) * -extent;
                let mut end   = (against - along) * -extent;
                let step = against * scale;
                for i in -slices..=slices {
                    self.draw_line_3D(start, end, if i.unsigned_abs() % minor_mod == 0 { major_color } else { minor_color });
                    start += step;
                    end   += step;
                }
            }
        }

        if let Some(color) = o_color {
            self.draw_line_3D(x_axis.0 * -extent, x_axis.0 * extent, color);
            self.draw_line_3D(y_axis.0 * -extent, y_axis.0 * extent, color);
            self.draw_line_3D(z_axis.0 * -extent, z_axis.0 * extent, color);
        }
    }
}

impl<D: RaylibDraw3D> DebugDraw3D for D {}
