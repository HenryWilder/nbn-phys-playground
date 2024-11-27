use std::{collections::VecDeque, ops::RangeInclusive};
use raylib::prelude::*;
use crate::extras::*;

pub trait Simulation {
    fn tick(&mut self, rl: &mut RaylibHandle, scale: f32);
    fn draw(&self, d: &mut impl RaylibDraw3D);
}

pub struct Path3D {
    pub points: VecDeque<Vector3>,
}

impl Path3D {
    pub fn new() -> Self {
        Self {
            points: VecDeque::new(),
        }
    }

    fn evaluate_formula(mut f: impl FnMut(f32) -> Vector3, domain: RangeInclusive<f32>, samples: usize) -> impl Iterator<Item = Vector3> {
        let (t_min, t_max) = (*domain.start(), *domain.end());
        let domain_size = t_max - t_min;
        let step = domain_size / samples as f32;
        (0..samples).map(move |i| {
            let t = t_min + i as f32 * step;
            assert!(t_min <= t && t <= t_max);
            f(t)
        })
    }

    pub fn from_fn(f: impl FnMut(f32) -> Vector3, domain: RangeInclusive<f32>, samples: usize) -> Self {
        Self {
            points: Self::evaluate_formula(f, domain, samples)
                .collect()
        }
    }

    pub fn from_piecewise<const N: usize>(pieces: [(&mut dyn FnMut(f32) -> Vector3, RangeInclusive<f32>, usize); N]) -> Self {
        Self {
            points: pieces
                .into_iter()
                .flat_map(|(f, domain, samples)| Self::evaluate_formula(f, domain, samples))
                .collect(),
        }
    }

    pub fn draw(&self, d: &mut impl RaylibDraw3D, color: Color) {
        for i in 1..self.points.len() {
            d.draw_line_3D(self.points[i - 1], self.points[i], color);
        }
        for pt in &self.points {
            d.draw_debug_sphere(*pt, 0.125, color);
        }
    }
}

pub struct Trail {
    path: Path3D,
    freq: f32,
    max_points: Option<usize>,
}

impl Trail {
    pub fn new(freq: f32, max_points: Option<usize>) -> Self {
        Self {
            path: Path3D::new(),
            freq,
            max_points,
        }
    }

    pub fn try_push_point(&mut self, point: Vector3) {
        let back = self.path.points.back();
        if back.is_none() || back.is_some_and(|back| back.distance_to(point) >= self.freq) {
            if let Some(max_points) = self.max_points {
                if self.path.points.len() >= max_points {
                    assert_eq!(self.path.points.len(), max_points, "points should always pop before exceeding max");
                    self.path.points.pop_front();
                }
            }
            self.path.points.push_back(point);
        }
    }

    pub fn draw(&self, d: &mut impl RaylibDraw3D, color: Color) {
        self.path.draw(d, color);
    }
}
