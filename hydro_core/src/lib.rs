//! hydro_core: 1D water hammer solver (Method of Characteristics).
//! Layout: reservoir (x = 0) -> penstock -> valve (x = L).
//! H is head above the valve datum (metres of water). Slope is ignored for now.

pub const G: f64 = 9.81;

// ---------------------------------------------------------------- physics

/// Joukowsky head rise (m) for an instant velocity change.
/// a: wave speed (m/s), v: velocity change (m/s).
pub fn joukowsky(a: f64, v: f64) -> f64 {
    a * v / G
}

/// Wave speed (m/s) in an elastic pipe.
/// k: water bulk modulus (Pa), rho: water density (kg/m3), d: diameter (m),
/// e: wall thickness (m), young: pipe Young's modulus (Pa),
/// c1: anchoring factor (1.0 is a good starting value).
pub fn wave_speed(k: f64, rho: f64, d: f64, e: f64, young: f64, c1: f64) -> f64 {
    (k / rho).sqrt() / (1.0 + (k * d) / (young * e) * c1).sqrt()
}

/// Darcy friction factor, Swamee-Jain (explicit approximation of Colebrook).
/// eps: roughness (m), d: diameter (m), re: Reynolds number.
pub fn swamee_jain(eps: f64, d: f64, re: f64) -> f64 {
    let term = eps / (3.7 * d) + 5.74 / re.powf(0.9);
    0.25 / term.log10().powi(2)
}

// ---------------------------------------------------------------- inputs / outputs

pub struct Params {
    pub length: f64,        // L, m
    pub diameter: f64,      // D, m
    pub wall: f64,          // e, m
    pub young: f64,         // pipe E, Pa  Measures how stiff the material is
    pub bulk: f64,          // water K, Pa
    pub rho: f64,           // kg/m3
    pub roughness: f64,     // eps, m
    pub c1: f64,            // anchoring factor
    pub h_reservoir: f64,   // m above valve
    pub q0: f64,            // steady flow, m3/s
    pub n: usize,           // number of reaches
    pub t_end: f64,         // s
    pub t_close: f64,       // valve closure time, s (0 = instant)
    pub friction: bool,     // false = f = 0 (for the Joukowsky check)
}

pub struct Results {
    pub a: f64,
    pub dt: f64,
    pub time: Vec<f64>,
    pub h_valve: Vec<f64>,
    pub q_valve: Vec<f64>,
    pub h_max: Vec<f64>, // envelope per node
    pub h_min: Vec<f64>,
}

/// Valve opening: 1 = open, 0 = closed. Linear for now.
fn tau(t: f64, t_close: f64) -> f64 {
    if t_close <= 0.0 {
        if t > 0.0 { 0.0 } else { 1.0 }
    } else {
        (1.0 - t / t_close).max(0.0)
    }
}

// ---------------------------------------------------------------- solver

pub fn simulate(p: &Params) -> Results {
    let n = p.n;
    let area = std::f64::consts::PI * p.diameter * p.diameter / 4.0;
    let a = wave_speed(p.bulk, p.rho, p.diameter, p.wall, p.young, p.c1);

    // Grid: pick N first, then dt so that dt = dx / a holds exactly.
    let dx = p.length / n as f64;
    let dt = dx / a;
    let steps = (p.t_end / dt).ceil() as usize;

    // Friction factor from the steady velocity.
    let v0 = p.q0 / area;
    let f = if p.friction {
        let re = v0 * p.diameter / 1.0e-6; // water ~20 C, nu = 1e-6 m2/s
        swamee_jain(p.roughness, p.diameter, re)
    } else {
        0.0
    };

    let b = a / (G * area);                                       // impedance
    let r = f * dx / (2.0 * G * p.diameter * area * area);        // resistance

    // Steady state: head falls linearly along the pipe because of friction.
    let mut h_old: Vec<f64> = (0..=n).map(|i| p.h_reservoir - i as f64 * r * p.q0 * p.q0).collect();
    let mut q_old: Vec<f64> = vec![p.q0; n + 1];
    let mut h_new = h_old.clone();
    let mut q_new = q_old.clone();

    // Valve constant: Q = cv * tau * sqrt(H)
    let cv = p.q0 / h_old[n].sqrt();

    let mut res = Results {
        a,
        dt,
        time: vec![0.0],
        h_valve: vec![h_old[n]],
        q_valve: vec![q_old[n]],
        h_max: h_old.clone(),
        h_min: h_old.clone(),
    };

    for step in 1..=steps {
        let t = step as f64 * dt;

        // Interior nodes: C+ from i-1, C- from i+1 (both read the OLD arrays).
        for i in 1..n {
            let cp = h_old[i - 1] + b * q_old[i - 1] - r * q_old[i - 1] * q_old[i - 1].abs();
            let cm = h_old[i + 1] - b * q_old[i + 1] + r * q_old[i + 1] * q_old[i + 1].abs();
            q_new[i] = (cp - cm) / (2.0 * b);
            h_new[i] = (cp + cm) / 2.0;
        }

        // Upstream boundary: reservoir at constant head, uses C- from node 1.
        {
            let cm = h_old[1] - b * q_old[1] + r * q_old[1] * q_old[1].abs();
            h_new[0] = p.h_reservoir;
            q_new[0] = (h_new[0] - cm) / b;
        }

        // Downstream boundary: valve, uses C+ from node n-1.
        {
            let cp = h_old[n - 1] + b * q_old[n - 1] - r * q_old[n - 1] * q_old[n - 1].abs();
            let k = cv * tau(t, p.t_close);
            let k2 = k * k;
            // Solve Q = k*sqrt(cp - B*Q)  ->  Q^2 + B*k2*Q - k2*cp = 0
            let q = if k2 == 0.0 {
                0.0
            } else {
                (-b * k2 + (b * b * k2 * k2 + 4.0 * k2 * cp.max(0.0)).sqrt()) / 2.0
            };
            q_new[n] = q;
            h_new[n] = cp - b * q;
        }

        // Swap buffers: new becomes old for the next step.
        std::mem::swap(&mut h_old, &mut h_new);
        std::mem::swap(&mut q_old, &mut q_new);

        for i in 0..=n {
            res.h_max[i] = res.h_max[i].max(h_old[i]);
            res.h_min[i] = res.h_min[i].min(h_old[i]);
        }
        res.time.push(t);
        res.h_valve.push(h_old[n]);
        res.q_valve.push(q_old[n]);
    }
    res
}

// ---------------------------------------------------------------- tests

#[cfg(test)]
mod tests {
    use super::*;

    fn steel_params() -> Params {
        Params {
            length: 1000.0,
            diameter: 1.0,
            wall: 0.02,
            young: 200.0e9,
            bulk: 2.2e9,
            rho: 1000.0,
            roughness: 0.000045,
            c1: 1.0,
            h_reservoir: 100.0,
            q0: 3.0,
            n: 50,
            t_end: 10.0,
            t_close: 0.0,
            friction: false,
        }
    }

    #[test]
    fn joukowsky_basic() {
        let h = joukowsky(1000.0, 4.0);
        assert!((h - 407.75).abs() < 0.01);
    }

    #[test]
    fn wave_speed_steel() {
        let a = wave_speed(2.2e9, 1000.0, 1.0, 0.02, 200.0e9, 1.0);
        assert!((a - 1191.0).abs() < 2.0, "a = {a}");
    }

    #[test]
    fn instant_closure_matches_joukowsky() {
        let p = steel_params();
        let res = simulate(&p);
        let area = std::f64::consts::PI * p.diameter * p.diameter / 4.0;
        let expected = joukowsky(res.a, p.q0 / area);
        let rise = res.h_max[p.n] - p.h_reservoir;
        let err = (rise - expected).abs() / expected * 100.0;
        assert!(err < 0.5, "MOC {rise:.3} vs Joukowsky {expected:.3} ({err:.3}%)");
    }

    #[test]
    fn period_is_4l_over_a() {
        let p = steel_params();
        let res = simulate(&p);
        let n = p.n;
        // Frictionless instant closure: valve head is high for 2L/a, low for 2L/a.
        let h0 = p.h_reservoir;
        assert!(res.h_valve[n] > h0 + 10.0);     // t = L/a   -> high
        assert!(res.h_valve[3 * n] < h0 - 10.0); // t = 3L/a  -> low
        assert!(res.h_valve[5 * n] > h0 + 10.0); // t = 5L/a  -> high again
    }

    #[test]
    fn friction_damps_the_waves() {
        let mut p = steel_params();
        p.friction = true;
        p.t_end = 20.0;
        let res = simulate(&p);
        let n = p.n;
        let first = res.h_valve[1..=2 * n].iter().cloned().fold(f64::MIN, f64::max);
        let late = res.h_valve[8 * n..].iter().cloned().fold(f64::MIN, f64::max);
        assert!(late < first, "first {first:.2}, late {late:.2}");
    }

    #[test]
    fn slow_closure_cuts_the_peak() {
        let mut p = steel_params();
        let fast = simulate(&p).h_max[p.n];
        p.t_close = 5.0 * 2.0 * p.length / wave_speed(p.bulk, p.rho, p.diameter, p.wall, p.young, p.c1);
        let slow = simulate(&p).h_max[p.n];
        assert!(slow < fast - 30.0, "fast {fast:.1}, slow {slow:.1}");
    }
}