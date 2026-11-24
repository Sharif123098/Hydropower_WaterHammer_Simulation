use hydro_core::{joukowsky, simulate, Params};

fn main() {
    let p = Params {
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
    };

    let res = simulate(&p);
    let area = std::f64::consts::PI * p.diameter * p.diameter / 4.0;
    let v0 = p.q0 / area;

    println!("wave speed a        = {:.1} m/s", res.a);
    println!("time step dt        = {:.5} s", res.dt);
    println!("steady velocity v0  = {:.3} m/s", v0);
    println!("Joukowsky rise      = {:.1} m", joukowsky(res.a, v0));
    println!("MOC max valve head  = {:.1} m", res.h_max[p.n]);
    println!("MOC min valve head  = {:.1} m", res.h_min[p.n]);
    println!("period 4L/a         = {:.3} s", 4.0 * p.length / res.a);

    println!("\n  t (s)   H valve (m)   Q valve (m3/s)");
    for i in (0..res.time.len()).step_by(20) {
        println!("{:7.3}   {:10.1}   {:10.3}", res.time[i], res.h_valve[i], res.q_valve[i]);
    }
}