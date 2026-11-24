const G: f64 = 9.81;

/// Joukowsky pressure rise from an instant velocity change.
/// a: wave speed (m/s), v: velocity change (m/s). Returns head rise in metres of water.
pub fn joukowsky(a: f64, v: f64) -> f64 {
    a * v / G
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn joukowsky_basic() {
        let h = joukowsky(1000.0, 4.0);
        assert!((h - 407.75).abs() < 0.01);
    }
}