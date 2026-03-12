/// Prior distributions for host galaxy association.

/// Uniform prior on fractional offset up to a maximum value.
///
/// P(d_FR) = 1/max_offset for d_FR ∈ [0, max_offset], else 0.
pub fn offset_prior(fractional_offset: f64, max_offset: f64) -> f64 {
    if fractional_offset >= 0.0 && fractional_offset <= max_offset {
        1.0 / max_offset
    } else {
        0.0
    }
}

/// Prior probability that the true host is outside the search radius.
///
/// Conservative default: small but non-zero.
pub fn p_outside(n_candidates: usize) -> f64 {
    if n_candidates == 0 {
        0.5
    } else {
        0.01
    }
}

/// Prior probability that the true host is too faint to be in the catalog.
///
/// Depends on survey depth; for deep surveys (e.g., DECaLS) this is small.
pub fn p_unobserved() -> f64 {
    0.01
}

/// Prior probability that the transient is hostless (e.g., intracluster).
pub fn p_hostless() -> f64 {
    0.005
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_offset_prior_in_range() {
        assert_relative_eq!(offset_prior(3.0, 10.0), 0.1);
    }

    #[test]
    fn test_offset_prior_at_boundary() {
        assert_relative_eq!(offset_prior(10.0, 10.0), 0.1);
        assert_relative_eq!(offset_prior(0.0, 10.0), 0.1);
    }

    #[test]
    fn test_offset_prior_out_of_range() {
        assert_relative_eq!(offset_prior(11.0, 10.0), 0.0);
        assert_relative_eq!(offset_prior(-1.0, 10.0), 0.0);
    }

    #[test]
    fn test_p_outside_no_candidates() {
        assert_relative_eq!(p_outside(0), 0.5);
    }

    #[test]
    fn test_p_outside_with_candidates() {
        assert_relative_eq!(p_outside(5), 0.01);
    }
}
