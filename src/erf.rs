//! Error function erf and complementary error function erfc to near f64 precision (R4.1).
//!
//! Accuracy: relative error below 1e-13 wherever the value is above 1e-300; measured at most
//! 3e-15 against a 400-digit reference (tests/erf.rs). No fitted coefficients are used:
//!
//! - |x| < 1: the series erf(x) = (2/√π) e^(−x²) Σ_{n≥0} 2ⁿ x^(2n+1) / (1·3·5···(2n+1)). All its
//!   terms are positive, so there is no cancellation. erfc = 1 − erf, which is at least 0.157
//!   here, so the subtraction loses no significant digits.
//! - x ≥ 1: erfc by the continued fraction
//!   erfc(x) = e^(−x²)/√π · 1/(x + (1/2)/(x + 1/(x + (3/2)/(x + 2/(x + …))))),
//!   evaluated by the modified Lentz method (at most about 190 terms, at x = 1); erf = 1 − erfc.
//! - Negative arguments by symmetry: erf(−x) = −erf(x), erfc(−x) = 2 − erfc(x).
//!
//! The factor e^(−x²) is evaluated as e^(−h)·e^(−l), where x² = h + l is split exactly with a
//! fused multiply-add. A plain `(-x * x).exp()` would carry the rounding error of x² (relative
//! 1.1e-16, absolute 7e-14 at x = 26) into the result as a relative error of the same size.
//!
//! Units: the argument and the results are dimensionless. A NaN argument gives NaN.

const FRAC_2_SQRT_PI: f64 = std::f64::consts::FRAC_2_SQRT_PI;

/// 1/√π. Units: none.
const FRAC_1_SQRT_PI: f64 = 0.5 * std::f64::consts::FRAC_2_SQRT_PI;

/// Above this argument erfc(x) is smaller than the smallest subnormal f64 (erfc(27.3) ≈ 1e-325).
const ERFC_UNDERFLOW: f64 = 28.0;

/// e^(−x²) for finite x with |x| ≤ ERFC_UNDERFLOW, without the rounding error of x².
/// Units: none.
fn exp_neg_square(x: f64) -> f64 {
    let hi = x * x;
    let lo = x.mul_add(x, -hi);
    (-hi).exp() * (-lo).exp()
}

/// erf(x) for 0 ≤ x < 1 by the positive series. Units: none.
fn erf_series(x: f64) -> f64 {
    let x2 = x * x;
    let mut term = x;
    let mut sum = 0.0;
    let mut n = 0.0_f64;
    // The term ratio is 2x²/(2n+3) ≤ 2/3, so at most about 90 terms reach 1e-17 relative.
    for _ in 0..200 {
        sum += term;
        n += 1.0;
        term *= 2.0 * x2 / (2.0 * n + 1.0);
        if term <= sum * 1e-17 {
            break;
        }
    }
    FRAC_2_SQRT_PI * exp_neg_square(x) * sum
}

/// erfc(x) for x ≥ 1 by the continued fraction (modified Lentz method). Units: none.
fn erfc_continued_fraction(x: f64) -> f64 {
    if x > ERFC_UNDERFLOW {
        return 0.0;
    }
    // f = x + a1/(x + a2/(x + ...)) with a_n = n/2; erfc(x) = e^(-x^2)/sqrt(pi) / f.
    const TINY: f64 = 1e-300;
    let mut f = x;
    let mut c = x;
    let mut d = 0.0_f64;
    let mut n = 1.0_f64;
    for _ in 0..1000 {
        let a = 0.5 * n;
        d = x + a * d;
        if d == 0.0 {
            d = TINY;
        }
        c = x + a / c;
        if c == 0.0 {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() <= 1e-16 {
            break;
        }
        n += 1.0;
    }
    exp_neg_square(x) * FRAC_1_SQRT_PI / f
}

/// The error function erf(x) = (2/√π) ∫₀ˣ e^(−t²) dt, to near f64 precision (relative error
/// below 1e-13; see the module documentation).
/// Units: x dimensionless; result dimensionless, in [−1, 1]. NaN gives NaN.
pub fn erf(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let ax = x.abs();
    let value = if ax < 1.0 {
        erf_series(ax)
    } else {
        1.0 - erfc_continued_fraction(ax)
    };
    if x < 0.0 { -value } else { value }
}

/// The complementary error function erfc(x) = 1 − erf(x), to near f64 precision: relative error
/// below 1e-13 wherever erfc(x) > 1e-300 (see the module documentation).
/// Units: x dimensionless; result dimensionless, in [0, 2]. NaN gives NaN.
pub fn erfc(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 1.0 && x > -1.0 {
        return 1.0 - erf_series(x.abs()) * x.signum();
    }
    if x >= 1.0 {
        erfc_continued_fraction(x)
    } else {
        2.0 - erfc_continued_fraction(-x)
    }
}
