//! Physical quantities with explicit units.
//!
//! Values are stored in SI base units (metres, radians). Parsing accepts a
//! number followed by an optional unit; a bare number takes the caller's
//! default unit. Full expressions arrive with variables in Phase 8.

use serde::{Deserialize, Serialize};

/// Physical dimension of a quantity.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Dimension {
    /// Length (metres).
    Length,
    /// Plane angle (radians).
    Angle,
    /// Pure number.
    Dimensionless,
}

/// A value with a dimension, stored in SI units.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Quantity {
    /// Value in SI units (metres or radians).
    pub si: f64,
    /// Dimension.
    pub dim: Dimension,
}

/// A length unit.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum LengthUnit {
    /// Millimetre.
    #[default]
    Millimetre,
    /// Centimetre.
    Centimetre,
    /// Metre.
    Metre,
    /// Inch.
    Inch,
    /// Foot.
    Foot,
}

impl LengthUnit {
    /// Metres per unit.
    pub fn metres(self) -> f64 {
        match self {
            LengthUnit::Millimetre => 1e-3,
            LengthUnit::Centimetre => 1e-2,
            LengthUnit::Metre => 1.0,
            LengthUnit::Inch => 0.0254,
            LengthUnit::Foot => 0.3048,
        }
    }

    /// Unit symbol.
    pub fn symbol(self) -> &'static str {
        match self {
            LengthUnit::Millimetre => "mm",
            LengthUnit::Centimetre => "cm",
            LengthUnit::Metre => "m",
            LengthUnit::Inch => "in",
            LengthUnit::Foot => "ft",
        }
    }

    /// Formats a length in metres in this unit.
    pub fn format(self, metres: f64) -> String {
        format!("{} {}", trim(metres / self.metres()), self.symbol())
    }
}

fn trim(x: f64) -> String {
    let s = format!("{x:.6}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" {
        "0".into()
    } else {
        s.to_string()
    }
}

/// Errors parsing a quantity.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum UnitError {
    /// The numeric part is missing or malformed.
    #[error("not a number: {0:?}")]
    BadNumber(String),
    /// The unit is not recognised.
    #[error("unknown unit: {0:?}")]
    UnknownUnit(String),
    /// The unit has the wrong dimension for this field.
    #[error("expected a {expected:?} but got a {got:?}")]
    WrongDimension {
        /// Expected dimension.
        expected: Dimension,
        /// Parsed dimension.
        got: Dimension,
    },
}

fn unit_factor(u: &str) -> Option<(f64, Dimension)> {
    Some(match u {
        "mm" => (1e-3, Dimension::Length),
        "cm" => (1e-2, Dimension::Length),
        "m" => (1.0, Dimension::Length),
        "um" | "µm" => (1e-6, Dimension::Length),
        "in" | "\"" => (0.0254, Dimension::Length),
        "ft" | "'" => (0.3048, Dimension::Length),
        "deg" | "°" => (std::f64::consts::PI / 180.0, Dimension::Angle),
        "rad" => (1.0, Dimension::Angle),
        _ => return None,
    })
}

impl Quantity {
    /// A length in metres.
    pub fn metres(si: f64) -> Self {
        Self {
            si,
            dim: Dimension::Length,
        }
    }

    /// An angle in radians.
    pub fn radians(si: f64) -> Self {
        Self {
            si,
            dim: Dimension::Angle,
        }
    }

    /// Parses `"<number> [unit]"`. A bare number is interpreted with
    /// `default_factor` (SI per unit) and `expected` dimension. The result
    /// must have dimension `expected`.
    pub fn parse(text: &str, expected: Dimension, default_factor: f64) -> Result<Self, UnitError> {
        let t = text.trim();
        // The longest prefix that parses as a number is the value.
        let (value, unit) = (1..=t.len())
            .rev()
            .filter(|&i| t.is_char_boundary(i))
            .find_map(|i| t[..i].trim().parse::<f64>().ok().map(|v| (v, &t[i..])))
            .ok_or_else(|| UnitError::BadNumber(text.to_string()))?;
        if !value.is_finite() {
            return Err(UnitError::BadNumber(text.to_string()));
        }
        let unit = unit.trim();
        let (factor, dim) = if unit.is_empty() {
            (default_factor, expected)
        } else {
            unit_factor(unit).ok_or_else(|| UnitError::UnknownUnit(unit.to_string()))?
        };
        if dim != expected {
            return Err(UnitError::WrongDimension { expected, got: dim });
        }
        Ok(Self {
            si: value * factor,
            dim,
        })
    }

    /// Parses a length; bare numbers are in `default`.
    pub fn parse_length(text: &str, default: LengthUnit) -> Result<f64, UnitError> {
        Self::parse(text, Dimension::Length, default.metres()).map(|q| q.si)
    }

    /// Parses an angle; bare numbers are degrees.
    pub fn parse_angle(text: &str) -> Result<f64, UnitError> {
        Self::parse(text, Dimension::Angle, std::f64::consts::PI / 180.0).map(|q| q.si)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_lengths() {
        let mm = LengthUnit::Millimetre;
        assert!((Quantity::parse_length("25", mm).unwrap() - 0.025).abs() < 1e-15);
        assert!((Quantity::parse_length("25 mm", mm).unwrap() - 0.025).abs() < 1e-15);
        assert!((Quantity::parse_length("1in", mm).unwrap() - 0.0254).abs() < 1e-15);
        assert!((Quantity::parse_length("-2.5 cm", mm).unwrap() + 0.025).abs() < 1e-15);
        assert!((Quantity::parse_length("1e3 um", mm).unwrap() - 1e-3).abs() < 1e-15);
        assert!((Quantity::parse_length("1.5e-1m", mm).unwrap() - 0.15).abs() < 1e-15);
    }

    #[test]
    fn parses_angles_and_rejects_bad_input() {
        assert!((Quantity::parse_angle("90").unwrap() - std::f64::consts::FRAC_PI_2).abs() < 1e-15);
        assert!((Quantity::parse_angle("1 rad").unwrap() - 1.0).abs() < 1e-15);
        assert!(matches!(
            Quantity::parse_angle("10 mm"),
            Err(UnitError::WrongDimension { .. })
        ));
        assert!(matches!(
            Quantity::parse_length("abc", LengthUnit::Millimetre),
            Err(UnitError::BadNumber(_))
        ));
        assert!(matches!(
            Quantity::parse_length("3 furlongs", LengthUnit::Millimetre),
            Err(UnitError::UnknownUnit(_))
        ));
    }

    #[test]
    fn formats() {
        assert_eq!(LengthUnit::Millimetre.format(0.025), "25 mm");
        assert_eq!(LengthUnit::Inch.format(0.0254), "1 in");
    }
}
