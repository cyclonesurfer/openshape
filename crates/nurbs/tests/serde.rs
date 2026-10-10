//! Serialization round trips and validation on deserialization.

use glam::DVec3;
use openshape_nurbs::{NurbsCurve, NurbsSurface};

#[test]
fn curve_round_trip() {
    let c = NurbsCurve::circular_arc(DVec3::new(1., 2., 3.), DVec3::X, DVec3::Y, 2.0, 0.1, 4.0)
        .unwrap();
    let json = serde_json::to_string(&c).unwrap();
    let back: NurbsCurve = serde_json::from_str(&json).unwrap();
    assert_eq!(back, c);
}

#[test]
fn surface_round_trip() {
    let profile = NurbsCurve::line(DVec3::new(1., 0., 0.), DVec3::new(2., 0., 1.)).unwrap();
    let s = NurbsSurface::revolution(&profile, DVec3::ZERO, DVec3::Z, 2.0).unwrap();
    let json = serde_json::to_string(&s).unwrap();
    let back: NurbsSurface = serde_json::from_str(&json).unwrap();
    assert_eq!(back, s);
}

#[test]
fn invalid_curve_data_is_rejected() {
    let ok =
        r#"{"degree":1,"knots":[0.0,0.0,1.0,1.0],"points":[[0.0,0.0,0.0,1.0],[1.0,0.0,0.0,1.0]]}"#;
    assert!(serde_json::from_str::<NurbsCurve>(ok).is_ok());
    let bad = [
        // Negative weight.
        r#"{"degree":1,"knots":[0.0,0.0,1.0,1.0],"points":[[0.0,0.0,0.0,1.0],[1.0,0.0,0.0,-1.0]]}"#,
        // Wrong knot count.
        r#"{"degree":1,"knots":[0.0,0.0,1.0],"points":[[0.0,0.0,0.0,1.0],[1.0,0.0,0.0,1.0]]}"#,
        // Decreasing knots.
        r#"{"degree":1,"knots":[0.0,1.0,0.5,1.0],"points":[[0.0,0.0,0.0,1.0],[1.0,0.0,0.0,1.0]]}"#,
        // Degree zero.
        r#"{"degree":0,"knots":[0.0,1.0],"points":[[0.0,0.0,0.0,1.0]]}"#,
    ];
    for json in bad {
        assert!(serde_json::from_str::<NurbsCurve>(json).is_err(), "{json}");
    }
}

#[test]
fn invalid_surface_data_is_rejected() {
    let profile = NurbsCurve::line(DVec3::new(1., 0., 0.), DVec3::new(2., 0., 1.)).unwrap();
    let s = NurbsSurface::extrusion(&profile, DVec3::Y).unwrap();
    let mut value = serde_json::to_value(&s).unwrap();
    value["count_v"] = serde_json::json!(3);
    assert!(serde_json::from_value::<NurbsSurface>(value).is_err());
    let mut value = serde_json::to_value(&s).unwrap();
    value["points"][0][3] = serde_json::json!(0.0);
    assert!(serde_json::from_value::<NurbsSurface>(value).is_err());
}
