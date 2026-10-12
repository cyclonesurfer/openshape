//! Feature parameter dialogs.
//!
//! Each dialog is a list of typed fields parsed with units. In Phase 9 these
//! are generated from the custom-feature language's parameter definitions;
//! for now the built-in feature templates are listed here.

use opendatum_core_model::{
    units::{LengthUnit, Quantity},
    CurveDef, DVec2, DVec3, Feature, FeatureEntry, FeatureId, PlaneRef, ProfileDef, ProfileSource,
    RevolveAxis,
};

/// Type of a dialog field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldKind {
    /// A length (bare numbers are millimetres).
    Length,
    /// An angle (bare numbers are degrees).
    Angle,
    /// A checkbox.
    Flag,
}

/// One editable field.
#[derive(Clone, Debug)]
pub struct Field {
    /// Label.
    pub label: &'static str,
    /// Kind.
    pub kind: FieldKind,
    /// Text being edited (for Flag: "true"/"false").
    pub text: String,
}

/// The feature templates the UI can create.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Template {
    /// Box primitive.
    Box,
    /// Extruded circle.
    Cylinder,
    /// Extruded rectangle with a round hole.
    Plate,
    /// Revolved circle (torus or partial torus).
    Ring,
    /// Revolved spline profile.
    Vase,
}

impl Template {
    /// All templates, in menu order.
    pub const ALL: [Template; 5] = [
        Template::Box,
        Template::Cylinder,
        Template::Plate,
        Template::Ring,
        Template::Vase,
    ];

    /// Menu label.
    pub fn label(self) -> &'static str {
        match self {
            Template::Box => "Box",
            Template::Cylinder => "Cylinder (extrude)",
            Template::Plate => "Plate with hole (extrude)",
            Template::Ring => "Ring (revolve)",
            Template::Vase => "Vase (spline revolve)",
        }
    }
}

/// An open feature dialog.
#[derive(Clone, Debug)]
pub struct FeatureDialog {
    /// Template.
    pub template: Template,
    /// Fields.
    pub fields: Vec<Field>,
    /// The feature being edited, or `None` when creating.
    pub editing: Option<FeatureId>,
    /// Last validation error.
    pub error: Option<String>,
}

fn f(label: &'static str, kind: FieldKind, text: &str) -> Field {
    Field {
        label,
        kind,
        text: text.to_string(),
    }
}

fn mm(v: f64) -> String {
    LengthUnit::Millimetre.format(v)
}

fn deg(v: f64) -> String {
    format!("{} deg", (v.to_degrees() * 1e6).round() / 1e6)
}

impl FeatureDialog {
    /// A dialog for a new feature with default values.
    pub fn new(template: Template) -> Self {
        use FieldKind::*;
        let fields = match template {
            Template::Box => vec![
                f("Length (X)", Length, "40 mm"),
                f("Width (Y)", Length, "30 mm"),
                f("Height (Z)", Length, "20 mm"),
            ],
            Template::Cylinder => vec![
                f("Diameter", Length, "20 mm"),
                f("Height", Length, "30 mm"),
                f("Centre X", Length, "60 mm"),
                f("Centre Y", Length, "15 mm"),
                f("Symmetric", Flag, "false"),
            ],
            Template::Plate => vec![
                f("Width", Length, "60 mm"),
                f("Depth", Length, "40 mm"),
                f("Thickness", Length, "5 mm"),
                f("Hole diameter", Length, "16 mm"),
            ],
            Template::Ring => vec![
                f("Ring radius", Length, "30 mm"),
                f("Tube radius", Length, "6 mm"),
                f("Angle", Angle, "360 deg"),
            ],
            Template::Vase => vec![
                f("Height", Length, "60 mm"),
                f("Base radius", Length, "15 mm"),
            ],
        };
        Self {
            template,
            fields,
            editing: None,
            error: None,
        }
    }

    /// A dialog prefilled from an existing feature, if it matches a template.
    pub fn edit(entry: &FeatureEntry) -> Option<Self> {
        let mut d = match &entry.feature {
            Feature::Box { size, .. } => {
                let mut d = Self::new(Template::Box);
                for (i, v) in [size.x, size.y, size.z].into_iter().enumerate() {
                    d.fields[i].text = mm(v);
                }
                d
            }
            Feature::Extrude {
                profile: ProfileSource::Inline(profile),
                distance,
                symmetric,
            } => match profile.curves.as_slice() {
                [CurveDef::Circle { center, radius }] => {
                    let mut d = Self::new(Template::Cylinder);
                    d.fields[0].text = mm(2.0 * radius);
                    d.fields[1].text = mm(*distance);
                    d.fields[2].text = mm(center.x);
                    d.fields[3].text = mm(center.y);
                    d.fields[4].text = symmetric.to_string();
                    d
                }
                [CurveDef::Line { .. }, _, _, CurveDef::Line { end, .. }, CurveDef::Circle { radius, .. }] =>
                {
                    let mut d = Self::new(Template::Plate);
                    d.fields[0].text = mm(-2.0 * end.x);
                    d.fields[1].text = mm(-2.0 * end.y);
                    d.fields[2].text = mm(*distance);
                    d.fields[3].text = mm(2.0 * radius);
                    d
                }
                _ => return None,
            },
            Feature::Revolve {
                profile: ProfileSource::Inline(profile),
                angle,
                ..
            } => match profile.curves.as_slice() {
                [CurveDef::Circle { center, radius }] => {
                    let mut d = Self::new(Template::Ring);
                    d.fields[0].text = mm(center.x);
                    d.fields[1].text = mm(*radius);
                    d.fields[2].text = deg(*angle);
                    d
                }
                [CurveDef::Spline { through }, ..] => {
                    let mut d = Self::new(Template::Vase);
                    let h = through.last().map(|p| p.y).unwrap_or(0.06);
                    d.fields[0].text = mm(h);
                    d.fields[1].text = mm(through[0].x);
                    d
                }
                _ => return None,
            },
            _ => return None,
        };
        d.editing = Some(entry.id);
        Some(d)
    }

    fn length(&self, i: usize) -> Result<f64, String> {
        Quantity::parse_length(&self.fields[i].text, LengthUnit::Millimetre)
            .map_err(|e| format!("{}: {e}", self.fields[i].label))
    }

    fn angle(&self, i: usize) -> Result<f64, String> {
        Quantity::parse_angle(&self.fields[i].text)
            .map_err(|e| format!("{}: {e}", self.fields[i].label))
    }

    fn flag(&self, i: usize) -> bool {
        self.fields[i].text == "true"
    }

    /// Builds the feature definition from the fields.
    pub fn build(&self) -> Result<Feature, String> {
        Ok(match self.template {
            Template::Box => Feature::Box {
                corner: DVec3::ZERO,
                size: DVec3::new(self.length(0)?, self.length(1)?, self.length(2)?),
            },
            Template::Cylinder => Feature::Extrude {
                profile: ProfileSource::Inline(ProfileDef {
                    plane: PlaneRef::Top,
                    curves: vec![CurveDef::Circle {
                        center: DVec2::new(self.length(2)?, self.length(3)?),
                        radius: self.length(0)? / 2.0,
                    }],
                }),
                distance: self.length(1)?,
                symmetric: self.flag(4),
            },
            Template::Plate => {
                let (w, d) = (self.length(0)? / 2.0, self.length(1)? / 2.0);
                let corners = [
                    DVec2::new(-w, -d),
                    DVec2::new(w, -d),
                    DVec2::new(w, d),
                    DVec2::new(-w, d),
                ];
                let mut curves: Vec<CurveDef> = (0..4)
                    .map(|i| CurveDef::Line {
                        start: corners[i],
                        end: corners[(i + 1) % 4],
                    })
                    .collect();
                curves.push(CurveDef::Circle {
                    center: DVec2::ZERO,
                    radius: self.length(3)? / 2.0,
                });
                Feature::Extrude {
                    profile: ProfileSource::Inline(ProfileDef {
                        plane: PlaneRef::Top,
                        curves,
                    }),
                    distance: self.length(2)?,
                    symmetric: false,
                }
            }
            Template::Ring => Feature::Revolve {
                profile: ProfileSource::Inline(ProfileDef {
                    plane: PlaneRef::Front,
                    curves: vec![CurveDef::Circle {
                        center: DVec2::new(self.length(0)?, 0.0),
                        radius: self.length(1)?,
                    }],
                }),
                axis: RevolveAxis::Inline {
                    origin: DVec2::ZERO,
                    dir: DVec2::Y,
                },
                angle: self.angle(2)?,
            },
            Template::Vase => {
                let (h, r) = (self.length(0)?, self.length(1)?);
                let through = vec![
                    DVec2::new(r, 0.0),
                    DVec2::new(1.6 * r, 0.3 * h),
                    DVec2::new(0.8 * r, 0.7 * h),
                    DVec2::new(1.1 * r, h),
                ];
                Feature::Revolve {
                    profile: ProfileSource::Inline(ProfileDef {
                        plane: PlaneRef::Front,
                        curves: vec![
                            CurveDef::Spline { through },
                            CurveDef::Line {
                                start: DVec2::new(1.1 * r, h),
                                end: DVec2::new(0.0, h),
                            },
                            CurveDef::Line {
                                start: DVec2::new(0.0, h),
                                end: DVec2::ZERO,
                            },
                            CurveDef::Line {
                                start: DVec2::ZERO,
                                end: DVec2::new(r, 0.0),
                            },
                        ],
                    }),
                    axis: RevolveAxis::Inline {
                        origin: DVec2::ZERO,
                        dir: DVec2::Y,
                    },
                    angle: std::f64::consts::TAU,
                }
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use opendatum_core_model::{OpId, ReplicaId};

    #[test]
    fn every_template_builds_and_round_trips_through_edit() {
        for t in Template::ALL {
            let d = FeatureDialog::new(t);
            let feature = d.build().unwrap();
            let entry = FeatureEntry {
                id: FeatureId(OpId {
                    lamport: 1,
                    replica: ReplicaId(1),
                }),
                name: "x".into(),
                suppressed: false,
                feature: feature.clone(),
            };
            let e = FeatureDialog::edit(&entry).expect("template is editable");
            assert_eq!(e.template, t);
            let rebuilt = e.build().unwrap();
            assert_eq!(rebuilt, feature, "{t:?}");
        }
    }

    #[test]
    fn bad_units_are_reported() {
        let mut d = FeatureDialog::new(Template::Box);
        d.fields[0].text = "10 deg".into();
        assert!(d.build().unwrap_err().contains("Length (X)"));
    }
}
