# File formats

All format code is pure Rust and lives in `crates/io`. Every importer gets a
fuzz target and a real-world corpus before it is ticked in the roadmap.

## Status

| Format | Import | Export | Phase | Notes |
|---|---|---|---|---|
| OpenShape native (`.osdoc`) | ✔ | ✔ | 0/5 | op-log CBOR; spec in NATIVE_FORMAT.md |
| STL ASCII/binary | – | ✔ | 1 | export in mm; import pending (mesh bodies) |
| STEP AP203/214/242 | – | – | 4 | ISO 10303-21 text; PMI where feasible |
| IGES | – | – | 13 | public spec (US PRO/IPO-100) |
| 3MF | – | – | 13 | open spec (3MF Consortium) |
| OBJ, PLY, glTF/GLB, VRML, AMF | – | – | 13 | open specs |
| DXF | – | – | 7 | Autodesk publishes the DXF reference |
| DWG | – | – | 7 | via the Open Design Alliance's publicly documented "Open Design Specification for .dwg files"; clean-room only |
| SVG, PDF | – | – | 7 | export |
| IDF / IDX (ECAD) | – | – | 12 | IDF 3.0/4.0 spec is public; IDX (ProSTEP iViP EDMD) is public |
| Parasolid text `.x_t` | – | – | 13 | only if a clean-room implementation from public documentation is feasible; whether Siemens' published XT format reference permits this must be assessed before any work starts |

## Proprietary formats

These are not implemented by reverse engineering. The legal routes are:

| Format | Status | Legal route |
|---|---|---|
| SolidWorks (`.sldprt`, `.sldasm`) | not supported | export STEP from SolidWorks; or a licensed SDK (Dassault/Spatial 3D InterOp) under a separate commercial agreement |
| CATIA V5 (`.CATPart`, `.CATProduct`) | not supported | export STEP; licensed SDK (Spatial/Dassault) |
| Siemens NX (`.prt`) | not supported | export STEP or Parasolid text (`.x_t`) |
| Autodesk Inventor (`.ipt`, `.iam`) | not supported | export STEP |
| PTC Creo (`.prt`, `.asm`) | not supported | export STEP |
| Parasolid binary (`.x_b`) | not supported | export `.x_t` (text) or STEP |

When a user opens one of these, the app will say which export to use.
