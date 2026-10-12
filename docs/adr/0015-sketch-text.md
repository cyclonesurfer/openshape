# ADR 0015 — Text in sketches

Status: accepted (Phase 2)

## Context

The sketcher needs text: a line of characters whose outlines are sketch
curves, so that regions of the letters can be extruded (embossed or
engraved lettering, part numbers). Publicly documented sketchers place text
in a box whose height sets the size and whose corners can be constrained;
the details below are ours (clean-room rule: public standards, public
documentation of behaviour and our own reasoning).

Five questions, and a sixth found while building the tool:

1. **What the document stores**: the string and font, or the geometry.
2. **How text takes part in the solver.**
3. **How a line of text is laid out** (shaping, kerning, size).
4. **Which fonts**, under which licence, from where.
5. **How regions of letters are named** so an extrude keeps its letters
   when the text moves.
6. **How the tool places text and how other tools treat it.**

## Options

### 1. Stored form

* **A. String and font name; lay out at every regeneration.** Smallest
  documents, but the geometry then depends on the font bytes and the
  shaping code of the version that replays the document: a font or shaper
  update would move or reshape letters in existing documents without any
  edit, which the project rules forbid (a model never changes silently).
* **B. String, font name and the laid-out outline**, normalised to the box.
  Replaying a document reproduces its letters exactly in any version;
  layout runs only when someone types. Costs a few kilobytes per text in
  the op log.

### 2. Solver

* **A. Text is a rigid block with a position, angle and size** as new
  unknowns. Needs new unknown kinds and leaves nothing for the user to
  constrain except the block.
* **B. Text owns four corner points** (bottom-left, bottom-right,
  top-right, top-left) held as a rectangle of the text's aspect ratio by
  four linear equations; the tool draws the box as four construction lines
  through the corners. Everything existing (coincidence, dimensions,
  horizontal, drag, transform, mirror, patterns) works on the corners and
  box lines unchanged. 4 freedoms: position, angle, height.

### 3. Layout

* **A. Character-to-glyph map and advance widths only.** Simple, but no
  kerning ("AV" gaps) and nothing beyond simple scripts.
* **B. A full OpenType shaper** (harfrust, the HarfBuzz project's
  MIT-licensed Rust port of its shaping algorithm): kerning, ligatures and
  contextual forms from the font's own tables. Metrics and outlines are
  read with skrifa (MIT OR Apache-2.0, the fontations project). The first
  draft used rustybuzz and ttf-parser; both were declared unmaintained in
  2026 (RUSTSEC-2026-0206 and RUSTSEC-2026-0192, which name harfrust and
  skrifa as their successors), and `cargo deny` refuses unmaintained
  direct dependencies, so the layout moved before anything was committed.

Box height: the font's line height (ascender − descender, the font's own
metrics), or the capital height (descenders and accents would leave the
box).

### 4. Fonts

* **A. System fonts** found at run time: documents would depend on what is
  installed (and option 1B only protects existing letters).
* **B. Fonts bundled in the program**, under a licence that allows it.

### 5. Region names

Region keys (design note 0006 §5) name boundary pieces by the sketch entity
that owns them. All of a text's curves share one entity, so pieces and
vertices could only be told apart by position, and an extrude would lose
its letters when the text moves.

* **A. Name each piece by its entity and its index among the entity's
  curves** (0 for every other entity, so existing keys are unchanged).
* **B. One sketch entity per glyph**: many entities per text, and edits of
  the string would have to add and remove entities.

## Decision

1. **B.** `Geometry::Text` stores the string, the font name, a `mirrored`
   flag, the aspect ratio and the outline: closed contours of pieces in box
   units (u from 0 to the aspect ratio, v from 0 to 1). A piece is a line
   (degree 1) or a clamped uniform B-spline of degree 2 or 3; consecutive
   pieces share their end points exactly, so the contours close exactly
   after any affine map. Editing the string or font lays it out again and
   stores the new outline.
2. **B.** With `u = tl − bl` and `R` the quarter turn clockwise
   (anticlockwise when mirrored): `br − bl = k·R(u)` and
   `tr − tl = br − bl`, where `k` is the aspect ratio. Linear, so they
   never slow the solver; the box lines are ordinary construction lines.
   The tool adds a horizontal constraint to the bottom line, leaving
   position and height free (3 freedoms).
3. **B**, with the box as the **line height**: the box encloses the text,
   so constraints on the box mean what they look like, and the height
   scales like a font size. One line only (a line break is refused); each
   run is shaped in its script's direction, without bidirectional
   reordering of mixed runs. A character the font lacks is refused by
   name, never drawn as a placeholder.
4. **B.** The Liberation fonts 2.1.5 (Sans, Serif and Mono, regular and
   bold; Red Hat and Google digitized data) under the SIL Open Font
   Licence 1.1, which allows bundling and embedding with any software if
   the copyright notice and licence travel with it; the files are
   unmodified, their licence is shipped next to them, NOTICE names them,
   and the `opendatum-text` crate's licence is
   `GPL-3.0-or-later AND OFL-1.1` so cargo-deny checks it. Our code stays
   GPL; the fonts stay OFL. (The UI toolkit's own fonts are already OFL
   and UFL, ADR 0009.)
5. **A.** `EdgeKey` and `VertexKey` gain a part index (serialised only when
   non-zero), so keys of sketches without text are byte-for-byte what they
   were and every text piece has a name that survives moving, turning and
   scaling the box.

### 6. The tool and editing (decided while building it)

* **Placing**: the Text tool takes the characters and font from the
  sketch dialog, then two clicks: a corner of the box, then the height.
  The box grows towards the pointer (a click below and to the left of the
  first makes the first the top-right corner) and always reads left to
  right with a horizontal bottom; its width follows from the aspect
  ratio. The corner at the first click is attached to what it snapped to.
  A string the font cannot draw is refused at the first click, by
  character name and code point.
* **Editing**: selecting a single text shows its characters and font in
  the dialog; Apply (or Enter in the field) lays them out again in one
  undo step. The box keeps its bottom-left corner, height and direction
  and takes the new width (its right-hand corners are moved, the solver
  then honours any constraint on them). Mirror images and pattern copies
  of the text, and texts they in turn were copied to, take the same
  characters and font, so derived texts never disagree with their source
  (the solver keeps an image a copy of its source only when their aspect
  ratios match, and refuses the sketch otherwise).
* **Other tools**: mirror, patterns and transforms carry a text's box
  lines with it; deleting a text deletes its box and corners unless
  something else uses them, while deleting a box line leaves the text.
  Points never snap onto letters; constraints and dimensions refuse the
  text entity itself ("select its corners or box lines instead"); trim,
  extend, split and fillet refuse letters, and box corners for fillets;
  letters cut other curves for trim and extend but never hold a cut end.

## Consequences

* Letters regenerate identically forever; a font or shaper upgrade changes
  only text laid out after it.
* Text is single-line, unstyled per character, in six bundled faces. More
  fonts, multi-line text, mixed-direction reordering and system fonts are
  later work.
* Trim, extend, split, fillet and offset refuse text; mirror, patterns and
  transforms copy it (a mirror image is mirrored text).
* Revisit the box height if users mostly dimension capital height; a
  read-out of the capital height in the dialog would bridge the two.
