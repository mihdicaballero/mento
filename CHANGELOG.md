# Changelog

All notable changes to mento are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). From 1.0.0
the public API is stable: breaking changes need a major release. Before 1.0.0 it could
change between minor releases, which is what ADR-0003 documented at the time.

Releases before 0.5.0 were not tracked in a changelog; the entries below were reconstructed
from the release history and are summaries rather than complete lists.

## [Unreleased]

### Added

- **A design keeps its alternatives.** `beam.flexure_design.bottom.options` and
  `.top.options` are tuples of `RebarOption` — `layers` (the same `RebarLayer` objects the
  applied bars are read as), `A_s` and the `functional` the search ranks by — and
  `beam.shear_design.options` a tuple of `StirrupOption` (`n_stirrups`, `d_b`, `s_l`, `s_w`,
  `A_v`, `functional`). They come best first, and `options[0]` is the layout the section
  carries: the search's ranked table used to be dropped after its first row, and is now
  kept with the mechanical cover the design finished on. The stirrup alternatives are one
  layout per other bar diameter the code offers, lighter and heavier alike, each sized at
  its own depth; `functional` says what each one adds in steel, the excess of `A_v` over
  its own requirement plus one per extra stirrup. An alternative, longitudinal or
  transverse, is offered only if the finished section passes with it, and carries the
  `DCR` it passes at. Which layout is built is unchanged: fewest stirrups first, least
  steel among those. How many options are kept is at most
  `BeamSettings(design_options=3)`. A check reports none, a footing offers none, and a
  face changed by hand after the design drops its own.

- **Detailing limits are reported as data.** `beam.warnings` and `node.warnings` are tuples
  of `DesignWarning` with a stable `code` (`As_below_min`, `As_above_max`,
  `clear_spacing_below_min`, `bar_spacing_below_min`, `bar_spacing_exceeds_max`,
  `bars_do_not_fit`, `As_below_required`, `stirrups_required`, `Av_below_min`,
  `stirrup_spacing_exceeds_max`, `stirrup_spacing_exceeds_compression_support`,
  `stirrup_diameter_below_compression_support`, `shear_exceeds_section_limit`), a
  `message` in the language of `mento.set_language`, the `values` it quotes as quantities,
  the `face` and the `combinations` it occurs under. They are the limit rows the detailed
  reports mark with ❌, which until now were only text. A warning does not change a DCR.
  There is no code for the stirrup diameter alone: neither code states a minimum for a
  stirrup placed for shear only; the 10 mm (6 mm under CIRSOC) the design starts from is
  the bottom of the catalogue. The minimum §9.7.6.4.2 does state is for the stirrups that
  brace compression bars, which has its own code.

- **A wall's mesh and shear results are readable as data.** `wall.mesh` is a `WallMesh`
  with `.horizontal` and `.vertical` `MeshDirection`s (`d_b`, `s`, `rho`, `n_curtains`,
  `A_s` per unit length); `wall.shear_checks` holds one `WallShearCheck` per combination
  (`V_u`, `V_capacity` = ØVn, `V_max` = ØVn,max, `rho_t`, `rho_t_req`, `rho_l`, `rho_l_min`,
  `s_h_max`, `s_v_max`, `DCR`), and `wall.shear_design` a `WallShearDesign` with the mesh and
  the envelope. `wall.shear_check_results(forces)` returns the per-combination results
  without building the report, and `wall.warnings` reports the mesh limits missed with the
  new codes `mesh_ratio_below_min` and `mesh_spacing_exceeds_max`, plus
  `shear_exceeds_section_limit` against ØVn,max.

- **The minimum a face has to meet, `A_s_min_eff`.** `FlexureFaceCheck` and
  `FlexureFaceDesign` carry it next to `A_s_min`: under ACI 318-19 and CIRSOC 201-25 it is
  the smaller of `A_s_min` and 4/3·`A_s_calc`, the relief of §9.6.1.3; under EN 1992-1-1,
  and on a member on the ground, it is `A_s_min` itself. A designed face can sit below
  `A_s_min` and comply; it cannot sit below `A_s_min_eff`.

- **The maximum a face has to meet, `A_s_max_eff`.** `FlexureFaceCheck` and
  `FlexureFaceDesign` carry it next to `A_s_max`: under ACI 318-19 and CIRSOC 201-25 it is
  the tension steel the face can take and stay tension-controlled with the compression
  steel of the other face, `A_s_max + A_s'·f_s'/f_y` (§9.3.3.1, Table 21.2.2); under
  EN 1992-1-1 it is `A_s_max` itself. The detailed report's *Max.* column prints it, and a
  doubly reinforced face is checked against it instead of being excused (`✅ D.R.` now
  means "past `A_s_max`, within `A_s_max_eff`").

- **The release workflow publishes a test count.** After uploading to PyPI it attaches
  `stats.json` (`{"tests": N}`) to the GitHub Release and sends a `mento-release`
  `repository_dispatch` to `mihdicaballero/mento-web`. `N` counts the tests marked
  `published_example`: the 31 whose expected numbers come from a document outside mento —
  the Calcpad sheets of the ACI and EN beam cases, The Concrete Centre's Eurocode 2 guide,
  and rows 27–49 of the ETABS/spreadsheet cross-check of the flexure suite. Each says where
  in that document the number is, in a `Source:` paragraph of its docstring, and
  `tests/test_published_examples.py` fails when one does not. 21 tests that pinned mento's
  own output — the 4/3 rule of §9.6.1.3 and the geometric floor where ETABS applies
  A_s,min, re-baselined EN checks, the slab tests whose Calcpad sheet still carries the
  beam minimum, and eight shear tests that cite no source — are not marked. With nothing
  marked the workflow publishes 0 and warns; it never publishes the size of the whole
  suite. The dispatch needs the repository secret `MENTO_WEB_DISPATCH_TOKEN`, a token with
  `Contents: read and write` on mento-web; without it the step warns and the release goes
  on.

- **The crack-control cap of Table 24.3.2 on the bars nearest the tension face.**
  ACI 318-19 / CIRSOC 201-25 §7.7.2.2 (slabs) and §9.7.2.2 (beams) send them to
  Table 24.3.2, s ≤ min(380·(280/f_s) − 2.5·c_c, 300·(280/f_s)), with f_s = (2/3)·f_y as
  §24.3.2.1 permits: 300 mm with Grade 420 and 25 mm of cover. Slabs take it beside
  §7.7.2.3 in the design and the check (an ACI 12 cm slab under 7.5 kN·m goes from Ø10/34
  to Ø10/30; a footing with 50 mm of cover from 300 to 255 mm); beams get two
  `Maximum spacing` rows in the flexure report and `bar_spacing_exceeds_max` on the
  tension face (a 60×50 beam with 2Ø25 is 505 mm past a 292.5 mm cap), and the bar search
  lays a beam out within the cap (a 60×50 beam under 150 kN·m goes from 2Ø20, 510 mm
  apart, to 3Ø20 at 255 mm). EN 1992-1-1 is unchanged: its slab cap stays §9.3.1.1(3) and
  its crack control needs the service stress mento does not have. New registry hook
  `max_bar_spacing_tension`.

- **Registry hooks for the alternatives and the compression bars.**
  `DesignCode.flexure_admissible` (tension-controlled under ACI 318-19 / CIRSOC 201-25,
  within the 4 % of §9.2.1.1(3) under EN 1992-1-1: the rule each flexure design already
  applied, now data of the registry), `DesignCode.stirrup_compression_support` and
  `DesignCode.min_stirrup_for_compression_bar`; `CompressionSupport` in
  `mento.codes.check_state`; `DCR` on `RebarOption` and `StirrupOption`.

### Fixed

- **A clear spacing equal to its limit is no longer lost to rounding.** The effective width
  came out of the unit conversions a hair short (12 cm − 2·(25 + 8) mm = 53.99999999999999
  mm), so two Ø12 bars sat 29.999999999999993 mm apart against the 30 mm vibrator limit and
  the selector dropped them: a 12×30 beam with 40 kN·m was designed with 4Ø10 = 3.14 cm²
  (DCR 1.53) instead of the 4Ø12 = 4.52 cm² that fit (DCR 1.13).

- **An over-reinforced section is no longer credited with strength it does not have.**
  Under ACI 318-19 and CIRSOC 201-25, a face past `A_s_max` had its tension steel cut back
  to `A_s_max + A_s'·f_s'/f_y` and kept at φ = 0.90, whatever the strain it reached. That
  overstated every section that is not tension-controlled: a CIRSOC 25×40 with 3Ø32 over
  4Ø16 at −206.90 kN·m read φMn = 213.4 kN·m (DCR 0.97) against the 196.2 it has
  (ε_t = 0.00252, φ = 0.685, DCR 1.055). The capacity now comes from strain
  compatibility with every bar at the stress its strain gives it and φ from ε_t
  (Table 21.2.2); nothing is capped. A tension-controlled doubly reinforced section reads
  exactly what it did. In a sweep of 960 ACI designs the old check had passed 32 whose real
  DCR reached 1.07.

- **A section that is not tension-controlled is reported as such.** §9.3.3.1 does not
  allow it in a beam, but a doubly reinforced section was excused from the maximum
  altogether. The check now holds the face in tension to `A_s_max_eff` and warns
  `As_above_max` past it; 195 of the 960 designs of that sweep had passed without it.

- **A flexure design passes its own check.** The design accepts a layout only if it
  carries the moment and keeps within the code's limits — tension-controlled under
  ACI 318-19 / CIRSOC 201-25, within the 4 % of §9.2.1.1(3) on both faces under
  EN 1992-1-1 — and otherwise says it found none. To get there:
  - where no bars land between A_s,req and A_s,max (a narrow web jumps from 4Ø12 to 2Ø20)
    it takes the smallest layout past A_s,req instead of one below it;
  - it sizes the compression steel for the tension steel **placed**, not the area asked
    for, which the bars round up;
  - it chooses compression steel by the depth each candidate sits at, not by area alone:
    two layers of thin bars sit deeper, reach less stress and need more — a 15×30 section
    at 83 kN·m ran away to 2Ø25 + 2Ø25 on top and failed, where 2Ø25 in one layer works;
  - when its iteration ends on a layout that fails, it tries every pair of the layouts
    visited on the two faces, the faces being coupled, instead of repairing one face
    against whatever the other held.

  In the 960-design sweep every ACI design now either passes its check or warns
  `As_below_required`, and a brute-force search finds no valid layout for any of the
  ones that warn. `tests/test_flexure_design_properties.py` holds the property.

- **The vibrator size only spaces the top bars in a design.** The check and the warnings
  already held the bottom face to 25 mm (1 in.) and the bar diameter, but the bar search
  applied the 30 mm of `vibrator_size` to both faces, so it discarded bottom layouts the
  check would pass: a 15 cm web could not take three bars a layer. The bottom of a 15×30
  beam at 40 kN·m is now 2Ø12 + 1Ø10 over 2Ø10 + 1Ø10 = 5.40 cm² (DCR 0.92). The stirrup
  the flexure design assumes stays the `stirrup_diameter_ini` of the settings, even when
  the shear design later settles on a thinner one.

- **A design that cannot reach the steel it needs says so.** When no layout that fits the
  width carries the moment (tension-controlled, under ACI 318-19 / CIRSOC 201-25), the
  design leaves the closest it found — 4Ø12 in a 12 cm web asked for 5.18 cm² — and only
  the DCR used to show it. It now warns `As_below_required` on the face that fell short,
  quoting `A_s` and `A_s_req`, for as long as the face carries what the design left.

- **A moment no tension steel alone can carry is designed doubly reinforced.** Under ACI
  318-19 / CIRSOC 201-25, when the equation for singly reinforced steel had no solution
  (a negative discriminant) the requirement was set to A_s,max and the compression-steel
  branch never ran: a 15×30 beam at −80 kN·m reported A_s,req = 5.6 cm² on top and none
  below. It now reports the couple, 10.94 cm² above and 8.48 cm² below, and A_s,req no
  longer drops as the moment grows. Where the compression bars sit too close to the
  neutral axis to help (f_s' − 0.85·f'c ≤ 0, a shallow section) it no longer asks for a
  negative area of them — −39 cm² — but for none.

- **Under ACI 318-19 / CIRSOC 201-25, `As_above_max` is only read on the face the
  combination puts in tension.** The bars a negative moment asks for on the bottom are
  compression steel, and a combination with no moment pulls neither face, yet both were
  held to A_s,max and warned. The detailed report skips the same check; it still prints the
  limit. EN 1992-1-1 keeps both faces: its 4 % caps "tension or compression
  reinforcement".

- **The flexure design loop no longer stops when only one face has settled.** It took a
  repeated layout on either face as a limit cycle, which is also what a face that has
  converged does while the other is still moving; it now waits for the pair of layouts to
  repeat.

- **`As_below_min` no longer fires on a face the 4/3 relief covers, and fires when it does
  not.** The warning read the flag that says the *requirement* adopted 4/3·A_s_calc, so it
  stayed silent on bars checked by hand between A_s_calc and 4/3·A_s_calc. It now compares
  the steel provided against `A_s_min_eff`, and quotes that minimum by its name,
  `A_s,min,eff`; `values` carries `A_s_min` as the clause writes it and `A_s_min_eff`, the
  minimum the face has to meet after the 4/3 relief of §9.6.1.3, which is the one it is
  compared against. The detailed flexure report marks a face the relief covers
  `✅ 9.6.1.3` instead of the bare article number.

- **`ShearWall` no longer answers as the beam it inherits from.** `shear_design`,
  `shear_checks` and `shear_check_results` described stirrups the wall does not have and
  DCRs that were not its own (7.37 and 15.6 for a wall whose check gives 0.58); they now
  return the wall's results. `reinforcement`, `flexure_design`, `flexure_checks` and
  `flexure_check_results` raise `mento.shear_wall.NotABeamError` pointing to `wall.mesh`
  — an `AttributeError` that is also a `NotImplementedError`, so
  `hasattr(wall, "reinforcement")` is False and `getattr(wall, "reinforcement", None)`
  takes its default, and a loop over beams and walls can ask for the member.

- **The stirrup spacing limit halves at 0.33√f'c·bw·d, not 0.083.** ACI 318-19 /
  CIRSOC 201-25 Table 9.7.6.2.2 halves the limits once the nominal `Vs,req = (Vu − φVc)/φ`
  passes `0.33·√f'c·bw·d` (4·√f'c·bw·d in psi), with no λ. 1.2.0 compared `Vu − φVc`
  against `0.083·λ·√f'c·bw·d`, the threshold of §9.6.3.1, so an ordinary beam was held to
  d/4: a 20×60 beam of f'c = 25 MPa under Vu = 120 kN now gets d/2 in the check and the
  design alike. (#161)

- **`design()` gives the same result every time.** It read the stirrup diameter the
  previous run had left on the section, so a second run on the same beam could detail
  1eØ10/27 after 1eØ10/28. A design now starts from the same state however it is called
  — the stirrup the settings assume and the placeholder bars — and every stirrup
  diameter is sized against the demand and the spacing limits of the section it would
  make, so the applied layout is read at the depth the finished beam has: the first
  design used to pick 28 cm at the depth of the 8 mm starter stirrup, past the 27.95 cm
  limit of the beam once its Ø10 was placed.

- **The stirrup design is sized at the depth of its own diameter.** Every bar the code
  offers is sized against the `A_v,req` and the Table 9.7.6.2.2 row read with that bar on
  the section, so the applied layout passes its own check by construction. The design
  used to repeat itself with the diameter it chose and, when two diameters kept trading
  places (8 → 10 → 8), applied a row sized at the other's depth: a CIRSOC 30×40 with
  Vu = 180 kN got 1eØ10/15 against a limit of 8.9 cm, and a 25×50 with Vu = 300 kN got
  1eØ10/10 with A_v = 15.71 cm²/m against 15.78 required. The outer loop is gone. An
  alternative's `functional` is measured against its own demand.

- **Stirrup alternatives are built and checked on the finished section.**
  `shear_design.options[1:]` are one layout per other bar diameter, lighter and heavier
  alike; each is built on the section, checked for shear and flexure under every
  combination, and kept only if it passes. `StirrupOption.DCR` says at what ratio;
  `options[0]` is the applied layout with its own. They were described as "the same cage
  in a heavier bar" and never built: 215 of 1152 failed when they were (an ACI 25×50 under
  350 kN offered 1eØ16/11 at DCR 1.007).

- **Longitudinal alternatives are verified on the finished beam.**
  `flexure_design.<face>.options[1:]` are rebuilt after the design ends — the stirrups it
  chose, the other face as applied — and kept only if the beam carries both moments with
  them within the code's reinforcement limits (tension-controlled under ACI 318-19 /
  CIRSOC 201-25, 4 % under EN 1992-1-1) and the bars fit; `RebarOption.DCR` says at what
  ratio. They came from the last Picard iteration with the 8 mm starter stirrup, ranked
  by area, and some failed when placed (ACI 20×60, Mu = 80 kN·m: a 3.93 cm² alternative
  at DCR 1.023; EN 20×60 at 400 kN·m: 2Ø32 on top took the bottom face to 1.001).

- **A footing offers no alternatives.** Its mat is chosen as a whole; the per-face rows
  the mat replaced were listed after it and were not alternatives to it (509 of 540
  footings in a sweep, some failing: an EN 1 m × 15 cm footing offered Ø12/15 at 1.065).

- **Stirrups that brace compression reinforcement.** A section whose flexure relies on
  compression steel now holds its stirrups to ACI 318-19 / CIRSOC 201-25 §9.7.6.4: spacing
  at most the least of 16 d_b of the compression bar, 48 d_b of the stirrup and the least
  dimension of the beam (§9.7.6.4.3), and a stirrup no thinner than §9.7.6.4.2 (ACI: No. 10
  up to a No. 32 bar, No. 13 above) / CIRSOC Tabla 9.7.6.4.2 (6 to 12 mm by bar) require.
  The shear design applies both; the check reports
  `stirrup_spacing_exceeds_compression_support` and
  `stirrup_diameter_below_compression_support`, read off the section and not a
  combination. mento's own designs broke both (ACI 20×50, Mu = 260 kN·m: 1eØ10/21 against
  200 mm; CIRSOC 20×40, Mu = 200 kN·m: Ø6 where the table asks Ø8). EN 1992-1-1
  §9.2.1.2(3) (15·φ) is not applied yet.

- **A slab strip is counted as bars per metre.** `Ø10/12` on a metre is 8.33 bars,
  6.54 cm²/m, not the 9 bars (7.07 cm²) that covered the strip. The design rounds the
  spacing down to the whole centimetre, so the strip never carries less steel than the
  search chose; a CIRSOC 100×25 strip under 80 kN·m used to be reported at DCR 0.995 with a
  real DCR of 1.071 per metre, and a 100×30 strip designed to its 5.40 cm² minimum
  carried 5.24. The footing's mat search counted the same way and is read the same.
  `RebarLayer.n` and the `n_bars` properties are floats now (still whole on beams); a
  section drawing shows the whole bars that cover the strip.

- **The slab minimum belongs to the face in tension.** ACI 318-19 §7.6.1.1 / CIRSOC 201-25
  §7.6.1 is a flexural minimum (R7.6.1.1 / C 7.6.1): under a combination with no moment
  both faces report `A_s,min = 0`, as a beam does, instead of 0.0018·Ag on each and an
  `As_below_min` on a face nothing pulls. The design still places the studio's 1.8‰ on a
  face with no moment.

- **ρl,min of a wall reads the horizontal mesh it carries.** ACI 318-19 / CIRSOC 201-25
  §11.6.2(a) put the plain ρt in Eq. (11.6.2) and cap ρl at the ρt required for strength by
  §11.5.4.3; mento fed the equation the required ratio, under which the cap could never
  bind. `WallShearCheck.rho_l_min`, `WallShearDesign.rho_l_min`, the report and the design
  now use max(0.0025, min(0.0025 + 0.5·(2.5 − hw/lw)·(ρt − 0.0025), ρt,req)) with the ρt
  provided: a 25×400 wall with Ø12/15 E.F. under 2000 kN needs ρl,min = 0.00337 (was
  0.00321), and its design gives Ø10/17 E.F. vertical instead of Ø12/27.
  `min_vertical_reinforcement_ratio(hw_lw, rho_t, rho_t_req)` takes the provided ratio as
  its second argument and the required one as the cap.

- **`ShearWallSummary.check()` fails a wall that misses a limit under any combination.**
  The status came from the pass flag of the last combination checked; a wall that missed
  ρl,min under the governing one and met it under the last came out ✅. It is now DCR ≤ 1
  under every combination and no `wall.warnings`, so a mesh spacing past §11.7 fails the
  status as well.

- **`wall.shear_design` never pairs one mesh with another's DCR.** `WallShearCheck`
  carries the `mesh` it was checked with and the design is built from it; a mesh set by
  hand afterwards drops the results (`shear_checks` and `warnings` empty, `shear_design`
  raises `DesignNotRunError`) until the next check.

- **An imperial wall reports in kip.** `WallShearCheck.V_u`, `V_capacity` and `V_max`, the
  wall's attributes and the detailed strength table came in kN (the table labelled them
  kip).

- **`hw` is the height of the entire wall, or of the segment considered** (ACI 318-19 /
  CIRSOC 201-25 Chapter 2), not the storey height: the class docstring and the user guides
  said otherwise, and hw/lw sets αc and ρl,min. Wall clause numbers in the docs follow
  ACI 318-19 (§11.5.4.3 for αc and Vn, §11.5.4.2 for Vn,max, §11.6.2(b) for ρt,min,
  §11.7.2.1 / §11.7.3.1 for the spacing).

- **`shear_exceeds_section_limit` reads the limit of the section however it is
  reinforced:** under ACI 318-19 / CIRSOC 201-25 the Eq. (22.5.1.2) limit with the V_c of
  Table 22.5.5.1 for a section carrying A_v,min (V_c rises from row (c) to rows (a)/(b)
  once it does), under EN 1992-1-1 V_Rd,max of Eq. (6.9) at θ = 45°. A section that is
  only short of stirrups is no longer told to enlarge the section: a 20×60 with 2Ø12 and
  no stirrups under 320 kN was, at φV_max = 309 kN, while 2eØ10/10 on it carries the load
  at DCR 0.92 (its limit is 354 kN). `ShearCheckState` / `ENShearCheckState` carry it as
  `section_shear_limit`; `phi_V_max`, `V_Rd_max` and the report rows are unchanged.

- **An EN section with no stirrups is asked for the shear reinforcement the demand
  needs:** the minimum of §9.2.2 while V_Ed ≤ V_Rd,c (§6.2.1(3)–(4)), the truss of §6.2.3
  at the angle the demand fixes past it (§6.2.1(5)). `stirrups_required` quoted A_v,min
  whatever the shear: 2.40 cm²/m for a 30×50 under 150 kN, which needs 3.28.

- **The clear spacing between bars follows the stirrup.** `set_transverse_rebar`
  recomputes it, so `clear_spacing_below_min` no longer depends on whether the stirrups
  were set before or after the bars, nor waits for a reporting check.

- **`bars_do_not_fit` is cleared for a face given bars by hand** (beam and slab setters);
  it used to outlive the design that raised it.

- **`clear_spacing_below_min` is not reported for a face whose layers hold one bar
  each:** there is no pair to measure.

- **`mesh_spacing_exceeds_max` says the limit is the one mento applies** (§11.7.3.1 /
  §11.7.2.1 with lw/5 and lw/3 taken always, a conservative choice), not the clause's.

- **A combination with no label is named `#n` by its position** in `combinations`; an
  empty tuple now means a limit of the section alone, as documented.

- **Warning messages print as many significant figures as it takes to tell a value from
  its limit** (`13 cm exceeds the maximum 12.95 cm`, not `13 cm … 13 cm`).

- **The slab minimum is documented as what it is.** `OneWaySlab` no longer says its
  minimum is sized with the beam rule of §9.6.1.2 "as a known open point";
  `Section.support` cites §13.3.2.1 → §7.6.1.1 instead of a "§9.6.1.1(b)" that does not
  exist; and the footing theory page and docstrings quote §8.6.1.1 as printed —
  "0.0018Ag, or as defined in 8.6.1.2" — noting that the minimum over the punching
  critical section of 8.6.1.2 is not implemented.

### Changed

- **A one-way slab takes the slab minimum, 0.0018·Ag, under ACI 318-19 and CIRSOC 201-25.**
  `OneWaySlab` used the beam minimum of §9.6.1.2, ρmin·b·d, and the 4/3 relief of
  §9.6.1.3 that goes with it. A one-way slab is designed under Chapter 7, and its minimum
  is §7.6.1.1 (CIRSOC §7.6.1): 0.0018·b·h on the gross section, which §9.6.1.3 does not
  relieve — the same clause a `Footing` already reached through §13.3.2.1. A 100×20 strip
  with ADN 420 goes from As,min = 5.63 cm² to 3.60 cm², and a face whose moment asks for
  less than that is now designed to the minimum itself rather than to 4/3·A_s,calc. The
  Calcpad sheet "ACI 318-19 Slab Flexure 01 - Metric" was written with the beam minimum
  and needs the same update. EN 1992-1-1 is unchanged.

- **The slab and footing minimum is 0.0018·Ag for every steel grade.** mento scaled it as
  0.0018·420/f_y with a 0.0014 floor, the Table 7.6.1.1 / 24.4.3.2 of ACI 318-14. ACI 318-19
  withdrew that reduction (R24.4.3.2) and CIRSOC 201-25 prints the flat ratio as well. With
  ADN 420 or Grade 60 nothing changes; a B500S footing goes from 1.51‰ to 1.8‰, and for
  f_y below 420 MPa the minimum goes **down**: the old scaling 0.0018·420/f_y gave 2.7‰ at
  f_y = 280 MPa (Grade 40) and 2.16‰ at 350 MPa, where the flat ratio is 1.8‰ (a 1 m × 200 mm
  strip with f_y = 280 MPa goes from 5.40 cm² to 3.60 cm²; ACI 318-14's Table 24.4.3.2 had
  0.0020 for those grades, ACI 318-19 prints 0.0018 for all).

- **mento runs on pint 0.26.** The cap added after 0.26 broke CI is lifted and the
  dependency is `pint>=0.24` again. pint 0.26 types every arithmetic result as
  `PlainQuantity`, the base class of the registry's `Quantity`, so annotating with
  `pint.Quantity` rejected the result of any calculation under strict mypy. mento now
  annotates with `mento.units.Quantity`: the base class for the type checker and the
  registry class at runtime, so `isinstance` checks and user-side construction are
  unchanged. Under 0.26 the pretty multiplication sign in printed quantities is `⋅`
  (U+22C5) rather than `·`; the test that pinned the old sign now formats the expected
  value with pint itself, as the imperial one already did.

## [1.2.0] - 2026-09-09

How a Word report looks is now something a caller decides. The tables carried Word's
own `Light Shading`, whose colours were unreachable and whose fill Word overrode with a
theme; mento writes its own style definition instead, and `set_table_style` chooses it.

### Added

- **The section over a report's first three tables is called "Section Data".** It was
  called "Materials", which named only the first of the three tables under it — the others
  are the geometry and the design forces. The shear section that follows the limit checks
  is "Strength Checks" rather than "Design checks", because it is the strength that is
  checked there whether the run was a check or a design. Both are renamed in the per-beam
  reports, the shear wall report and the summaries, and both are translated. (#157)

- **The report headings are numbered and coloured.** The title is a `Heading 1` and the
  sections under it are `Heading 2`, numbered `1`, `1.1`, `1.2` … The numbers are Word's
  own rather than text mento wrote: the heading styles are attached to a multilevel list
  defined once in the document, so a section moved, deleted or inserted in Word renumbers
  the rest. The `Heading 1` is `#0A3E81`; the sub-headings, the running text and the tables
  are all `#323232`, so the one colour that appears reads as a heading rather than as
  decoration — `TableStyle.text_color` follows, and its default is no longer `1A1A1A`.
  Word's built-in heading styles are named in terms of the document theme twice over —
  `w:themeColor` beside `w:color`, `w:asciiTheme` beside `w:ascii` — and Word resolves the
  theme side first, the same trap as the theme fill in a table style. The colour survives
  because python-docx replaces the whole `w:color` element; the font needed the theme names
  removed by hand, which showed up in a place nobody looks: a heading's number is drawn in
  the paragraph mark's font, and the mark kept the theme's Calibri while the heading text,
  set run by run, was Lato. (#157)

- **The detailed annexes close on one page.** `flexure_results_detailed_doc()` and
  `shear_results_detailed_doc()` ran to two pages, and the title sat below the top margin
  rather than on it. Both came from python-docx's template defaults, which mento never
  overrode: 1.15-line spacing and a 10 pt gap after every paragraph. The extra leading of a
  multiple-spaced line goes *above* the text, and a heading's paragraph mark keeps the
  heading style's 14 pt whatever size its runs are, so a 10 pt title was set at the foot of
  a 14 pt line. `Normal` now carries single spacing and no trailing gap, a heading's mark
  is sized to its text, the space around a heading belongs to the heading, and the
  paragraph that keeps two tables from merging is 3 pt rather than a full line of body
  text. With the document text at 8.5 pt and `TableStyle.cell_padding_pt` — a new field —
  defaulting to 0.8 pt, a detailed annex produces
  one-page flexure and shear annexes under both ACI 318-19 and EN 1992-1-1, measured in
  Word. They stay one page whatever is checked: each report is built from the governing
  combination rather than from every one of them, so its tables have a fixed shape — 42
  rows for a flexure report, 40 for a shear one. `TableStyle(cell_padding_pt=1.4)` restores
  Word's own padding. (#157)

- **The Word tables have a style of their own, and it can be set.** Every table was tagged
  with Word's built-in `Light Shading`, whose grey and rule weights live in python-docx's
  template where no argument reaches them, whose fill names a colour and a *theme* colour
  at once — Word resolves the theme, so the grey on the page was never the grey declared —
  and which bolds the first column, which the builder then undid cell by cell. mento now
  writes its own definition into the document's `styles.xml`, once per document however
  many tables point at it, and `mento.set_table_style(TableStyle(...))` chooses it for the
  rest of the session the way `set_language` chooses the language. The knobs are the band
  fill and size, the header fill, colour and weight, the text and border colours, and the
  rules above the header, under it, at the foot, between rows and between columns. Because
  the look is a rule in the document rather than paint on the cells that existed when it
  was written, a row added in Word afterwards is banded like the rest. Colours are six hex
  digits and a `#` raises rather than being dropped silently by Word; thicknesses are
  points, clamped to the 0.25–12 pt Word will draw, and a thickness of zero means no rule.
  The default look is unchanged in kind — banded grey, bold header, ruled top, header and
  foot — and the green and red of a verdict column still outrank it. (#157)

### Added

- **`set_language` now reaches the summaries.** It translated every detailed report but
  stopped at `BeamSummary` and `ShearWallSummary`, which printed English tables and wrote
  English Word documents whatever the language was set to — the user guide said as much.
  `check()`, `flexure_results()` and `shear_results()` now return their DataFrames in the
  active language, and `results_detailed_doc()` writes a translated document: it was
  building its `DocumentBuilder` without passing the language along, and naming its
  headings with f-strings, so the interpolated text could never match a catalog key. Only
  the columns holding words are translated — `Beam`, `Label`, `Level`, `Position` and its
  `Top`/`Bottom`, `Status`. `b`, `As,bot`, `Av` and `DCRv` are the notation of the design
  code and stay identical in every language, as everywhere else in the catalog. English
  output is unchanged. (#155)

- **The sign convention is documented.** `N_x` is positive in compression, which is what
  feeds σ_Nu into v_c under ACI 318-19 §22.5.5.1 and σ_cp under EN 1992-1-1 §6.2.2(1), so
  a tension force entered positive is unconservative and silent; `M_y` is positive
  sagging and selects the tension face; `V_z` is a magnitude, and the design routine sizes
  stirrups from the largest required `A_v`, so a negative one reads as a smaller demand.
  The coordinate system page states all three with the clause each comes from, and the
  forces page states them where the components are introduced. No behaviour change — this
  is what the code already did, written down. (#155)

### Fixed

- **A translated summary table lost its pass/fail shading.** The Word builder looked its
  verdict column up by the English name and raised `KeyError` on a frame that had already
  been translated. It resolves the name against the document's language now, so the
  shading finds its column whichever spelling arrives. (#155)

## [1.1.0] - 2026-09-02

A footing is a section mento now designs as it is built, and the results a check returns
carry the resistance their ratio was formed from. Both came from mako, which had to
correct the engine's output for the first and divide it back out for the second.

### Added

- **`Footing`: a one-way slab bearing on the ground.** Everything about it is a
  `OneWaySlab`; what changes is the minimum longitudinal reinforcement, and that change is
  the design codes' rather than the class's. `Section.support` (`"free"` / `"soil"`) is
  a ClassVar the codes read, so the clause lives in `codes/`: ACI 318-19 §9.6.1.1(b)
  grants the exemption and §13.3.1.2 substitutes the shrinkage and temperature steel of
  §24.4.3.2 on the gross section (CIRSOC 201-25 shares the clause); EN 1992-1-1 takes the
  halved geometric minimum of a foundation on every face and, on a face that is bending,
  the larger of that and the crack-control minimum of §7.3.2(2). A footing is detailed
  between 100 and 300 mm, both faces set out at one spacing (or the top at twice the
  bottom), and a thin one warns. mento does no geotechnical calculation: bearing,
  settlement, sliding, overturning, plan size and thickness stay with the caller, and
  a face the moment does not put in tension is left unreinforced, because whether a
  footing carries steel there is the consumer's call. A theory page and a user guide
  page open on exactly that. (#150, #151)

- **A designed slab reads as a spacing, and its bars may only sit so far apart.** The
  rebar search answers in groups of bars, and `OneWaySlab` applied that answer the beam's
  way, so a designed strip came out as `2Ø12 + 4Ø12` with the spacing that actually
  details a slab left at zero. Each layer is now spread over the strip when the design is
  applied and the spacing is what the slab stores — `str(slab.reinforcement.bottom)` is
  `Ø12 mm/17 cm` — with `RebarLayer.s` carrying it and reading `None` on a beam. Both
  codes also cap the centre-to-centre spacing (ACI 318-19 §7.7.2.3, EN 1992-1-1 §9.3.1.1)
  through a registry hook, so a lightly loaded strip no longer covers its area with a bar
  every half metre. Areas and DCRs are unchanged. (#149)

- **The frozen results carry the capacity the DCR was formed from.** `FlexureFaceCheck`
  and `FlexureFaceDesign` gain `M_capacity`; `ShearCheck` and `ShearDesign` gain
  `V_capacity`. One neutral name for every code: `ØMn` and `ØVn` (capped by `ØVmax`)
  under ACI 318-19 and CIRSOC 201-25, `MRd` and `VRd` under EN 1992-1-1. A caller
  printing the resistance next to the demand no longer divides the demand by a rounded
  ratio, gets a real number for a face reinforced against no demand, and needs none of
  the `_phi_M_n_*`, `_M_Rd_*`, `_phi_V_n` or `_V_Rd` attributes of the compatibility
  layer. On a design the capacity is the governing combination's, so `demand / DCR`
  gives it back there too. Brief: `docs/architecture/capacity-in-design-results.md`.
  (#152)

- **`A_s_calc`: the steel the moment alone asks for, next to `A_s_req`.** `A_s_req` on
  `FlexureFaceCheck` and `FlexureFaceDesign` has the minimum folded in — `max(mechanical,
  minimum)`, or the 4/3 rule under ACI — which is right for choosing bars and wrong for
  anchoring them: a development length scaled by `A_s,nec / A_s,prov` needs the
  mechanical area, and a footing mat governed by its minimum carries little of the
  stress the minimum is sized for. Both design codes now return that area before the
  `max` — `reinforcement_for_moment` or `A_s1_lim + A_s2` under EN, `A_s_calc` under
  ACI — and it reaches the results as `A_s_calc`, enveloped over the combinations like
  `A_s_req`. Zero with no moment; equal to `A_s_req` once the moment governs, compression
  steel included. `A_s_req` does not change meaning.

### Fixed

- **`flexure_check_results()` returned `A_s_req`, `A_s_min` and `A_s_max` as bare
  floats** in the design code's canonical units (mm² or in²), where `check_flexure()`
  returned them as quantities in cm² or in². Both entry points now read the check state
  through the same path and return quantities. (#152)

- **The EN crack-control minimum governs the thin footings, not the thick ones.** Two
  docstrings had it inverted; the calculation and the theory page were right all along.
  Written on `A_ct = b·h/2`, its ratio goes with `k/2`, and `k` decays from 1.00 to 0.65
  between 300 and 800 mm, so the geometric minimum takes over around h = 640 mm for C25
  with B500S. The trend is now pinned by tests. (#151)

## [1.0.1] - 2026-08-31

### Changed

- **`BeamSummary.check()` returns the columns the report prints, and only those.** It
  carried the section, the reinforcement, the required areas, the demands, the capacities
  and the three DCRs, and the Word report then selected a shorter set on its way to the
  page — so a notebook and the document disagreed about what the summary is. Both now
  show `Beam`, `b`, `h`, `As,top`, `As,bot`, `Av`, the governing demands, the three DCRs
  and `Ok?`: the shape the shear-wall summary already had, ending on the DCRs and the
  verdict.

  Gone from `check(capacity_check=False)`: `cc`, `As,req,top`, `As,req,bot`, `Av,req`,
  `Av,real`, and the capacity columns (`ØMn,top`/`MRd,top`, `ØMn,bot`/`MRd,bot`,
  `ØVn`/`VRd`). None of it is lost — `flexure_results()` and `shear_results()` report the
  required areas and the capacities per combination, which is where they mean something,
  and `check(capacity_check=True)` still returns the capacities.

  This removes columns from a public DataFrame, which under the policy above would call
  for a major release. It ships as a patch deliberately: 1.0.0 was a day old and the
  change is narrow. Recorded here rather than passed over quietly.

## [1.0.0] - 2026-08-31

The API is stable from this release on (ADR-0003): breaking changes need a major version
from here. What made 1.0 worth declaring is the architecture work below — checks that
return values instead of writing them to the section, so a caller can run one over many
sections and trust the answer.

### Added

- **`section.reinforcement`** — what a section carries right now, as a frozen
  `SectionReinforcement(bottom, top, transverse)`. It never raises, so a layout that was
  just set with `set_longitudinal_rebar_bot` can be read back without a check having run.
  The design-result objects were previously the only way to ask, and they are gated behind
  `DesignNotRunError`.
- **`beam.shear_check_results(forces)` and `beam.flexure_check_results(forces)`** — one
  frozen result per load combination, and no report built at all. Three to five times
  faster than the reporting path, and the numbers are identical; a test asserts that for
  both design codes.
- **`beam.shear_checks` and `beam.flexure_checks`** — every combination checked so far, so
  a caller that wants to envelope them differently no longer has to reach into the element.

### Changed

- **A check no longer writes to the section.** `shear_check_results` and
  `flexure_check_results` leave the element untouched — measured across both design codes,
  both checks and with and without stirrups: zero attributes changed, zero created. The
  reporting entry points still copy their results back, because the report tables read
  them off the element; that compatibility layer is deprecated and goes in a later release.
- **The shear check keeps the assumed stirrup diameter, as the flexure check already did.**
  On a section with no stirrups configured, `check_shear` used to drop the diameter the
  settings assume and recompute the effective depths on the section, so a `check_flexure`
  run afterwards reported a different DCR purely because of call order — 0.593 against
  0.582 on the same beam. Both checks now read the same effective depth. If you check a
  section that has no transverse reinforcement, say so with
  `set_transverse_rebar(n_stirrups=0, d_b=0, s_l=0)` and the depths follow.
- **`BeamSummary.check()` names its verdict column `Ok?`**, not `Status`.
- Adding a design code no longer means editing an element. A code declares itself in
  `mento/codes/<code>/code.py` and is found by walking that directory; nothing under
  `beam.py`, `shear_wall.py`, `rebar.py` or the summaries names a code any more.
- The Word reports read better on the page: tables are sized to their content rather than
  stretched across the line, forces show one decimal and DCRs two, the pass/fail cells are
  shaded green or red, and the summaries drop the capacity ticks that repeat what the DCR
  column beside them already says.

### Performance

- **A check is roughly ten times faster.** Shear and flexure on one section went from
  2.35 ms to 0.23 ms, so 20,000 sections take 4.6 s instead of 47 s. The equations run on
  plain floats and a section publishes its geometry and materials converted once
  (ADR-0005); pint stays at the boundary, where the inputs and the results are.
- Flexural design is four times faster: 306 ms to 72 ms, from the reinforcement search no
  longer building a `Quantity` per candidate.

### Fixed

- The report tables no longer write to the section they describe.


### Removed

- **Support for Python 3.10 and 3.11.** `requires-python` is now `>=3.12`, so pip refuses
  to install mento on the older interpreters instead of installing it and failing later.
  The classifiers and the conda recipe's `python_min` follow. Stay on 0.5.2 if you are
  pinned to 3.10 or 3.11.
- Ubuntu from the test matrix. Tests run on windows-latest against Python 3.12 and 3.13,
  down from eight jobs to two. mento is pure Python and nothing in it is
  platform-specific, but Linux is no longer verified on every pull request. The lint and
  docs jobs still run on ubuntu-latest, and the PyPI distributions are still built there.

### Fixed

- **Shear design ignored the spacing limit across the width of the section.** The stirrups
  were sized from the required area alone, so a wide beam came back with a single
  two-legged stirrup whose legs sat far further apart than ACI 318-19 Table 9.7.6.2.2 or
  EN 1992-1-1 9.2.2(8) allow — 44 cm against a 28 cm limit on a 50 cm section. The check
  reported the violation, but the design would not avoid it, so `design_shear` handed back
  a layout its own detailed report then marked as not compliant. The design now starts from
  the fewest legs the width admits and adds stirrups rather than only tightening the
  longitudinal spacing. Over a sweep of 168 width and demand combinations across the three
  codes, 95 designs were in violation and none are now. Closes
  [#94](https://github.com/mihdicaballero/mento/issues/94).
- The spacing across the width was computed with whatever stirrup diameter the previous
  pass had left on the beam instead of the one being tried, so the value stored for each
  candidate was off by the difference between the two diameters.

## [0.5.2] - 2026-08-25

### Added

- Spanish detailed reports. `mento.set_language("es")` switches
  `flexure_results_detailed()`, `shear_results_detailed()` and their `_doc()` counterparts
  to Spanish, for beams, one-way slabs and shear walls, in the console and in the generated
  Word documents. English remains the default, and `mento.get_language()` and
  `mento.available_languages()` report the current and the available choices. Variable
  names, units, the design code designation and the generated file names are not
  translated. A label with no translation is written in English rather than raising. See
  [Report language](https://mento-docs.readthedocs.io/en/latest/user_guide/language.html).
  Closes [#79](https://github.com/mihdicaballero/mento/issues/79) and
  [#126](https://github.com/mihdicaballero/mento/issues/126).
- A DOI. Releases are archived on Zenodo, and
  [10.5281/zenodo.21956634](https://doi.org/10.5281/zenodo.21956634) always resolves to the
  latest one. It is in `CITATION.cff`, in the README badge and in the citing guide.
- A [Theory section](https://mento-docs.readthedocs.io/en/latest/theory/index.html) in the
  documentation, with one page per element and design code, so the tool can be audited
  against the codes it implements. Each page names what is out of scope, the places where
  mento takes a position the code leaves open, and a table mapping every check to the test
  that pins it and the external source it was verified against.

### Changed

- Flexure output is labelled with the symbols of the active design code. An EN 1992-2004
  result was reported with ACI symbols — `M_u` and `\phi M_n` in the markdown line, and
  `Mu,top` / `Mu,bot` heading the design forces of the detailed table — where the Eurocode
  writes `M_Ed` and `M_Rd`. Only the limiting-case branch was affected, which is the one
  taken when no explicit force is passed. ACI 318-19 and CIRSOC 201-25 are unchanged.
- The conda recipe is pinned to 0.5.1 and its checksum.

### Fixed

- **EN 1992-2004 flexural design produced layouts that its own check then rejected.** Over
  a sweep of 288 combinations, 29 designs came back with DCR > 1. Several independent
  defects: the lever arm applied `lambda` twice, under-sizing the tension steel by about
  2%; the ductility limits mixed the neutral-axis and block-depth conventions, leaving
  `M_lim` about 20% too high; the redistribution limit ignored `k_3`/`k_4` above C50/60;
  `M_Rd` counted compression steel the section did not need, which made capacity
  non-monotonic; and the top face was sized with the bottom face's effective depth.
  Cross-checked against the Concise Eurocode 2 closed form and plain equilibrium. Flexural
  design is now driven by one shared strategy for both codes.
- **The ACI size effect factor was not capped at 1.0.** ACI 318-19 Eq. 22.5.5.1.3 writes
  `lambda_s = sqrt(2/(1 + d/10in)) <= 1.0`. Without the cap, sections with `d` below about
  250 mm got a factor above 1, inflating `V_c` instead of reducing it — the opposite of
  what the provision is for, and on the unsafe side. It only applies in the
  `A_v < A_v_min` branch, so in practice this moves slabs; beams under test are unaffected.
- Checking a section with no reinforcement for ACI shear raised `ZeroDivisionError`
  instead of reporting an insufficient section. With `rho_w` at zero, Table 22.5.5.1 puts
  `phi*V_n` at exactly zero and the DCR division blew up. It now reports `DCR = inf`,
  matching the guard already in the wall module. EN is unaffected, since `V_Rd,c` carries
  the `v_min` floor of 6.2.2(2).
- A one-way slab is no longer reported as a beam. `OneWaySlab` inherited the titles of
  `RectangularBeam`, so its detailed reports were headed `BEAM FLEXURE DETAILED RESULTS`
  and its Word file named `Beam S1 flexure check ACI 318-19.docx`. Slabs now use their own
  wording and file name; beams are unchanged.
- The forces table of the detailed flexure report is headed `Design forces`, matching the
  shear report and the Word output. It was `Design_forces` in the console output only.
- The Word reports spell the `Limit checks` heading the same way everywhere. Flexure used
  `Limit Checks` and shear used `Limit checks`, in both the element and the summary
  documents.

## [0.5.1] - 2026-08-15

### Added

- Public design results API: `beam.flexure_design` and `beam.shear_design` return plain,
  frozen data objects with the reinforcement a check or design produced, so results no
  longer have to be read from private attributes such as `_A_s_bot` or `_stirrup_s_l`.
  Reading either before running a check raises `DesignNotRunError`. Required areas and
  DCRs are the envelope over every load combination checked, so they describe the
  combination that governs each face. See
  [Design results](https://mento-docs.readthedocs.io/en/latest/user_guide/design_results.html).
- A [citing guide](https://mento-docs.readthedocs.io/en/latest/getting_started/citing.html)
  in the documentation, covering which version to cite and a BibTeX entry.
- A conda recipe under `conda-recipe/`, kept in step with `pyproject.toml`, ready to be
  submitted to conda-forge. The submission and update procedure, and the one-time Zenodo
  setup that gives releases a DOI, are documented in CONTRIBUTING.
- Two example notebooks in Spanish, design and check of a rectangular beam under
  CIRSOC 201-2025.

### Changed

- `mento.summary` was renamed to `mento.beam_summary`, matching `mento.shear_wall_summary`.
  The old module still works and emits a `DeprecationWarning`; `from mento import
  BeamSummary` is unaffected.

### Fixed

- **`pip install mento` on Google Colab.** The IPython requirement was `>=8.0`, while Colab
  pins ipython to 7.34.0, so installing mento there could not resolve without upgrading
  IPython out from under the running session. The floor is now `>=7.34`; mento only uses
  `IPython.display.Markdown` and `display`, which both bounds cover.

## [0.5.0] - 2026-08-01

Infrastructure release. No changes to the design calculations, and no changes to the public
API beyond one addition.

### Added

- `mento.__version__` exposes the installed package version.
- Contributing guide, code of conduct, security policy and citation metadata at the
  repository root, plus issue and pull request templates.
- Python 3.13 is tested and supported.
- Optional dependency groups: `mento[test]`, `mento[dev]` and `mento[docs]`.
- Continuous integration now runs mypy, `ruff format --check` and a documentation build in
  addition to the test suite.
- Releases are published to PyPI automatically from a GitHub Release, using trusted
  publishing.
- A disclaimer of professional responsibility in the README and the documentation.

### Changed

- Dependencies are declared in `pyproject.toml` instead of being read from
  `requirements.txt`, with lower bounds rather than exact pins.
- **numpy is no longer pinned to 1.26.4.** The pin blocked numpy 2.x and caused install
  conflicts in environments with other scientific packages. The test suite passes on both
  numpy 1.26 and numpy 2.5.
- The ruff pre-commit hook was updated from v0.6.2 to v0.15.22, so the hook, CI and a local
  `ruff format` all agree. Eleven files that had drifted were reformatted.
- mypy configuration moved from `mypy.ini` into `pyproject.toml` as a single source. Modules
  that do not pass strict checking yet are listed explicitly, so the modules that are clean
  cannot regress.
- Read the Docs builds on Python 3.12 and installs the `docs` extra.

### Fixed

- `RectangularSection._ax` was annotated as `plt.Axes`, which is not a valid type; it is now
  `matplotlib.axes.Axes`.
- The `LICENSE` file still contained the unfilled `<year>` and `<author>` placeholders.
- The test workflow checked for `requirements.txt` but installed `requirements_dev.txt`.

### Removed

- Dead `[tool.black]` and `[tool.pylint]` configuration, an empty `.github/workflow` file,
  and unused documentation dependencies (dask, xarray, sparse, mip and others that mento
  never imported).

## [0.4.1] - 2026-05-17

### Added

- Shear wall summary, with results across multiple walls and a Word export.

## [0.4.0] - 2026-05-15

### Fixed

- Minimum reinforcement calculation for ACI 318-19.

## [0.3.6] - 2026-05-15

### Added

- Shear wall check and design for ACI 318-19 and CIRSOC 201-25.

## [0.3.5] - 2026-01-04

## [0.3.4] - 2025-12-20

### Added

- One way slab check and design.

## [0.3.0] - 2025-11-09

### Added

- Flexure check and design for EN 1992-2004.

## [0.2.8] - 2025-10-21

## [0.2.7] - 2025-08-17

## [0.2.6] - 2025-07-27

## [0.2.5] - 2025-03-23

First public release on PyPI: rectangular concrete beam check and design for flexure and
shear under ACI 318-19 and CIRSOC 201-25, unit aware calculations, results as pandas
DataFrames, and Word calculation reports.

[Unreleased]: https://github.com/mihdicaballero/mento/compare/v1.2.0...HEAD
[1.2.0]: https://github.com/mihdicaballero/mento/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/mihdicaballero/mento/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/mihdicaballero/mento/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/mihdicaballero/mento/compare/v0.5.2...v1.0.0
[0.5.2]: https://github.com/mihdicaballero/mento/compare/v0.5.1...v0.5.2
[0.5.1]: https://github.com/mihdicaballero/mento/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/mihdicaballero/mento/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com/mihdicaballero/mento/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/mihdicaballero/mento/compare/v0.3.6...v0.4.0
[0.3.6]: https://github.com/mihdicaballero/mento/compare/v0.3.5...v0.3.6
[0.3.5]: https://github.com/mihdicaballero/mento/compare/v0.3.4...v0.3.5
[0.3.4]: https://github.com/mihdicaballero/mento/compare/v0.3.0...v0.3.4
[0.3.0]: https://github.com/mihdicaballero/mento/compare/v0.2.8...v0.3.0
[0.2.8]: https://github.com/mihdicaballero/mento/compare/v0.2.7...v0.2.8
[0.2.7]: https://github.com/mihdicaballero/mento/compare/v0.2.6...v0.2.7
[0.2.6]: https://github.com/mihdicaballero/mento/compare/v0.2.5...v0.2.6
[0.2.5]: https://github.com/mihdicaballero/mento/releases/tag/v0.2.5
