# Brief — the frozen results should carry the capacity

Status: proposed (2026-09-02). Raised from mako, which hit it building the footing
reinforcement tables.

## The gap in one sentence

A checked or designed section publishes **a ratio but not the resistance the ratio was
formed from**, so a caller that wants to print `ØMn` / `MRd` / `ØVn` / `VRd` next to the
demand has nowhere public to read it.

## What the public results carry today

`mento/design_results.py`:

```python
@dataclass(frozen=True)
class FlexureFaceCheck:
    A_s_req: Optional[Quantity]
    A_s_min: Optional[Quantity]
    A_s_max: Optional[Quantity]
    DCR: float

@dataclass(frozen=True)
class FlexureFaceDesign:
    layers: Tuple[RebarLayer, ...]
    A_s: Quantity
    A_s_req: Quantity
    A_s_min: Quantity
    A_s_max: Quantity
    DCR: float

@dataclass(frozen=True)
class ShearCheck:
    label: str
    A_v_req: Optional[Quantity]
    A_v_min: Optional[Quantity]
    DCR: float

@dataclass(frozen=True)
class ShearDesign:
    n_stirrups: int
    ...
    A_v_req: Quantity
    A_v_min: Quantity
    DCR: float
```

Areas, ratios, bars. No resistance anywhere.

## The value already exists — it just does not leave the check

It is computed and held in the check state, one file away from the frozen result:

- `mento/codes/check_state.py` → `FlexureCheckState.phi_M_n_bot` / `.phi_M_n_top`
- `mento/codes/check_state.py` → `ShearCheckState.phi_V_n`
- `mento/codes/ACI_318_19_beam.py:818` → `st.DCR_bot = st.M_u_bot / st.phi_M_n_bot`
- `mento/codes/EN_1992_2004_beam.py:677` → `st.DCR_bot = round(st.M_Ed_bot / st.M_Rd_bot, 3)`

So the ratio the public result *does* publish is formed from a number the public result
*does not*. From the outside the only route is the compatibility layer —
`section._phi_M_n_bot`, `section._M_Rd_bot`, `section._phi_V_n`, `section._V_Rd` — which
`mento/beam.py` declares under `TYPE_CHECKING` with the comment:

> They leave with the compatibility layer. A new design code should not add to this list —
> its report tables should read the check state instead.

A consumer reading those is coupling to something mento has already announced it will
delete, so mako does not.

## What mako does instead, and why it is not good enough

mako divides the demand back out:

```python
capacity = demand / DCR   # exact by definition, when DCR > 0
```

Two problems, both visible in office output:

1. **It inherits the ratio's rounding.** EN rounds `DCR` to 3 decimals, so two faces
   carrying identical reinforcement come back with different resistances — a real footing
   printed `MRd` 43.3 on one face and 43.4 on the other.
2. **It says nothing when `DCR = 0`.** A footing's top mat usually carries minimum
   reinforcement against no demand. It has a real flexural capacity; nothing published
   says what it is, so the cell goes blank.

Neither is worth working around downstream: computing a section capacity in mako would be
a second implementation of mento's own equations, diverging the first time either changes.

## Proposal

Add one field to each of the four frozen results.

```python
@dataclass(frozen=True)
class FlexureFaceCheck:
    ...
    M_capacity: Optional[Quantity]   # ØM_n under ACI/CIRSOC, M_Rd under EN

@dataclass(frozen=True)
class FlexureFaceDesign:
    ...
    M_capacity: Quantity

@dataclass(frozen=True)
class ShearCheck:
    ...
    V_capacity: Optional[Quantity]   # ØV_n under ACI/CIRSOC, V_Rd under EN

@dataclass(frozen=True)
class ShearDesign:
    ...
    V_capacity: Quantity
```

Notes on the shape:

- **One neutral name per quantity, not one per code.** The code already owns the
  *symbol* — `DesignCode.summary_columns` publishes `ØMn,bot` against `MRd,bot` — so the
  presentation layer keeps naming things per code and the data structure does not have to.
  This is the same split ADR-0004 draws.
- **The envelope helpers need a rule.** `envelope_flexure_face` and `envelope_shear`
  currently take the worst of each quantity independently. A capacity is a property of the
  section, not of the combination, so every combination should report the same one; taking
  the minimum present is the safe reading if that ever stops holding.
- **`Optional` on the check results only**, matching the fields beside them: a field is
  `None` when the design code did not set it for that combination.
- Nothing in `codes/` needs new physics. `phi_M_n_bot`, `phi_M_n_top` and `phi_V_n` are
  already on the state; this is plumbing them into the frozen result the way `A_s_req`
  already is.

One wrinkle to decide on the way past: `EN_1992_2004_beam.py` writes `st.M_Rd_bot` and
`st.M_Ed_bot`, which are **not** fields of `FlexureCheckState` — the dataclass is neither
frozen nor slotted, so they attach dynamically. Either EN should fill the declared
`phi_M_n_*` fields, or the state should declare the EN names too. Silent dynamic
attributes on a state object are how two codes end up disagreeing about what a field
means.

## Acceptance

- `beam.check_flexure(...)` / `check_shear(...)` and `slab.design(...)` /
  `Footing.design(...)` all return results whose capacity field is populated for both ACI
  318-19 and EN 1992-2004.
- For any combination with `DCR > 0`, `demand / DCR == capacity` to within the code's own
  rounding — a regression test worth writing, since it is the invariant mako is
  currently leaning on.
- A face carrying only minimum reinforcement reports a real capacity rather than nothing.
- No consumer needs `section._phi_M_n_bot`, `_M_Rd_bot`, `_phi_V_n` or `_V_Rd` any more —
  which is one more thing the compatibility layer can drop on its way out.

## Who is waiting on it

mako's `get_footing_rebar_summary()`, one row per footing type, direction and face:

| Tipo | Dir. | Cara | MEd (kNm) | MRd (kNm) | VEd (kN) | VRd (kN) | Ø/s | DCR |
|---|---|---|---|---|---|---|---|---|
| 6 | x | inferior | 131.3 | 157.3 | 129.5 | 261.4 | Ø10/15 | 0.835 |
| 3 | y | inferior | 2.0 | 38.6 | 0.0 | *(vacío)* | Ø10/24 | 0.051 |

The `MRd` column is the division above; the blank `VRd` is the `DCR = 0` case.
