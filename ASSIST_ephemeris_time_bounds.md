# ASSIST failure mode: integrating outside the JPL ephemeris time span

**One-line summary.** ASSIST propagates a test particle using the Sun, planets,
Moon, and large asteroids whose positions come from JPL SPK ephemeris files.
Those files cover a *finite* interval of time. Certain orbits, integrated far
enough in time, drive the integrator's evaluation epoch toward (and past) the
edge of that interval — where the ephemeris can no longer be interpolated. The
practical symptoms are a long integrator "grind" near the boundary and,
ultimately, undefined behavior if the boundary is crossed.

This note describes the failure mode and the conditions that trigger it, so it
can be investigated directly in the ASSIST source. It was written from the
*downstream* (layup) side; the exact ASSIST-internal behavior at the boundary
is flagged below as the thing to verify.

---

## Background: where the time bound comes from

ASSIST does not integrate the planets — it reads their states from JPL
ephemerides and integrates the massless (or test) particle in that
time-dependent force field. Two files are involved in the configuration this
was observed with:

- **Planetary ephemeris** — `linux_p1550p2650.440` (JPL DE440). The filename
  encodes the span (~years 1550–2650).
- **Small-body perturbers** — `sb441-n16.bsp` (16 largest asteroids, from the
  DE441/sb441 small-body set).

Every call that advances the integrator (`assist_integrate_or_interpolate(t)`)
needs *all* of these bodies' states at the requested epoch `t`. The usable
interval is the **intersection** of the two files' spans — i.e. the *narrower*
of the planetary and small-body coverage binds.

> **Verified in ASSIST (2026-06-18).** Measured directly with the new
> `assist_ephem_time_bounds()` API (see below) on the configuration above
> (`linux_p1550p2650.440` + `sb441-n16.bsp`): the binding interval is
> **JD 2287184.5 … 2688976.5 = MJD −112816 … 288976 (≈ years 1550 … 2650)**.
> The **planetary file binds**; the small-body file is equal-or-wider. This does
> *not* match layup's hardcoded guard (MJD −163545 … 236455, ≈ 1411 … 2506): the
> real floor is ~year **1550**, ~140 yr *later* than layup's 1411, so layup's
> lower guard is currently **too permissive** (it would let an integration step
> to ~1411 and hit a real coverage error before its own guard fires), while its
> upper guard is over-conservative (2506 vs. the true 2650). Callers should use
> `assist_ephem_time_bounds()` instead of hardcoding.

---

## The failure mode

When the integrator's epoch `t` approaches the edge of the ephemeris coverage:

1. **Grind near the boundary.** Reaching the boundary often coincides with the
   physically hard part of the trajectory (deep perihelion passage, a close
   planetary encounter, or a near-parabolic/hyperbolic limb). IAS15's adaptive
   step size collapses there, so the wall-clock cost per unit of integrated
   time balloons — the run appears to "hang" while it is in fact taking a very
   large number of tiny steps. This is slow but *bounded* and *completes* if the
   orbit stays inside the span.

2. **Crossing the boundary.** If `t` is pushed outside the covered interval, the
   ephemeris for at least one body can no longer be interpolated. **Verified
   behavior (2026-06-18):** this is a *clean, guarded error*, not undefined
   behavior. Both backends bounds-check `jd_ref + jd_rel` against the target's
   `beg`/`end` *before* any array indexing and return `ASSIST_ERROR_COVERAGE`
   ("The requested time is outside the coverage provided by the ephemeris
   file.") — see `spk.c:416`, `spk.c:563`, `ascii_ephem.c:294`. No silent
   extrapolation, no out-of-range read. In the force routine this becomes
   `reb_simulation_error(...)` + `sim->status = REB_STATUS_GENERIC_ERROR`
   (`forces.c:77`), which halts the integration. As of the fixes below,
   `assist_integrate_or_interpolate` also short-circuits on that status instead
   of interpolating from the incomplete final step.

---

## Conditions under which it is encountered

The common thread is **integrating far from the observation epoch**, which
distant and weakly-bound orbits force:

- **Long-period / near-parabolic comets.** Computing an "original" or "future"
  orbit means integrating until the object is far from the planetary region
  (e.g. heliocentric distance ≈ 250 AU, where the barycentric orbit is taken as
  the original element set). **Measured on the CODE long-period-comet catalogue
  (369 objects): each integration runs ~300 years *backward* in time** (296–306 yr
  across the sampled comets) to reach the inbound 250 AU crossing, finishing
  around the **early-to-mid 1500s**. The ephemeris floor here is only ~440 yr back
  (MJD −163545 ≈ **year 1411**), so a typical comet clears it by just **~70–170
  years** — a thin margin. Objects with earlier epochs, smaller perihelia, or
  more hyperbolic `e` need a bit more than ~300 yr and **cross the floor before
  reaching 250 AU**: that is the boundary case. (The full span of the files above
  is ~years 1411–2506 vs a ~2025 epoch — ~500–600 yr each way.)

- **Hyperbolic / interstellar objects** (e ≥ 1): the integration never "turns
  around," so propagating far enough in time always eventually exits the span.

- **Any long-baseline propagation** (orbit determination or prediction over many
  centuries), independent of object type.

- **Step-size and direction interact:** larger requested steps overshoot the
  boundary; near perihelion the adaptive step shrinks and the grind dominates.

Conversely, ordinary main-belt / NEO / TNO fits over their observed arcs stay
comfortably inside the span and integrate in well under a second.

---

## Worked example (downstream, for reference)

In **layup**, `comet._apply_comet` computes long-period-comet original orbits
(inverse original semimajor axis, `inv_ao`) by stepping a comet out to 250 AU:

```python
# src/layup/comet.py  (paraphrased)
ASSIST_TIMEFRAME_MAX_MJD = 236455     # ≈ year 2506
ASSIST_TIMEFRAME_MIN_MJD = -163545    # ≈ year 1411

while of.d > 250 and oi.d > of.d and sim.t < ASSIST_TIMEFRAME_MAX_MJD:
    oi, of, sim = _assist_integrate(sim, ex, dt, ephem, include_assist=True)
...
if sim.t >= ASSIST_TIMEFRAME_MAX_MJD or sim.t <= ASSIST_TIMEFRAME_MIN_MJD:
    # comet left the ASSIST-covered interval before reaching 250 AU;
    # fall back to a pure-REBOUND (no ASSIST perturbations) continuation
    rebound_only.append(comet)
```

So layup **guards the loop with hardcoded MJD bounds** and switches to a
pure-REBOUND continuation once a comet would step outside the ASSIST interval —
i.e. it works around the failure mode rather than hitting it. The guard call is
`ex.integrate_or_interpolate(t)` → `clibassist.assist_integrate_or_interpolate`.

**Observed cost:** the CODE-catalogue LPC regression (`test_apply_comet`) takes
**~246 s** to complete on a 2024 Apple-silicon Mac — entirely inside ASSIST's
integrator, dominated by the near-boundary / deep-perihelion grind described
above. It *passes*; it is just slow. (It was mistaken for a crash only because a
120 s test timeout killed the worker mid-integration.)

---

## Investigation outcome (2026-06-18)

The original questions, with what the source actually does:

1. **Exact covered interval / which file binds / API to query it.** *Resolved.*
   The per-target `beg`/`end` (Julian Days, TDB; `spk.c:_jul`) define the span;
   the binding interval is their intersection across all required bodies. For
   the configuration above the **planetary file binds at ≈ years 1550–2650**
   (see the verified block at the top). A new public accessor —
   **`assist_ephem_time_bounds(ephem, &t_beg, &t_end)`** (Python:
   `Ephem.time_bounds()`) — returns this interval in the same time convention as
   the `t` passed to the integrator (relative to `jd_ref`), so callers can guard
   precisely instead of hardcoding.
2. **What it does outside the span.** *Resolved.* Explicit `ASSIST_ERROR_COVERAGE`
   → `REB_STATUS_GENERIC_ERROR`, bounds-checked before any indexing. No clamp,
   no silent extrapolation, no out-of-range read. (See "Crossing the boundary".)
3. **Clear, fast out-of-bounds signal.** *Addressed.* The signal already existed
   internally but `assist_integrate_or_interpolate` swallowed it — it proceeded
   to interpolate from an incomplete step without checking `sim->status`. It now
   short-circuits on `REB_STATUS_GENERIC_ERROR` and leaves the status set, so a
   caller can detect-and-defer by inspecting `sim.status` after the call.
4. **Is the near-boundary grind avoidable?** *Open.* Still no adaptive-step floor
   or max-steps cap in this path; the grind is physical (IAS15 step collapse at
   deep perihelion) and bounded but slow. Not addressed by the changes above.
   One subtlety surfaced: IAS15 evaluates forces at *trial substeps* that can
   overshoot the requested `t`, so a request sitting exactly at the boundary can
   trip `ASSIST_ERROR_COVERAGE` even when the endpoint is in-bounds — keep guard
   bounds slightly inside the reported span.
5. **Frame/timescale consistency.** *Resolved.* All comparisons are JD in the
   ephemeris timescale (TDB): `target->beg/end` are JD via `_jul`, compared
   against `jd_ref + jd_rel`. MJD↔JD and TDB↔UTC conversions are the caller's
   responsibility; ASSIST works in JD/TDB throughout. `assist_ephem_time_bounds`
   reports relative to `jd_ref` (add `jd_ref` for absolute JD).

### Changes made (this commit)

- **`assist_ephem_time_bounds()`** — new C API + `assist.h` declaration +
  `Ephem.time_bounds()` Python wrapper. Returns the binding coverage interval.
- **`assist_integrate_or_interpolate`** — checks `sim->status` after
  `reb_simulation_integrate` and skips interpolation on a coverage/integration
  failure, leaving the error status for the caller. Function stays `void` (no
  ABI/binding change).

Both verified against `linux_p1550p2650.440` + `sb441-n16.bsp`: bounds report
≈1550–2650; an in-bounds step yields `status==0`, an out-of-bounds step yields
`status==1` with no spurious interpolation.

---

*Provenance: distilled from a layup CI investigation (2026-06-18), then verified
and fixed in the ASSIST repo the same day. The grind in
`assist_integrate_or_interpolate` surfaced while running long-period-comet
original-orbit tests; root cause traced to finite JPL ephemeris time spans. The
boundary behavior turned out to be a clean guarded error (`ASSIST_ERROR_COVERAGE`),
now also surfaced through `sim->status` and the new `assist_ephem_time_bounds()`
API.*
