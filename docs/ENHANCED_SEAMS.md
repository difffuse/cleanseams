# Enhanced seam policy

`enhanced_v1` is an explicit, experimental policy. `native_baseline` remains the
default. The enhancement pass follows the baseline result and evaluates bounded
proposals against constraints and a fixed run-origin quality envelope.

The default enabled stages include endpoint sliding, relief, rings, cavities,
shell proposals and spur removal. Preference-aware region proposals are disabled
by default (`enhanced.regions.enabled=false`); enable them per run. Symmetry defaults
to asymmetric operation, and physical scale is used for distance policies.
`include/cleanseams/options.hpp` and `option_registry()` are the authority for defaults,
limits and UI controls.

Locked edges, protected seams, explicit pins and painted required/forbidden edges
remain hard constraints. Candidates are transactions: a refused candidate cannot
partially alter the retained state. Candidate and solve budgets, deadlines and
cancellation bound the search. Reports include rejection reasons and work counts;
an enabled stage does not promise an accepted improvement.

The native objective is `seam_objective.v2`. Acceptance combines its ordering with
quality/constraint checks and a practical-effect threshold. Strict overlap checks
account for stored float32 UVs. Existing invalid layouts must be repaired within the
policy limits or reported as invalid; a successful process exit does not mean that
every quality target passed.

Select the policy through native options, `cleanseams_cli run --behavior enhanced_v1`,
or the standalone Blender seam-policy control. Painted attributes are transported
by the host and interpreted by the native library. No Python computational fallback
is used. Review seam placement and UVs visually before adopting an enhanced result.
