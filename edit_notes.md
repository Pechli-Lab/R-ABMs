- # ABM Tutorial — Edit Notes
  *April 8th 2026*

  ## Bug fixes
  1. `Probs()` was subsetting the population before computing probabilities, causing row indices to drift — SQ individuals were getting wrong infection probabilities. Rewrote using full-population logical masks.
  2. Redundant `v_names_states` re-declaration inside `Probs()` shadowed the global. Removed.
  3. `Probs()` was using global `n_i` to size the output matrix instead of `length(M_t)` — silently breaks if called outside the main sim.

  ## Performance
  1. Replaced string states (`"S"`, `"I"`, etc.) with integer constants (`STATE_S = 0L` through `STATE_R = 4L`). `m_M` is now an integer matrix.
  2. `home_exposure()` was a for-loop over families. Now uses `tabulate()` + vectorized indexing — one pass, no repeated scans.
  3. `work_exposure()` was a growing-vector loop with repeated full-table scans. Now pre-builds a `split()`-based lookup, `lapply` over infected workers, single `unlist()` at the end.
  4. `create_abmpop()` contact sampling replaced `rowwise() %>% mutate()` with `vapply()`.
  5. `v_Ts` update in `MicroSim()` swapped `dplyr::if_else` for base R logical indexing.

  ## Code quality
  1. Renamed `expl` → `work_exposure_count` and `expl_home` → `home_exposure_count` everywhere including the field table in Section 04.1.
  2. Added `debug = FALSE` parameter to `MicroSim()` — transition checks only run when asked.
  3. Added inline comments for: `infect_prob()` formula, SQ home-only risk in `Probs()`, quarantine clock reset logic, `v_Ts` purpose, day-17 recovery assumption.

  ## To discuss
  1. `m_M` is now integers not characters — any external scripts reading it expecting strings need updating.
  2. `debug = FALSE` by default — ok for performance but maybe keep `TRUE` in the tutorial version?
