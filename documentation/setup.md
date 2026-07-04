# Setup and running

`base-edit-safety` has no pip dependencies of its own. You install the prediction tools yourself, under
their own licenses; thin adapters call them (`src/base_edit_safety/adapters/`), and this repo never
bundles their source or weights.

## Requirements
- **Python 3.10+.**
- **BE-DICT** (edit-outcome predictor, MIT) for Step 2. Install its checkout and its own Python env, then
  set `BE_SAFETY_BEDICT_REPO` (holds `trained_models/`) and `BE_SAFETY_BEDICT_PYTHON`.
- **CRISPOR** (guide design + off-target search) for Step 2. Install its checkout and an indexed genome,
  then set `BE_SAFETY_CRISPOR_REPO`, `BE_SAFETY_CRISPOR_PYTHON` (optionally `BE_SAFETY_CRISPOR_GENOME`,
  `BE_SAFETY_CRISPOR_GENOMEDIR`).
- **AlphaGenome API key** for the optional Step 3 consequence scores (non-commercial). Set
  `ALPHAGENOME_API_KEY`. The key is read from the environment only; it is never logged, printed, or
  written to disk.
- **Network** on the first run, to fetch the reference sequence from UCSC (cached under `data/inputs/` after).

## Install and run
```
pip install -e . pytest          # or use a virtualenv
pytest -q                        # unit tests always run; integration tests skip if tools unset
python -m base_edit_safety.pipeline.spine                # Step 1: one-guide spine (demo)
python -m base_edit_safety.pipeline.enumerate_guides     # Steps 2+4: both promoters, full annotation
python -m base_edit_safety.pipeline.annotate_consequence # Step 3 (optional): AlphaGenome consequence
python -m base_edit_safety.pipeline.rank_designs         # Steps 5+6: ranking + credibility
python -m base_edit_safety.pipeline.design_records       # Main output: HBG design records
```
Outputs are written to `data/outputs/`.

## Full validation (with external tools)
The unit tests always run. The two integration tests (`test_bedict_integration`, `test_step2_integration`)
skip until both external tools are installed and pointed to:
```
# BE-DICT (edit-outcome predictor)
export BE_SAFETY_BEDICT_REPO=/path/to/bedict-checkout      # contains trained_models/
export BE_SAFETY_BEDICT_PYTHON=/path/to/bedict-env/bin/python
# CRISPOR (guide design + off-target search)
export BE_SAFETY_CRISPOR_REPO=/path/to/crispor-checkout    # contains crispor.py
export BE_SAFETY_CRISPOR_PYTHON=/path/to/crispor-env/bin/python
export BE_SAFETY_CRISPOR_GENOME=hg38                       # optional (default hg38)
export BE_SAFETY_CRISPOR_GENOMEDIR=/path/to/genomes        # optional (default <repo>/genomes)

pytest -q    # now also runs the integration tests
```
Re-running `enumerate_guides` then `annotate_consequence` then `design_records` with the tools configured
regenerates the real outputs after a code change. See [example-run.md](example-run.md) for what one run
produces.

The 12-gene clonal-hematopoiesis panel is built into `annotation.py`; per-mutation driver scores are out
of scope (if wired in later, the source URLs are in a comment there: boostDM-CH CC-BY-NC, AlphaMissense
CC-BY-4.0).
