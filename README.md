# base-edit-safety

A research tool for choosing base-editing guides at the **HBG1/HBG2 fetal-hemoglobin promoters**, the
target used to switch fetal hemoglobin back on in sickle cell disease and beta-thalassemia.

You give it candidate guide + editor designs. For each one, it writes a record you can read before
doing any lab work: what the edit does, what could go wrong, and which lab tests would confirm it. It
does not reduce a design to a single "safety" number.

Scope is narrow on purpose. The HBG1/HBG2 coordinates and the in-scope edits are fixed in the code
(see `docs/step0-coordinates.md` and `io/step0_anchor.py`), not something to generalize yet in this first version of the tool. This helps you pick which designs to test first. It is not a general editing platform and not a clinical tool.

## What each record contains
For every candidate guide + editor:
- the intended edit, and the predicted on-target and nearby bystander edits;
- where else in the genome the guide could match (off-targets), split by kind (Cas9-style sequence
  match, base-editor-specific, guide-independent, and RNA), so a sequence-match guess is never treated
  as a measured off-target;
- whether the guide sits in the near-identical HBG1/HBG2 region, which can cause large deletions and
  needs extra checking;
- whether any off-target falls in a clonal-hematopoiesis driver gene (a documentation flag, not a risk
  estimate);
- how strong the evidence is: predicted only, similar to a known edit, or actually measured;
- how the edit relates to known natural fetal-hemoglobin variants and the clinical HBG programs;
- the lab tests you need to run to confirm it.

There is one safety score, but it is a rough tiebreaker for ordering, not the answer. Read the evidence
and the required tests first.

## Who this is for
A researcher deciding which guide + editor to take into the lab. You have an HBG1/HBG2 promoter edit to
install and many candidate guides, but lab work (editing CD34+ cells, amplicon sequencing, off-target
sequencing) is slow and expensive, so you can only test a few at a time. This tool writes down the case for each candidate first, so you can see the predicted safety profiles of such candidate before moving forward. It does not necessarily replace lab work.

## Status
Research prototype for HBG1/HBG2 promoter base-editing. Outputs are predicted records for lab
confirmation (for example by amplicon NGS or CHANGE-seq-BE). Every candidate is marked with how strong its evidence is (`predicted_only`, `analog_supported`, or `empirically_measured`). Every predicted number is labelled model-predicted, _not_ measured.

## Pipeline
1. Target window + edit positions (`io/step0_anchor.py`), taken from `docs/step0-coordinates.md`.
2. Enumerate guides over both promoters and predict each guide's edits (`pipeline/enumerate_guides.py`).
3. Optional: add regulatory-consequence scores with AlphaGenome (`pipeline/annotate_consequence.py`).
   Without this step each record's consequence stays `deferred`. Coding-missense scoring is out of scope.
4. Flag off-targets that fall in a clonal-hematopoiesis driver gene (`annotation.py`).
5. Score and rank designs by the safety rubric (`scoring/rubric.py`, run via `pipeline/rank_designs.py`).
6. Check the ranking differs from naive sorts, and add per-design uncertainty (`scoring/credibility.py`).
7. Write the design records: the main output (`design_record.py`, run via `pipeline/design_records.py`).

## External tools
You install the prediction tools yourself, under their own licenses; thin adapters call them
(`src/base_edit_safety/adapters/`). They are not bundled here. Point to them with environment variables:
- Edit-outcome predictor: `BE_SAFETY_BEDICT_REPO`, `BE_SAFETY_BEDICT_PYTHON`
- Guide design / off-target search: `BE_SAFETY_CRISPOR_REPO`, `BE_SAFETY_CRISPOR_PYTHON`
  (optionally `BE_SAFETY_CRISPOR_GENOME`, `BE_SAFETY_CRISPOR_GENOMEDIR`)

The 12-gene clonal-hematopoiesis panel is built into `annotation.py`. Per-mutation driver scores are out
of scope; if you wire them in later, the source URLs are noted in a comment there (boostDM-CH, CC-BY-NC;
AlphaMissense, CC-BY-4.0).

## Running
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

### Full validation (with external tools)
The unit tests always run. The two integration tests (`test_bedict_integration`,
`test_step2_integration`) skip until both external tools are installed and pointed to:
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
Reference sequence is fetched from UCSC on the first run (needs network) and cached under `data/inputs/`.

## Related
Off-target search here uses the reference genome only (CRISPOR, up to 4 mismatches, no bulges, no
population variants). For off-target search that accounts for population and ancestry variants on the
same HBG1/HBG2 guides (CRISPRme with gnomAD and 1000 Genomes), see the companion project
`ancestry-aware-offtarget` which is still ongoing under my Github account.

## License
Apache-2.0. Some external tools are non-commercial or copyleft; using them means resolving those
upstream terms.
