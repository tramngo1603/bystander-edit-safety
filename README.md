# Safety triage for HBG1/HBG2 fetal-hemoglobin base-editing guides

A research tool for choosing base-editing guides at the HBG1/HBG2 fetal-hemoglobin promoters, the target
used to switch fetal hemoglobin back on in sickle cell disease and beta-thalassemia.

You give it candidate guide + editor designs. For each one, it writes a record you can read before doing
any lab work: what the edit does, what could go wrong, and which lab tests would confirm it. It does not
reduce a design to a single "safety" number.

Scope is narrow on purpose. The HBG1/HBG2 coordinates and the in-scope edits are fixed in the code (see
[documentation/step0-coordinates.md](documentation/step0-coordinates.md) and `io/step0_anchor.py`), not something to
generalize yet in this first version of the tool. This helps you pick which designs to test first. It is
not a general editing platform and not a clinical tool.

## What each record contains
For every candidate guide + editor:
- the intended edit, and the predicted on-target and nearby bystander edits;
- where else in the genome the guide could match (off-targets), split by kind (Cas9-style sequence match,
  base-editor-specific, guide-independent, and RNA), so a sequence-match guess is never treated as a
  measured off-target;
- whether the guide sits in the near-identical HBG1/HBG2 region, which can cause large deletions and needs
  extra checking;
- whether any off-target falls in a clonal-hematopoiesis driver gene (a documentation flag, not a risk
  estimate);
- how strong the evidence is: predicted only, similar to a known edit, or actually measured;
- how the edit relates to known natural fetal-hemoglobin variants and the clinical HBG programs;
- the lab tests you need to run to confirm it.

## Who this is for
A researcher deciding which guide + editor might be flagged for safety before testing further. You have an
HBG1/HBG2 promoter edit to install and many candidate guides, but lab work (editing CD34+ cells, amplicon
sequencing, off-target sequencing) is slow and expensive, so you can only test a few at a time. This tool
writes down the case for each candidate first, so you can see the predicted safety profiles of such
candidate before moving forward with lab work.

## Status
Research prototype for HBG1/HBG2 promoter base-editing. Outputs are predicted records for lab confirmation
(for example by amplicon NGS or CHANGE-seq-BE). Every candidate is marked with how strong its evidence is
(`predicted_only`, `analog_supported`, or `empirically_measured`). Every predicted number is labelled
model-predicted, _not_ measured.

## Quick start
```
pip install -e . pytest    # or use a virtualenv
pytest -q                  # unit tests; integration tests skip until the external tools are set up
python -m base_edit_safety.pipeline.design_records    # main output: HBG design records -> data/outputs/
```
Installing the external prediction tools and running the full pipeline is covered in
[documentation/setup.md](documentation/setup.md).

## Documentation
- [Setup and running](documentation/setup.md) - requirements, external tools, the AlphaGenome API key, and full validation.
- [Pipeline and scope](documentation/pipeline.md) - the seven steps, the off-target scope, and related work.
- [Example run](documentation/example-run.md) - a real HBG1/HBG2 run with its input, output, and findings.
- [Coordinates](documentation/step0-coordinates.md) - the HBG1/HBG2 windows and the in-scope edits.

## License
Apache-2.0. Some external tools are non-commercial or copyleft; using them means resolving those upstream
terms.
