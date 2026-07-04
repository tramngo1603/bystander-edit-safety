# Pipeline and scope

## Steps
1. Target window + edit positions (`io/step0_anchor.py`), taken from
   [step0-coordinates.md](step0-coordinates.md).
2. Enumerate guides over both promoters and predict each guide's edits (`pipeline/enumerate_guides.py`).
3. Optional: add regulatory-consequence scores with AlphaGenome (`pipeline/annotate_consequence.py`).
   Without this step each record's consequence stays `deferred`. Coding-missense scoring is out of scope.
4. Flag off-targets that fall in a clonal-hematopoiesis driver gene (`annotation.py`).
5. Score and rank designs by the safety rubric (`scoring/rubric.py`, run via `pipeline/rank_designs.py`).
6. Check the ranking differs from naive sorts, and add per-design uncertainty (`scoring/credibility.py`).
7. Write the design records: the main output (`design_record.py`, run via `pipeline/design_records.py`).

Installing the external tools and running these steps is covered in [setup.md](setup.md).

## Off-target scope and related work
Off-target search here uses the reference genome only (CRISPOR, up to 4 mismatches, no bulges, no
population variants). For off-target search that accounts for population and ancestry variants on the
same HBG1/HBG2 guides (CRISPRme with gnomAD and 1000 Genomes), see the companion project
`ancestry-aware-offtarget`, still ongoing under my GitHub account.
