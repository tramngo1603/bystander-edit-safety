# Example run (HBG1/HBG2)

One full run with all tools configured (CRISPOR + BE-DICT + AlphaGenome). Exact numbers depend on tool
versions and will vary.

**Input:** the two promoter windows (HBG1 `chr11:5,249,937-5,250,067`, HBG2 `chr11:5,254,861-5,254,991`),
CRISPOR-enumerated NGG guides, and the ABE8e and BE4max editors.

**Pipeline:** enumerate guides, predict edits (BE-DICT), flag CH-driver off-targets, score regulatory
consequence (AlphaGenome), rank, write design records.

**Output:** 158 edit records (36 on-target anchor edits, 122 bystanders) reduced to 30 candidate designs.
The only CH-driver gene hit by any off-target was DNMT3A.

**Findings:**
- Lead candidate: `HBG1/HBG2_60forw_ABE8e`, installing the -124/-123 pair, about 49% predicted editing,
  no CH off-target (flagged tier PREFERRED).
- A second -124/-123 pair (`61forw_ABE8e`) carries the DNMT3A off-target, so it needs off-target
  sequencing before use.
- The other 26 designs fall below the evidenced-activity efficacy floor (currently 16%, the
  preclinical-activity anchor; the separate 30% clinical-precedent floor is the higher bar the PREFERRED
  lead clears).
- All 30 designs map to a known natural fetal-hemoglobin edit (`analog_supported`: 22 with
  published-HSPC-strength support, 8 matching a clinical program), and all 30 are flagged for HBG1/HBG2
  paralog confirmation (long-read or paralog-resolving assays), since every one edits the duplicated
  promoter.

The full per-design records (evidence, off-target taxonomy, and the required-assay checklist) are written
to `data/outputs/design_records_HBG1_HBG2.json` (gitignored, not committed).
