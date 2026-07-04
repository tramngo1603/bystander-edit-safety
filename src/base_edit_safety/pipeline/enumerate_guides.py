"""Step 2: enumerate guides across both promoter windows and predict each guide's edit spectrum.

For each HBG1 and HBG2 proximal-promoter window this does the following. It enumerates NGG guides
with the external guide-design tool, which runs a genome-wide off-target search. It predicts every
editable guide's outcome spectrum with the external bystander predictor. It reduces these to
per-position EditRecords, classified against the locked in-scope anchor set. It attaches the
off-target summary, marks paralog co-targeting, marks the consequence as deferred (the Step 3 seam),
and writes the combined result to data/outputs.
"""
from __future__ import annotations

from collections import defaultdict

from ..adapters.bedict import BedictBystanderPredictor
from ..adapters.consequence import annotate_all
from ..adapters.crispor import CrisporGuideDesigner
from ..annotation import annotate_records as annotate_ch, ch_driver_offtargets
from ..io import outputs
from ..io import reference
from ..io import step0_anchor as anchor
from . import environment
from .model import EditRecord, spectrum_to_records

# Step 0 proximal-promoter modelling windows (-210..-80 vs cap), GRCh38. These HBG1/HBG2 coordinates
# are the intentional, hard-coded product scope of this prototype (see documentation/step0-coordinates.md),
# not a placeholder to be generalized.
WINDOWS = {
    "HBG1": ("chr11", 5_249_937, 5_250_067),
    "HBG2": ("chr11", 5_254_861, 5_254_991),
}
EDITOR_NAME = {"ABE": "ABE8e", "CBE": "BE4max"}


def _anchor_coords_by_gene_editor() -> dict[tuple[str, str], set[int]]:
    """Genomic coordinates of in-scope anchor edits, keyed by (gene, editor_class)."""
    coords: dict[tuple[str, str], set[int]] = defaultdict(set)
    for pos in anchor.IN_SCOPE:
        coords[("HBG1", pos.editor_class)].add(pos.hbg1_pos)
        coords[("HBG2", pos.editor_class)].add(pos.hbg2_pos)
    return coords


def run_enumeration(write: bool = True) -> list[EditRecord]:
    """Enumerate guides over both promoters and return all annotated EditRecords."""
    anchor.check_homology_invariant()  # G8 hard gate before any downstream work
    designer = CrisporGuideDesigner()
    predictor = BedictBystanderPredictor()
    anchor_coords = _anchor_coords_by_gene_editor()

    # Enumerate candidates per gene and record which protospacers occur in each promoter (co-targeting).
    candidates_by_gene = {}
    protospacers_by_gene: dict[str, set[str]] = {}
    ambiguous: list[str] = []   # CRISPOR targets that mapped to >1 window site (dropped, but reported)
    for gene, (chrom, start, end) in WINDOWS.items():
        cands = designer.design(chrom, start, end, gene)
        ambiguous.extend(getattr(designer, "last_ambiguous", []))
        candidates_by_gene[gene] = cands
        protospacers_by_gene[gene] = {c.protospacer for c in cands}

    all_records: list[EditRecord] = []
    for gene, cands in candidates_by_gene.items():
        other_gene = "HBG2" if gene == "HBG1" else "HBG1"
        # Build editor-specific guides for every candidate that can be edited.
        editor_guides = []
        candidate_by_guide_id = {}
        for cand in cands:
            for editor_class in ("ABE", "CBE"):
                if cand.has_editable_base(editor_class):
                    guide = cand.for_editor(editor_class, EDITOR_NAME[editor_class])
                    editor_guides.append(guide)
                    candidate_by_guide_id[guide.guide_id] = cand

        spectra = predictor.predict_many(editor_guides)

        for guide in editor_guides:
            cand = candidate_by_guide_id[guide.guide_id]
            spectrum = spectra[guide.guide_id]
            on_coords = anchor_coords.get((gene, guide.editor_class), set())
            records = spectrum_to_records(spectrum, on_target_coords=on_coords)
            co_targets = guide.protospacer in protospacers_by_gene[other_gene]
            ch_offtargets = ch_driver_offtargets(cand.offtarget_sites)  # CH-driver off-target hits
            for record in records:
                record.provenance = {
                    **record.provenance,
                    "offtarget": cand.offtarget,
                    "co_targets_paralog": co_targets,
                    "co_target_gene": other_gene if co_targets else None,
                    "ch_driver_offtargets": ch_offtargets,
                }
            all_records.extend(records)

    all_records = annotate_all(all_records)   # deferred consequence (Step 3 seam)
    all_records = annotate_ch(all_records)    # Step 4: per-edit CH-driver flag (static)

    if write:
        on_target_n = sum(1 for r in all_records if r.on_target_or_bystander == "on_target")
        ch_genes_hit = sorted({
            h["ch_driver_gene"]
            for r in all_records
            for h in r.provenance.get("ch_driver_offtargets", [])
        })
        repro = environment.reproducibility(
            tools=[designer.tool_metadata(), predictor.tool_metadata()],
            editors=EDITOR_NAME,
            ensemble_runs=predictor.runs,
            reference_windows={g: reference.window_provenance(c, s, e)
                               for g, (c, s, e) in WINDOWS.items()},
        )
        outputs.write_records(
            all_records,
            filename="enumeration_HBG1_HBG2.json",
            meta={
                "step": "2 (locus->guides->spine) + 4 (CH-driver annotation)",
                "windows": WINDOWS,
                "guide_design_tool": designer.name,
                "predictor": predictor.name,
                "editors": EDITOR_NAME,
                "record_count": len(all_records),
                "on_target_records": on_target_n,
                "bystander_records": len(all_records) - on_target_n,
                "ch_driver_set": "boostDM-CH 12-gene panel (CC-BY-NC)",
                "ch_driver_offtarget_genes_hit": ch_genes_hit,
                "ambiguous_candidate_count": len(ambiguous),   # guides dropped for multi-site mapping
                "ambiguous_candidates": ambiguous,             # the search space was incomplete by this many
                "reproducibility": repro,
            },
        )
    return all_records


if __name__ == "__main__":
    recs = run_enumeration()
    on = [r for r in recs if r.on_target_or_bystander == "on_target"]
    # CH-driver off-target summary (the safety signal): guides with off-targets in CH genes.
    ch_hits = {}
    for r in recs:
        for h in r.provenance.get("ch_driver_offtargets", []):
            gid = r.provenance.get("guide_id")
            ch_hits.setdefault((r.gene, gid, h["ch_driver_gene"]),
                               (h["chrom"], h["start"], h["mismatch_count"]))
    print(f"total records: {len(recs)}  |  on-target (anchor) records: {len(on)}  "
          f"|  edits in CH-driver genes: {sum(1 for r in recs if r.ch_driver_flag)}")
    print("CH-driver OFF-TARGET hits (guide -> CH gene):")
    for (gene, gid, ch_gene), (chrom, start, mm) in sorted(ch_hits.items()):
        print(f"  {gid} ({gene}) -> off-target in {ch_gene} at {chrom}:{start} ({mm} mismatches)")
    print(f"\ntotal records: {len(recs)}  |  on-target (anchor) records: {len(on)}")
    print("on-target anchor edits found (guide -> anchor):")
    seen = set()
    for r in sorted(on, key=lambda r: (r.gene, r.promoter_offset)):
        key = (r.gene, r.promoter_offset, r.provenance.get("guide_id"))
        if key in seen:
            continue
        seen.add(key)
        ot = r.provenance.get("offtarget", {})
        print(
            f"  {r.gene} {-abs(r.promoter_offset):>4} ({r.chrom}:{r.pos_1based}) "
            f"{r.editor_class} rate={r.predicted_rate} via {r.provenance.get('guide_id')} "
            f"| co-target={r.provenance.get('co_targets_paralog')} "
            f"| offtargets={ot.get('offtarget_count')} mitSpec={ot.get('mit_specificity')}"
        )
