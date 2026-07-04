"""Step 1 spine: one anchor guide -> predicted edit spectrum -> annotated EditRecords.

Builds one reference-checked NGG-PAM guide that places an in-scope anchor position inside the editing
window, predicts its edit spectrum with the external bystander predictor, reduces it to per-position
EditRecords (on-target + bystander), marks the consequence as deferred (the Step 3 seam), and writes
the result under data/outputs. It needs the external predictor environment (see config.py). Step 2
(enumerate_guides) does the full guide search across the promoter windows.
"""
from __future__ import annotations

from ..adapters.bedict import BedictBystanderPredictor
from ..adapters.consequence import annotate_all
from ..io import outputs
from ..io import step0_anchor as anchor
from ..io.reference import fetch_window
from . import environment
from .model import (
    EDIT_SOURCE_BASE, EDIT_TARGET_BASE, EDIT_WINDOW, EditRecord, Guide,
    reverse_complement, spectrum_to_records,
)

PROTO_LEN = 20


def _gene_tss(gene: str) -> int:
    return anchor.HBG1_TSS if gene == "HBG1" else anchor.HBG2_TSS


def build_anchor_guide(gene: str, promoter_offset: int, editor_class: str,
                       editor_name: str | None = None) -> Guide:
    """Return one valid guide whose protospacer puts `promoter_offset` in the editing window."""
    tss = _gene_tss(gene)
    target = tss + promoter_offset
    edit_from, edit_to = EDIT_SOURCE_BASE[editor_class], EDIT_TARGET_BASE[editor_class]
    # Fetch enough context to cover any 20 nt protospacer + 3 nt PAM around the target on either strand.
    win = fetch_window("chr11", target - 30, target + 30)
    plus = win.plus_sequence
    lo = win.start_1based

    def plus_sub(a: int, b: int) -> str:  # inclusive genomic [a, b]
        return plus[a - lo:b - lo + 1]

    candidates = []
    # Plus-strand protospacers: 5' end at low coord g; PAM immediately 3' (higher coords).
    for g in range(lo, win.end_1based - PROTO_LEN - 2):
        proto = plus_sub(g, g + PROTO_LEN - 1)
        pam = plus_sub(g + PROTO_LEN, g + PROTO_LEN + 2)
        if len(proto) != PROTO_LEN or pam[1:] != "GG":
            continue
        tpos = target - g + 1
        if tpos in EDIT_WINDOW and proto[tpos - 1] == edit_from:
            candidates.append(("+", proto, pam, g, tpos))
    # Minus-strand protospacers: 5' end at high coord g; PAM immediately 3' (lower coords).
    for g in range(lo + PROTO_LEN + 2, win.end_1based + 1):
        proto = reverse_complement(plus_sub(g - PROTO_LEN + 1, g))
        pam = reverse_complement(plus_sub(g - PROTO_LEN - 2, g - PROTO_LEN))
        if len(proto) != PROTO_LEN or pam[1:] != "GG":
            continue
        tpos = g - target + 1
        if tpos in EDIT_WINDOW and proto[tpos - 1] == edit_from:
            candidates.append(("-", proto, pam, g, tpos))

    if not candidates:
        raise ValueError(
            f"No NGG-PAM guide places {gene} offset {promoter_offset} in the editing window."
        )
    # Prefer the placement nearest window center (position 6).
    strand, proto, pam, proto_5p, tpos = min(candidates, key=lambda c: abs(c[4] - 6))
    guide_id = f"{gene}_{'m' if promoter_offset >= 0 else 'p'}{abs(promoter_offset)}_{editor_class}"
    return Guide(
        guide_id=guide_id,
        protospacer=proto,
        pam=pam,
        strand=strand,
        genome_build=anchor.GENOME_BUILD,
        chrom=anchor.CHROM,
        gene=gene,
        gene_strand=anchor.GENE_STRAND,
        tss_1based=tss,
        editor_class=editor_class,
        editor_name=editor_name or "",
        edit_from=edit_from,
        edit_to=edit_to,
        proto_5p_genomic=proto_5p,
    )


def run_anchor_spine(
    gene: str = "HBG1",
    promoter_offset: int = 113,
    editor_class: str = "ABE",
    editor_name: str = "ABE8e",
    write: bool = True,
) -> list[EditRecord]:
    """Run the spine for one anchor guide and return the annotated EditRecords."""
    guide = build_anchor_guide(gene, promoter_offset, editor_class, editor_name)
    predictor = BedictBystanderPredictor()
    spectrum = predictor.predict(guide)
    records = spectrum_to_records(spectrum, on_target_offset=promoter_offset)
    records = annotate_all(records)  # deferred consequence (Step 3 seam)
    if write:
        repro = environment.reproducibility(
            tools=[predictor.tool_metadata()],
            editors={editor_class: editor_name},
            ensemble_runs=predictor.runs,
            reference_windows={},
        )
        outputs.write_records(
            records,
            filename=f"spine_{guide.guide_id}.json",
            meta={
                "step": "1 (spine on one guide)",
                "gene": gene,
                "anchor_offset": -abs(promoter_offset),
                "guide_id": guide.guide_id,
                "protospacer": guide.protospacer,
                "pam": guide.pam,
                "strand": guide.strand,
                "editor": editor_name,
                "predictor": spectrum.model_name,
                "reproducibility": repro,
            },
        )
    return records


if __name__ == "__main__":
    recs = run_anchor_spine()
    for r in recs:
        role = r.on_target_or_bystander
        print(
            f"{r.gene} offset {r.promoter_offset:>4} ({r.chrom}:{r.pos_1based}) "
            f"genomic {r.genomic_plus_ref}>{r.genomic_plus_alt} / sense "
            f"{r.promoter_sense_ref}>{r.promoter_sense_alt}  rate={r.predicted_rate}  [{role}]"
        )
