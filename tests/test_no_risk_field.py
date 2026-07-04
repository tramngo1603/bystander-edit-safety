"""Guardrail G3: the data model and serialized output carry NO clonal/lifetime-risk field.

The clonal-hematopoiesis layer is annotation only. Any field name implying a risk-over-time,
clonal-expansion, or leukemia/MDS prediction is forbidden in the schema and outputs.
"""
import dataclasses

from base_edit_safety.pipeline.model import EditRecord
from base_edit_safety.io.outputs import records_to_dicts

FORBIDDEN_SUBSTRINGS = (
    "clonal", "expansion", "leukemia", "leukaemia", "mds",
    "lifetime", "risk", "trajectory", "incidence",
)


def _record():
    return EditRecord(
        genome_build="GRCh38", chrom="chr11", pos_1based=5_249_970,
        genomic_plus_ref="T", genomic_plus_alt="C",
        promoter_sense_ref="A", promoter_sense_alt="G",
        gene="HBG1", gene_strand="-", promoter_offset=113,
        editor_class="ABE", on_target_or_bystander="on_target",
    )


def test_dataclass_fields_have_no_risk_field():
    names = {f.name.lower() for f in dataclasses.fields(EditRecord)}
    for field_name in names:
        assert not any(bad in field_name for bad in FORBIDDEN_SUBSTRINGS), field_name


def test_serialized_record_keys_have_no_risk_field():
    keys = set()
    for d in records_to_dicts([_record()]):
        keys.update(k.lower() for k in d.keys())
    for key in keys:
        assert not any(bad in key for bad in FORBIDDEN_SUBSTRINGS), key


def test_ch_layer_is_flag_and_gene_only():
    field_names = {f.name for f in dataclasses.fields(EditRecord)}
    assert "ch_driver_flag" in field_names
    assert "ch_driver_gene" in field_names
    # nothing beyond a boolean flag and a gene symbol for the CH layer
    ch_fields = {n for n in field_names if n.startswith("ch_")}
    assert ch_fields == {"ch_driver_flag", "ch_driver_gene"}
