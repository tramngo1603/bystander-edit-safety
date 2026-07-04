"""Hard build-stop: HBG1/HBG2 homology invariant. Failure must stop the build."""
HBG1_TSS_GRCH38 = 5_249_857
HBG2_TSS_GRCH38 = 5_254_781
OFFSET = 4924

# (promoter_offset, hbg1_pos, hbg2_pos) for in-scope positions; from docs/step0-coordinates.md
IN_SCOPE = [
    (198, 5_250_055, 5_254_979),
    (175, 5_250_032, 5_254_956),
    (124, 5_249_981, 5_254_905),
    (123, 5_249_980, 5_254_904),
    (117, 5_249_974, 5_254_898),
    (114, 5_249_971, 5_254_895),   # corrected: HBG2 -114 = 5,254,895
    (113, 5_249_970, 5_254_894),
]

def test_tss_separation():
    assert HBG2_TSS_GRCH38 - HBG1_TSS_GRCH38 == OFFSET

def test_paired_positions_differ_by_offset():
    for off, hbg1, hbg2 in IN_SCOPE:
        assert hbg2 - hbg1 == OFFSET, f"offset {off}: {hbg2-hbg1} != {OFFSET}"

def test_offset_arithmetic():
    for off, hbg1, hbg2 in IN_SCOPE:
        assert hbg1 == HBG1_TSS_GRCH38 + off
        assert hbg2 == HBG2_TSS_GRCH38 + off
