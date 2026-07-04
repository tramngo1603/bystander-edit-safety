# Test fixtures

These are **SYNTHETIC, hand-authored** fixtures, not captured real-tool outputs. They exist to test
the CRISPOR/BE-DICT parsers, validation, and strand-aware mapping deterministically without installing
the external tools. They are representative in shape (columns, JSON structure, value ranges), but the
numbers are made up.

They are **not** regression fixtures. Recorded real-tool outputs (a captured CRISPOR run over the
HBG1/HBG2 windows and a captured BE-DICT ensemble response) are a follow-up: they should live here
alongside these, clearly named as real captures, and back a separate regression test once available.
