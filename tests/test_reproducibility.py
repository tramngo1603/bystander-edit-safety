"""Reproducibility metadata helpers: checksums, cache provenance, tool + run environment."""
from base_edit_safety.io import reference
from base_edit_safety.pipeline import environment


def test_sequence_checksum_stable_and_sensitive():
    c = reference.sequence_checksum("ACGTACGT")
    assert c == reference.sequence_checksum("ACGTACGT") and len(c) == 16
    assert reference.sequence_checksum("ACGTACGA") != c


def test_window_provenance_reports_uncached():
    prov = reference.window_provenance("chr11", 999_000, 999_130)  # not in the cache
    assert prov["cached"] is False
    assert prov["cache_key"] == "hg38_chr11_999000_999130.json"


def test_tool_environment_drops_missing_and_keeps_set_fields():
    env = environment.tool_environment("BE-DICT", repo_root=None, python_executable="/env/py")
    assert env["name"] == "BE-DICT" and env["python_executable"] == "/env/py"
    assert "repo_root" not in env and "genome_index" not in env  # None fields dropped


def test_git_commit_never_raises():
    assert environment.git_commit(None) is None
    assert environment.git_commit("/no/such/path/xyz") is None


def test_tool_environment_records_failed_version_lookup_visibly():
    env = environment.tool_environment("CRISPOR", repo_root="/no/such/repo/xyz", python_executable="/x/py")
    assert env["repo_root"] == "/no/such/repo/xyz"                # exact path kept
    assert env["python_executable"] == "/x/py"
    assert env["commit"] is None and env["commit_lookup"] == "unavailable"  # failure visible, not dropped


def test_reproducibility_block_shape():
    block = environment.reproducibility(
        tools=[{"name": "CRISPOR"}], editors={"ABE": "ABE8e"}, ensemble_runs=5,
        reference_windows={"HBG1": {"cached": False}},
    )
    assert block["genome_build"] == "GRCh38"
    assert block["ensemble_runs"] == 5
    assert block["tools"][0]["name"] == "CRISPOR"
    assert block["reference_windows"]["HBG1"]["cached"] is False
