import io
import json
from contextlib import redirect_stdout
from types import SimpleNamespace

from scrnaseq_tools.commands import scan as scan_command
from scrnaseq_tools.commands import summary as summary_command


def test_scan_lists_supported_files(tmp_path):
    (tmp_path / "a.h5ad").write_text("x", encoding="utf-8")
    (tmp_path / "b.csv").write_text("x", encoding="utf-8")
    (tmp_path / "ignore.md").write_text("x", encoding="utf-8")

    args = SimpleNamespace(path=str(tmp_path), recursive=False, json=False)
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        rc = scan_command.run(args)

    output = buffer.getvalue()
    assert rc == 0
    assert str(tmp_path / "a.h5ad") in output
    assert str(tmp_path / "b.csv") in output
    assert "ignore.md" not in output


def test_scan_json_output(tmp_path):
    nested = tmp_path / "nested"
    nested.mkdir()
    (nested / "matrix.mtx").write_text("x", encoding="utf-8")

    args = SimpleNamespace(path=str(tmp_path), recursive=True, json=True)
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        rc = scan_command.run(args)

    payload = json.loads(buffer.getvalue())
    assert rc == 0
    assert payload["recursive"] is True
    assert payload["count"] == 1
    assert str(nested / "matrix.mtx") in payload["files"]


def test_summary_json_for_generic_file(tmp_path):
    target = tmp_path / "data.bin"
    target.write_bytes(b"abc")

    args = SimpleNamespace(path=str(target), json=True)
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        rc = summary_command.run(args)

    payload = json.loads(buffer.getvalue())
    assert rc == 0
    assert payload["type"] == "file"
    assert payload["path"] == str(target)
    assert payload["size_bytes"] == 3


def test_summary_missing_file_returns_error(tmp_path):
    target = tmp_path / "missing.h5ad"
    args = SimpleNamespace(path=str(target), json=False)
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        rc = summary_command.run(args)

    output = buffer.getvalue()
    assert rc == 1
    assert f"not found: {target}" in output
