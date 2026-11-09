"""Договор о схеме FASTA: на синтетике и на реальных файлах из data/."""
import re, pathlib
import pytest
from seqstats import read_fasta, AA20

ROOT = pathlib.Path(__file__).resolve().parents[1]
FILES = sorted((ROOT / "data").glob("*.fasta"))
ACC = re.compile(r"^>(sp|tr)\|[A-Z0-9]{6,10}\||^>[A-Z]{1,2}_?\d+")


def violations(header, seq, kind="protein"):
    """Список кодов нарушений договора (пустой — договор выполнен)."""
    v = []
    if not header or not ACC.match(header): v.append("no_accession")
    if not seq: v.append("empty"); return v
    alphabet = set(AA20) if kind == "protein" else set("ACGT")
    if not set(seq) <= alphabet: v.append("bad_alphabet")
    return v


def test_contract_ok(tmp_path):
    f = tmp_path / "x.fasta"; f.write_text(">sp|P37840|SYUA_HUMAN x\nMDVF\nMKGL\n")
    h, s = read_fasta(f); assert violations(h, s) == [] and len(s) == 8


@pytest.mark.parametrize("seq,code", [("MDVX", "bad_alphabet"), ("", "empty")])
def test_contract_catches(tmp_path, seq, code):
    f = tmp_path / "x.fasta"; f.write_text(">sp|P37840|SYUA_HUMAN x\n" + seq + "\n")
    assert code in violations(*read_fasta(f))


def test_contract_catches_missing_accession(tmp_path):
    f = tmp_path / "x.fasta"; f.write_text(">unknown\nMDVF\n")
    assert "no_accession" in violations(*read_fasta(f))


@pytest.mark.skipif(not FILES, reason="нет data/*.fasta")
@pytest.mark.parametrize("path", FILES, ids=lambda p: p.name)
def test_contract_real_data(path):
    h, s = read_fasta(path)
    kind = "dna" if "dna" in path.name.lower() or "nuc" in path.name.lower() else "protein"
    assert violations(h, s, kind) == [], path.name


@pytest.mark.skipif(not (ROOT / "data/SNCA_P37840.fasta").exists(), reason="нет SNCA")
def test_snca_length_140():
    assert len(read_fasta(ROOT / "data/SNCA_P37840.fasta")[1]) == 140
#код написан ИИ-агентом
