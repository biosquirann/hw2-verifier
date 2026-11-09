"""Проверка на реальных данных проекта: белки из data/ids.csv."""
import csv, pathlib
import pytest
import seqstats as S

ROOT = pathlib.Path(__file__).resolve().parents[1]
ROWS = list(csv.DictReader((ROOT / "data" / "ids.csv").open()))
EXPECTED = {r["gene"]: int(r["length"]) for r in csv.DictReader((ROOT / "checks/expected_lengths.csv").open())}


def load(row):
    path = ROOT / "data" / f"{row['gene']}_{row['uniprot_id']}.fasta"
    assert path.exists(), f"нет файла {path.name}"
    return S.read_fasta(path)


@pytest.mark.parametrize("row", ROWS, ids=lambda r: r["gene"])
def test_length_matches_uniprot(row):
    _, seq = load(row)
    assert len(seq) == EXPECTED[row["gene"]]


@pytest.mark.parametrize("row", ROWS, ids=lambda r: r["gene"])
def test_header_and_alphabet(row):
    header, seq = load(row)
    assert f"|{row['uniprot_id']}|" in header, header
    assert f"GN={row['gene']} " in header, header
    assert set(seq) <= set(S.AA20), f"недопустимые символы: {set(seq) - set(S.AA20)}"


@pytest.mark.parametrize("row", ROWS, ids=lambda r: r["gene"])
def test_composition_matches_biopython(row):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    _, seq = load(row)
    mine = S.aa_composition(seq)
    assert sum(mine.values()) == len(seq)
    theirs = {k: v for k, v in ProteinAnalysis(seq).count_amino_acids().items() if v}
    assert mine == theirs


def snca_row():
    return next(r for r in ROWS if r["gene"] == "SNCA")


def test_snca_matches_diary_table():
    import json
    _, seq = load(snca_row())
    expected = json.loads((ROOT / "checks/expected_snca_diary.json").read_text())
    assert S.aa_composition(seq) == expected


def test_analyze_snca_script_agrees_with_seqstats():
    import analyze_snca
    _, seq = load(snca_row())
    result = analyze_snca.analyze(analyze_snca.parse_fasta(analyze_snca.FASTA_PATH))
    assert result["length"] == len(seq)
    mine = S.aa_composition(seq)
    assert set(result["composition"]) == set(mine)
    for aa, n in mine.items():
        assert result["composition"][aa] == round(n / len(seq), 4)
