"""Известные значения и проверка обеих цепей: ловят ошибки, которые свойства P4-P6 пропускают."""
from hypothesis import given, strategies as st
import seqstats as S

ORF = "ATGGCTAAATAA"  # ATG GCT AAA TAA


def test_revcomp_known_values():
    assert S.revcomp("AAC") == "GTT"
    assert S.revcomp("ATGC") == "GCAT"


def test_normalize_known_value():
    assert S.normalize(" a1c\n gt\t2") == "ACGT"


@given(st.text(alphabet="ACGT", max_size=100))
def test_normalize_keeps_clean_input(s):
    assert S.normalize(s) == s


def test_orf_forward_found():
    assert any(o[0] == "+" and o[3] == ORF for o in S.find_orfs("CC" + ORF + "CC", 9))


def test_orf_minus_strand_found():
    seq = "CC" + S.revcomp(ORF) + "CC"
    assert any(o[0] == "-" and o[3] == ORF for o in S.find_orfs(seq, 9))


def test_orf_min_len_boundary():
    orf = "ATGAAATAA"  # ровно 9 нт
    seq = "CC" + orf + "CC"
    assert any(o[3] == orf for o in S.find_orfs(seq, 9))
    assert not any(o[3] == orf for o in S.find_orfs(seq, 10))
