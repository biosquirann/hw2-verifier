"""Свойства расчетов + проверка проверки: каждое свойство обязано ловить свой мутант."""
import random #для написания кода использовала ИИ-агента и правки
import pytest
from collections import Counter
from hypothesis import given, strategies as st
import seqstats as S
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

prot = st.text(alphabet=S.AA20, min_size=1, max_size=300)
dna = st.text(alphabet="ACGT", min_size=3, max_size=600)
CODONS = ["ATG", "TAA", "TAG", "TGA", "GCT", "AAA", "GGC", "CTG"]
# ДНК из малого набора кодонов: старт и три стопа в одной рамке встречаются часто
# (на равномерной ДНК ORF с внутренним TGA hypothesis почти не порождает)
codon_dna = st.lists(st.sampled_from(CODONS), min_size=1, max_size=100).map("".join)
dirty = st.text(alphabet="acgtACGT \n\t0123456789", max_size=200)
bodies = st.lists(st.integers(0, 3), min_size=3, max_size=30)
STOP = {"TAA", "TAG", "TGA"}


#внесенные ошибки, по одному на свойство 
def comp_drop_last(seq): return dict(Counter(seq[:-1]))
def comp_first_half(seq): return dict(Counter(seq[:len(seq) // 2]))
def gc_only_g(seq): return seq.count("G") / len(seq)
def revcomp_c_fixed(seq):
    return "".join({"A": "T", "T": "A", "G": "C", "C": "C"}[b] for b in reversed(seq))
def normalize_swapcase(seq): return seq.strip().swapcase()
def orfs_no_tga(seq, min_len):            # забыт стоп-кодон TGA
    res = []
    for i in range(len(seq) - 2):
        if seq[i:i + 3] == "ATG":
            for j in range(i, len(seq) - 2, 3):
                if seq[j:j + 3] in {"TAA", "TAG"}:
                    if j + 3 - i >= min_len: res.append(("+", i, j + 3, seq[i:j + 3]))
                    break
    return res
def orfs_frame0_only(seq, min_len): return [o for o in S.find_orfs(seq, min_len) if o[1] % 3 == 0]


#свойства как функции от реализации
def p1(comp, p):                      # сумма состава = длине
    assert sum(comp(p).values()) == len(p)
def p2(comp, p, k):                   # состав не зависит от порядка
    q = list(p); random.Random(k).shuffle(q)
    assert comp(p) == comp("".join(q))
def p3(gc, d):                        # физические границы доли
    assert 0 <= gc(d) <= 1 and gc("GC" * 3) == 1 and gc("AT" * 3) == 0
def p3b(gc, rc, d):                   # GC не зависит от цепи
    assert gc(rc(d)) == gc(d)
def p4(rc, d):                        # обратный комплемент — инволюция
    assert rc(rc(d)) == d
def p5(norm, s):                      # идемпотентность конвейера
    assert norm(norm(s)) == norm(s)
def p6(orfs, d):                      # ORF корректны
    for _, _, _, o in orfs(d, 9):
        assert o.startswith("ATG") and len(o) % 3 == 0 and o[-3:] in STOP
        assert not any(o[i:i + 3] in STOP for i in range(0, len(o) - 3, 3))
def p7(orfs, body_idx, pad):          # внедрённый ORF находится при любом сдвиге
    codons = ["GCT", "AAA", "GGC", "CTG"]
    orf = "ATG" + "".join(codons[i] for i in body_idx) + "TAA"
    assert any(o[3] == orf for o in orfs("C" * pad + orf + "CCCC", 9))


# ─── тесты свойств на настоящих функциях ──────────────────────────────────────
@given(prot)
def test_p1_composition_sum(p): p1(S.aa_composition, p)

@given(prot, st.integers(0, 10**6))
def test_p2_composition_order_invariant(p, k): p2(S.aa_composition, p, k)

@given(dna)
def test_p3_gc_bounds(d): p3(S.gc_content, d)

@given(dna)
def test_p3b_gc_strand_invariant(d): p3b(S.gc_content, S.revcomp, d)

@given(dna)
def test_p4_revcomp_involution(d): p4(S.revcomp, d)

@given(dirty)
def test_p5_normalize_idempotent(s): p5(S.normalize, s)

@given(codon_dna)
def test_p6_orfs_wellformed(d): p6(S.find_orfs, d)

@given(bodies, st.integers(0, 7))
def test_p7_planted_orf_found(b, pad): p7(S.find_orfs, b, pad)

@given(prot, dna)
def test_oracle_biopython(p, d):
    assert S.aa_composition(p) == {k: v for k, v in ProteinAnalysis(p).count_amino_acids().items() if v}
    assert abs(S.gc_content(d) - gc_fraction(d)) < 1e-12


# ─── проверка проверки: на мутанте свойство обязано упасть ────────────────────
def test_mutant_p1_caught():
    @given(prot)
    def run(p): p1(comp_drop_last, p)
    with pytest.raises(AssertionError): run()

def test_mutant_p2_caught():
    @given(prot, st.integers(0, 10**6))
    def run(p, k): p2(comp_first_half, p, k)
    with pytest.raises(AssertionError): run()

def test_mutant_p3_caught():
    @given(dna)
    def run(d): p3b(gc_only_g, S.revcomp, d)
    with pytest.raises(AssertionError): run()

def test_mutant_p4_caught():
    @given(dna)
    def run(d): p4(revcomp_c_fixed, d)
    with pytest.raises(AssertionError): run()

def test_mutant_p5_caught():
    @given(dirty)
    def run(s): p5(normalize_swapcase, s)
    with pytest.raises(AssertionError): run()

def test_mutant_p6_caught():
    @given(codon_dna)
    def run(d): p6(orfs_no_tga, d)
    with pytest.raises(AssertionError): run()

def test_mutant_p7_caught():
    @given(bodies, st.integers(0, 7))
    def run(b, pad): p7(orfs_frame0_only, b, pad)
    with pytest.raises(AssertionError): run()
