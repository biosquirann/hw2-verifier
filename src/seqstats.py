"""Подсчеты по последовательностям: нормализация, состав, GC, обратный комплемент, ORF."""
import re
from collections import Counter

AA20 = "ACDEFGHIKLMNPQRSTVWY"
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
STOPS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line
            elif line:
                seq.append(line)
    return header, normalize("".join(seq))


def normalize(seq):
    """Конвейер очистки: убрать пробелы, цифры и переводы строк, привести к верхнему регистру."""
    return re.sub(r"[\s\d]+", "", seq).upper()


def aa_composition(seq):
    return dict(Counter(seq))


def gc_content(seq):
    return seq.count("G") / len(seq)


def revcomp(seq):
    return "".join(COMP[b] for b in reversed(seq))


def _orfs_forward(seq, min_len):
    out = []
    for i in range(len(seq) - 2):
        if seq[i:i + 3] != "ATG":
            continue
        for j in range(i, len(seq) - 2, 3):
            if seq[j:j + 3] in STOPS:
                if j + 3 - i >= min_len:
                    out.append((i, j + 3))
                break
    return out

#код написан ии агентом + мои правки
def find_orfs(seq, min_len=30):
    """ORF от ATG до стопа включительно в обеих цепях: (цепь, начало, конец, орф)."""
    res = [("+", s, e, seq[s:e]) for s, e in _orfs_forward(seq, min_len)]
    rc = revcomp(seq)
    res += [("-", s, e, rc[s:e]) for s, e in _orfs_forward(rc, min_len)]
    return res
