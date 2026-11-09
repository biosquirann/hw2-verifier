"""Вносит по одной ошибке в src/seqstats.py и проверяет, что checks/ падает. Файл всегда восстанавливается."""
import subprocess, sys
P = "src/seqstats.py"
src = open(P).read()
MUT = {
 "revcomp без разворота": ('for b in reversed(seq))', 'for b in seq)'),
 "normalize без upper": ('.upper()\n\n\ndef aa_comp', '\n\n\ndef aa_comp'),
 "normalize не убирает цифры": (r're.sub(r"[\s\d]+", "", seq)', r're.sub(r"[\s]+", "", seq)'),
 "ORF: min_len игнорируется": ('if j + 3 - i >= min_len:', 'if True:'),
 "ORF: min_len сдвинут на 1": ('if j + 3 - i >= min_len:', 'if j + 3 - i > min_len:'),
 "ORF: нет минус-цепи": ('    res += [("-", s, e, rc[s:e]) for s, e in _orfs_forward(rc, min_len)]\n', '\n'),
 "ORF: неверные координаты минус-цепи": ('rc[s:e]) for s, e in _orfs_forward(rc, min_len)]', 'rc[s:e]) for s, e in _orfs_forward(rc[::-1], min_len)]'),
}
try:
    for name, (a, b) in MUT.items():
        assert a in src, f"не найден фрагмент для мутации: {name}"
        open(P, "w").write(src.replace(a, b, 1))
        r = subprocess.run([sys.executable, "-m", "pytest", "checks", "-q", "-x", "-p", "no:cacheprovider"], capture_output=True)
        print(f"{name}: {'ПОЙМАН' if r.returncode else 'ПРОПУЩЕН'}")
finally:
    open(P, "w").write(src)
