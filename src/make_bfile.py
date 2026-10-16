"""Write bfile/b400741.txt from results/orbit_count.txt, after checking it
against the published terms and against results/line_count.txt.

usage: python3 src/make_bfile.py
"""
PUBLISHED = [1, 1, 1, 1, 2, 3, 6, 7, 16, 25, 57, 68, 210, 293, 541, 1017,
             2570, 2778, 9198, 11310, 27530, 47906]  # A400741 data, 2026-10-16


def load(path):
    values = {}
    for line in open(path):
        n, value = map(int, line.split()[:2])
        if values.get(n, value) != value:
            raise SystemExit(f"{path}: two different values for n = {n}")
        values[n] = value
    return values


orbit = load("results/orbit_count.txt")
line = load("results/line_count.txt")
last = max(orbit)
assert sorted(orbit) == list(range(1, last + 1)), "gap in orbit_count results"
for n, value in enumerate(PUBLISHED, 1):
    assert orbit[n] == value, f"n = {n} disagrees with OEIS"
for n, value in line.items():
    assert orbit[n] == value, f"n = {n}: orbit_count and line_count disagree"
print(f"published terms 1..{len(PUBLISHED)} match; "
      f"line_count agrees on n = 1..{max(line)}; writing 1..{last}")

with open("bfile/b400741.txt", "w") as out:
    for n in range(1, last + 1):
        out.write(f"{n} {orbit[n]}\n")
