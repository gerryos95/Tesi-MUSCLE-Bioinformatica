"""
Script di verifica per il Capitolo 3 della tesi (MUSCLE, Fase 1 - Draft Progressive Alignment).

Ricalcola da zero, e stampa, ogni numero citato nel capitolo:
  - distanza k-mer (estratto a mano + sequenze complete)
  - clustering UPGMA
  - allineamento progressivo (Needleman-Wunsch, BLOSUM62)
  - identita' percentuali di controllo

Dipendenze:  pip install biopython --break-system-packages

Sequenze reali (fonti):
  Umano   - UniProt P99999 (Homo sapiens)
  Cavallo - forma matura, Equus caballus (cfr. BMRB entry 274)
  Lievito - UniProt P00044, iso-1-citocromo c (Saccharomyces cerevisiae)
"""

from itertools import combinations
from Bio.Align import substitution_matrices, PairwiseAligner

SEQS = {
    "Umano":   "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE",
    "Cavallo": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "Lievito": "TEFKAGSAKKGATLFKTRCLQCHTVEKGGPHKVGPNLHGIFGRHSGQAEGYSYTDANIKKNVLWDENNMSEYLTNPKKYIPGTKMAFGGLKKEKDRNDLITYLKKACE",
}


# ----------------------------------------------------------------------
# 1. DISTANZA K-MER
# ----------------------------------------------------------------------
def kmer_set(seq, k=3):
    return set(seq[i:i + k] for i in range(len(seq) - k + 1))


def kmer_distance(a, b, k=3):
    ka, kb = kmer_set(a, k), kmer_set(b, k)
    common = ka & kb
    frac_common = len(common) / min(len(ka), len(kb))
    return 1 - frac_common


def print_section(title):
    print("\n" + "=" * len(title))
    print(title)
    print("=" * len(title))


print_section("1. ESTRATTO A MANO (primi 12 residui, Umano vs Cavallo)")
h12, e12 = SEQS["Umano"][:12], SEQS["Cavallo"][:12]
print("Umano  [1:12] =", h12)
print("Cavallo[1:12] =", e12)
kh = [h12[i:i + 3] for i in range(len(h12) - 2)]
ke = [e12[i:i + 3] for i in range(len(e12) - 2)]
common = set(kh) & set(ke)
print("3-meri comuni:", sorted(common), f"({len(common)} su {min(len(kh), len(ke))})")
print(f"Distanza estratto = 1 - {len(common)}/{min(len(kh), len(ke))} = {kmer_distance(h12, e12):.4f}")

print_section("2. MATRICE DI DISTANZA K-MER (sequenze complete, k=3)")
names = list(SEQS.keys())
D = {}
for a, b in combinations(names, 2):
    d = kmer_distance(SEQS[a], SEQS[b])
    D[(a, b)] = d
    print(f"D({a}, {b}) = {d:.4f}")


def get_dist(x, y):
    return D.get((x, y), D.get((y, x)))


print_section("3. CLUSTERING UPGMA")
closest = min(D, key=D.get)
print(f"Coppia piu' vicina: {closest}  distanza={D[closest]:.4f}")
h1 = D[closest] / 2
print(f"Altezza fusione 1 = {D[closest]:.4f}/2 = {h1:.4f}")
third = [n for n in names if n not in closest][0]
avg = (get_dist(closest[0], third) + get_dist(closest[1], third)) / 2
print(f"Distanza media cluster-{third} = ({get_dist(closest[0], third):.4f} + {get_dist(closest[1], third):.4f})/2 = {avg:.4f}")
h_root = avg / 2
print(f"Altezza radice = {avg:.4f}/2 = {h_root:.4f}")
print(f"Newick: (({closest[0]}:{h1:.4f}, {closest[1]}:{h1:.4f}):{h_root - h1:.4f}, {third}:{h_root:.4f});")


print_section("4. ALLINEAMENTO PROGRESSIVO (Needleman-Wunsch, BLOSUM62)")
aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.mode = "global"

a_name, b_name = closest
aln1 = aligner.align(SEQS[a_name], SEQS[b_name])[0]
human_aligned, horse_aligned = str(aln1[0]), str(aln1[1])
print(f"\nStep 1.3: {a_name} x {b_name}  (score={aln1.score:.1f})")
print(a_name.ljust(9), human_aligned)
print(b_name.ljust(9), horse_aligned)

# mappa esplicita colonna-step1 <-> posizione senza-gap del rappresentante (b_name),
# con gestione delle colonne "orfane" (dove b_name ha gia' un gap allo step 1)
col_map, orphans_before, pending = [], {}, []
for i, ch in enumerate(horse_aligned):
    if ch == "-":
        pending.append(i)
    else:
        if pending:
            orphans_before[len(col_map)] = pending
            pending = []
        col_map.append(i)
trailing_orphans = pending
rep_nogap = "".join(horse_aligned[i] for i in col_map)
assert rep_nogap == SEQS[b_name]

aln2 = aligner.align(rep_nogap, SEQS[third])[0]
rep_realigned, third_aligned = str(aln2[0]), str(aln2[1])
print(f"\nStep 2.3: profilo({a_name},{b_name}) x {third}  (score={aln2.score:.1f})")
print("profilo".ljust(9), rep_realigned)
print(third.ljust(9), third_aligned)

final = {a_name: [], b_name: [], third: []}


def emit_orphans(cols):
    for i in cols:
        final[a_name].append(human_aligned[i])
        final[b_name].append("-")
        final[third].append("-")


j = 0
third_slot = []  # segnaposto per i caratteri di 'third', riempiti alla fine
for col in rep_realigned:
    if col == "-":
        final[a_name].append("-")
        final[b_name].append("-")
        third_slot.append(None)
    else:
        if j in orphans_before:
            emit_orphans(orphans_before[j])
            third_slot.extend([None] * 0)  # gli orfani hanno gia' '-' fissato sopra
        step1_col = col_map[j]
        final[a_name].append(human_aligned[step1_col])
        final[b_name].append(horse_aligned[step1_col])
        assert horse_aligned[step1_col] == col
        third_slot.append(None)
        j += 1
if j in orphans_before:
    emit_orphans(orphans_before[j])
emit_orphans(trailing_orphans)

# ricostruisco correttamente la riga di 'third' includendo i gap gia' fissati dagli orfani
final[third] = []
ti = 0
tchars = list(third_aligned)
ptr_orphan_positions = set()
# ricostruzione piu' semplice e robusta: rifaccio il merge in un unico passaggio lineare
final = {a_name: [], b_name: [], third: []}
j = 0
tchars_iter = iter(third_aligned)
for col in rep_realigned:
    if col != "-" and j in orphans_before:
        for i in orphans_before[j]:
            final[a_name].append(human_aligned[i])
            final[b_name].append("-")
            final[third].append("-")
    if col == "-":
        final[a_name].append("-")
        final[b_name].append("-")
        final[third].append(next(tchars_iter))
    else:
        step1_col = col_map[j]
        final[a_name].append(human_aligned[step1_col])
        final[b_name].append(horse_aligned[step1_col])
        final[third].append(next(tchars_iter))
        j += 1
if j in orphans_before:
    for i in orphans_before[j]:
        final[a_name].append(human_aligned[i])
        final[b_name].append("-")
        final[third].append("-")
for i in trailing_orphans:
    final[a_name].append(human_aligned[i])
    final[b_name].append("-")
    final[third].append("-")

L = len(final[a_name])
assert L == len(final[b_name]) == len(final[third])
print(f"\n=== MSA1 finale (lunghezza {L}) ===")
for name in [a_name, b_name, third]:
    print(f"{name.ljust(9)}{''.join(final[name])}")

print_section("5. IDENTITA' PERCENTUALE DI CONTROLLO (allineamento globale a coppie)")
for a, b in combinations(names, 2):
    aln = aligner.align(SEQS[a], SEQS[b])[0]
    a1, a2 = aln[0], aln[1]
    matches = sum(1 for x, y in zip(a1, a2) if x == y and x != "-")
    print(f"{a}-{b}: {matches}/{len(a1)} = {matches/len(a1)*100:.1f}% identita'")

print_section("6. COLONNE USATE NELL'ESEMPIO LOG-EXPECTATION (indici 1-based)")
for idx in (7, 17):
    col = tuple(final[name][idx - 1] for name in [a_name, b_name, third])
    print(f"Colonna {idx}: {a_name}={col[0]}  {b_name}={col[1]}  {third}={col[2]}")