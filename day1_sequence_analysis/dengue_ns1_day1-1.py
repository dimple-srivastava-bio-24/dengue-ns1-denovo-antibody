"""
=============================================================
De Novo Antibody Development Against Dengue NS1
DAY 1: Sequence QC, Physicochemical Analysis & Conservation
Author: Dimple Srivastava
=============================================================
"""

from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────
# 1. SEQUENCES
# ─────────────────────────────────────────────

DEN1_AA = (
    "DSGCVINWKGRELKCGSGIFVTNEVHTWTEQYKFQADSPKRLSAAIGKAWEEGVCGIRSATR"
    "LENIMWKQISNELNHILLENDMKFTVVVGDVSGILAQGKKMIRPQPMEHKYSWKSWGKAKII"
    "GADVQNTTFIIDGPNTPECPDNQRAWNIWEVEDYGFGIFTTNIWLKLRDSYTQVCDHRLMSA"
    "AIKDSKAVHADMGYWIESEKNETWKLARASFIEVKTCIWPKSHTLWSNGVLESEMIIPKIYGG"
    "PISQHNYRPGYFTQTAGPWHLGKLELDFDLCEGTTVVVDEHCGNRGPSLRTTTVTGKTIHEW"
    "CCRSCTLPPLRFKGEDGCWYGMEIRPVKEKEENLVKSMVSA"
)

DEN2_AA = (
    "DSGCVVSWKNKELKCGSGIFITDNVHTWTEQYKFQPESPSKLASAIQKAHEEGICGIRSV"
    "TRLENLMWKQITPELNHILSENEVKLTIMTGDIKGIMQAGKRSLRPQPTELKYSWKTWGK"
    "AKMLSTESHNQTFLIDGPETAECPNTNRAWNSLEVEDYGFGVFTTNIWLKLKEKQDVFCDS"
    "KLMSAAIKDNRAVHADMGYWIESALNDTWKIEKASFIEVKNCHWPKSHTLWSNGVLESEMII"
    "PKNLAGPVSQHNYRPGYHTQITGPWHLGKLEMDFDFCDGTTVVVTEDCGNRGPSLRTTTASG"
    "KLITEWCCRSCTLPPLRYRGEDGCWYGMEIRPLKEKEENLVNSLVTA"
)

DEN1_DNA = (
    "GATAGTGGTTGCGTGATCAATTGGAAGGGTAGAGAACTGAAGTGTGGCTCCGGCATCTTCGT"
    "GACCAATGAGGTGCACACCTGGACCGAGCAGTACAAGTTCCAGGCCGATTCTCCCAAGCGGC"
    "TGTCTGCCGCCATCGGCAAGGCCTGGGAGGAAGGTGTGTGCGGCATCCGCTCTGCTACCAGA"
    "CTCGAAAATATCATGTGGAAGCAGATCTCCAACGAGCTGAACCACATCCTGCTGGAGAACGA"
    "CATGAAGTTTACCGTGGTCGTTGGAGATGTGTCTGGCATCCTGGCTCAGGGCAAGAAGATGA"
    "TCCGGCCTCAGCCTATGGAACACAAGTACTCCTGGAAAAGCTGGGGCAAAGCTAAGATCATT"
    "GGCGCCGACGTGCAGAACACCACCTTCATCATCGACGGCCCTAACACCCCTGAGTGTCCTGA"
    "CAACCAGCGGGCTTGGAACATCTGGGAAGTGGAAGACTACGGCTTCGGCATTTTTACTACCA"
    "ACATCTGGCTGAAGCTGCGGGACTCCTACACCCAAGTGTGTGATCATAGACTGATGTCCGCC"
    "GCTATCAAGGACTCCAAGGCCGTGCACGCCGACATGGGCTACTGGATCGAGTCTGAGAAGAA"
    "CGAGACATGGAAGCTGGCCAGAGCCAGCTTCATCGAGGTCAAGACCTGCATCTGGCCAAAGT"
    "CTCATACACTCTGGTCCAACGGCGTGCTGGAATCCGAGATGATCATCCCCAAGATCTACGGC"
    "GGCCCCATCTCTCAGCACAACTACCGGCCTGGCTACTTCACCCAGACCGCTGGACCTTGGCA"
    "CCTGGGAAAGCTGGAACTGGACTTCGACCTGTGCGAGGGCACCACCGTGGTGGTGGACGAGC"
    "ACTGCGGCAACAGAGGCCCAAGCCTGAGAACCACCACAGTGACAGGAAAAACCATCCACGAG"
    "TGGTGCTGCAGATCCTGCACCCTGCCTCCTCTGCGGTTCAAGGGCGAGGATGGCTGCTGGT"
    "ATGGCATGGAGATCAGACCTGTGAAAGAAAAAGAGGAAAACCTGGTGAAGTCCATGGTGTCC"
    "GCT"
)

DEN2_DNA = (
    "GATAGTGGTTGCGTGGTGAGTTGGAAGAATAAAGAACTGAAGTGCGGCTCTGGCATCTTCAT"
    "CACCGACAACGTGCATACATGGACCGAGCAGTACAAGTTCCAACCTGAGTCTCCTTCTAAGC"
    "TGGCCTCTGCTATCCAGAAGGCTCACGAGGAAGGCATCTGCGGAATCCGCTCCGTGACCAGA"
    "CTGGAAAACCTGATGTGGAAGCAGATCACACCCGAGCTGAATCACATTCTGTCCGAGAACGA"
    "GGTGAAGCTGACCATCATGACCGGCGACATCAAGGGCATCATGCAGGCCGGCAAGCGGTCTC"
    "TGAGACCTCAGCCCACCGAGCTGAAGTACTCCTGGAAGACCTGGGGCAAGGCCAAGATGCTG"
    "TCTACCGAGAGCCATAACCAGACCTTTCTGATCGACGGCCCTGAAACCGCCGAGTGTCCTAA"
    "CACCAATAGAGCCTGGAACAGCCTGGAAGTGGAAGATTACGGCTTCGGCGTGTTCACCACCA"
    "ACATCTGGCTGAAGCTGAAAGAGAAGCAGGACGTGTTCTGCGACTCCAAACTGATGTCCGCC"
    "GCTATCAAGGACAACCGGGCCGTGCACGCCGACATGGGCTACTGGATCGAGAGCGCTCTGAA"
    "CGACACCTGGAAGATCGAGAAGGCCAGCTTTATCGAGGTCAAGAACTGCCACTGGCCCAAGT"
    "CCCACACCCTGTGGTCCAACGGTGTGCTGGAATCTGAAATGATCATCCCTAAGAACCTGGCT"
    "GGCCCTGTGTCCCAGCACAACTACCGGCCCGGCTACCACACCCAGATCACCGGGCCTTGGCA"
    "CCTGGGCAAACTGGAGATGGACTTCGACTTCTGCGATGGCACCACCGTGGTGGTCACAGAGG"
    "ATTGCGGCAACAGAGGACCTTCCCTGCGGACCACCACAGCTTCTGGAAAGCTGATCACCGAG"
    "TGGTGCTGTAGATCCTGTACACTGCCTCCACTGAGATACAGAGGCGAGGATGGCTGCTGGTA"
    "TGGCATGGAAATCCGGCCACTGAAGGAGAAAGAAGAGAATTTGGTGAACTCCCTGGTGACCG"
    "CT"
)

print("=" * 60)
print("  DENGUE NS1 — DE NOVO ANTIBODY PROJECT | DAY 1")
print("=" * 60)

# ─────────────────────────────────────────────
# 2. SEQUENCE QC
# ─────────────────────────────────────────────

print("\n[1] SEQUENCE QC")
print("-" * 40)

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def qc_aa(name, seq):
    invalid = [c for c in seq if c not in VALID_AA]
    print(f"\n  {name}")
    print(f"    Length      : {len(seq)} aa")
    print(f"    Invalid aa  : {set(invalid) if invalid else 'None'}")
    print(f"    Starts with : {seq[:5]}...")
    print(f"    Ends with   : ...{seq[-5:]}")
    return len(invalid) == 0

def qc_dna(name, seq):
    valid_nt = set("ATGC")
    invalid = [c for c in seq.upper() if c not in valid_nt]
    translated = str(Seq(seq).translate(to_stop=True))
    print(f"\n  {name}")
    print(f"    Length      : {len(seq)} nt  ({len(seq)//3} codons)")
    print(f"    Invalid nt  : {set(invalid) if invalid else 'None'}")
    print(f"    Translated  : {translated[:20]}...")
    return translated

d1_ok = qc_aa("Den1 NS1 (AA)", DEN1_AA)
d2_ok = qc_aa("Den2 NS1 (AA)", DEN2_AA)
d1_trans = qc_dna("Den1 NS1 (DNA)", DEN1_DNA)
d2_trans = qc_dna("Den2 NS1 (DNA)", DEN2_DNA)

# DNA→AA consistency check
print("\n  DNA→AA Consistency Check:")
print(f"    Den1: {'✓ MATCH' if d1_trans[:20] == DEN1_AA[:20] else '✗ MISMATCH — check reading frame'}")
print(f"    Den2: {'✓ MATCH' if d2_trans[:20] == DEN2_AA[:20] else '✗ MISMATCH — check reading frame'}")

# ─────────────────────────────────────────────
# 3. PHYSICOCHEMICAL ANALYSIS
# ─────────────────────────────────────────────

print("\n\n[2] PHYSICOCHEMICAL ANALYSIS")
print("-" * 40)

def physchem(name, seq):
    pa = ProteinAnalysis(seq)
    mw      = pa.molecular_weight()
    pi      = pa.isoelectric_point()
    instab  = pa.instability_index()
    gravy   = pa.gravy()
    aroma   = pa.aromaticity()
    helix, turn, sheet = pa.secondary_structure_fraction()
    aa_comp = pa.amino_acids_percent

    print(f"\n  {name}")
    print(f"    Molecular Weight   : {mw:,.1f} Da  ({mw/1000:.1f} kDa)")
    print(f"    Isoelectric Point  : {pi:.2f}")
    print(f"    Instability Index  : {instab:.2f}  ({'Stable' if instab < 40 else 'Unstable'})")
    print(f"    GRAVY Score        : {gravy:.4f}  ({'Hydrophilic' if gravy < 0 else 'Hydrophobic'})")
    print(f"    Aromaticity        : {aroma:.4f}")
    print(f"    2° Structure Est.  : Helix {helix*100:.1f}%  Turn {turn*100:.1f}%  Sheet {sheet*100:.1f}%")

    return {
        "name": name, "mw": mw, "pi": pi, "instab": instab,
        "gravy": gravy, "aroma": aroma,
        "helix": helix, "turn": turn, "sheet": sheet,
        "aa_comp": aa_comp
    }

d1_prop = physchem("Den1 NS1", DEN1_AA)
d2_prop = physchem("Den2 NS1", DEN2_AA)

# ─────────────────────────────────────────────
# 4. PAIRWISE ALIGNMENT & CONSERVATION
# ─────────────────────────────────────────────

print("\n\n[3] PAIRWISE ALIGNMENT & CONSERVATION")
print("-" * 40)

from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")
alignments = pairwise2.align.globalds(
    DEN1_AA, DEN2_AA,
    blosum62,
    -10, -0.5
)

best = alignments[0]
aln1, aln2 = best.seqA, best.seqB

# Identity & similarity
matches = sum(1 for a, b in zip(aln1, aln2) if a == b and a != '-')
aligned_len = sum(1 for a, b in zip(aln1, aln2) if a != '-' and b != '-')
identity = matches / aligned_len * 100

# Conservation per position
conserved_pos = [i for i, (a, b) in enumerate(zip(aln1, aln2)) if a == b and a != '-']
variable_pos  = [i for i, (a, b) in enumerate(zip(aln1, aln2)) if a != b and a != '-' and b != '-']

print(f"\n  Alignment length   : {len(aln1)} positions")
print(f"  Identical residues : {matches}")
print(f"  % Identity         : {identity:.1f}%")
print(f"  Conserved positions: {len(conserved_pos)}")
print(f"  Variable positions : {len(variable_pos)}")

# Print alignment excerpt (first 120 chars)
print(f"\n  Alignment preview (first 120 aa):")
print(f"  Den1: {aln1[:120]}")
match_str = ''.join(['|' if a==b else ' ' for a,b in zip(aln1[:120], aln2[:120])])
print(f"        {match_str}")
print(f"  Den2: {aln2[:120]}")

# ─────────────────────────────────────────────
# 5. KNOWN DOMAIN BOUNDARIES (NS1 literature)
# ─────────────────────────────────────────────

print("\n\n[4] NS1 DOMAIN IDENTIFICATION")
print("-" * 40)

domains = {
    "β-roll (dimerization)"  : (1,  29,  "Mediates NS1 dimerization; partially buried"),
    "Wing subdomain"         : (30, 180, "Surface-exposed; contains immunodominant loops"),
    "β-ladder connector"     : (181, 223,"Links wing to β-ladder; partially accessible"),
    "β-ladder"               : (224, 352,"Forms hexamer interface; secreted NS1 outer face"),
}

print("\n  Domain                  |  Residues  |  Structural Role")
print("  " + "-"*75)
for dom, (s, e, func) in domains.items():
    print(f"  {dom:<24}|  {s:>3}–{e:<3}   |  {func}")
 
print("\n  Note: Surface accessibility analysis in Day 2 will determine")
print("  which conserved candidate regions fall in accessible domains.")

# ─────────────────────────────────────────────
# 6. CONSERVATION WINDOW ANALYSIS
# ─────────────────────────────────────────────

print("\n\n[5] CONSERVATION WINDOW ANALYSIS (window=10)")
print("-" * 40)

window = 10
conservation_scores = []
positions = []

min_len = min(len(DEN1_AA), len(DEN2_AA))
for i in range(0, min_len - window + 1):
    seg1 = DEN1_AA[i:i+window]
    seg2 = DEN2_AA[i:i+window]
    score = sum(1 for a, b in zip(seg1, seg2) if a == b) / window
    conservation_scores.append(score)
    positions.append(i + window//2)

# Top conserved windows — report correct start:end residue ranges
top_idx = sorted(range(len(conservation_scores)),
                 key=lambda x: conservation_scores[x], reverse=True)[:5]
print("\n  Top 5 most conserved 10-aa windows (candidate conserved regions):")
print(f"  {'Residues (1-based)':>20}  {'Den1 seq':>12}  {'Den2 seq':>12}  {'Score':>6}")
print("  " + "-"*60)
for idx in sorted(top_idx):
    start_i = idx           # 0-based index into sequence
    end_i   = idx + window - 1
    start_r = start_i + 1  # 1-based residue number
    end_r   = end_i + 1
    seg1 = DEN1_AA[start_i:start_i+window]
    seg2 = DEN2_AA[start_i:start_i+window]
    note = "← highly conserved candidate" if conservation_scores[idx] == 1.0 else ""
    print(f"  aa {start_r:>3}–{end_r:<3}  {' '*9}  {seg1:>12}  {seg2:>12}  {conservation_scores[idx]*100:>5.0f}%  {note}")
print("\n  NOTE: conservation alone does not confirm epitope status.")
print("  These regions require structural (surface accessibility) and")
print("  experimental/literature evidence before classification as epitopes.")

# ─────────────────────────────────────────────
# 7. FIGURES
# ─────────────────────────────────────────────

print("\n\n[6] GENERATING FIGURES...")
print("-" * 40)

fig = plt.figure(figsize=(18, 14))
fig.patch.set_facecolor('#0D1B2A')
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

TEAL   = '#00B4D8'
AMBER  = '#F4A261'
GREEN  = '#57CC99'
RED    = '#E63946'
LIGHT  = '#E0E0E0'
BG     = '#0D1B2A'
PANEL  = '#162032'

def style_ax(ax, title):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=LIGHT, labelsize=8)
    ax.xaxis.label.set_color(LIGHT)
    ax.yaxis.label.set_color(LIGHT)
    ax.title.set_color(TEAL)
    ax.set_title(title, fontsize=10, fontweight='bold', pad=8)
    for spine in ax.spines.values():
        spine.set_edgecolor('#2A3A4A')

# ── Plot 1: Conservation score along sequence ──
ax1 = fig.add_subplot(gs[0, :])
ax1.fill_between(positions, conservation_scores, alpha=0.3, color=TEAL)
ax1.plot(positions, conservation_scores, color=TEAL, linewidth=1.2)
ax1.axhline(0.8, color=AMBER, linewidth=1, linestyle='--', label='80% conserved')
ax1.axhline(1.0, color=GREEN, linewidth=1, linestyle='--', label='100% conserved')

# Domain shading
dom_colors = ['#1a3a5c', '#1a4a3a', '#3a1a4a', '#4a3a1a']
for (dom, (s, e, _func)), col in zip(domains.items(), dom_colors):
    ax1.axvspan(s, min(e, min_len-window), alpha=0.15, color=col,
                label=dom if s < min_len-window else "")

ax1.set_xlabel("Residue Position")
ax1.set_ylabel("Conservation Score")
ax1.set_xlim(0, min_len)
ax1.set_ylim(0, 1.05)
ax1.legend(fontsize=7, loc='lower right',
           facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax1, "Den1 vs Den2 NS1 — Conservation Profile (Window = 10 aa)")

# ── Plot 2: AA composition comparison ──
ax2 = fig.add_subplot(gs[1, 0])
aa_list = sorted(VALID_AA)
# amino_acids_percent already gives fractional proportions (sum=1.0), multiply by 100 for %
pa_d1 = ProteinAnalysis(DEN1_AA)
pa_d2 = ProteinAnalysis(DEN2_AA)
d1_comp = [pa_d1.amino_acids_percent.get(aa, 0) for aa in aa_list]
d2_comp = [pa_d2.amino_acids_percent.get(aa, 0) for aa in aa_list]
x = np.arange(len(aa_list))
ax2.bar(x - 0.2, d1_comp, 0.4, label='Den1', color=TEAL, alpha=0.85)
ax2.bar(x + 0.2, d2_comp, 0.4, label='Den2', color=AMBER, alpha=0.85)
ax2.set_xticks(x)
ax2.set_xticklabels(list(aa_list), fontsize=6)
ax2.set_ylabel("Percentage (%)")
ax2.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax2, "Amino Acid Composition")

# ── Plot 3: Physicochemical comparison radar-like bar ──
ax3 = fig.add_subplot(gs[1, 1])
props = ['MW\n(kDa)', 'pI', 'Instab.\nIndex', 'GRAVY\n(×10)', 'Aroma.\n(×100)']
d1_vals = [d1_prop['mw']/1000, d1_prop['pi'], d1_prop['instab'],
           d1_prop['gravy']*10, d1_prop['aroma']*100]
d2_vals = [d2_prop['mw']/1000, d2_prop['pi'], d2_prop['instab'],
           d2_prop['gravy']*10, d2_prop['aroma']*100]
x3 = np.arange(len(props))
ax3.bar(x3 - 0.2, d1_vals, 0.4, label='Den1', color=TEAL, alpha=0.85)
ax3.bar(x3 + 0.2, d2_vals, 0.4, label='Den2', color=AMBER, alpha=0.85)
ax3.set_xticks(x3)
ax3.set_xticklabels(props, fontsize=7)
ax3.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax3, "Physicochemical Properties")

# ── Plot 4: Secondary structure ──
ax4 = fig.add_subplot(gs[1, 2])
cats = ['α-Helix', 'β-Turn', 'β-Sheet']
d1_ss = [d1_prop['helix']*100, d1_prop['turn']*100, d1_prop['sheet']*100]
d2_ss = [d2_prop['helix']*100, d2_prop['turn']*100, d2_prop['sheet']*100]
x4 = np.arange(3)
ax4.bar(x4 - 0.2, d1_ss, 0.4, color=TEAL, alpha=0.85, label='Den1')
ax4.bar(x4 + 0.2, d2_ss, 0.4, color=AMBER, alpha=0.85, label='Den2')
ax4.set_xticks(x4)
ax4.set_xticklabels(cats, fontsize=8)
ax4.set_ylabel("Fraction (%)")
ax4.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax4, "Predicted 2° Structure Fractions")

# ── Plot 5: Heatmap of conservation in domains ──
ax5 = fig.add_subplot(gs[2, :2])
# Per-residue identity array (no gaps)
min_len2 = min(len(DEN1_AA), len(DEN2_AA))
per_res = np.array([1 if DEN1_AA[i] == DEN2_AA[i] else 0 for i in range(min_len2)])
# Reshape for heatmap (pad to multiple of 50)
cols = 50
rows = (min_len2 // cols)
heatmap_data = per_res[:rows*cols].reshape(rows, cols)
sns.heatmap(heatmap_data, ax=ax5, cmap='YlOrRd_r', linewidths=0,
            xticklabels=False, yticklabels=[f"{i*cols+1}" for i in range(rows)],
            cbar_kws={'label': 'Conserved (1) / Variable (0)'})
ax5.set_xlabel("Position within row (50 aa each)")
ax5.set_ylabel("Residue Block Start")
style_ax(ax5, "Per-Residue Conservation Heatmap (Den1 vs Den2)")
ax5.collections[0].colorbar.ax.yaxis.label.set_color(LIGHT)
ax5.collections[0].colorbar.ax.tick_params(colors=LIGHT)

# ── Plot 6: Top conserved windows ──
ax6 = fig.add_subplot(gs[2, 2])
top10_idx = sorted(range(len(conservation_scores)),
                   key=lambda x: conservation_scores[x], reverse=True)[:10]
top10_pos  = [positions[i] for i in top10_idx]
top10_scr  = [conservation_scores[i]*100 for i in top10_idx]
sorted_pairs = sorted(zip(top10_pos, top10_scr))
ax6.barh([f"aa {p}" for p,_ in sorted_pairs],
         [s for _,s in sorted_pairs],
         color=GREEN, alpha=0.85)
ax6.set_xlabel("Conservation (%)")
ax6.axvline(80, color=AMBER, linestyle='--', linewidth=1)
style_ax(ax6, "Top 10 Conserved Windows")

# Title
fig.suptitle(
    "Dengue NS1 — Day 1: Sequence QC, Physicochemical & Conservation Analysis\n"
    "De Novo Antibody Development Project | Dimple Srivastava",
    fontsize=13, color=LIGHT, fontweight='bold', y=0.98
)

plt.savefig("dengue_ns1_day1_analysis.png",
            dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("  ✓ Figure saved.")

# ─────────────────────────────────────────────
# 8. SUMMARY REPORT
# ─────────────────────────────────────────────

print("\n\n" + "=" * 60)
print("  DAY 1 SUMMARY REPORT")
print("=" * 60)

print(f"""
  Sequences processed   : Den1 NS1 ({len(DEN1_AA)} aa) | Den2 NS1 ({len(DEN2_AA)} aa)
  Sequence identity     : {identity:.1f}%
  Conserved positions   : {len(conserved_pos)} / {aligned_len}
  Variable positions    : {len(variable_pos)}

  Den1 NS1:
    MW = {d1_prop['mw']/1000:.1f} kDa | pI = {d1_prop['pi']:.2f} | 
    Instability = {d1_prop['instab']:.1f} ({'Stable' if d1_prop['instab']<40 else 'Unstable'}) | 
    GRAVY = {d1_prop['gravy']:.3f}

  Den2 NS1:
    MW = {d2_prop['mw']/1000:.1f} kDa | pI = {d2_prop['pi']:.2f} | 
    Instability = {d2_prop['instab']:.1f} ({'Stable' if d2_prop['instab']<40 else 'Unstable'}) | 
    GRAVY = {d2_prop['gravy']:.3f}

  ► DAY 1 CONCLUSION:
    Den1NS1 and Den2NS1 are both 352-aa proteins with identical coding
    lengths of 1056 nt. Pairwise comparison shows 72.7% amino-acid
    sequence identity, with 256 conserved and 96 variable positions.
    Several highly conserved sequence regions were identified as candidate
    regions for subsequent structural and epitope-oriented analysis.
    These regions should NOT yet be considered confirmed antibody epitopes
    — surface accessibility and literature evidence are required (Day 2).

  ► NEXT STEP (Day 2):
    1. Fetch PDB structures: 4O6B (Den1 NS1) and 4OIG (Den2 NS1)
    2. Map conserved candidate regions onto 3D structure in PyMOL
    3. Assess surface accessibility (DSSP / FreeSASA / server)
    4. Submit to BepiPred 3.0 + ElliPro for B-cell epitope prediction
    5. Cross-reference with NS1 literature for experimental epitope data
    6. Select final target region based on: conservation + accessibility + evidence
""")

print("  Output: dengue_ns1_day1_analysis.png")
print("=" * 60)
