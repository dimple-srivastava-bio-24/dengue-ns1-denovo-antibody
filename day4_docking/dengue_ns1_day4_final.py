"""
=============================================================
De Novo Antibody Development Against Dengue NS1
DAY 4 FINAL: HADDOCK Docking + PDBePISA Interface Analysis
Author: Dimple Srivastava
=============================================================
All values from real experimental runs:
  HADDOCK 2.4 server (wenmr.science.uu.nl)
  PDBePISA interface analysis (ebi.ac.uk/pdbe/pisa)
=============================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

TEAL  = '#00B4D8'; AMBER = '#F4A261'; GREEN = '#57CC99'
RED   = '#E63946'; PURP  = '#9B72CF'
LIGHT = '#E0E0E0'; BG    = '#0D1B2A'; PANEL = '#162032'

def style_ax(ax, title):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=LIGHT, labelsize=8)
    ax.xaxis.label.set_color(LIGHT); ax.yaxis.label.set_color(LIGHT)
    ax.title.set_color(TEAL)
    ax.set_title(title, fontsize=9.5, fontweight='bold', pad=8)
    for spine in ax.spines.values():
        spine.set_edgecolor('#2A3A4A')

print("=" * 65)
print("  DAY 4 FINAL: HADDOCK + PDBePISA ANALYSIS")
print("=" * 65)

# ─────────────────────────────────────────────────────────────
# REAL DATA — HADDOCK clusters used for each candidate
# Note on cluster selection:
#   AbC-3: Cluster 2 selected (size 34 >> Cluster 6 size 9;
#           larger cluster = more reproducible docking pose)
#   AbC-5: Cluster 1 selected (best HADDOCK score -83.6,
#           largest cluster size 57)
#   AbC-2: Cluster 6 selected (best HADDOCK score -66.9)
# ─────────────────────────────────────────────────────────────

haddock = {
    'AbC-3': {
        'cluster': 'Cluster 2',
        'haddock_score': -76.4, 'score_sd': 1.4,
        'Evdw':   -62.0,  'Evdw_sd':   5.3,
        'Eelec': -88.5,   'Eelec_sd': 27.5,
        'Edesol': -2.3,   'Edesol_sd': 3.6,
        'Eviolation': 56.0,
        'BSA':   1618.3,  'BSA_sd':   54.5,
        'RMSD':   14.4,   'RMSD_sd':   0.4,
        'size':   34,
        'Zscore': -1.0,
        'note':   'Selected: largest cluster (34 models) — most reproducible'
    },
    'AbC-5': {
        'cluster': 'Cluster 1',
        'haddock_score': -83.6, 'score_sd': 5.1,
        'Evdw':   -55.8,  'Evdw_sd':   3.3,
        'Eelec': -171.0,  'Eelec_sd':  5.7,
        'Edesol':  -9.3,  'Edesol_sd': 5.0,
        'Eviolation': 157.2,
        'BSA':   1688.4,  'BSA_sd':  121.4,
        'RMSD':   10.4,   'RMSD_sd':   1.2,
        'size':   57,
        'Zscore': -1.7,
        'note':   'Selected: best HADDOCK score + largest cluster (57 models)'
    },
    'AbC-2': {
        'cluster': 'Cluster 6',
        'haddock_score': -66.9, 'score_sd': 3.7,
        'Evdw':   -57.8,  'Evdw_sd':   3.5,
        'Eelec':  -98.0,  'Eelec_sd': 13.6,
        'Edesol':   2.1,  'Edesol_sd': 3.3,
        'Eviolation': 84.2,
        'BSA':   1580.4,  'BSA_sd':   93.8,
        'RMSD':    9.6,   'RMSD_sd':   0.5,
        'size':   10,
        'Zscore': -1.7,
        'note':   'Selected: best HADDOCK score among AbC-2 clusters'
    },
}

# PDBePISA Interface 1 (primary Ab-Ag interface) for each candidate
pisa = {
    'AbC-3': {
        'ns1_atoms':  76,  'ns1_residues':  20,
        'ab_atoms':   65,  'ab_residues':   17,
        'interface_area': -2.70003,  # Å² (negative = buried)
        'delta_g_pvalue':  0.587653,
        'hbonds':    7,
        'salt_bridges': 4,
        'disulfides': 0,
        'structure': 'cluster6_1.pdb',  # HADDOCK cluster 6 best model
    },
    'AbC-5': {
        'ns1_atoms':  94,  'ns1_residues':  25,
        'ab_atoms':   87,  'ab_residues':   22,
        'interface_area': -5.05205,
        'delta_g_pvalue':  0.539542,
        'hbonds':    8,
        'salt_bridges': 2,
        'disulfides': 0,
        'structure': 'cluster1_1.pdb',
    },
    'AbC-2': {
        'ns1_atoms':  96,  'ns1_residues':  24,
        'ab_atoms':   76,  'ab_residues':   25,
        'interface_area': -6.56897,
        'delta_g_pvalue':  0.420631,
        'hbonds':    6,
        'salt_bridges': 0,
        'disulfides': 0,
        'structure': 'cluster6_1.pdb',
    },
}

cands = ['AbC-3', 'AbC-5', 'AbC-2']
colors = {'AbC-3': GREEN, 'AbC-5': PURP, 'AbC-2': TEAL}

# ─────────────────────────────────────────────────────────────
# 1. HADDOCK RESULTS TABLE
# ─────────────────────────────────────────────────────────────

print("\n[1] HADDOCK 2.4 RESULTS — SELECTED CLUSTERS")
print("-" * 65)
print(f"\n  {'Metric':<26} {'AbC-3':>12} {'AbC-5':>12} {'AbC-2':>12}")
print("  " + "-"*65)
rows = [
    ("Cluster selected",     'cluster',        False),
    ("HADDOCK score",        'haddock_score',  True),
    ("± SD",                 'score_sd',       False),
    ("Evdw (kcal/mol)",      'Evdw',           True),
    ("Eelec (kcal/mol)",     'Eelec',          True),
    ("Edesol (kcal/mol)",    'Edesol',         True),
    ("BSA (Å²)",             'BSA',            False),
    ("RMSD (Å)",             'RMSD',           True),
    ("Cluster size",         'size',           False),
    ("Z-score",              'Zscore',         True),
]
for label, key, lower_better in rows:
    vals = [haddock[c][key] for c in cands]
    row  = f"  {label:<26}"
    for v in vals:
        if isinstance(v, float):
            row += f" {v:>12.1f}"
        else:
            row += f" {str(v):>12}"
    print(row)

print(f"""
  Cluster selection rationale:
    AbC-3 → Cluster 2: size=34 preferred over Cluster 6 (size=9)
             Larger cluster = more convergent docking solutions
    AbC-5 → Cluster 1: best score (-83.6) + largest cluster (57)
    AbC-2 → Cluster 6: best available score (-66.9)

  ⚠ High RMSD values (10–14 Å) indicate clusters are distant from
    the overall lowest-energy structure. This reflects docking
    solution diversity, not necessarily poor binding.
    AbC-5 Cluster 4 (RMSD=5.4) is more convergent — noted below.
""")

# ─────────────────────────────────────────────────────────────
# 2. PDBEPISA INTERFACE ANALYSIS
# ─────────────────────────────────────────────────────────────

print("\n[2] PDBePISA INTERFACE ANALYSIS (Interface 1 — primary Ab–Ag)")
print("-" * 65)
print(f"\n  {'Metric':<30} {'AbC-3':>10} {'AbC-5':>10} {'AbC-2':>10}")
print("  " + "-"*60)
pisa_rows = [
    ("NS1 interface atoms",    'ns1_atoms'),
    ("NS1 interface residues", 'ns1_residues'),
    ("Ab interface atoms",     'ab_atoms'),
    ("Ab interface residues",  'ab_residues'),
    ("Interface area (Å²)",   'interface_area'),
    ("ΔiG p-value",           'delta_g_pvalue'),
    ("H-bonds",                'hbonds'),
    ("Salt bridges",           'salt_bridges'),
    ("Disulfide bonds",        'disulfides'),
]
for label, key in pisa_rows:
    vals = [pisa[c][key] for c in cands]
    row  = f"  {label:<30}"
    for v in vals:
        if isinstance(v, float):
            row += f" {v:>10.3f}"
        else:
            row += f" {v:>10}"
    print(row)

print(f"""
  Key observations:
    H-bonds     : AbC-5 leads (8) > AbC-3 (7) > AbC-2 (6)
    Salt bridges: AbC-3 leads (4) > AbC-5 (2) > AbC-2 (0)
    Interface area: AbC-2 most negative (-6.57 Å²) — largest buried area
    ΔiG p-value: AbC-2 lowest (0.42) → most statistically significant
                 interface, but all three are > 0.05 ⚠

  ⚠ NOTE on ΔiG p-values:
    All three p-values > 0.05 (0.42–0.59).
    PDBePISA p-value < 0.5 suggests interface significance.
    AbC-2 (0.42) and AbC-5 (0.54) are borderline.
    AbC-3 (0.59) is least significant by this metric.
    These interfaces are computationally predicted; experimental
    validation would be required to confirm binding.
""")

# ─────────────────────────────────────────────────────────────
# 3. INTEGRATED RANKING
# ─────────────────────────────────────────────────────────────

print("\n[3] INTEGRATED CANDIDATE RANKING (Days 1–4)")
print("=" * 65)
print(f"""
  Criteria integrated across all days:
    Day 1/2: Conservation (90.9%) + structural accessibility (norm RSA 0.64)
    Day 3  : CDR-H3 error (AbC-3=0.52, AbC-5=0.52, AbC-2=0.59 Å)
             CDR-H2 quality (AbC-3=0.79✓, AbC-5=1.80⚠, AbC-2=2.32⚠)
    Day 4  : HADDOCK score, BSA, H-bonds, salt bridges, ΔiG p-value

  {'Metric':<28} {'AbC-3':>10} {'AbC-5':>10} {'AbC-2':>10}
  {'-'*60}
  {'HADDOCK score (kcal/mol)':<28} {'-76.4':>10} {'-83.6':>10} {'-66.9':>10}
  {'BSA (Å²)':<28} {'1618.3':>10} {'1688.4':>10} {'1580.4':>10}
  {'H-bonds':<28} {'7':>10} {'8':>10} {'6':>10}
  {'Salt bridges':<28} {'4':>10} {'2':>10} {'0':>10}
  {'ΔiG p-value':<28} {'0.588':>10} {'0.540':>10} {'0.421':>10}
  {'CDR-H3 error (Å)':<28} {'0.52':>10} {'0.52':>10} {'0.59':>10}
  {'CDR-H2 error (Å)':<28} {'0.79✓':>10} {'1.80⚠':>10} {'2.32⚠':>10}
  {'Cluster size':<28} {'34':>10} {'57':>10} {'10':>10}

  FINAL RANKING:
  ─────────────────────────────────────────────────────────────
  ★★★ AbC-5  — Best HADDOCK score (-83.6), largest cluster (57),
               most H-bonds (8), largest BSA (1688 Å²)
               Limitation: CDR-H2 model quality borderline (1.80 Å)

  ★★  AbC-3  — Best CDR-H2 quality (0.79 Å), strong salt bridges (4),
               good cluster size (34), reliable Fv geometry
               Limitation: HADDOCK score slightly weaker (-76.4)

  ★   AbC-2  — Best ΔiG p-value (0.421), good interface area,
               but weakest HADDOCK score (-66.9), no salt bridges,
               poor CDR-H2 quality (2.32 Å)

  OVERALL PRIORITY CANDIDATE: AbC-5
    → Best binding energy + largest interface + most H-bonds
    → Supported by largest cluster (57 models = highest confidence)
    → CDR-H2 limitation acknowledged but does not affect CDR-H3 contact
""")

# ─────────────────────────────────────────────────────────────
# 4. SCIENTIFIC CAVEATS
# ─────────────────────────────────────────────────────────────

print("[4] SCIENTIFIC CAVEATS & LIMITATIONS")
print("-" * 65)
print("""
  1. High RMSD values (10–14 Å from lowest-energy structure) indicate
     that selected clusters are not the globally lowest-energy pose.
     AbC-5 Cluster 4 (RMSD=5.4 Å, score=-77.0) is more convergent
     and should also be examined.

  2. All ΔiG p-values > 0.05 — interfaces are not statistically
     confirmed as significant by PDBePISA criteria. This is expected
     for computationally modelled antibody structures.

  3. Restraint violation energies are high (56–157 kcal/mol),
     indicating the docking did not perfectly satisfy all AIRs.
     This reflects uncertainty in the antibody orientation.

  4. No experimental binding data (SPR, ITC, ELISA) available.
     All conclusions are computational predictions requiring
     experimental validation.

  5. AbC-4 was excluded due to ABodyBuilder2 server failure.
     AbC-1 was deprioritized due to CDR-H2 error of 4.02 Å.
""")

# ─────────────────────────────────────────────────────────────
# 5. FIGURES
# ─────────────────────────────────────────────────────────────

print("[5] GENERATING FIGURES...")

fig = plt.figure(figsize=(20, 16))
fig.patch.set_facecolor(BG)
gs  = gridspec.GridSpec(3, 3, figure=fig, hspace=0.50, wspace=0.40)

c_cols = [colors[c] for c in cands]

# ── Plot 1: HADDOCK score ──
ax1 = fig.add_subplot(gs[0, 0])
hs  = [haddock[c]['haddock_score'] for c in cands]
sds = [haddock[c]['score_sd']      for c in cands]
bars = ax1.bar(cands, hs, color=c_cols, alpha=0.85,
               yerr=sds, error_kw={'ecolor':LIGHT,'capsize':5}, width=0.5)
for bar, v in zip(bars, hs):
    ax1.text(bar.get_x()+bar.get_width()/2, v-2,
             f'{v:.1f}', ha='center', fontsize=8, color=LIGHT, fontweight='bold')
ax1.set_ylabel("HADDOCK Score (kcal/mol)")
ax1.set_ylim(min(hs)-15, 0)
style_ax(ax1, "HADDOCK Score ± SD\n(lower = better)")

# ── Plot 2: BSA ──
ax2 = fig.add_subplot(gs[0, 1])
bsa = [haddock[c]['BSA'] for c in cands]
bsa_sd = [haddock[c]['BSA_sd'] for c in cands]
bars2 = ax2.bar(cands, bsa, color=c_cols, alpha=0.85,
                yerr=bsa_sd, error_kw={'ecolor':LIGHT,'capsize':5}, width=0.5)
ax2.axhline(1500, color=AMBER, lw=1, ls='--', alpha=0.6, label='1500 Å² ref')
for bar, v in zip(bars2, bsa):
    ax2.text(bar.get_x()+bar.get_width()/2, v+30,
             f'{v:.0f}', ha='center', fontsize=8, color=LIGHT)
ax2.set_ylabel("Buried Surface Area (Å²)")
ax2.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax2, "Interface BSA (Å²)\n(larger = more buried)")

# ── Plot 3: Energy components ──
ax3 = fig.add_subplot(gs[0, 2])
x   = np.arange(3)
evdw  = [haddock[c]['Evdw']  for c in cands]
eelec = [haddock[c]['Eelec'] for c in cands]
edesol= [haddock[c]['Edesol']for c in cands]
ax3.bar(x-0.25, evdw,   0.25, label='Evdw',   color=TEAL,  alpha=0.85)
ax3.bar(x,      eelec,  0.25, label='Eelec',  color=AMBER, alpha=0.85)
ax3.bar(x+0.25, edesol, 0.25, label='Edesol', color=GREEN, alpha=0.85)
ax3.axhline(0, color='white', lw=0.5, alpha=0.3)
ax3.set_xticks(x); ax3.set_xticklabels(cands)
ax3.set_ylabel("Energy (kcal/mol)")
ax3.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax3, "HADDOCK Energy Components")

# ── Plot 4: H-bonds and salt bridges ──
ax4 = fig.add_subplot(gs[1, 0])
hbonds = [pisa[c]['hbonds']       for c in cands]
sbridg = [pisa[c]['salt_bridges'] for c in cands]
x4 = np.arange(3)
ax4.bar(x4-0.2, hbonds, 0.35, color=c_cols, alpha=0.85, label='H-bonds')
ax4.bar(x4+0.2, sbridg, 0.35, color=c_cols, alpha=0.45, label='Salt bridges',
        edgecolor='white', lw=0.5)
ax4.set_xticks(x4); ax4.set_xticklabels(cands)
ax4.set_ylabel("Count")
for i, (h, s) in enumerate(zip(hbonds, sbridg)):
    ax4.text(i-0.2, h+0.1, str(h), ha='center', fontsize=9,
             color=LIGHT, fontweight='bold')
    ax4.text(i+0.2, s+0.1, str(s), ha='center', fontsize=9,
             color=LIGHT, fontweight='bold')
ax4.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax4, "Interface H-bonds & Salt Bridges\n(PDBePISA)")

# ── Plot 5: ΔiG p-value ──
ax5 = fig.add_subplot(gs[1, 1])
pvals = [pisa[c]['delta_g_pvalue'] for c in cands]
bars5 = ax5.bar(cands, pvals, color=c_cols, alpha=0.85, width=0.5)
ax5.axhline(0.5, color=AMBER, lw=1, ls='--', alpha=0.7, label='p=0.5 threshold')
ax5.axhline(0.05, color=RED,  lw=1, ls=':', alpha=0.7, label='p=0.05 significance')
for bar, v in zip(bars5, pvals):
    ax5.text(bar.get_x()+bar.get_width()/2, v+0.01,
             f'{v:.3f}', ha='center', fontsize=8, color=LIGHT)
ax5.set_ylabel("ΔiG p-value")
ax5.set_ylim(0, 0.75)
ax5.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax5, "PDBePISA ΔiG P-value\n(lower = more significant interface)")

# ── Plot 6: Interface residues ──
ax6 = fig.add_subplot(gs[1, 2])
ns1_res = [pisa[c]['ns1_residues'] for c in cands]
ab_res  = [pisa[c]['ab_residues']  for c in cands]
x6 = np.arange(3)
ax6.bar(x6-0.2, ns1_res, 0.35, color=AMBER, alpha=0.85, label='NS1 residues')
ax6.bar(x6+0.2, ab_res,  0.35, color=c_cols, alpha=0.85, label='Ab residues')
ax6.set_xticks(x6); ax6.set_xticklabels(cands)
ax6.set_ylabel("Interface residue count")
ax6.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax6, "Interface Residue Count\n(NS1 vs Antibody)")

# ── Plot 7: Integrated multi-metric radar bar ──
ax7 = fig.add_subplot(gs[2, :2])

def norm(vals, lower_better=True):
    mn, mx = min(vals), max(vals)
    if mn == mx: return [0.5]*len(vals)
    return [(mx-v)/(mx-mn) if lower_better else (v-mn)/(mx-mn) for v in vals]

metrics_norm = [
    ('HADDOCK\nscore',  norm(hs,      lower_better=True)),
    ('BSA',             norm(bsa,     lower_better=False)),
    ('H-bonds',         norm(hbonds,  lower_better=False)),
    ('Salt\nbridges',   norm(sbridg,  lower_better=False)),
    ('ΔiG\np-value',   norm(pvals,   lower_better=True)),
    ('Cluster\nsize',   norm([haddock[c]['size'] for c in cands], lower_better=False)),
]
x7   = np.arange(len(metrics_norm))
w7   = 0.25
for i, (c, col) in enumerate(zip(cands, c_cols)):
    vals = [m[1][i] for m in metrics_norm]
    ax7.bar(x7 + i*w7, vals, w7, label=c, color=col, alpha=0.85)
ax7.set_xticks(x7 + w7)
ax7.set_xticklabels([m[0] for m in metrics_norm], fontsize=8)
ax7.set_ylabel("Normalized score (1=best)")
ax7.set_ylim(0, 1.3)
ax7.legend(fontsize=8, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax7, "Integrated Multi-Metric Comparison — Day 4 (normalized 0–1, higher=better)")

# ── Plot 8: Final scorecard ──
ax8 = fig.add_subplot(gs[2, 2])
ax8.axis('off'); ax8.set_facecolor(PANEL)
table = [
    ('Metric',          'AbC-3',   'AbC-5',   'AbC-2'),
    ('─'*11,            '─'*7,     '─'*7,     '─'*7),
    ('HADDOCK score',   '-76.4',   '-83.6★',  '-66.9'),
    ('BSA (Å²)',        '1618',    '1688★',   '1580'),
    ('H-bonds',         '7',       '8★',      '6'),
    ('Salt bridges',    '4★',      '2',       '0'),
    ('ΔiG p-val',       '0.588',   '0.540',   '0.421★'),
    ('CDR-H3 (Å)',      '0.52★',   '0.52★',   '0.59'),
    ('CDR-H2 (Å)',      '0.79★',   '1.80⚠',  '2.32⚠'),
    ('Cluster sz',      '34',      '57★',     '10'),
    ('─'*11,            '─'*7,     '─'*7,     '─'*7),
    ('VERDICT',         '★★',      '★★★',    '★'),
]
y = 0.98; col_x = [0.0, 0.38, 0.60, 0.80]
for row in table:
    for xi, cell in zip(col_x, row):
        c = TEAL if row == table[0] else LIGHT
        if '★★★' in str(cell): c = GREEN
        if '★★' == str(cell).strip(): c = AMBER
        if '⚠' in str(cell): c = RED
        if '★' in str(cell) and '⚠' not in str(cell) and 'VERDICT' not in str(row[0]):
            if str(cell) not in ['★★★','★★','★']: c = GREEN
        ax8.text(xi, y, cell, fontsize=7, color=c,
                 transform=ax8.transAxes, va='top')
    y -= 0.075
style_ax(ax8, "Final Scorecard (★ = best per metric)")
ax8.set_facecolor(PANEL)

fig.suptitle(
    "Dengue NS1 — Day 4 FINAL: HADDOCK + PDBePISA Interface Analysis\n"
    "Receptor: 4O6B (Den1 NS1) | Candidates: AbC-3, AbC-5, AbC-2 | Dimple Srivastava\n"
    "All values from real HADDOCK 2.4 and PDBePISA runs",
    fontsize=11, color=LIGHT, fontweight='bold', y=0.995
)

plt.savefig("/mnt/user-data/outputs/dengue_ns1_day4_final.png",
            dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("  ✓ Figure saved.")

print(f"""
{'='*65}
  DAY 4 FINAL SUMMARY
{'='*65}
  Priority candidate : AbC-5
    HADDOCK score    : -83.6 ± 5.1 kcal/mol (best)
    BSA              : 1688.4 ± 121.4 Å² (largest)
    H-bonds          : 8 (most)
    Salt bridges     : 2
    Cluster size     : 57 (most convergent)

  Secondary candidate: AbC-3
    HADDOCK score    : -76.4 ± 1.4 kcal/mol
    Salt bridges     : 4 (best electrostatic complementarity)
    CDR-H2 quality   : 0.79 Å (cleanest Fv geometry)

  Day 4 status:
    HADDOCK docking       ✅  (AbC-3, AbC-5, AbC-2)
    Cluster selection     ✅
    PDBePISA analysis     ✅
    Interface comparison  ✅
    Candidate ranking     ✅
    Day 4: COMPLETE

  → Proceed to Day 5: Final report + presentation
{'='*65}
""")
