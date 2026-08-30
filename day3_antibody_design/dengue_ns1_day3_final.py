"""
=============================================================
De Novo Antibody Development Against Dengue NS1
DAY 3 — FINAL: Antibody Structure Validation Summary
Author: Dimple Srivastava
=============================================================
Sources:
  ABodyBuilder2 modelling scores — Job ID: 20260818_0701560
  PyMOL structural validation — Cα counts, CDR-H SASA
  CDR-H region: residues 95–102 (Chothia numbering)
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
print("  DAY 3 FINAL: ANTIBODY STRUCTURE VALIDATION SUMMARY")
print("=" * 65)

# ─────────────────────────────────────────────────────────────
# DATA
# ─────────────────────────────────────────────────────────────

candidates = ['AbC-3', 'AbC-2', 'AbC-5']
colors = {'AbC-3': GREEN, 'AbC-2': TEAL, 'AbC-5': PURP}

ab2_scores = {
    'AbC-3': {'FW_H':0.75,'CDR_H1':1.71,'CDR_H2':0.79,'CDR_H3':0.52,
              'FW_L':0.51,'CDR_L1':1.90,'CDR_L2':0.33,'CDR_L3':0.67},
    'AbC-2': {'FW_H':0.88,'CDR_H1':1.43,'CDR_H2':2.32,'CDR_H3':0.59,
              'FW_L':0.56,'CDR_L1':0.43,'CDR_L2':0.40,'CDR_L3':0.60},
    'AbC-5': {'FW_H':0.95,'CDR_H1':1.25,'CDR_H2':1.80,'CDR_H3':0.52,
              'FW_L':0.65,'CDR_L1':0.73,'CDR_L2':0.37,'CDR_L3':0.81},
}

pymol = {
    'AbC-3': {'ca_H':114, 'ca_L':98, 'cdr_sasa':882.545,  'hchain_sasa':12779.681},
    'AbC-2': {'ca_H':114, 'ca_L':98, 'cdr_sasa':881.385,  'hchain_sasa':12986.202},
    'AbC-5': {'ca_H':117, 'ca_L':98, 'cdr_sasa':810.845,  'hchain_sasa':10866.554},
}

# Derived values
for c in candidates:
    d = pymol[c]
    d['pct_hchain']  = d['cdr_sasa'] / d['hchain_sasa'] * 100
    d['per_res']     = d['cdr_sasa'] / 8
    d['hchain_mean'] = d['hchain_sasa'] / d['ca_H']
    d['norm_rsa']    = d['per_res'] / d['hchain_mean']

# ─────────────────────────────────────────────────────────────
# 1. PYMOL VALIDATION SUMMARY
# ─────────────────────────────────────────────────────────────

print("\n[1] PyMOL STRUCTURAL VALIDATION")
print("    CDR-H region: residues 95–102 (Chothia numbering, 8 residues)")
print("-" * 65)

print(f"\n  {'Parameter':<35} {'AbC-3':>10} {'AbC-2':>10} {'AbC-5':>10}")
print("  " + "-"*65)

rows = [
    ("H-chain Cα atoms",         "ca_H",        "",     ""),
    ("L-chain Cα atoms",         "ca_L",        "",     ""),
    ("CDR-H SASA (Å²)",          "cdr_sasa",    "",     ""),
    ("H-chain total SASA (Å²)",  "hchain_sasa", "",     ""),
    ("CDR-H % of H-chain",       "pct_hchain",  "%",    ""),
    ("CDR-H per-residue (Å²/res)","per_res",    "",     ""),
    ("H-chain mean per-res",     "hchain_mean", "",     ""),
    ("Normalized RSA (CDR/mean)","norm_rsa",    "",     ""),
]

for label, key, suffix, _ in rows:
    vals = [pymol[c][key] for c in candidates]
    fmt  = ".1f" if key in ('per_res','hchain_mean','hchain_sasa','cdr_sasa') else ".2f"
    if key in ('ca_H','ca_L'): fmt = "d"
    row  = f"  {label:<35}"
    for v in vals:
        if fmt == "d": row += f" {int(v):>10}"
        elif suffix == "%": row += f" {v:>9.2f}%"
        else: row += f" {v:>10{fmt}}"
    print(row)

print(f"""
  Coordinate gaps noted (non-CDR framework regions):
    AbC-3 H: 10, 58–63, 66–69, 72–73  |  L: 29–37, 57–64, 66–68, 73, 81–82
    AbC-2 H: 10, 33, 58–63, 66–72     |  L: 24–26, 29–36, 57–64, 66–68, 73
    AbC-5 H: 10, 33, 58–63, 66–71     |  L: 25–26, 29–36, 39, 57–64, 73

  ⚠ These gaps are in non-CDR framework regions — standard ABodyBuilder2
    modeling omissions. CDR-H3 (95–102) confirmed fully modeled in all 3.
  AbC-5 contains PDB insertion codes (111A, 112A, 112B) — valid numbering.
""")

# ─────────────────────────────────────────────────────────────
# 2. CDR-H SASA INTERPRETATION
# ─────────────────────────────────────────────────────────────

print("[2] CDR-H SASA INTERPRETATION")
print("-" * 65)
print(f"""
  Normalized RSA values (CDR-H per-residue ÷ H-chain per-residue mean):
    AbC-3: 110.3 / 112.1 = 0.98  → CDR-H exposure ≈ whole-chain mean
    AbC-2: 110.2 / 113.9 = 0.97  → CDR-H exposure ≈ whole-chain mean
    AbC-5: 101.4 /  92.9 = 1.09  → CDR-H slightly above chain mean

  All three CDR-H regions are well-exposed (norm RSA ~1.0):
    → The antigen-contact loop is solvent-accessible in all candidates
    → This is the expected and desired result for a functional paratope
    → AbC-5 CDR-H is marginally more exposed than AbC-3/AbC-2

  Note: CDR-H SASA measured on the isolated antibody Fv structure.
  Upon antigen binding, a portion of this surface will become buried
  at the interface — quantified by HADDOCK + PDBePISA in Day 4.
""")

# ─────────────────────────────────────────────────────────────
# 3. COMBINED RANKING
# ─────────────────────────────────────────────────────────────

print("[3] COMBINED DAY 3 RANKING")
print("    Criteria: ABodyBuilder2 CDR-H3 score + CDR-H2 quality + CDR-H SASA")
print("-" * 65)

print(f"""
  {'Candidate':<10} {'CDR-H3 (Å)':>11} {'CDR-H2 (Å)':>11} {'Norm RSA':>10} {'Verdict'}
  {'-'*65}
  AbC-3      0.52 ✓        0.79 ✓       0.98        ★★★ PRIORITY — proceed to docking
  AbC-5      0.52 ✓        1.80 ⚠       1.09        ★★  Strong — proceed to docking
  AbC-2      0.59 ✓        2.32 ⚠       0.97        ★★  Proceed with caution (CDR-H2)
  AbC-1      0.57 ✓        4.02 ✗        —          ★   Deprioritized — poor H2
  AbC-4       —             —            —          ✗   EXCLUDED — server failure

  FORWARD TO DAY 4 DOCKING: AbC-3, AbC-5, AbC-2 (in priority order)
""")

print("[4] DAY 3 CONCLUSION")
print("=" * 65)
print("""
  Three antibody Fv structures were successfully produced by ABodyBuilder2
  and validated in PyMOL. All three (AbC-3, AbC-2, AbC-5) contain both H
  and L chains with CDR-H3 (residues 95–102) fully modeled and
  solvent-accessible (normalized RSA ~1.0).

  Coordinate gaps in framework regions are consistent with standard
  ABodyBuilder2 modeling behavior and do not affect the CDR-H3 loop
  used for antigen contact.

  ABodyBuilder2 modelling scores confirm CDR-H3 prediction error ≤ 0.59 Å
  for all three candidates — high confidence loop geometry.

  AbC-3 is ranked as the priority candidate:
    • Best CDR-H3 error (0.52 Å)
    • Only candidate with CDR-H2 < 1.0 Å (0.79 Å) — clean Fv geometry
    • CDR-H normalized RSA = 0.98 — well-exposed paratope

  Day 3 status:
    Antibody candidate design       ✅
    ABodyBuilder2 structure pred.   ✅  (4/5 successful)
    PyMOL chain verification        ✅
    Cα atom count verification      ✅
    Coordinate gap documentation    ✅
    CDR-H3 coordinate confirmation  ✅
    CDR-H SASA quantification       ✅
    Candidate ranking               ✅
    Day 3 structural validation: COMPLETE

  → Proceed to Day 4: HADDOCK docking of AbC-3/AbC-5/AbC-2 vs NS1 (4O6B)
""")

# ─────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────

print("[5] GENERATING FIGURES...")

fig = plt.figure(figsize=(20, 14))
fig.patch.set_facecolor(BG)
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.50, wspace=0.40)

# ── Plot 1: ABodyBuilder2 heatmap ──
ax1 = fig.add_subplot(gs[0, :2])
regions = ['FW_H','CDR_H1','CDR_H2','CDR_H3','FW_L','CDR_L1','CDR_L2','CDR_L3']
matrix  = np.array([[ab2_scores[c][r] for r in regions] for c in candidates])
im = ax1.imshow(matrix, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=4.5)
ax1.set_xticks(range(len(regions)))
ax1.set_xticklabels([r.replace('_','-') for r in regions], fontsize=8, color=LIGHT)
ax1.set_yticks(range(len(candidates)))
ax1.set_yticklabels(candidates, fontsize=9, color=LIGHT)
for i in range(len(candidates)):
    for j in range(len(regions)):
        val = matrix[i,j]
        col = 'black' if val < 2.5 else 'white'
        ax1.text(j, i, f'{val:.2f}', ha='center', va='center',
                 fontsize=8.5, color=col, fontweight='bold')
cbar = plt.colorbar(im, ax=ax1, pad=0.02)
cbar.set_label('Prediction Error (Å)', color=LIGHT, fontsize=8)
cbar.ax.tick_params(colors=LIGHT)
style_ax(ax1, "ABodyBuilder2 Modelling Scores — Lower = Better (green) | ★ = CDR-H3")

# ── Plot 2: CDR-H3 vs CDR-H2 scatter ──
ax2 = fig.add_subplot(gs[0, 2])
for c in candidates:
    h3 = ab2_scores[c]['CDR_H3']
    h2 = ab2_scores[c]['CDR_H2']
    ax2.scatter(h2, h3, s=180, color=colors[c], zorder=5,
                edgecolors='white', lw=0.8)
    ax2.annotate(c, (h2, h3), textcoords="offset points",
                 xytext=(6, 4), fontsize=8, color=colors[c])
ax2.axvline(2.0, color=RED,   lw=1, ls='--', alpha=0.6, label='Poor threshold (2.0 Å)')
ax2.axhline(1.0, color=AMBER, lw=1, ls='--', alpha=0.6, label='Good threshold (1.0 Å)')
ax2.set_xlabel("CDR-H2 Prediction Error (Å)")
ax2.set_ylabel("CDR-H3 Prediction Error (Å)")
ax2.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
# Best quadrant label
ax2.text(0.3, 0.4, "Best\nquadrant", fontsize=7, color=GREEN,
         alpha=0.5, transform=ax2.transAxes)
style_ax(ax2, "CDR-H3 vs CDR-H2 Error\n(bottom-left = best)")

# ── Plot 3: CDR-H SASA comparison ──
ax3 = fig.add_subplot(gs[1, 0])
sasa_vals = [pymol[c]['cdr_sasa'] for c in candidates]
bars = ax3.bar(candidates, sasa_vals,
               color=[colors[c] for c in candidates], alpha=0.85, width=0.5)
for bar, v in zip(bars, sasa_vals):
    ax3.text(bar.get_x()+bar.get_width()/2, v+8,
             f'{v:.1f}', ha='center', fontsize=8, color=LIGHT, fontweight='bold')
ax3.set_ylabel("CDR-H (95–102) SASA (Å²)")
ax3.set_ylim(700, 950)
style_ax(ax3, "CDR-H Absolute SASA (Å²)")

# ── Plot 4: Normalized RSA ──
ax4 = fig.add_subplot(gs[1, 1])
norm_vals = [pymol[c]['norm_rsa'] for c in candidates]
bars4 = ax4.bar(candidates, norm_vals,
                color=[colors[c] for c in candidates], alpha=0.85, width=0.5)
ax4.axhline(1.0, color='white', lw=1, ls='--', alpha=0.5,
            label='H-chain mean (1.0)')
for bar, v in zip(bars4, norm_vals):
    ax4.text(bar.get_x()+bar.get_width()/2, v+0.01,
             f'{v:.2f}', ha='center', fontsize=9, color=LIGHT, fontweight='bold')
ax4.set_ylabel("Normalized RSA (CDR ÷ H-chain mean)")
ax4.set_ylim(0, 1.4)
ax4.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax4, "CDR-H Normalized RSA\n(~1.0 = well-exposed paratope)")

# ── Plot 5: Final scorecard ──
ax5 = fig.add_subplot(gs[1, 2])
ax5.axis('off')
ax5.set_facecolor(PANEL)

table_data = [
    ('Metric',          'AbC-3',    'AbC-2',    'AbC-5'),
    ('─'*12,            '─'*8,      '─'*8,      '─'*8),
    ('H-chain Cα',      '114',      '114',      '117'),
    ('L-chain Cα',      '98',       '98',       '98'),
    ('CDR-H3 err (Å)',  '0.52 ✓',  '0.59 ✓',  '0.52 ✓'),
    ('CDR-H2 err (Å)',  '0.79 ✓',  '2.32 ⚠',  '1.80 ⚠'),
    ('CDR SASA (Å²)',   '882.5',    '881.4',    '810.8'),
    ('Norm. RSA',       '0.98',     '0.97',     '1.09'),
    ('Verdict',         '★★★',     '★★ ⚠H2',  '★★'),
]
y = 0.97; col_x = [0.01, 0.35, 0.58, 0.78]
for row in table_data:
    for xi, cell in zip(col_x, row):
        c = TEAL if row == table_data[0] else LIGHT
        if '★★★' in str(cell): c = GREEN
        if '⚠' in str(cell):   c = AMBER
        if '✓' in str(cell) and '★' not in str(cell): c = GREEN
        ax5.text(xi, y, cell, fontsize=7.5, color=c,
                 transform=ax5.transAxes, va='top')
    y -= 0.10
style_ax(ax5, "Day 3 Validation Scorecard")

fig.suptitle(
    "Dengue NS1 — Day 3 FINAL: Antibody Structure Validation\n"
    "ABodyBuilder2 scores + PyMOL CDR-H SASA | Dimple Srivastava\n"
    "Forward to docking: AbC-3 (priority) → AbC-5 → AbC-2",
    fontsize=11, color=LIGHT, fontweight='bold', y=0.995
)

plt.savefig("/mnt/user-data/outputs/dengue_ns1_day3_final.png",
            dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("  ✓ Figure saved.")
print("=" * 65)
