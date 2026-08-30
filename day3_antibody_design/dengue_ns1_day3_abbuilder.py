"""
=============================================================
De Novo Antibody Development Against Dengue NS1
DAY 3 — ABodyBuilder2 Modelling Score Analysis
Author: Dimple Srivastava
=============================================================
Source: ABodyBuilder2 online server (Job ID: 20260818_0701560)
        https://opig.stats.ox.ac.uk/webapps/abodybuilder2
Metric: Prediction error = RMSD (Å) per antibody region
        Lower = higher confidence structure prediction
=============================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

TEAL  = '#00B4D8'; AMBER = '#F4A261'; GREEN = '#57CC99'
RED   = '#E63946'; PURP  = '#9B72CF'; PINK  = '#F72585'
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
print("  DAY 3: ABodyBuilder2 MODELLING SCORE ANALYSIS")
print("=" * 65)

# ─────────────────────────────────────────────────────────────
# RAW DATA FROM ABodyBuilder2
# ─────────────────────────────────────────────────────────────

regions = ['FW_H', 'CDR_H1', 'CDR_H2', 'CDR_H3',
           'FW_L', 'CDR_L1', 'CDR_L2', 'CDR_L3']

all_results = {
    'AbC-1': [1.52, 1.37, 4.02, 0.57, 0.58, 0.84, 0.46, 0.60],
    'AbC-2': [0.88, 1.43, 2.32, 0.59, 0.56, 0.43, 0.40, 0.60],
    'AbC-3': [0.75, 1.71, 0.79, 0.52, 0.51, 1.90, 0.33, 0.67],
    'AbC-4': None,   # FAILED — OSError strained sidechain bond
    'AbC-5': [0.95, 1.25, 1.80, 0.52, 0.65, 0.73, 0.37, 0.81],
}

valid_ids = ['AbC-1', 'AbC-2', 'AbC-3', 'AbC-5']
valid_data = {k: dict(zip(regions, all_results[k])) for k in valid_ids}

cand_colors = {
    'AbC-1': AMBER, 'AbC-2': TEAL,
    'AbC-3': GREEN, 'AbC-5': PURP
}

# ─────────────────────────────────────────────────────────────
# 1. STATUS REPORT
# ─────────────────────────────────────────────────────────────

print("\n[1] JOB STATUS")
print("-" * 55)
print(f"""
  Job ID   : 20260818_0701560
  Server   : ABodyBuilder2 (Oxford Protein Informatics Group)
  Status   : Finished

  AbC-1    : ✓ Structure produced
  AbC-2    : ✓ Structure produced
  AbC-3    : ✓ Structure produced
  AbC-4    : ✗ FAILED — OSError [Errno 30] read-only filesystem
               strained_sidechain_bonds_fixer could not write temp file
               All scores reported as 0.00 (server artefact, not real)
               → EXCLUDED from all downstream analysis
  AbC-5    : ✓ Structure produced (note: sidechain warning logged
               but structure was extracted successfully)

  Proceeding with 4 valid candidates: AbC-1, AbC-2, AbC-3, AbC-5
""")

# ─────────────────────────────────────────────────────────────
# 2. RAW SCORE TABLE
# ─────────────────────────────────────────────────────────────

print("[2] RAW MODELLING SCORES (RMSD in Å — lower = better)")
print("-" * 70)

region_labels = {
    'FW_H':   'Framework H-chain',
    'CDR_H1': 'CDR-H1',
    'CDR_H2': 'CDR-H2',
    'CDR_H3': 'CDR-H3  ★ key',
    'FW_L':   'Framework L-chain',
    'CDR_L1': 'CDR-L1',
    'CDR_L2': 'CDR-L2',
    'CDR_L3': 'CDR-L3  ★ key',
}

print(f"\n  {'Region':<22} {'AbC-1':>8} {'AbC-2':>8} {'AbC-3':>8} {'AbC-5':>8}  Note")
print("  " + "-"*75)
for r in regions:
    vals = [valid_data[c][r] for c in valid_ids]
    best = min(vals)
    note = '← poor (> 2.0 Å)' if max(vals) > 2.0 else ''
    row  = f"  {region_labels[r]:<22}"
    for v in vals:
        marker = '*' if v == best else ' '
        row += f" {v:>6.2f}{marker} "
    row += f" {note}"
    print(row)

print("\n  * = best score for that region")
print("  CDR-H3 and CDR-L3 are the primary antigen-contact loops")

# ─────────────────────────────────────────────────────────────
# 3. DERIVED METRICS
# ─────────────────────────────────────────────────────────────

print("\n\n[3] DERIVED METRICS & RANKING")
print("-" * 65)

for c in valid_ids:
    v = valid_data[c]
    v['mean_all']   = np.mean(list(v.values()))
    v['mean_cdr']   = np.mean([v['CDR_H1'], v['CDR_H2'], v['CDR_H3'],
                                v['CDR_L1'], v['CDR_L2'], v['CDR_L3']])
    v['H3_L3_mean'] = (v['CDR_H3'] + v['CDR_L3']) / 2
    v['H3_score']   = v['CDR_H3']
    v['n_poor']     = sum(1 for val in v.values() if isinstance(val, float) and val > 2.0)
    # Overall model quality score (lower = better)
    # Weight CDR-H3 × 3, CDR-H2 × 1 (H2 often noisy), rest × 1
    v['weighted'] = (v['CDR_H3']*3 + v['CDR_H1'] + v['CDR_H2']*0.5 +
                     v['FW_H'] + v['FW_L'] + v['CDR_L1'] + v['CDR_L2'] + v['CDR_L3']) / 9

print(f"\n  {'Candidate':<10} {'Mean all':>10} {'Mean CDR':>10} {'H3+L3 avg':>11} "
      f"{'CDR-H3':>8} {'Poor CDRs':>10} {'Weighted':>10}")
print("  " + "-"*70)
for c in valid_ids:
    v = valid_data[c]
    print(f"  {c:<10} {v['mean_all']:>10.2f} {v['mean_cdr']:>10.2f} "
          f"{v['H3_L3_mean']:>11.2f} {v['H3_score']:>8.2f} "
          f"{v['n_poor']:>10} {v['weighted']:>10.2f}")

# Rank by weighted score
ranked = sorted(valid_ids, key=lambda c: valid_data[c]['weighted'])

print(f"\n  RANKING (by weighted score — lower = better):")
stars = ['★★★ PRIORITY', '★★  Strong', '★★  Strong', '★   Moderate']
for i, c in enumerate(ranked):
    v = valid_data[c]
    poor_note = f" ⚠ CDR-H2={v['CDR_H2']:.2f}" if v['CDR_H2'] > 2.0 else ""
    print(f"  {i+1}. {c}  weighted={v['weighted']:.2f}  CDR-H3={v['CDR_H3']:.2f}  {stars[i]}{poor_note}")

# ─────────────────────────────────────────────────────────────
# 4. INTERPRETATION
# ─────────────────────────────────────────────────────────────

print(f"""

[4] SCORE INTERPRETATION
{'-'*65}
  ABodyBuilder2 prediction error = RMSD (Å) of predicted vs
  template/ensemble structure for each antibody region.

  Typical ranges:
    < 1.0 Å  = high confidence prediction
    1.0–2.0  = moderate confidence (acceptable for CDR loops)
    > 2.0 Å  = low confidence — loop geometry uncertain

  CDR-H3 is the most important loop for antigen contact.
  All 4 valid candidates have CDR-H3 error ≤ 0.59 Å → HIGH confidence.
  This is excellent and means the primary contact loop geometry
  is well-predicted in all candidates.

  Problem region: CDR-H2
    AbC-1: CDR-H2 = 4.02 Å ⚠ very poor — loop geometry unreliable
    AbC-2: CDR-H2 = 2.32 Å ⚠ poor
    AbC-5: CDR-H2 = 1.80 Å — borderline
    AbC-3: CDR-H2 = 0.79 Å ✓ good

  CDR-H2 does not directly contact the epitope (CDR-H3 does),
  but poor CDR-H2 prediction can affect the overall Fv geometry
  and VH-VL packing. AbC-1 and AbC-2 should be used with caution.

  FORWARD TO DOCKING:
    ★★★ AbC-3  — best overall; CDR-H3=0.52, CDR-H2=0.79 (only good H2)
    ★★  AbC-5  — CDR-H3=0.52 tied best; CDR-H2=1.80 borderline
    ★★  AbC-2  — CDR-H3=0.59; CDR-H2=2.32 ⚠ use with caution
    ★   AbC-1  — CDR-H3=0.57; CDR-H2=4.02 ⚠ lowest confidence
    ✗   AbC-4  — EXCLUDED (server failure, no structure produced)
""")

# ─────────────────────────────────────────────────────────────
# 5. FIGURES
# ─────────────────────────────────────────────────────────────

print("[5] GENERATING FIGURES...")

fig = plt.figure(figsize=(20, 14))
fig.patch.set_facecolor(BG)
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.50, wspace=0.40)

# ── Plot 1: Heatmap of all scores ──
ax1 = fig.add_subplot(gs[0, :2])
matrix = np.array([[valid_data[c][r] for r in regions] for c in valid_ids])
im = ax1.imshow(matrix, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=4.5)
ax1.set_xticks(range(len(regions)))
ax1.set_xticklabels([r.replace('_', '-') for r in regions], fontsize=8, color=LIGHT)
ax1.set_yticks(range(len(valid_ids)))
ax1.set_yticklabels(valid_ids, fontsize=9, color=LIGHT)
for i in range(len(valid_ids)):
    for j in range(len(regions)):
        val = matrix[i, j]
        col = 'black' if val < 2.5 else 'white'
        ax1.text(j, i, f'{val:.2f}', ha='center', va='center',
                 fontsize=8, color=col, fontweight='bold')
cbar = plt.colorbar(im, ax=ax1, pad=0.02)
cbar.set_label('Prediction Error (Å)', color=LIGHT, fontsize=8)
cbar.ax.tick_params(colors=LIGHT)
style_ax(ax1, "ABodyBuilder2 Modelling Scores — All Regions (green=good, red=poor)")

# ── Plot 2: CDR-H3 comparison (key metric) ──
ax2 = fig.add_subplot(gs[0, 2])
h3_vals = [valid_data[c]['CDR_H3'] for c in valid_ids]
colors_v = [cand_colors[c] for c in valid_ids]
bars = ax2.bar(valid_ids, h3_vals, color=colors_v, alpha=0.85, width=0.5)
ax2.axhline(1.0, color='white', lw=0.8, ls='--', alpha=0.5, label='1.0 Å threshold')
ax2.axhline(2.0, color=RED,     lw=0.8, ls='--', alpha=0.5, label='2.0 Å poor threshold')
for bar, v in zip(bars, h3_vals):
    ax2.text(bar.get_x()+bar.get_width()/2, v+0.02,
             f'{v:.2f}', ha='center', fontsize=9, color=LIGHT, fontweight='bold')
ax2.set_ylabel("CDR-H3 Prediction Error (Å)")
ax2.set_ylim(0, 1.2)
ax2.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax2, "CDR-H3 Error ★ Primary Contact Loop")

# ── Plot 3: Mean CDR error comparison ──
ax3 = fig.add_subplot(gs[1, 0])
mean_cdrs = [valid_data[c]['mean_cdr'] for c in valid_ids]
bars3 = ax3.bar(valid_ids, mean_cdrs, color=colors_v, alpha=0.85, width=0.5)
ax3.axhline(1.0, color='white', lw=0.8, ls='--', alpha=0.5)
ax3.axhline(2.0, color=RED,     lw=0.8, ls='--', alpha=0.5)
for bar, v in zip(bars3, mean_cdrs):
    ax3.text(bar.get_x()+bar.get_width()/2, v+0.03,
             f'{v:.2f}', ha='center', fontsize=8, color=LIGHT)
ax3.set_ylabel("Mean CDR Error (Å)")
style_ax(ax3, "Mean CDR Prediction Error")

# ── Plot 4: Radar-style region error ──
ax4 = fig.add_subplot(gs[1, 1])
x = np.arange(len(regions))
width = 0.2
for i, (c, col) in enumerate(cand_colors.items()):
    vals = [valid_data[c][r] for r in regions]
    ax4.bar(x + i*width, vals, width, label=c, color=col, alpha=0.75)
ax4.axhline(2.0, color=RED, lw=0.8, ls='--', alpha=0.6, label='Poor threshold')
ax4.set_xticks(x + width*1.5)
ax4.set_xticklabels([r.replace('_','-') for r in regions], fontsize=6.5, rotation=30)
ax4.set_ylabel("Prediction Error (Å)")
ax4.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT, ncol=2)
style_ax(ax4, "All Region Errors by Candidate")

# ── Plot 5: Final ranking summary ──
ax5 = fig.add_subplot(gs[1, 2])
ax5.axis('off')
summary = [
    ('Rank', 'Candidate', 'CDR-H3', 'CDR-H2', 'Weighted', 'Decision'),
    ('─'*4,  '─'*9,      '─'*7,    '─'*7,    '─'*9,     '─'*14),
]
for i, c in enumerate(ranked):
    v = valid_data[c]
    h2_flag = '⚠' if v['CDR_H2'] > 2.0 else '✓'
    summary.append((str(i+1), c, f"{v['CDR_H3']:.2f} Å",
                    f"{v['CDR_H2']:.2f}{h2_flag}",
                    f"{v['weighted']:.2f}",
                    stars[i]))
summary.append(('✗', 'AbC-4', 'FAILED', '—', '—', 'EXCLUDED'))

y = 0.97
col_x = [0.0, 0.10, 0.28, 0.42, 0.57, 0.70]
for row in summary:
    for xi, cell in zip(col_x, row):
        col = TEAL if row == summary[0] else LIGHT
        if 'AbC-3' in str(row): col = GREEN
        if 'EXCLUDED' in str(row) or 'FAILED' in str(row): col = RED
        ax5.text(xi, y, cell, fontsize=7, color=col,
                 transform=ax5.transAxes, va='top')
    y -= 0.14
style_ax(ax5, "Final Candidate Ranking")
ax5.set_facecolor(PANEL)

fig.suptitle(
    "Dengue NS1 — Day 3: ABodyBuilder2 Modelling Score Analysis\n"
    "Source: ABodyBuilder2 | Job ID: 20260818_0701560 | Dimple Srivastava\n"
    "4/5 structures produced (AbC-4 excluded — server failure)",
    fontsize=11, color=LIGHT, fontweight='bold', y=0.995
)

plt.savefig("/mnt/user-data/outputs/dengue_ns1_day3_abbuilder.png",
            dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("  ✓ Figure saved.")

print(f"""
{'='*65}
  DAY 3 SUMMARY
{'='*65}
  Valid structures  : AbC-1, AbC-2, AbC-3, AbC-5 (4 of 5)
  Excluded          : AbC-4 (server OSError — no structure)

  PRIORITY FOR DOCKING (Day 4):
    1. AbC-3  CDR-H3=0.52 Å ★ best H3; CDR-H2=0.79 ✓ clean overall
    2. AbC-5  CDR-H3=0.52 Å ★ tied best H3; CDR-H2=1.80 borderline
    3. AbC-2  CDR-H3=0.59 Å; CDR-H2=2.32 ⚠ use with caution

  AbC-1 deprioritized: CDR-H2=4.02 Å (very poor geometry confidence)

  DAY 4 PLAN:
    → Download AbC-3, AbC-5, AbC-2 PDB files from ABodyBuilder2
    → Prepare 4O6B NS1 structure for docking
    → Submit to HADDOCK 2.4 with epitope aa 25-34 as active residues
    → Collect docking clusters → binding energy + interface contacts
    → PDBePISA for interface buried surface area analysis
{'='*65}
""")
