"""
=============================================================
De Novo Antibody Development Against Dengue NS1
DAY 2 — FINAL (PyMOL-validated)
Structural Mapping, DSSP & Quantitative SASA Analysis
Author: Dimple Srivastava
=============================================================
SASA source: PyMOL get_area on 4O6B (Cand A) and 4OIG (Cand B)
DSSP source: PyMOL DSSP plugin on both structures
This replaces all estimated/proxy accessibility values from Day 2 v2.
=============================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

TEAL  = '#00B4D8'; AMBER = '#F4A261'; GREEN = '#57CC99'
RED   = '#E63946'; PURP  = '#9B72CF'; LIGHT = '#E0E0E0'
BG    = '#0D1B2A'; PANEL = '#162032'

def style_ax(ax, title):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=LIGHT, labelsize=8)
    ax.xaxis.label.set_color(LIGHT); ax.yaxis.label.set_color(LIGHT)
    ax.title.set_color(TEAL)
    ax.set_title(title, fontsize=9.5, fontweight='bold', pad=8)
    for spine in ax.spines.values():
        spine.set_edgecolor('#2A3A4A')

print("=" * 65)
print("  DENGUE NS1 — DAY 2 FINAL (PyMOL-validated)")
print("=" * 65)

# ─────────────────────────────────────────────────────────────
# REAL PyMOL DATA (from your analysis)
# ─────────────────────────────────────────────────────────────

data = {
    "A": {
        "pdb":          "4O6B",
        "chain":        "A",
        "seq_region":   "aa 24–34",
        "pdb_region":   "PDB 25–35",
        "sequence":     "EVHTWTEQYKF",
        "n_residues":   11,
        "dssp":         {24:'L',25:'L',26:'L',27:'L',28:'L',
                         29:'L',30:'L',31:'L',32:'L',33:'S',34:'S'},
        "region_sasa":  1310.156,
        "struct_sasa":  65127.566,
        "visual":       "Partially surface-accessible; some residues slightly embedded",
        "pdb_offset":   "+1 (seq aa 24 = PDB residue 25)",
        "conservation": "90.9% (E24→N in Den2; core VHTWTEQYKF 100%)",
        "literature":   "Akey et al. 2014 — confirmed mAb contact residues",
    },
    "B": {
        "pdb":          "4OIG",
        "chain":        "A",
        "seq_region":   "aa 193–207",
        "pdb_region":   "PDB 193–207",
        "sequence":     "AVHADMGYWIESEKN",
        "n_residues":   15,
        "dssp":         {**{i:'E' for i in range(193,197)},
                         **{i:'L' for i in range(197,200)},
                         **{i:'E' for i in range(200,207)},
                         207:'L'},
        "region_sasa":  1448.054,
        "struct_sasa":  69512.312,
        "visual":       "Partially exposed; portions facing outward, some partially embedded",
        "pdb_offset":   "None (seq and PDB numbering identical)",
        "conservation": "87% Den1/Den2 (AVHADMGYWIES fully conserved)",
        "literature":   "Chen et al. 2015 — cross-reactive Den1/Den2 epitope",
    }
}

# ─────────────────────────────────────────────────────────────
# 1. STRUCTURE VERIFICATION
# ─────────────────────────────────────────────────────────────

print("\n[1] STRUCTURE & CHAIN VERIFICATION")
print("-" * 55)
for cid, d in data.items():
    print(f"\n  Candidate {cid}:")
    print(f"    PDB structure  : {d['pdb']}  (chain {d['chain']})")
    print(f"    Sequence region: {d['seq_region']}  →  {d['sequence']}")
    print(f"    PDB mapping    : {d['pdb_region']}")
    print(f"    Numbering note : {d['pdb_offset']}")

# ─────────────────────────────────────────────────────────────
# 2. DSSP RESULTS
# ─────────────────────────────────────────────────────────────

print("\n\n[2] DSSP SECONDARY STRUCTURE RESULTS")
print("    Source: PyMOL DSSP plugin on 4O6B / 4OIG")
print("-" * 55)

ss_labels = {'L': 'Loop/Irregular', 'S': 'Strand', 'E': 'Strand', 'H': 'Helix'}

for cid, d in data.items():
    print(f"\n  Candidate {cid} ({d['pdb']}, {d['seq_region']}):")
    print(f"  {'Residue':<10} {'DSSP':<6} {'SS Type'}")
    print("  " + "-"*35)
    for res, ss in d['dssp'].items():
        print(f"  {res:<10} {ss:<6} {ss_labels.get(ss, ss)}")

print(f"""
  Summary:
    Candidate A: aa 24–32 loop/irregular | aa 33–34 strand
    Candidate B: aa 193–196 strand | aa 197–199 loop | aa 200–206 strand | aa 207 loop

  Note: DSSP characterizes secondary structure only.
  It was NOT used as evidence of solvent exposure (SASA measured separately).
""")

# ─────────────────────────────────────────────────────────────
# 3. QUANTITATIVE SASA (PyMOL get_area)
# ─────────────────────────────────────────────────────────────

print("\n[3] QUANTITATIVE SASA — PyMOL get_area")
print("-" * 55)

for cid, d in data.items():
    n   = d['n_residues']
    reg = d['region_sasa']
    tot = d['struct_sasa']
    per_res      = reg / n
    struct_mean  = tot / 352          # NS1 monomer = 352 aa
    norm_rsa     = per_res / struct_mean
    pct_struct   = reg / tot * 100
    d['per_res']     = per_res
    d['norm_rsa']    = norm_rsa
    d['pct_struct']  = pct_struct
    d['struct_mean'] = struct_mean

    print(f"\n  Candidate {cid} ({d['pdb']}, {d['seq_region']}):")
    print(f"    Region SASA (absolute)    : {reg:>10.3f} Å²")
    print(f"    Whole-structure SASA      : {tot:>10.3f} Å²")
    print(f"    Residues in region        : {n:>10}")
    print(f"    Per-residue SASA          : {per_res:>10.1f} Å²/residue")
    print(f"    Whole-structure mean/res  : {struct_mean:>10.1f} Å²/residue")
    print(f"    Normalized RSA            : {norm_rsa:>10.2f}  (region ÷ structure mean)")
    print(f"    % of whole-structure SASA : {pct_struct:>10.2f}%")

# ─────────────────────────────────────────────────────────────
# 4. NORMALIZED COMPARISON & INTERPRETATION
# ─────────────────────────────────────────────────────────────

print("\n\n[4] NORMALIZED SASA COMPARISON & INTERPRETATION")
print("-" * 65)

a = data['A']; b = data['B']
print(f"""
  {'Parameter':<35} {'Candidate A':>14} {'Candidate B':>14}
  {'-'*65}
  {'PDB structure':<35} {'4O6B':>14} {'4OIG':>14}
  {'Sequence region':<35} {'aa 24–34':>14} {'aa 193–207':>14}
  {'Residues':<35} {a['n_residues']:>14} {b['n_residues']:>14}
  {'Absolute SASA (Å²)':<35} {a['region_sasa']:>14.3f} {b['region_sasa']:>14.3f}
  {'Per-residue SASA (Å²/res)':<35} {a['per_res']:>14.1f} {b['per_res']:>14.1f}
  {'Whole-structure mean (Å²/res)':<35} {a['struct_mean']:>14.1f} {b['struct_mean']:>14.1f}
  {'Normalized RSA':<35} {a['norm_rsa']:>14.2f} {b['norm_rsa']:>14.2f}
  {'% of whole-structure SASA':<35} {a['pct_struct']:>13.2f}% {b['pct_struct']:>13.2f}%
  {'DSSP':<35} {'Mostly loop':>14} {'Mixed strand/loop':>14}
  {'Visual (PyMOL surface)':<35} {'Partial':>14} {'Partial':>14}

  ⚠ Raw SASA comparison is NOT valid here:
    Candidate B (15 res, 1448 Å²) vs Candidate A (11 res, 1310 Å²)
    Larger region → larger absolute SASA by default.

  Normalized per-residue comparison:
    Candidate A: {a['per_res']:.1f} Å²/res  (norm RSA = {a['norm_rsa']:.2f})
    Candidate B: {b['per_res']:.1f} Å²/res  (norm RSA = {b['norm_rsa']:.2f})

  Both normalized RSA values < 1.0:
    → Both regions are BELOW the whole-structure per-residue mean.
    → This is consistent with visual observation: \"partially accessible,
      portions slightly embedded.\"
    → Neither region is fully buried (norm RSA would be ~0.1–0.2)
      nor maximally exposed (norm RSA would be ~1.5–2.0).

  Candidate A is relatively more exposed per residue (0.64 vs 0.49).
  This is consistent with its predominantly loop DSSP character —
  loop regions are generally more solvent-accessible than strands.

  Biological relevance:
    Partially exposed loop regions are frequently targeted by
    neutralizing antibodies — they are accessible but conformationally
    constrained, which aids binding specificity.
    Both candidates remain structurally justified targets.
""")

# ─────────────────────────────────────────────────────────────
# 5. FINAL CANDIDATE TABLE
# ─────────────────────────────────────────────────────────────

print("\n[5] FINAL DAY 2 CANDIDATE TABLE")
print("=" * 65)
print(f"""
  {'Parameter':<28} {'Candidate A':<25} {'Candidate B'}
  {'-'*75}
  {'Structure':<28} {'4O6B (chain A)':<25} {'4OIG (chain A)'}
  {'Sequence region':<28} {'aa 24–34':<25} {'aa 193–207'}
  {'PDB mapping':<28} {'PDB 25–35 (+1 offset)':<25} {'PDB 193–207 (direct)'}
  {'Sequence':<28} {'EVHTWTEQYKF':<25} {'AVHADMGYWIESEKN'}
  {'Length':<28} {'11 aa':<25} {'15 aa'}
  {'Conservation':<28} {'90.9% (core 100%)':<25} {'87%'}
  {'DSSP':<28} {'Mostly loop (33–34: S)':<25} {'Strand/loop mixed'}
  {'Absolute SASA':<28} {'1310.156 Å²':<25} {'1448.054 Å²'}
  {'Per-residue SASA':<28} {'119.1 Å²/res':<25} {'96.5 Å²/res'}
  {'Normalized RSA':<28} {'0.64':<25} {'0.49'}
  {'Visual':<28} {'Partially exposed':<25} {'Partially exposed'}
  {'Literature':<28} {'Akey 2014 ✓':<25} {'Chen 2015 ✓'}
  {'Decision':<28} {'★★★ PRIMARY TARGET':<25} {'★★  Secondary target'}
""")

# ─────────────────────────────────────────────────────────────
# 6. DAY 2 CONCLUSION
# ─────────────────────────────────────────────────────────────

print("[6] DAY 2 CONCLUSION")
print("=" * 65)
print("""
  The β-roll loop region encompassing aa 24–34 (EVHTWTEQYKF) represents
  the leading conserved candidate region based on:
    • 90.9% sequence conservation (core VHTWTEQYKF: 100%)
    • PDB mapping confirmed on 4O6B chain A (PDB residues 25–35)
    • DSSP: predominantly loop/irregular — favourable for Ab access
    • Per-residue SASA: 119.1 Å²/res (norm RSA = 0.64)
    • Literature: confirmed mAb contact residues (Akey et al. 2014)

  The connector region aa 193–207 (AVHADMGYWIESEKN) is the secondary
  candidate:
    • 87% conservation; confirmed cross-reactive Den1/Den2 epitope
    • DSSP: mixed strand/loop; per-residue SASA 96.5 Å²/res (norm 0.49)
    • PDB mapping direct on 4OIG (no numbering offset)

  Both regions are partially surface-accessible (normalized RSA < 1.0),
  consistent with PyMOL molecular-surface visualization showing
  partial embedding. This is structurally appropriate for antibody
  targeting — conformationally constrained loop regions provide
  binding specificity.

  Day 2 status:
    Sequence conservation     ✅  (Day 1)
    PDB structure mapping     ✅  (4O6B / 4OIG)
    Chain verification        ✅
    3D visualization          ✅
    DSSP secondary structure  ✅
    Quantitative SASA         ✅  (PyMOL get_area)
    Normalized comparison     ✅
    Literature cross-ref      ✅
    Day 2 structural validation: COMPLETE

  → Proceed to Day 3: de novo antibody CDR design targeting aa 25–34
""")

# ─────────────────────────────────────────────────────────────
# 7. FIGURE
# ─────────────────────────────────────────────────────────────

print("[7] GENERATING FIGURE...")

fig = plt.figure(figsize=(18, 12))
fig.patch.set_facecolor(BG)
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.50, wspace=0.40)

# ── Plot 1: Absolute SASA comparison ──
ax1 = fig.add_subplot(gs[0, 0])
labels = ['Cand A\n(4O6B, 11 res)', 'Cand B\n(4OIG, 15 res)']
sasas  = [1310.156, 1448.054]
bars = ax1.bar(labels, sasas, color=[RED, AMBER], alpha=0.85, width=0.5)
for bar, v in zip(bars, sasas):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 15,
             f'{v:.1f} Å²', ha='center', fontsize=8, color=LIGHT)
ax1.set_ylabel("Absolute SASA (Å²)")
ax1.set_ylim(0, 1700)
style_ax(ax1, "Absolute SASA (⚠ not directly comparable)")

# ── Plot 2: Per-residue SASA — fair comparison ──
ax2 = fig.add_subplot(gs[0, 1])
per_res_vals = [119.1, 96.5]
struct_means = [185.0, 197.5]
x = np.arange(2)
ax2.bar(x - 0.2, per_res_vals, 0.35, color=[RED, AMBER], alpha=0.85, label='Region (Å²/res)')
ax2.bar(x + 0.2, struct_means, 0.35, color=[TEAL, TEAL], alpha=0.45, label='Struct mean (Å²/res)')
ax2.set_xticks(x); ax2.set_xticklabels(['Cand A\n(4O6B)', 'Cand B\n(4OIG)'])
ax2.set_ylabel("SASA per Residue (Å²/res)")
ax2.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax2, "Per-Residue SASA vs Whole-Structure Mean")

# ── Plot 3: Normalized RSA ──
ax3 = fig.add_subplot(gs[0, 2])
norm_vals = [0.64, 0.49]
bars3 = ax3.bar(['Cand A\n(4O6B)', 'Cand B\n(4OIG)'],
                norm_vals, color=[RED, AMBER], alpha=0.85, width=0.5)
ax3.axhline(1.0, color='white', lw=1.0, ls='--', alpha=0.6,
            label='Whole-structure mean (1.0)')
ax3.axhline(0.2, color=PURP, lw=0.8, ls=':', alpha=0.5, label='~Buried threshold')
for bar, v in zip(bars3, norm_vals):
    ax3.text(bar.get_x() + bar.get_width()/2, v + 0.02,
             f'{v:.2f}', ha='center', fontsize=9, color=LIGHT, fontweight='bold')
ax3.set_ylabel("Normalized RSA (region ÷ struct mean)")
ax3.set_ylim(0, 1.3)
ax3.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax3, "Normalized RSA — Fair Accessibility Comparison")

# ── Plot 4: DSSP breakdown Candidate A ──
ax4 = fig.add_subplot(gs[1, 0])
dssp_a = {'Loop/Irreg': 9, 'Strand': 2}
ax4.pie(dssp_a.values(), labels=dssp_a.keys(),
        colors=[TEAL, AMBER], autopct='%1.0f%%',
        startangle=90, textprops={'color': LIGHT, 'fontsize': 9})
style_ax(ax4, "Cand A DSSP (4O6B, aa 24–34)")

# ── Plot 5: DSSP breakdown Candidate B ──
ax5 = fig.add_subplot(gs[1, 1])
dssp_b = {'Strand': 10, 'Loop/Irreg': 5}
ax5.pie(dssp_b.values(), labels=dssp_b.keys(),
        colors=[AMBER, TEAL], autopct='%1.0f%%',
        startangle=90, textprops={'color': LIGHT, 'fontsize': 9})
style_ax(ax5, "Cand B DSSP (4OIG, aa 193–207)")

# ── Plot 6: Summary scorecard ──
ax6 = fig.add_subplot(gs[1, 2])
ax6.axis('off')
metrics = [
    ('Conservation',      '90.9%\n(core 100%)',   '87%'),
    ('DSSP',              'Loop (82%)',            'Strand/loop'),
    ('Norm. RSA',         '0.64',                  '0.49'),
    ('Per-res SASA',      '119.1 Å²',              '96.5 Å²'),
    ('Literature',        'Akey 2014 ✓',           'Chen 2015 ✓'),
    ('Decision',          '★★★ PRIMARY',           '★★  Secondary'),
]
y_start = 0.95
ax6.text(0.0,  y_start + 0.05, 'Metric',      fontsize=8, color=TEAL,  fontweight='bold', transform=ax6.transAxes)
ax6.text(0.42, y_start + 0.05, 'Cand A',      fontsize=8, color=RED,   fontweight='bold', transform=ax6.transAxes)
ax6.text(0.72, y_start + 0.05, 'Cand B',      fontsize=8, color=AMBER, fontweight='bold', transform=ax6.transAxes)
ax6.axhline(y_start + 0.02, color='#2A3A4A', lw=0.8, xmin=0, xmax=1)
for i, (metric, va, vb) in enumerate(metrics):
    y = y_start - (i + 1) * 0.14
    ax6.text(0.0,  y, metric, fontsize=7.5, color=LIGHT, transform=ax6.transAxes)
    ax6.text(0.42, y, va,     fontsize=7.5, color=RED,   transform=ax6.transAxes)
    ax6.text(0.72, y, vb,     fontsize=7.5, color=AMBER, transform=ax6.transAxes)
style_ax(ax6, "Candidate Summary Scorecard")
ax6.set_facecolor(PANEL)

fig.suptitle(
    "Dengue NS1 — Day 2 FINAL: PyMOL-Validated Structural Analysis\n"
    "SASA source: PyMOL get_area | DSSP source: PyMOL plugin | "
    "Structures: 4O6B (Den1) + 4OIG (Den2) | Dimple Srivastava",
    fontsize=11, color=LIGHT, fontweight='bold', y=0.995
)

plt.savefig("/mnt/user-data/outputs/dengue_ns1_day2_final.png",
            dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("  ✓ Figure saved.")
print("\n  Outputs: dengue_ns1_day2_final.png")
print("=" * 65)
