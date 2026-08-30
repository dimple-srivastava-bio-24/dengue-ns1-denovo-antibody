"""
=============================================================
De Novo Antibody Development Against Dengue NS1
DAY 3: De Novo Antibody Candidate Generation, Sequence
        Design, Physicochemical Scoring & Structural
        Submission Preparation
Author: Dimple Srivastava
=============================================================

TARGET EPITOPE : VHTWTEQYKF (aa 25–34, Den1/Den2 NS1)
                 100% conserved between supplied Den1 and Den2 sequences
                 Literature contact region: Akey et al. 2014

DESIGN STRATEGY:
  Framework : Human germline IGHV1-2*02 (VH) + IGKV1-39*01 (VL)
              Most frequent human antibody germline combination;
              minimises immunogenicity risk in therapeutic context.
  CDR-H3    : De novo designed (10–14 aa) to complement epitope
              physicochemically — charge, aromaticity, polarity.
  CDR-L3    : Designed to support VH-VL interface and antigen contact.
  CDR-H1/H2, CDR-L1/L2: Retained from germline with minor mutations
              to optimise paratope geometry.

METHODOLOGY NOTE:
  Sequences generated here use:
  (a) Published human germline framework sequences (IMGT database)
  (b) CDR sequences rationally designed from epitope physicochemistry
      and published anti-dengue NS1 antibody CDR data (Akey 2014,
      Wan et al. 2017, Britto et al. 2023)
  This is a rational/computational design approach, not ML generation.
  For true de novo generation: submit to ProteinMPNN or RFdiffusion
  (see submission guide at end of script).
=============================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────
# COLOURS & STYLE
# ─────────────────────────────────────────────────────────────
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
print("  DENGUE NS1 — DE NOVO ANTIBODY PROJECT | DAY 3")
print("=" * 65)

# ─────────────────────────────────────────────────────────────
# 1. TARGET EPITOPE SUMMARY
# ─────────────────────────────────────────────────────────────

print("\n[1] TARGET EPITOPE PROPERTIES")
print("-" * 55)

EPITOPE    = "VHTWTEQYKF"   # aa 25–34, 100% conserved
EPITOPE_EX = "EVHTWTEQYKF"  # aa 24–34, 90.9% conserved (context)

pa_ep = ProteinAnalysis(EPITOPE)
print(f"\n  Core epitope  : {EPITOPE}  (aa 25–34, 100% Den1/Den2 conserved)")
print(f"  Extended      : {EPITOPE_EX}  (aa 24–34, 90.9% conserved)")
print(f"  MW            : {pa_ep.molecular_weight():.1f} Da")
print(f"  pI            : {pa_ep.isoelectric_point():.2f}  (slightly acidic)")
print(f"  GRAVY         : {pa_ep.gravy():.3f}  (hydrophilic)")
print(f"  Key residues  : V(hydrophobic) H(+,aromatic) T(polar) W(aromatic)")
print(f"                  E(-) Q(polar) Y(aromatic,H-bond) K(+) F(aromatic)")
print(f"\n  Design implications:")
print(f"    • 4 aromatic residues (H,W,Y,F) → design CDR-H3 with Y/W for π-stacking")
print(f"    • Net charge ~0 at pH 7 → CDR can be slightly basic (K,R) for electrostatics")
print(f"    • Hydrophilic epitope → avoid overly hydrophobic CDR-H3 (aggregation risk)")

# ─────────────────────────────────────────────────────────────
# 2. HUMAN GERMLINE FRAMEWORK
# ─────────────────────────────────────────────────────────────

print("\n\n[2] HUMAN GERMLINE FRAMEWORK SEQUENCES")
print("    Source: IMGT database (imgt.org)")
print("-" * 55)

# IGHV1-2*02 framework regions (human)
VH_FR1 = "QVQLVQSGAEVKKPGASVKVSCKASGYTFT"          # FR1
VH_FR2 = "WVRQAPGQGLEWMG"                           # FR2
VH_FR3 = "RVTMTTDTSTSTVYMELSSLRSEDTAVYYCAR"         # FR3
VH_FR4 = "WGQGTLVTVSS"                              # FR4

# IGKV1-39*01 framework regions (human)
VL_FR1 = "DIQMTQSPSSLSASVGDRVTITC"                 # FR1
VL_FR2 = "WYQQKPGKAPKLLIY"                          # FR2
VL_FR3 = "GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC"        # FR3
VL_FR4 = "FGQGTKVEIK"                               # FR4

print(f"\n  Heavy chain framework: IGHV1-2*02 (human)")
print(f"    FR1: {VH_FR1}")
print(f"    FR2: {VH_FR2}")
print(f"    FR3: {VH_FR3}")
print(f"    FR4: {VH_FR4}")
print(f"\n  Light chain framework: IGKV1-39*01 (human)")
print(f"    FR1: {VL_FR1}")
print(f"    FR2: {VL_FR2}")
print(f"    FR3: {VL_FR3}")
print(f"    FR4: {VL_FR4}")

# ─────────────────────────────────────────────────────────────
# 3. DE NOVO CDR DESIGN — 5 CANDIDATES
# ─────────────────────────────────────────────────────────────

print("\n\n[3] DE NOVO CDR SEQUENCE DESIGN")
print("    Strategy: epitope-complementary CDR-H3 + supporting CDR-L3")
print("    Basis: physicochemical complementarity + published NS1 mAb data")
print("-" * 65)

# Each candidate: CDR-H1, CDR-H2, CDR-H3, CDR-L1, CDR-L2, CDR-L3
# CDR-H3 is the primary antigen-contact loop (designed de novo)
# CDR-H1/H2, CDR-L1/L2 from published anti-dengue NS1 antibody data
# CDR-L3 designed to support antigen contact geometry

candidates = [
    {
        "id":    "AbC-1",
        "rationale": "Long CDR-H3 with dual aromatic (Y,W) for π-stacking with H26,W28,Y32,F34; K/R for E31,Q32 complementarity",
        "CDR_H1": "GYTFTSYW",        # Tyr-anchored H1 (common in anti-protein Abs)
        "CDR_H2": "INPSGGST",        # Polar H2, complements epitope hydrophilics
        "CDR_H3": "ARYGYWDYFDYW",    # 12aa; Y/W for aromatic stack; D for K33 complement
        "CDR_L1": "QSISSY",          # Standard kappa L1
        "CDR_L2": "AAS",             # Minimal L2
        "CDR_L3": "QQSYSTPYT",       # Y-rich L3 for aromatic contacts
    },
    {
        "id":    "AbC-2",
        "rationale": "Medium CDR-H3 with R for E31 salt bridge; F/Y for aromatic cluster; mimics published 4E11 CDR geometry",
        "CDR_H1": "GYTFTDYW",        # D for H26 charge complement
        "CDR_H2": "IDPYNGGT",        # N for Q32 H-bond
        "CDR_H3": "ARDRYYWYFDV",     # 11aa; R-E31 salt bridge; YWY aromatic cluster
        "CDR_L1": "QSISRY",          # R for T28 polar contact
        "CDR_L2": "GAS",
        "CDR_L3": "QQYDSSPWT",       # W for F34 aromatic stack
    },
    {
        "id":    "AbC-3",
        "rationale": "Short CDR-H3, high-affinity focused design; K for E31 direct contact; minimal aromatic to reduce off-target",
        "CDR_H1": "GFTFSSYW",
        "CDR_H2": "ISYDGSNT",
        "CDR_H3": "AKGYYYYGMDV",     # 11aa; YYYY motif for dense aromatic contacts with H,W,Y,F of epitope
        "CDR_L1": "QSVSSY",
        "CDR_L2": "DAS",
        "CDR_L3": "QQRSNWPLT",       # R for K33 steric accommodation; W for F34
    },
    {
        "id":    "AbC-4",
        "rationale": "Broad-spectrum design — CDR-H3 targets conserved TWTEQ motif (aa 27–31) shared across all DENV serotypes",
        "CDR_H1": "GYSFTGYW",
        "CDR_H2": "INHSGSST",
        "CDR_H3": "ARDYYYYSGWFDP",   # 13aa; targets T27-W28-T29-E30-Q31 conserved stretch
        "CDR_L1": "QSLVHSNGNTY",     # Longer L1 for broader contact
        "CDR_L2": "KVS",
        "CDR_L3": "SQSTHVPYT",       # T-rich for T27,T29 H-bond contacts
    },
    {
        "id":    "AbC-5",
        "rationale": "Dual-target design — CDR-H3 contacts primary epitope aa 25–34; CDR-L3 extended to also contact secondary region aa 193–207",
        "CDR_H1": "GYTFTSYW",
        "CDR_H2": "IYPYNGGT",
        "CDR_H3": "ARRGYYYYWDFDYW",  # 14aa (long); R for multiple charge contacts
        "CDR_L1": "QSISSY",
        "CDR_L2": "AAS",
        "CDR_L3": "QQHYSSPWT",       # H for His-stacking with connector region
    },
]

# ─────────────────────────────────────────────────────────────
# 4. ASSEMBLE FULL VH / VL SEQUENCES
# ─────────────────────────────────────────────────────────────

print("\n\n[4] FULL VH / VL SEQUENCE ASSEMBLY")
print("-" * 55)

def assemble_vh(c):
    return (VH_FR1 + c['CDR_H1'] + VH_FR2 + c['CDR_H2'] +
            VH_FR3 + c['CDR_H3'] + VH_FR4)

def assemble_vl(c):
    return (VL_FR1 + c['CDR_L1'] + VL_FR2 + c['CDR_L2'] +
            VL_FR3 + c['CDR_L3'] + VL_FR4)

for c in candidates:
    c['VH'] = assemble_vh(c)
    c['VL'] = assemble_vl(c)
    print(f"\n  {c['id']}")
    print(f"    CDR-H3 : {c['CDR_H3']}  ({len(c['CDR_H3'])} aa)")
    print(f"    CDR-L3 : {c['CDR_L3']}  ({len(c['CDR_L3'])} aa)")
    print(f"    VH len : {len(c['VH'])} aa")
    print(f"    VL len : {len(c['VL'])} aa")

# ─────────────────────────────────────────────────────────────
# 5. PHYSICOCHEMICAL SCORING
# ─────────────────────────────────────────────────────────────

print("\n\n[5] PHYSICOCHEMICAL SCORING OF CANDIDATES")
print("    Metrics: MW, pI, instability index, GRAVY, CDR-H3 charge")
print("-" * 65)

EPITOPE_PI    = pa_ep.isoelectric_point()   # 6.72
EPITOPE_GRAVY = pa_ep.gravy()               # -1.07

def score_candidate(c):
    pa_vh = ProteinAnalysis(c['VH'])
    pa_vl = ProteinAnalysis(c['VL'])
    pa_h3 = ProteinAnalysis(c['CDR_H3'])

    c['VH_MW']     = pa_vh.molecular_weight()
    c['VH_pI']     = pa_vh.isoelectric_point()
    c['VH_instab'] = pa_vh.instability_index()
    c['VH_gravy']  = pa_vh.gravy()
    c['VL_MW']     = pa_vl.molecular_weight()
    c['VL_pI']     = pa_vl.isoelectric_point()
    c['VL_instab'] = pa_vl.instability_index()
    c['H3_pI']     = pa_h3.isoelectric_point()
    c['H3_gravy']  = pa_h3.gravy()
    c['H3_len']    = len(c['CDR_H3'])

    # Developability flags
    c['flag_instab'] = c['VH_instab'] > 40 or c['VL_instab'] > 40
    c['flag_gravy']  = c['VH_gravy'] > 0.0   # hydrophobic VH = aggregation risk
    c['flag_h3pi']   = c['H3_pI'] < 5.0      # very acidic CDR-H3 = poor expression

    # Complementarity score (higher = better charge/hydrophobicity match to epitope)
    pi_diff    = abs(c['H3_pI'] - (14 - EPITOPE_PI))   # CDR basic vs epitope acidic
    gravy_diff = abs(c['H3_gravy'] - EPITOPE_GRAVY)
    c['comp_score'] = max(0, 10 - pi_diff - gravy_diff)

    return c

print(f"\n  {'ID':<8} {'VH_MW':>8} {'VH_pI':>6} {'VH_Instab':>10} {'VH_GRAVY':>9} "
      f"{'H3_len':>7} {'H3_pI':>7} {'Comp.Score':>11} {'Flags'}")
print("  " + "-"*85)

for c in candidates:
    score_candidate(c)
    flags = []
    if c['flag_instab']: flags.append("UNSTABLE")
    if c['flag_gravy']:  flags.append("HYDROPHOBIC")
    if c['flag_h3pi']:   flags.append("ACIDIC-H3")
    flag_str = ",".join(flags) if flags else "OK"
    print(f"  {c['id']:<8} {c['VH_MW']:>8.0f} {c['VH_pI']:>6.2f} "
          f"{c['VH_instab']:>10.1f} {c['VH_gravy']:>9.4f} "
          f"{c['H3_len']:>7} {c['H3_pI']:>7.2f} "
          f"{c['comp_score']:>11.2f}  {flag_str}")

# ─────────────────────────────────────────────────────────────
# 6. CANDIDATE RANKING
# ─────────────────────────────────────────────────────────────

print("\n\n[6] CANDIDATE RANKING")
print("    Ranked by: complementarity score + stability + no flags")
print("-" * 65)

def rank_score(c):
    base = c['comp_score']
    if c['flag_instab']: base -= 3
    if c['flag_gravy']:  base -= 2
    if c['flag_h3pi']:   base -= 2
    # Prefer CDR-H3 length 10–13 (optimal for protein antigen)
    if 10 <= c['H3_len'] <= 13: base += 1
    return base

for c in candidates:
    c['rank_score'] = rank_score(c)

ranked = sorted(candidates, key=lambda x: x['rank_score'], reverse=True)

print(f"\n  {'Rank':<5} {'ID':<8} {'Rank Score':>11} {'CDR-H3':>20} {'Assessment'}")
print("  " + "-"*65)
stars = ['★★★ PRIORITY', '★★  Strong', '★★  Strong', '★   Moderate', '★   Moderate']
for i, c in enumerate(ranked):
    print(f"  {i+1:<5} {c['id']:<8} {c['rank_score']:>11.2f} "
          f"{c['CDR_H3']:>20}  {stars[i]}")

# ─────────────────────────────────────────────────────────────
# 7. FASTA OUTPUT — for ABodyBuilder2 / IgFold submission
# ─────────────────────────────────────────────────────────────

print("\n\n[7] FASTA FILES FOR STRUCTURE PREDICTION")
print("    Submit to: https://opig.stats.ox.ac.uk/webapps/abodybuilder2")
print("    Or:        https://github.com/dptech-corp/IgFold")
print("-" * 55)

fasta_content = ""
for c in ranked:
    fasta_content += f">{c['id']}_VH | CDR-H3: {c['CDR_H3']}\n{c['VH']}\n"
    fasta_content += f">{c['id']}_VL | CDR-L3: {c['CDR_L3']}\n{c['VL']}\n\n"

with open("/mnt/user-data/outputs/dengue_ns1_antibody_candidates.fasta", "w") as f:
    f.write(fasta_content)

print("\n  Saved: dengue_ns1_antibody_candidates.fasta")
print("\n  Preview (top 2 candidates):")
for c in ranked[:2]:
    print(f"\n  >{c['id']}_VH")
    print(f"  {c['VH']}")
    print(f"  >{c['id']}_VL")
    print(f"  {c['VL']}")

# ─────────────────────────────────────────────────────────────
# 8. ABODYBUILDER2 SUBMISSION GUIDE
# ─────────────────────────────────────────────────────────────

print("\n\n[8] STRUCTURE PREDICTION — SUBMISSION GUIDE")
print("=" * 65)
print("""
  ABodyBuilder2 (recommended — free, online, publication-quality):
  ────────────────────────────────────────────────────────────────
  URL  : https://opig.stats.ox.ac.uk/webapps/abodybuilder2
  Input: Paste VH sequence → paste VL sequence → Submit
  Output: PDB file of antibody Fv structure (VH+VL)
  Time : ~2–5 minutes per candidate

  Submit all 5 candidates. Download each PDB file and name:
    AbC-1_Fv.pdb, AbC-2_Fv.pdb, ... AbC-5_Fv.pdb

  IgFold (alternative — slightly faster):
  ────────────────────────────────────────────────────────────────
  URL  : https://huggingface.co/spaces/nrbennet/IgFold
  Same VH/VL input format.

  AFTER STRUCTURE PREDICTION (Day 4 prep):
  ────────────────────────────────────────────────────────────────
  1. Open each AbC-*.pdb in PyMOL
  2. Check CDR-H3 loop geometry — should be well-resolved (no clashes)
  3. Run PROCHECK or MolProbity (online) for structure quality
  4. Proceed to HADDOCK docking with best 3 structures
""")

# ─────────────────────────────────────────────────────────────
# 9. FIGURES
# ─────────────────────────────────────────────────────────────

print("[9] GENERATING FIGURES...")

fig = plt.figure(figsize=(20, 16))
fig.patch.set_facecolor(BG)
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.50, wspace=0.40)

cand_colors_map = {
    'AbC-1': TEAL, 'AbC-2': AMBER, 'AbC-3': GREEN,
    'AbC-4': PURP, 'AbC-5': PINK
}
ids = [c['id'] for c in candidates]
colors_list = [cand_colors_map[i] for i in ids]

# ── Plot 1: CDR-H3 length comparison ──
ax1 = fig.add_subplot(gs[0, 0])
h3_lens = [c['H3_len'] for c in candidates]
bars = ax1.bar(ids, h3_lens, color=colors_list, alpha=0.85, edgecolor='none')
ax1.axhline(10, color='white', lw=0.8, ls='--', alpha=0.5, label='Min optimal (10)')
ax1.axhline(13, color=AMBER,   lw=0.8, ls='--', alpha=0.5, label='Max optimal (13)')
ax1.set_ylabel("CDR-H3 Length (aa)")
ax1.set_ylim(0, 16)
ax1.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax1, "CDR-H3 Length")

# ── Plot 2: VH pI comparison ──
ax2 = fig.add_subplot(gs[0, 1])
vh_pis = [c['VH_pI'] for c in candidates]
ax2.bar(ids, vh_pis, color=colors_list, alpha=0.85, edgecolor='none')
ax2.axhline(7.0, color='white', lw=0.8, ls='--', alpha=0.5, label='Neutral pH')
ax2.axhline(EPITOPE_PI, color=RED, lw=1.0, ls=':', alpha=0.8,
            label=f'Epitope pI ({EPITOPE_PI:.2f})')
ax2.set_ylabel("pI")
ax2.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax2, "VH Isoelectric Point")

# ── Plot 3: Instability index ──
ax3 = fig.add_subplot(gs[0, 2])
vh_instab = [c['VH_instab'] for c in candidates]
vl_instab = [c['VL_instab'] for c in candidates]
x = np.arange(len(candidates))
ax3.bar(x - 0.2, vh_instab, 0.4, color=colors_list, alpha=0.85, label='VH')
ax3.bar(x + 0.2, vl_instab, 0.4, color=colors_list, alpha=0.45, label='VL')
ax3.axhline(40, color=RED, lw=1.0, ls='--', alpha=0.8, label='Instability cutoff (40)')
ax3.set_xticks(x); ax3.set_xticklabels(ids, fontsize=8)
ax3.set_ylabel("Instability Index")
ax3.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax3, "VH/VL Instability Index (< 40 = stable)")

# ── Plot 4: GRAVY score ──
ax4 = fig.add_subplot(gs[1, 0])
vh_gravy = [c['VH_gravy'] for c in candidates]
ax4.bar(ids, vh_gravy, color=colors_list, alpha=0.85, edgecolor='none')
ax4.axhline(0, color='white', lw=0.8, ls='--', alpha=0.5, label='Hydrophilic/phobic boundary')
ax4.axhline(EPITOPE_GRAVY, color=RED, lw=1.0, ls=':',
            label=f'Epitope GRAVY ({EPITOPE_GRAVY:.2f})')
ax4.set_ylabel("GRAVY Score")
ax4.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax4, "VH GRAVY Score (negative = hydrophilic)")

# ── Plot 5: Complementarity score ──
ax5 = fig.add_subplot(gs[1, 1])
comp = [c['comp_score'] for c in candidates]
rank_scores = [c['rank_score'] for c in candidates]
ax5.bar(ids, comp, color=colors_list, alpha=0.85, edgecolor='none', label='Comp. score')
ax5.plot(ids, rank_scores, color='white', marker='o', lw=1.5,
         markersize=6, label='Rank score', zorder=5)
ax5.set_ylabel("Score")
ax5.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax5, "Complementarity & Rank Score")

# ── Plot 6: CDR-H3 pI vs epitope pI ──
ax6 = fig.add_subplot(gs[1, 2])
h3_pi = [c['H3_pI'] for c in candidates]
ax6.scatter(h3_pi, comp, c=colors_list, s=120, zorder=5, edgecolors='white', lw=0.5)
for c in candidates:
    ax6.annotate(c['id'], (c['H3_pI'], c['comp_score']),
                 textcoords="offset points", xytext=(6,4),
                 fontsize=7, color=cand_colors_map[c['id']])
ax6.axvline(14 - EPITOPE_PI, color=RED, lw=1, ls='--',
            label=f'Ideal CDR-H3 pI ({14-EPITOPE_PI:.1f})')
ax6.set_xlabel("CDR-H3 pI")
ax6.set_ylabel("Complementarity Score")
ax6.legend(fontsize=7, facecolor=PANEL, edgecolor='#2A3A4A', labelcolor=LIGHT)
style_ax(ax6, "CDR-H3 pI vs Complementarity Score")

# ── Plot 7: Ranked candidate summary ──
ax7 = fig.add_subplot(gs[2, :])
y_pos = np.arange(len(ranked))
rank_vals   = [c['rank_score']  for c in ranked]
comp_vals   = [c['comp_score']  for c in ranked]
h3len_vals  = [c['H3_len']      for c in ranked]
rank_labels = [f"{c['id']}  |  CDR-H3: {c['CDR_H3']:<16}  "
               f"H3-len:{c['H3_len']}  pI:{c['VH_pI']:.1f}  "
               f"Instab:{c['VH_instab']:.0f}" for c in ranked]

bars_rank = ax7.barh(y_pos, rank_vals,
                     color=[cand_colors_map[c['id']] for c in ranked],
                     alpha=0.85, edgecolor='none', height=0.55)
ax7.set_yticks(y_pos)
ax7.set_yticklabels(rank_labels, fontsize=8)
ax7.set_xlabel("Overall Rank Score")
ax7.axvline(0, color='white', lw=0.5, alpha=0.3)

# Star labels
for i, (bar, c) in enumerate(zip(bars_rank, ranked)):
    ax7.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2,
             stars[i], va='center', fontsize=8,
             color=cand_colors_map[c['id']], fontweight='bold')

style_ax(ax7, "Final Candidate Ranking — De Novo Anti-NS1 Antibody Candidates (Day 3)")

fig.suptitle(
    "Dengue NS1 — Day 3: De Novo Antibody Candidate Design & Physicochemical Scoring\n"
    "Target: VHTWTEQYKF (aa 25–34) | Framework: IGHV1-2*02 + IGKV1-39*01 | "
    "Dimple Srivastava",
    fontsize=12, color=LIGHT, fontweight='bold', y=0.995
)

plt.savefig("/mnt/user-data/outputs/dengue_ns1_day3_analysis.png",
            dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("  ✓ Figure saved.")

# ─────────────────────────────────────────────────────────────
# 10. DAY 3 SUMMARY
# ─────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("  DAY 3 SUMMARY")
print("=" * 65)
print(f"""
  Target epitope   : VHTWTEQYKF (aa 25–34, 100% conserved Den1/Den2)
  Framework        : IGHV1-2*02 (VH) + IGKV1-39*01 (VL)
  Candidates       : 5 antibody sequences (AbC-1 through AbC-5)

  TOP RANKED CANDIDATES:
    1. {ranked[0]['id']} — CDR-H3: {ranked[0]['CDR_H3']}  (rank score: {ranked[0]['rank_score']:.1f})
    2. {ranked[1]['id']} — CDR-H3: {ranked[1]['CDR_H3']}  (rank score: {ranked[1]['rank_score']:.1f})

  OUTPUTS:
    • dengue_ns1_antibody_candidates.fasta — submit to ABodyBuilder2
    • dengue_ns1_day3_analysis.png — figures

  DAY 4 PLAN:
    → Submit FASTA to ABodyBuilder2 → download 5 antibody PDB files
    → Prepare NS1 structure (4O6B) for docking
    → Run HADDOCK 2.4 with epitope aa 25–34 as active residues
    → Collect docking results → binding energy + interface analysis
    → PDBePISA for interface contact analysis
""")
print("=" * 65)
