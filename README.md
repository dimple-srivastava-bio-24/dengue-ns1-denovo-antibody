# De Novo Computational Antibody Development Against Dengue NS1

**Author:** Dimple Srivastava  
**Affiliation:** MSc Biotechnology, J.C. Bose University of Science & Technology (YMCA), Faridabad  
**Dissertation:** In silico identification of SARS-CoV-2 spike protein inhibitors — Dept. of Virology, Jamia Hamdard  

---

## Project Summary

A complete 4-day computational pipeline for de novo antibody development targeting the Dengue virus NS1 protein (DENV1 and DENV2). Starting from raw amino acid sequences, this project identifies a conserved surface-accessible epitope, designs five antibody candidates with human germline frameworks, predicts their 3D structures, and evaluates binding through protein-protein docking and interface analysis.

**Priority candidate identified: AbC-5**
- HADDOCK score: −83.6 ± 5.1 kcal/mol
- Buried Surface Area: 1688.4 ± 121.4 Å²
- Interface H-bonds: 8
- ABodyBuilder2 CDR-H3 prediction error: 0.52 Å

---

## Scientific Rationale

Dengue virus causes ~400 million infections annually across four serotypes (DENV1–4). Sequential infection with different serotypes can trigger severe dengue through antibody-dependent enhancement (ADE). NS1 is a secreted glycoprotein essential for viral replication and immune evasion — it activates complement pathways and contributes to vascular leakage in severe dengue.

A broad-spectrum antibody targeting conserved NS1 epitopes shared between DENV1 and DENV2 could offer cross-protective potential. This project identifies such a conserved, surface-accessible target region and develops computationally designed antibody candidates against it.

---

## Pipeline Overview

```
Den1 NS1 + Den2 NS1 sequences
          ↓
Day 1: Sequence QC & conservation analysis
          ↓
Day 2: PDB structure mapping + PyMOL DSSP + quantitative SASA
          ↓
Day 3: De novo antibody CDR design + ABodyBuilder2 structure prediction
          ↓
Day 4: HADDOCK 2.4 docking + PDBePISA interface analysis
          ↓
Day 5: Candidate ranking + final report
```

---

## Target Epitope

| Property | Value |
|---|---|
| Sequence | VHTWTEQYKF |
| Region | aa 25–34 (PDB 26–35 in 4O6B) |
| Conservation | 100% between DENV1 and DENV2 (supplied sequences) |
| Domain | β-roll loop, Den1 NS1 |
| Normalized RSA | 0.64 (PyMOL get_area on 4O6B) |
| DSSP | Mostly loop/irregular (aa 24–32); aa 33–34 strand |
| Literature support | Akey et al. 2014 — confirmed mAb contact residues |

---

## Day-by-Day Results

### Day 1 — Sequence QC & Conservation
- Den1 NS1: 352 aa, 1056 nt — stable (instability index 34.5)
- Den2 NS1: 352 aa, 1056 nt — borderline stable (41.6)
- Pairwise identity: **72.7%** (256 conserved / 96 variable positions)
- Multiple 100% conserved 10-aa windows identified
- Top candidate: aa 25–34 (VHTWTEQYKF) — 100% conserved

**Tools:** BioPython, Matplotlib, BLOSUM62 alignment

---

### Day 2 — Structural Mapping & Surface Accessibility
- PDB structures: **4O6B** (Den1 NS1) + **4OIG** (Den2 NS1)
- Loaded in PyMOL; chains verified; residue numbering offset confirmed (+1 for Candidate A)

| Metric | Candidate A (aa 24–34) | Candidate B (aa 193–207) |
|---|---|---|
| Structure | 4O6B chain A | 4OIG chain A |
| Absolute SASA | 1310.2 Å² | 1448.1 Å² |
| Per-residue SASA | 119.1 Å²/res | 96.5 Å²/res |
| Normalized RSA | 0.64 | 0.49 |
| DSSP | Mostly loop | Mixed strand/loop |
| Literature | Akey 2014 ✓ | Chen 2015 ✓ |

**Tools:** PyMOL (get_area, DSSP plugin), BioPython

---

### Day 3 — De Novo Antibody Design & Structure Prediction

**Framework:** Human germline IGHV1-2\*02 (VH) + IGKV1-39\*01 (VL)  
**CDR-H3 design strategy:** Epitope-guided — aromatic residues (Y/W) for π-stacking with H26/W28/Y32/F34; R/K for E30 charge complementarity

| Candidate | CDR-H3 | CDR-H3 error (Å) | CDR-H2 error (Å) | Status |
|---|---|---|---|---|
| AbC-1 | ARYGYWDYFDYW | 0.57 ✓ | 4.02 ✗ | Deprioritized |
| AbC-2 | ARDRYYWYFDV | 0.59 ✓ | 2.32 ⚠ | Proceed with caution |
| AbC-3 | AKGYYYYGMDV | 0.52 ✓ | 0.79 ✓ | ★★ Strong |
| AbC-4 | — | — | — | Failed (server error) |
| AbC-5 | ARRGYYYYWDFDYW | 0.52 ✓ | 1.80 ⚠ | ★★★ Priority |

**PyMOL CDR-H SASA validation (residues 95–102):**

| Candidate | CDR-H SASA | Normalized RSA |
|---|---|---|
| AbC-3 | 882.5 Å² | 0.98 |
| AbC-2 | 881.4 Å² | 0.97 |
| AbC-5 | 810.8 Å² | 1.09 |

**Tools:** ABodyBuilder2, PyMOL, BioPython

---

### Day 4 — Docking & Interface Analysis

**HADDOCK 2.4 setup:**
- Receptor: 4O6B chain A (Den1 NS1, 320 residues)
- Ligand: AbC-3/AbC-5/AbC-2 Fv structures (~212 residues)
- Active residues (NS1): PDB 26,27,28,29,30,31,32,33,34,35
- Active residues (antibody): CDR-H3 residues 95–102

**HADDOCK results (selected clusters):**

| Metric | AbC-3 (Cl.2) | AbC-5 (Cl.1) | AbC-2 (Cl.6) |
|---|---|---|---|
| HADDOCK score | −76.4 ± 1.4 | **−83.6 ± 5.1** ★ | −66.9 ± 3.7 |
| Eelec (kcal/mol) | −88.5 | **−171.0** ★ | −98.0 |
| Evdw (kcal/mol) | −62.0 | −55.8 | −57.8 |
| BSA (Å²) | 1618 ± 55 | **1688 ± 121** ★ | 1580 ± 94 |
| Cluster size | 34 | **57** ★ | 10 |
| Z-score | −1.0 | **−1.7** ★ | −1.7 ★ |

**PDBePISA interface analysis:**

| Metric | AbC-3 | AbC-5 | AbC-2 |
|---|---|---|---|
| H-bonds | 7 | **8** ★ | 6 |
| Salt bridges | **4** ★ | 2 | 0 |
| ΔiG p-value | 0.588 | 0.540 | 0.421 |

**Tools:** HADDOCK 2.4, PDBePISA, PyMOL

---

## Final Candidate Ranking

| Rank | Candidate | Key strengths | Limitation |
|---|---|---|---|
| ★★★ | **AbC-5** | Best HADDOCK score, BSA, H-bonds, cluster size | CDR-H2 borderline (1.80 Å) |
| ★★ | AbC-3 | Cleanest Fv geometry (CDR-H2=0.79 Å), 4 salt bridges | HADDOCK score slightly weaker |
| ★ | AbC-2 | Best ΔiG p-value (0.421) | No salt bridges, poor CDR-H2 |

---

## Scientific Caveats & Limitations

- SASA values are from PyMOL get_area (4O6B/4OIG) — not MD-derived
- All ΔiG p-values > 0.05 — interfaces not statistically confirmed by PDBePISA criteria
- High HADDOCK RMSD (10–14 Å) indicates docking solution diversity
- High restraint violation energies reflect incomplete AIR satisfaction
- AbC-4 excluded due to ABodyBuilder2 server failure (OSError)
- No experimental binding data — all results are computational predictions

---

## Tools & Resources Used

| Tool | Purpose | Access |
|---|---|---|
| BioPython | Sequence analysis, alignment | pip install biopython |
| PyMOL | Structure visualization, DSSP, SASA | pymol.org |
| ABodyBuilder2 | Antibody structure prediction | opig.stats.ox.ac.uk/webapps/abodybuilder2 |
| HADDOCK 2.4 | Protein-protein docking | wenmr.science.uu.nl/haddock2.4 |
| PDBePISA | Interface analysis | ebi.ac.uk/pdbe/pisa |
| RCSB PDB | Structure download (4O6B, 4OIG) | rcsb.org |
| Matplotlib/Seaborn | Data visualization | pip install matplotlib seaborn |

---

## Repository Structure

```
dengue-ns1-denovo-antibody/
├── README.md
├── requirements.txt
├── sequences/
│   ├── Den1_NS1_AA.fasta
│   ├── Den2_NS1_AA.fasta
│   └── antibody_candidates.fasta
├── day1_sequence_analysis/
│   ├── dengue_ns1_day1.py
│   └── dengue_ns1_day1_analysis.png
├── day2_structural_analysis/
│   ├── dengue_ns1_day2_final.py
│   └── dengue_ns1_day2_final.png
├── day3_antibody_design/
│   ├── dengue_ns1_day3.py
│   ├── dengue_ns1_day3_abbuilder.py
│   ├── dengue_ns1_day3_final.py
│   └── dengue_ns1_day3_final.png
├── day4_docking/
│   ├── dengue_ns1_day4_final.py
│   └── dengue_ns1_day4_final.png
└── results/
    └── DengueNS1_Antibody_Project_Day5.pptx
```

---

## References

1. Akey DL et al. (2014) Flavivirus NS1 structures reveal surfaces for associations with membranes and the immune system. *Science* 343:881–885.
2. Chen HR et al. (2015) Anti-dengue virus nonstructural protein 1 antibodies cause NO-mediated endothelial cell apoptosis. *PNAS* 112(19):6076–6081.
3. Edeling MA et al. (2014) Potent dengue virus neutralization by a therapeutic antibody with low monovalent affinity. *PLOS Pathogens* 10(9).
4. Flamand M et al. (1999) Dengue virus type 1 nonstructural glycoprotein NS1 is secreted from mammalian cells as a soluble hexamer. *J Virology* 73(7):6104–6110.
5. van Zundert GCP et al. (2016) The HADDOCK2.2 web server. *J Mol Biol* 428(4):720–725.
6. Tien MZ et al. (2013) Maximum allowed solvent accessibilities of residues in proteins. *PLOS ONE* 8(11).
7. Schneider C et al. (2022) Bi-specific antibodies designed to target NS1. *J Virology* (ABodyBuilder2 reference).

---

## Contact

**Dimple Srivastava**  
MSc Biotechnology | J.C. Bose University (YMCA), Faridabad  
Dissertation: In silico identification of SARS-CoV-2 spike protein inhibitors  
Dept. of Virology, Jamia Hamdard
