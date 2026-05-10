# 🧬 Zebrafish Fin Regeneration — Bulk RNA-seq Analysis

This repository contains a comprehensive bulk RNA-seq analysis workflow focused on **zebrafish caudal fin regeneration** across multiple regeneration stages (0dpa, 1dpa, and 4dpa).

The project investigates transcriptomic changes during regeneration with a particular emphasis on **RNA Polymerase I–related genes**, regenerative signaling pathways, and differential gene expression dynamics.

The dataset used in this project was obtained from **GEO Omnibus Series GSE231956**.

---

# 🔬 Project Overview

Zebrafish possess remarkable regenerative abilities, making them a valuable model organism for studying tissue regeneration and repair.

This project explores gene expression changes occurring during early regeneration stages following caudal fin amputation using bulk RNA sequencing approaches.

Key objectives include:

* Identification of differentially expressed genes (DEGs)
* Analysis of regeneration-associated pathways
* Investigation of RNA Polymerase I machinery during regeneration
* Functional enrichment and pathway visualization
* Generation of publication-style visualizations

---

# ⚙️ Analysis Workflow

The computational workflow includes:

1. Transcript quantification using **Salmon**
2. Transcript-to-gene mapping via **makeTxDbFromGFF**
3. Gene-level quantification using **tximport**
4. Differential expression analysis using **DESeq2**
5. Pairwise comparisons:

   * 1dpa vs 0dpa
   * 4dpa vs 0dpa
   * 4dpa vs 1dpa
6. Functional enrichment analysis:

   * Gene Ontology (GO)
   * KEGG pathway analysis
7. Visualization and interpretation of biological pathways

---

# 🧬 RNA Polymerase I–Focused Analysis

Special emphasis was placed on genes associated with RNA Polymerase I activity and ribosomal biogenesis, including:

* `rrn3`
* `taf1b`
* `polr1b`
* `polr1e`

Expression dynamics, differential expression status, and pathway associations of these genes were analyzed throughout the regeneration process.

---

# 📊 Key Visualizations

The project includes multiple downstream visualizations such as:

* MA plots
* Volcano plots
* Heatmaps of top differentially expressed genes
* Venn diagrams of overlapping DEGs
* GO enrichment dotplots
* KEGG pathway enrichment plots
* Gene expression visualization for RNA Polymerase I–related genes
* Enrichment network maps (emapplots)

---

# 📁 Repository Structure

```bash
bulk_rnaseq_zebrafish_regeneration/
├── README.md
├── environment.yml
├── samples.csv
├── data/
│   └── salmon_quant/
├── scripts/
│   ├── 01_tximport_deseq2.R
│   ├── 02_plots_volcano_heatmap.R
│   ├── 03_enrichment_analysis.R
│   └── 04_polI_analysis.R
├── results/
│   ├── plots/
│   └── tables/
```

---

# 🛠️ Tools & Technologies

* R
* DESeq2
* tximport
* Salmon
* clusterProfiler
* enrichplot
* ggplot2
* AnnotationDbi
* Bioconductor packages

---

# 📦 Installation & Environment Setup

Clone the repository and install dependencies using Conda:

```bash
conda env create -f environment.yml
conda activate rnaseq-zfish
```

---

# 📈 Sample Outputs

## Volcano Plot (4dpa vs 0dpa)

```markdown
results/plots/volcano_4dpa_vs_0dpa.png
```

## GO Enrichment Dotplot

```markdown
results/plots/GO_dotplot_4dpa_vs_0dpa.png
```

---

# 🚀 Future Improvements

Planned additions to this project include:

* GSEA (Gene Set Enrichment Analysis)
* Advanced pathway visualization
* Time-series transcriptomic analysis
* Integration with single-cell RNA-seq datasets
* Publication-ready figure refinement

---

# 👨‍🔬 Author

## Akshat Jaiswal

Biomedical Engineering | Computational Biology | Bioinformatics

Research Interests:

* Aging Biology
* Regeneration Biology
* RNA Biology
* Transcriptomics
* Computational Genomics

---

# 🌐 Connect & Collaborate

* 🔗 LinkedIn: https://www.linkedin.com/in/akshat-jaiswal-i-omics-flow-lab-9234b5384/
* 🌍 Portfolio: https://omicsflowlab.netlify.app/
* ✉️ Work Email: [omicsflow.lab@gmail.com](mailto:omicsflow.lab@gmail.com)

---

# 📜 License

MIT License

---

# 📚 References & Inspiration

* Love MI et al. (2014) — DESeq2, *Genome Biology*
* Yu G et al. (2012) — clusterProfiler, *OMICS*
* Marques IJ et al. — Zebrafish regeneration studies

---

If you find this repository useful, feel free to ⭐ the project, fork it, or contribute.

Open to collaborations, computational biology opportunities, and bioinformatics-related projects.
