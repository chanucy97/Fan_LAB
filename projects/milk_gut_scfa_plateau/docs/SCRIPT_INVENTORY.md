# Script Inventory

## Analysis

| Script | Purpose |
| --- | --- |
| `scripts/analysis/gut_microbes_lagged_sensitivity.R` | Refit nominally significant lagged milk-to-SCFA models after dropping one dyad at a time. |

## Figure redraws

| Figure | Script | Purpose |
| --- | --- | --- |
| Fig. 1 | `scripts/figures/fig1_cohort_design_overview/redraw_fig1_coverage_heatmap.R` | Build modular cohort coverage heatmap and source tables. |
| Fig. 1 | `scripts/figures/fig1_cohort_design_overview/assemble_fig1_complete.R` | Assemble Fig. 1 overview and coverage panels. |
| Fig. 2 | `scripts/figures/fig2_milk_bioactive_remodeling/redraw_fig2_milk_bioactive_remodeling.R` | Redraw milk LTF/HMO/LCFA remodeling figure. |
| Fig. 3 | `scripts/figures/fig3_infant_microbiome_maturation/redraw_fig3_infant_microbiome_maturation.R` | Redraw infant metagenomic maturation and Bifidobacterium panels. |
| Fig. 4 | `scripts/figures/fig4_scfa_maturation/redraw_fig4_scfa_maturation.R` | Redraw infant fecal SCFA maturation figure. |
| Fig. 5 | `scripts/figures/fig5_milk_scfa_lagged_candidates/redraw_fig5_milk_scfa_lagged_candidates.R` | Redraw lagged HMO/LTF-to-SCFA candidate association figure. |
| Fig. S1 | `scripts/figures/fig6_maternal_fecal_context/redraw_fig6_maternal_fecal_context.R` | Redraw maternal fecal metagenomic context figure. |
| Fig. S2 | `scripts/figures/fig7_clinical_modifier_layer/redraw_fig7_clinical_modifier_layer.R` | Redraw exploratory clinical modifier audit figure. |

## Tables and manuscript delivery

| Script | Purpose |
| --- | --- |
| `scripts/tables/build_manuscript_tables.mjs` | Build manuscript table package and figure/source-data index workbook. |
| `scripts/tables/build_gut_microbes_supplementary_tables.mjs` | Build supplementary tables workbook from reviewed source CSVs. |
| `scripts/manuscript_delivery/build_all_figures_pdf.py` | Assemble all final figures into a single PDF. |
| `scripts/manuscript_delivery/build_gut_microbes_manuscript_pdf.py` | Build manuscript review PDF with figures. |
| `scripts/manuscript_delivery/build_gut_microbes_manuscript_docx.py` | Build manuscript Word document with high-resolution figures. |
| `scripts/manuscript_delivery/build_gut_microbes_expert_review_pdf.py` | Build expert-review brief PDF. |
| `scripts/manuscript_delivery/build_expert_review_cn_pdf.py` | Build Chinese expert-review README PDF. |
| `scripts/manuscript_delivery/build_leadership_manuscript_pdf.py` | Build internal leadership-facing manuscript PDF. |
