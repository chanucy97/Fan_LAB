# Input And Output Boundary

## Required controlled inputs

The scripts expect derived, controlled analysis files under the workspace root. Important input families include:

- `direction1_publication_figures/source_data/all_milk_long_values.csv`
- `direction3_multiomics_blueprint/tables/direction3_top_milk_dynamic_features.csv`
- `direction3_multiomics_blueprint/tables/direction3_lagged_milk_infant_association_screen.csv`
- `direction3_multiomics_blueprint/tables/direction3_analysis_matrix_lagged_long.csv`
- `p0_multiomics_extension/tables/p0_hmo_ltf_to_scfa_lagged_model_highlights.csv`
- `p0_multiomics_extension/tables/p0_scfa_maturation_score.csv`
- `p1_clinical_modifiers/tables/p1_scfa_maturation_clinical_analysis_matrix.csv`
- `p1_clinical_modifiers/tables/p1_scfa_maturation_by_clinical_group_summary.csv`
- `p1_clinical_modifiers/tables/p1_scfa_maturation_clinical_lme_results.csv`
- `p2_bridge_analyses/tables/p2_milk_module_clinical_matrix.csv`
- `final_manuscript_planning_20260605/figures_redrawn/**/Fig*_source_*.csv`
- `final_manuscript_planning_20260605/gut_microbes_revision/gut_microbes_manuscript_expert_review.md`

These files may contain derived but still sensitive cohort-level or participant-level information. They are not included in the public repository.

## Main generated outputs

Typical outputs are written under:

- `final_manuscript_planning_20260605/figures_redrawn/`
- `final_manuscript_planning_20260605/compiled_figures/`
- `final_manuscript_planning_20260605/manuscript_tables/`
- `final_manuscript_planning_20260605/gut_microbes_revision/`

Generated outputs include figure PNG/PDF/SVG/TIFF files, source-data CSV exports, leave-one-dyad-out sensitivity tables, supplementary workbook files, and manuscript Word/PDF deliverables. These generated outputs are not uploaded here because they may contain sensitive source data, manuscript material, or large binary artifacts.
