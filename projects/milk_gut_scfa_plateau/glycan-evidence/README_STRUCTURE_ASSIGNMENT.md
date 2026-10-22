# Joint local-structure assignment — experimental dev9

`structure-assign` analyzes paired-fragment evidence against an explicit local candidate panel. It combines compatibility across the interval with direct and paired evidence at both internal boundaries. Expected profiles must have been generated with the same extraction and mapping contract as the sample.

```bash
glycan-evidence structure-assign --sam query_name_grouped.sam \
  --candidates candidates.json --profiles profiles.json \
  --config cli_config.json --out new_assignment_output
glycan-evidence verify --directory new_assignment_output
```

Install the `inference` extra. This command needs Python 3.12 or later with the pinned numerical runtime. It neither runs the mapper nor downloads references. `cli_config.json` has `extraction` settings and `fit` settings, as documented for `structure-infer`. The fit settings are `tv_slack`, `numerical_tolerance`, and optional `candidate_order`/`candidate_groups`.

Each accepted QNAME contributes one observed feature. The joint model allocates that feature's mass across candidates allowed by its observed alignments or generating profile. The allocated mass must agree with the candidate's expected feature frequencies within an explicitly recorded total-variation budget. Group weights receive minimum and maximum bounds. Boundary support receives a minimum and a best-fit point; it does not receive a full interval. Profile-only allocation edges never contribute observed boundary evidence.

`fit.json` contains:

- `joint_bestfit_tv` and `marginal_tv_at_joint_bestfit`, the latter evaluated at the joint optimum rather than separately optimized;
- `group_ranges` and `boundary_bounds`, including unrounded fractional-fragment lower bounds;
- `observed_support`, with distinct actual informative, possible and exclusive QNAME counts;
- `same_mask_no_profile` and `observed_mask_no_profile` allocation-only comparisons;
- unexplained fragments, profile-only source edges, numerical residual checks and explicit no-evidence output.

The public expert-input CLI records input/output hashes and refuses existing output directories. It does not establish that arbitrary input profiles share the sample's generation contract. The dated study runner additionally binds the candidate, profile, mapping and extraction receipts. Do not substitute a different profile dictionary while retaining a study's calibration claim.

The study uses a fixed slack of 0.02; it produces descriptive bounds, not confidence intervals. Candidate proportions are conditional local signal fractions, not organism abundance. Separate left and right support does not physically phase an entire interval or establish its carrier genome. An unsupported call remains unresolved. The same-mask ablation removes the frequency constraints but retains profile-permitted compatibility edges.

The 2026-10-15/16 comparison preserves v1 and records this v2 revision after examining v1 failures. All 32 simulation conditions are exposed internal development data; the 50 cohort libraries have no architecture ground truth. Main comparison includes expanded reference with possible shared-boundary support, an exclusive-boundary baseline and ordinary MEGAHIT on the same broadly recruited reads. Implementation and installation checks are separate from evidence of scientific improvement.
