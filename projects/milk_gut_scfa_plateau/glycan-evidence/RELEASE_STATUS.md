# Public source release: 22 October 2026

Version: **0.1.0.dev10**. This publication updates the dev7 source published on 12 October 2026 and retains the original scientific implementation, preset values and frozen identity rules.

## Included

The public core includes marker/locus evidence, recruitment and local assembly, anchor-interval readback, base and coding-sequence checks, fragment-compatible structure inference, joint assignment, assembly-assisted recovery and internal-region conditioning. Source code, resource JSON files, artificial test fixtures and command documentation are included. Companion [comparison](../methods/structure-comparison/) and [sequence/population](../methods/sequence-and-population/) modules document their input contracts and portability changes.

The source implementation was compared byte-for-byte with the accepted study dev10 wheel. The original study wheel includes a machine-specific provenance file and is not redistributed. A public wheel is rebuilt from the screened source; its hash differs from that historical artifact. Public contents and validation are recorded in the project release manifests.

## Interpretation

The development comparison reported 23/35 correct structural classes and zero additional class calls for the primary assembly-assisted rule, matching ordinary-assembly geometry. The software makes structure, read support, reference discordance and sequence checks explicit. These evidence records do not themselves establish biological function or universal diagnostic accuracy.

## Reproducing a study run

Supply reads, references, saved alignments and bound manifests according to each command's schema. Several research-bundle commands need additional artifacts; installing the package alone does not recreate all manuscript figures. Companion public-adapter files record original-source hashes and portability changes. Cohort-specific execution wrappers, participant tables and assemblies follow the study data-access route.
