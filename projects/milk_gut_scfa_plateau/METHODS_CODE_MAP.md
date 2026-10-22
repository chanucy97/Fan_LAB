# Manuscript methods to public code

Public snapshot: 22 October 2026; GlycanEvidence **0.1.0.dev10**.

| Method | Public entry point | Required input | Main output |
|---|---|---|---|
| Separate identity and glycan-locus evidence | [Core](glycan-evidence/), `run`, `call` | Paired FASTQ, verified marker panel, frozen rules | Marker/locus metrics and identity evidence |
| Local recruitment and anchor architecture | `assemble`, `structure` | Explicit recruitment reference, anchors, reads or contigs | Recruited assembly, interval geometry and readback |
| Fragment-compatible candidate structures | `structure-infer`, `structure-assign` | Candidate metadata, alignments, bound candidate profiles | Compatibility sets, allocations and ambiguity |
| Ordinary-assembly recovery | `structure-recover` | Bound, truth-free roster of reference calls and assembly evidence | Supported intervals, classes and reference discordance |
| Expanded-reference and assembly comparisons | [Comparison components](methods/structure-comparison/) | Fixed truth/condition ledger, prediction files, assemblies and saved history | Matched-condition accuracy and geometry summaries |
| CDS and local sequence consistency | `cds-replay`, `bases` | Bound sequences/HSPs or structure/readback receipts | Coding comparisons and per-position read support |
| Population-background comparisons | [Sequence/population components](methods/sequence-and-population/) | Existing inStrain profiles, reference configuration and comparison inputs | Coverage-gated shared-background summaries |
| Nutritional-enzyme sequence context | [Sequence/population components](methods/sequence-and-population/) | User-supplied assemblies, protein/CDS references and alignment results | Homologue/context and sequence-integrity records |
| Prespecified sequence-event readback | [Sequence/population components](methods/sequence-and-population/) | One explicitly selected candidate, original BAMs and bound configuration | Event/retained/unresolved fragment counts and parity checks |

## Version boundaries

Dev7 was published on 12 October 2026. Dev8 introduced fragment-compatible inference; dev9 added joint allocation; dev10 adds evidence from ordinary assemblies. The current publication preserves these accepted implementations. Internal-region conditioning remains a Python API. Historical research-bundle commands require the additional input bundles identified by their contracts.

The companion directories identify portable public adapters separately from the scripts executed in the original study. Each adapter records its original workspace-relative source and SHA-256, public SHA-256 and the portability change. The public comparison package does not include every cohort-specific roster-building wrapper. The public sequence package does not include participant inventories, individual assemblies or the selected study candidate's sequence configuration.

## Comparison interpretation

Use identical read inputs, the same condition ledger and the same geometry/read-support thresholds when comparing reference panels and assembly. Assembly-assisted evidence retains what assembly already recovered; it is not counted as an improvement over ordinary assembly. Previously inspected test cases remain exposed development data for later method changes. Report the evaluated source/condition denominators separately from cohort observations.

## Dependencies and data

Install the core according to [its README](glycan-evidence/README.md), then follow each companion directory's dependencies and test commands. Raw-read tools, inStrain profiles, reference databases and controlled cohort inputs are external inputs. [Data/code availability](docs/DATA_CODE_AVAILABILITY.md) describes access; the source release itself contains no participant-level dataset.
