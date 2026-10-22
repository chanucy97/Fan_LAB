"""Only synthetic reads/handles; no real BAM, covT, Rdic or phenotype IO."""
import copy
from dataclasses import replace
from types import SimpleNamespace as NS
import unittest
import core


LEFT = "ACGTACGTACGTACGTACGT"
RIGHT = "TGCATGCATGCATGCATGCA"
D = LEFT+RIGHT
R = LEFT+"ACGC"+RIGHT
N1 = "T"+R[1:]
N2 = R[:-1]+"T"
CONTEXTS = (core.SAMPLE_D, core.SAMPLE_R, core.REF_PREFIX+"r1", core.REF_PREFIX+"r2")
TEMPLATES = tuple(zip(CONTEXTS, (D,R,N1,N2)))
ASSAY = core.Assay("synthetic", "S", "r1", 100,144,((120,124),),TEMPLATES,20,30,10)
RULES = core.TechnicalRules(10,3,2,0.9,10)


def read(qname="q", seq=D, cigar="20M4D20M", **kwargs):
    values = dict(sample_id="S",qname=qname,reference="r1",start0=100,cigar=cigar,
                  sequence=seq,qualities=(30,)*len(seq),flag=67,mapq=10)
    values.update(kwargs)
    return core.Read(**values)


def call(r=None, assay=ASSAY, rdic=None):
    r = r or read()
    return core.classify_read(r,assay,{r.qname} if rdic is None else rdic)


class EvidenceTests(unittest.TestCase):
    def test_four_contexts_and_cigar(self):
        for label,seq in TEMPLATES:
            c=call(read(seq=seq,cigar="20M4D20M" if label==core.SAMPLE_D else "44M"))
            self.assertEqual(c["allele"],label)
            self.assertTrue(c["main_record_eligible"])

    def test_raw_quality_boundary_and_missing(self):
        self.assertEqual(call()["status"],"ALLELE")
        self.assertEqual(call(read(qualities=(29,)+(30,)*39))["status"],"ORIGINAL_BQ_BELOW_THRESHOLD")
        self.assertEqual(call(read(qualities=None))["status"],"ORIGINAL_QUALITIES_UNAVAILABLE")

    def test_retained_middle_base_quality_required(self):
        q=[30]*44;q[21]=29
        self.assertEqual(call(read(seq=R,cigar="44M",qualities=tuple(q)))["status"],"ORIGINAL_BQ_BELOW_THRESHOLD")

    def test_short_anchor_and_softclip(self):
        self.assertEqual(call(read(seq=D[1:],cigar="19M4D20M",start0=101))["status"],"NOT_SPANNING_BOTH_ANCHORS")
        self.assertEqual(call(read(seq="AAA"+D,cigar="3S20M4D20M"))["allele"],core.SAMPLE_D)

    def test_sequence_and_cigar_agree(self):
        self.assertEqual(call(read(cigar="19M4D21M"))["status"],"NOT_SPANNING_BOTH_ANCHORS")
        self.assertEqual(call(read(cigar="20M4N20M"))["status"],"REFERENCE_SKIP_UNRESOLVED")
        self.assertEqual(call(read(seq=LEFT+"T"+RIGHT,cigar="20M4D1I20M"))["status"],"OTHER_CIGAR_UNRESOLVED")
        self.assertEqual(call(read(seq="G"+D[1:]))["status"],"ANCHOR_SEQUENCE_MISMATCH")

    def test_original_flags_and_mapq(self):
        for mq in (9,255):
            self.assertFalse(call(read(mapq=mq))["main_record_eligible"])
        for flag in (65,67|256,67|512,67|2048):
            self.assertFalse(call(read(flag=flag))["main_record_eligible"])
        self.assertTrue(call(read(flag=67|1024))["main_record_eligible"])

    def test_reverse_direction_not_reverse_complemented_again(self):
        c=call(read(flag=83))
        self.assertEqual((c["allele"],c["direction"]),(core.SAMPLE_D,"reverse"))
        c=call(read(flag=83),replace(ASSAY,reference_to_candidate_orientation=-1))
        self.assertEqual(c["direction"],"forward")

    def test_mates_conflict(self):
        c=[call(),call(read(seq=R,cigar="44M",flag=131))]
        f=core.collapse_qnames(c,"main")
        self.assertEqual(f[0]["state"],"CONFLICT_EXCLUDED")

    def test_two_native_contexts_conflict(self):
        c=[call(read(seq=N1,cigar="44M")),call(read(seq=N2,cigar="44M",flag=131))]
        self.assertEqual(core.collapse_qnames(c,"main")[0]["state"],"CONFLICT_EXCLUDED")

    def test_unmatched_high_quality_retained_mate_conflicts(self):
        c=[call(),call(read(seq="G"+R[1:],cigar="44M",flag=131))]
        self.assertEqual(c[1]["status"],"ANCHOR_SEQUENCE_MISMATCH")
        self.assertEqual(c[1]["structural_class"],"RETAINED")
        self.assertEqual(core.collapse_qnames(c,"main")[0]["state"],"CONFLICT_EXCLUDED")
        c[1]=call(read(seq="G"+R[1:],cigar="44M",flag=131,qualities=(29,)*44))
        self.assertEqual(core.collapse_qnames(c,"main")[0]["state"],"ADJUDICABLE")

    def test_homologous_reference_starts_not_double_counted(self):
        a2=replace(ASSAY,reference="r2",start0=200,end0=244,gap_intervals0=((220,224),))
        c=[call(read(qname="one")),call(read(qname="two",reference="r2",start0=200),a2)]
        s=core.summarize_fragments(core.collapse_qnames(c,"main"),RULES,CONTEXTS)
        self.assertEqual(s["deletion_distinct_selected_read_starts"],1)
        self.assertEqual(s["deletion_distinct_reference_start_coordinates_descriptive"],2)

    def test_qualifying_mate1_preferred(self):
        c=[call(read(flag=147)),call(read(flag=67))]
        f=core.collapse_qnames(c,"main")[0]
        self.assertEqual(f["representative"]["mate"],1)
        self.assertEqual(f["qualifying_mates"],[1,2])
        c[1]=call(read(qualities=(29,)*40))
        self.assertEqual(core.collapse_qnames(c,"main")[0]["representative"]["mate"],2)

    def test_double_reference_qname_and_duplicate_fetch(self):
        a2=replace(ASSAY,reference="r2")
        c1=call();c2=call(read(reference="r2"),a2)
        f=core.collapse_qnames([c1,c1,c2],"main")
        self.assertEqual(len(f),1)
        self.assertEqual(f[0]["record_count_after_duplicate_fetch_removal"],2)
        self.assertTrue(f[0]["multiple_reference_mappings"])
        c3=copy.deepcopy(c1);c3["status"]="OTHER"
        with self.assertRaises(ValueError):core.collapse_qnames([c1,c3],"main")

    def test_rdic_diagnostics_never_enter_main(self):
        c=call(rdic=set())
        self.assertEqual(core.collapse_qnames([c],"main"),[])
        self.assertEqual(len(core.collapse_qnames([c],"diagnostic")),1)
        with self.assertRaises(ValueError): core.summarize_fragments(core.collapse_qnames([c],"diagnostic"),RULES,CONTEXTS)

    def test_ambiguous_context_excluded(self):
        a=replace(ASSAY,templates=TEMPLATES[:-1]+((CONTEXTS[-1],N1),))
        with self.assertRaises(ValueError): call(read(seq=N1,cigar="44M"),a)

    def test_identity_mismatch_rejected(self):
        with self.assertRaises(ValueError):call(read(sample_id="other"))
        with self.assertRaises(ValueError):call(read(reference="other"))
        with self.assertRaises(ValueError):call(read(cigar="40M",seq=D+"A"))

    def test_support_gates_and_context_denominators(self):
        calls=[]
        for i in range(10):
            pad=i%3
            calls.append(call(read(qname=str(i),start0=100-pad,seq="A"*pad+D,cigar=f"{20+pad}M4D20M",flag=83 if i<2 else 67)))
        calls += [call(read(qname="R",seq=R,cigar="44M"))]
        calls += [call(read(qname="N"+str(i),seq=N1,cigar="44M")) for i in range(10)]
        s=core.summarize_fragments(core.collapse_qnames(calls,"main"),RULES,CONTEXTS)
        self.assertEqual(s["technical_label"],"READ_SUPPORTED_ASSEMBLY_ALLELE_IN_SAMPLE_FLANK_CONTEXT")
        self.assertEqual(s["coexistence_label"],"COEXISTING_REFERENCE_LIKE_CONTEXT")
        self.assertEqual(s["sample_flank_deletion_fraction"],10/11)
        self.assertEqual(s["all_context_deletion_fraction_descriptive"],10/21)
        self.assertEqual(s["qnames_by_allele_context"][CONTEXTS[-1]],0)
        self.assertEqual(s["deletion_distinct_selected_read_starts"],3)
        s=core.summarize_fragments(core.collapse_qnames(calls[1:],"main"),RULES,CONTEXTS)
        self.assertEqual(s["technical_label"],"TECHNICAL_SUPPORT_CRITERIA_NOT_MET")

    def test_empty_has_four_explicit_zeros(self):
        s=core.summarize_fragments([],RULES,CONTEXTS)
        self.assertEqual(s["qnames_by_allele_context"],dict.fromkeys(CONTEXTS,0))
        self.assertIsNone(s["sample_flank_deletion_fraction"])


class CoverageTests(unittest.TestCase):
    def test_original_post_overlap_counter(self):
        def p(base,qname="accepted",delete=False,skip=False):
            return NS(is_del=delete,is_refskip=skip,alignment=NS(query_name=qname,query_sequence=base),query_position=0)
        col=NS(pileups=[p("A"),p("C"),p("T"),p("G"),p("N"),p("A","notRdic"),p("A",delete=True),p("A",skip=True)])
        self.assertEqual(core.original_column_counts(col,{"accepted"}),dict.fromkeys("ACTG",1))

    def test_exact_original_pileup_call(self):
        class Bam:
            def pileup(self,*args,**kw):self.args=args;self.kw=kw;return [NS(pos=4,pileups=[])]
        b=Bam();self.assertEqual(list(core.original_split_coverage(b,"r",0,9,set())),[(4,0)])
        self.assertEqual(b.args,("r",))
        self.assertEqual(b.kw,dict(truncate=True,max_depth=100000,stepper="nofilter",compute_baq=True,ignore_orphans=True,ignore_overlaps=True,min_base_quality=30,start=0,stop=10))

    def test_parity_sparse_zero_low_depth_and_unavailable(self):
        self.assertEqual(core.compare_covt_depth([(1,1),(2,4)],{1:1,2:4},0,4)["status"],"MATCH")
        self.assertEqual(core.compare_covt_depth([],{1:1},0,4)["status"],"MISMATCH_BLOCK_EVENT_INTERPRETATION")
        self.assertIsNone(core.compare_covt_depth([],None,0,4)["parity"])
        with self.assertRaises(ValueError):core.compare_covt_depth([(1,1),(1,1)],{},0,4)

    def test_actual_split_formula(self):
        self.assertEqual(list(core.original_splits(25000,10000)),[(0,8332),(8333,16665),(16666,24999)])
        self.assertEqual(list(core.original_splits(10000,10000)),[(0,4999),(5000,9999)])


if __name__=="__main__": unittest.main()
