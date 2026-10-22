"""Describe position-level QNAME ambiguity without editing a sequence."""

BASES = 'ACGT'


def assess_position(row):
    """Consume the frozen base-concordance table's genomic-oriented counts.

    Each QNAME has one A/C/G/T vote after consistent-mate collapse. The
    orientation groups F/R/both partition these votes; both contributes to
    each orientation check but is counted once in the total denominator.
    No flag is an independently validated variant or a sequence correction.
    """
    counts, strands = {}, {}
    for base in BASES:
        values = {}
        for key in ('F', 'R', 'both', 'total'):
            value = row.get(base + '_qname_' + key)
            if type(value) is not int or value < 0:
                raise ValueError('Nonnegative integer QNAME counts required')
            values[key] = value
        if values['total'] != values['F'] + values['R'] + values['both']:
            raise ValueError('QNAME orientation groups do not sum to total')
        counts[base], strands[base] = values['total'], values
    n = sum(counts.values())
    if type(row.get('acgt_qname_total')) is not int or row['acgt_qname_total'] != n:
        raise ValueError('ACGT denominator differs from QNAME totals')
    assembly = row.get('assembled_genomic_base')
    if not isinstance(assembly, str) or len(assembly) != 1 or assembly not in 'ACGTN':
        raise ValueError('Expected one genomic-oriented assembly base')
    maximum = max(counts.values())
    top = [b for b in BASES if counts[b] == maximum] if n else []
    alternatives = [b for b in BASES if b != assembly]
    assembly_n = counts.get(assembly, 0)
    unknown = assembly == 'N'
    exceeds = not unknown and any(counts[b] > assembly_n for b in alternatives)
    tied = bool(n and assembly in top and len(top) > 1)
    strong = []
    for base in alternatives:
        s = strands[base]
        if (not unknown and n and counts[base] >= 5 and counts[base] / n >= .8 and
                s['F'] + s['both'] >= 2 and s['R'] + s['both'] >= 2):
            strong.append(base)
    return dict(acgt_qnames=n, top_bases=top, assembly_base_unknown=unknown,
                assembly_fraction=assembly_n/n if n and not unknown else None,
                no_acgt_evidence=n == 0, fewer_than_five_qnames=n < 5,
                nonassembly_exceeds_assembly=exceeds, assembly_tied_for_top=tied,
                strong_nonassembly_bases=strong,
                strong_nonassembly_discordance=bool(strong),
                scope='base_evidence_flags_only; no_indel_audit_or_sequence_correction')
