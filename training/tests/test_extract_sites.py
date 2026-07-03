import random
import re

import pytest

from geneid_train.core.seq import revcomp
from geneid_train.prepare.base import Exon, GeneModel
from geneid_train.prepare.sites import (
    acceptor_windows,
    collect_sites,
    donor_windows,
    is_canonical_acceptor,
    is_canonical_donor,
    is_canonical_start,
    start_window,
)

from .conftest import ref_dir


def _synthetic_chrom():
    """A 200 bp chromosome with a 2-exon + gene: exons [41,70] and [121,150],
    intron [71,120] (GT..AG), start codon ATG at 41-43."""
    random.seed(1)
    c = [random.choice("ACGT") for _ in range(200)]

    def put(pos1, s):
        for i, ch in enumerate(s):
            c[pos1 - 1 + i] = ch

    put(41, "ATG")  # start codon
    put(71, "GT")  # donor: first intron bases
    put(119, "AG")  # acceptor: last intron bases
    return "".join(c)


def test_plus_strand_geometry():
    chrom = _synthetic_chrom()
    m = GeneModel("g+", "chr", "+", [Exon(41, 70), Exon(121, 150)])
    d = donor_windows(m, chrom)
    a = acceptor_windows(m, chrom)
    s = start_window(m, chrom)
    assert len(d) == 1 and len(a) == 1
    assert len(d[0]) == len(a[0]) == len(s) == 60
    # anchor at profile position 31; invariant motifs at their fixed offsets
    assert is_canonical_donor(d[0]) and d[0][31:33] == "GT"
    assert is_canonical_acceptor(a[0]) and a[0][28:30] == "AG"
    assert is_canonical_start(s) and s[30:33] == "ATG"
    # windows equal the exact genomic slices, computed independently of _window
    assert d[0] == chrom[39:99]  # donor anchor 70 (last exon base)
    assert a[0] == chrom[90:150]  # acceptor anchor 121 (first exon base)
    assert s == chrom[10:70]  # start anchor 41 (A of ATG)


def test_minus_strand_matches_plus():
    """Mirroring the genome and the gene onto the minus strand must yield the
    identical transcription-oriented windows."""
    chrom = _synthetic_chrom()
    n = len(chrom)
    mp = GeneModel("g+", "chr", "+", [Exon(41, 70), Exon(121, 150)])
    dp, ap, sp = donor_windows(mp, chrom), acceptor_windows(mp, chrom), start_window(mp, chrom)

    cm = revcomp(chrom)

    def mirror(a, b):
        return Exon(n + 1 - b, n + 1 - a)

    mm = GeneModel("g-", "chr", "-", sorted(
        [mirror(41, 70), mirror(121, 150)], key=lambda e: e.start
    ))
    dm, am, sm = donor_windows(mm, cm), acceptor_windows(mm, cm), start_window(mm, cm)
    assert dm == dp
    assert am == ap
    assert sm == sp


def test_windows_skipped_at_sequence_ends():
    chrom = _synthetic_chrom()
    # a gene whose start codon is at position 5 has no room for 30 bp upstream
    m = GeneModel("edge", "chr", "+", [Exon(5, 40)])
    assert start_window(m, chrom) is None


# ---- real-data validation: canonical fraction on Xerocrassa eval genes -------

REF = ref_dir()
GENOME = REF.parents[2] / "xgXerMont_curated.no_mt.scrubbed.fa.gz" if REF else None


def _load_geneid_gff(path, seqid_prefixes):
    """Parse a geneid-format gff (First/Internal/Terminal/Single, genome coords)
    into GeneModels, renaming SUPERn -> SUPER_n to match the genome headers."""
    tx: dict[tuple[str, str, str], list[tuple[int, int]]] = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if f[0] not in seqid_prefixes:
            continue
        seqid = re.sub(r"^SUPER", "SUPER_", f[0])
        key = (f[8].strip(), seqid, f[6])
        tx.setdefault(key, []).append((int(f[3]), int(f[4])))
    return [
        GeneModel(gene_id=t, seqid=s, strand=st, exons=[Exon(a, b) for a, b in sorted(ex)])
        for (t, s, st), ex in tx.items()
    ]


@pytest.mark.integration
@pytest.mark.skipif(
    REF is None or not (GENOME and GENOME.exists()),
    reason="set GENEID_TRAIN_REFDIR and provide the xgXerMont genome",
)
def test_real_extraction_canonical_fraction():
    from geneid_train.core.fasta import read_fasta_subset

    gff = REF / "Xerocrassa_montserratensis.eval.geneid.gff_sorted"
    models = _load_geneid_gff(gff, {"SUPER1", "SUPER2"})
    seqids = {m.seqid for m in models}
    genome = read_fasta_subset(str(GENOME), seqids)
    s = collect_sites(models, genome)

    ndon = len(s.donor) + len(s.noncanonical_donor)
    nacc = len(s.acceptor) + len(s.noncanonical_acceptor)
    assert ndon > 100 and nacc > 100  # a meaningful sample
    # real introns are overwhelmingly GT-AG; a 1-off geometry error collapses this
    assert len(s.donor) / ndon > 0.95
    assert len(s.acceptor) / nacc > 0.98
    # every annotated CDS starts with ATG under the standard code
    assert len(s.noncanonical_start) == 0 and len(s.start) == len(models)
    # geometry sanity on the extracted windows themselves
    assert all(len(w) == 60 for w in s.donor)
    assert all(w[31:33] == "GT" for w in s.donor)
