# geneid regression suite

A behavior-locking guardrail for refactoring work. Each case runs `geneid`
with a fixed parameter file, sequence, and flags, then compares the GFF
output (minus the volatile `# date` header) against a committed golden in
`golden/`. Any change that alters a covered code path shows up as a `FAIL`.

## Usage

```sh
test/regression/run.sh             # build (MEM=medium) and verify every case
test/regression/run.sh --no-build  # verify using the existing bin/geneid
test/regression/run.sh --bless     # (re)generate goldens from current output
test/regression/run.sh snake human # run only the named case(s)
```

Exit status is non-zero if any case fails, so it can gate a commit/PR.

## Cases

| Case | Param | Sequence | What it exercises |
|------|-------|----------|-------------------|
| `snake` | Hemorrhois (1 isochore) | rHemHip 1Mb fragment | Full output paths (`-3UDTA`): genes, proteins, cDNA, tDNA on a large realistic genome (~48 genes) |
| `longprot` | Hemorrhois (1 isochore) | rHemHip 155kb window | Very long transcript: a 153-exon / 23313-aa gene (the longest protein on SUPER_1) exercises the growable whole-protein and cDNA/tDNA buffers |
| `human` | human3iso (3 isochores) | example1.fa | Isochore selection + protein/cDNA/tDNA on the documented example |
| `morc_u12` | human3isoU12 | chr21 MORC3 locus | U12 intron prediction (`-3UnDTA`); the MORC3 gene has two real U12 introns |
| `rnaseq` | human.rnaseq | chr21 MORC3 locus | RNA-seq evidence: intron junctions (`-R`), expression coverage (`-S`), UTRs (`-u`) |
| `human_intron` | human.chr12.intron_length (weight 0.2) | Red-masked chr14 250kb slice | Soft intron-length penalty: two far-band introns (>27kb) are penalized away vs weight 0, pinning the convex-hinge penalty + far-band fast-DP path (`-3Un`) |
| `human_intron_multifrag` | human.chr12.intron_length (weight 0.2) | Red-masked chr14 550kb slice | Same penalty on a >500kb (LENGTHSi) slice → two fragments, pinning the near-band deque's cross-fragment maintenance in `BackupArrayD` |
| `human_intron_bb` | human.chr12.intron_length (weight 0.2) | Red-masked chr14 250kb slice + `-R` bigBed | bigBed `-R` evidence (per-split range query): the same two Intron junctions via bigBed span both far-band introns; exercises `bbOpen`→`bbQuery`→`AddEvidenceExon`, byte-identical to the GFF `-R` result |
| `human_o` | human3iso (3 isochores) | example1.fa + `-O` | Assemble-only (`-O`): forces example1's own 8 typed exons; guards the `-O` path's independent `nExons` accounting |
| `human_o_multilocus` | human3iso (3 isochores) | 2-record example1 + `-O` | Multi-locus `-O`: the same 8 exons under two sequences → two genes assembled (one per locus) |
| `morc_jo_u12` | human3isoU12 | whole chr21 (multi-split) + `-J -O -U` | Annotation-scoring (`-J`) + U12 typing over forced MORC3 CDS exons: guards `-U` allowed under `-O` and the multi-split `-J` classify (walks the printed `GOptim` chain → two U12 introns typed) |
| `morc_o_utr` | human.rnaseq | whole chr21 (multi-split) + `-O -u` | UTR assemble-only: forces MORC3 CDS + both UTR halves; guards `-u` under `-O` and the UTR-exon frame/remainder fix (the 3′ `UTR_Terminal_Half`, length not a multiple of 3, must not be dropped) |

The `morc_u12` and `rnaseq` cases mirror the documented MORC3 workflow from
the Current Protocols geneid chapter (Alioto et al.). They restrict
processing to the locus with `-j 36315000 -k 36380000`, which keeps each run
~1s while preserving native chr21 coordinates so the evidence GFFs align.

## Fixtures

`chr21.fa` is large and is **not** tracked; `run.sh` derives it on demand
from the tracked `samples/GRCh38.chr21.fa.gz` (the same `gunzip` step the
chapter describes). All other inputs (params, example1.fa, the snake
fragment + its param, the `longprot` 155 kb window cut from SUPER_1
(~199.14-199.30 Mb), the `human_intron` 250 kb Red-masked chr14 slice
(GRCh38 chr14:67,900,000-68,150,000) and the `human_intron_multifrag` 550 kb
slice (chr14:67,700,000-68,250,000) + their chr12-trained param, the
MORC evidence GFFs, and the `-O` annotation fixtures — `example1.geneid.gff`,
the 2-record `example1.2locus.fa` + `example1.2locus.geneid.gff`,
`MORC.CDS.geneid.gff` (17 typed CDS exons), and `MORC.UTR.geneid.gff`
(CDS + both UTR halves)) are tracked.

## Guarantee

At each refactoring step the goldens must stay byte-identical (modulo the
`# date` line) and the build must stay AddressSanitizer-clean. All
cases are deterministic across runs once the date line is stripped.
