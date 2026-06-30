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

The `morc_u12` and `rnaseq` cases mirror the documented MORC3 workflow from
the Current Protocols geneid chapter (Alioto et al.). They restrict
processing to the locus with `-j 36315000 -k 36380000`, which keeps each run
~1s while preserving native chr21 coordinates so the evidence GFFs align.

## Fixtures

`chr21.fa` is large and is **not** tracked; `run.sh` derives it on demand
from the tracked `samples/GRCh38.chr21.fa.gz` (the same `gunzip` step the
chapter describes). All other inputs (params, example1.fa, the snake
fragment + its param, the `longprot` 155 kb window cut from SUPER_1
(~199.14-199.30 Mb), and the MORC evidence GFFs) are tracked.

## Guarantee

At each refactoring step the goldens must stay byte-identical (modulo the
`# date` line) and the build must stay AddressSanitizer-clean. All four
cases are deterministic across runs once the date line is stripped.
