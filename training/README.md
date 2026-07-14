# geneid-train

Pure-Python training pipeline for the [geneid](../README.md) gene predictor,
with first-class U2/U12 intron support. Replaces the legacy Perl/AWK/C training
scripts. See [DESIGN.md](DESIGN.md) for the full plan and the frozen `.param`
format spec.

Status: **working end-to-end.** From a GFF3 annotation + genome it trains a
complete geneid `.param` — splice/start site profiles, the coding Markov model,
the gene model, and optional U12 and U2 branch-point profiles — and also
`optimize`s the exon weights, `evaluate`s predictions (SN/SP), and runs
`jackknife` cross-validation. See [DESIGN.md](DESIGN.md) for details.

## Setup

`geneid-train` is **pure Python with no third-party runtime dependencies** — just
the standard library. There is nothing to compile and no wheels to build, so
setup is only "get a recent Python and install."

**Requirement: Python ≥ 3.10.** Check with `python3 --version`. On macOS the
system `python3` is usually 3.9 — install a newer one (`brew install python@3.12`,
or use pyenv) and point the commands below at it with `PYTHON=python3.12`.

### One-command setup (recommended)

From this `training/` directory:

```bash
make venv                    # if python3 is already >= 3.10
make venv PYTHON=python3.12  # otherwise, point at a 3.10+ interpreter
```

That creates an isolated `./.venv`, upgrades pip, and installs `geneid-train`
(editable) plus the dev tools. **Nothing to activate** — the other targets
auto-detect `./.venv`:

```bash
make test        # fast unit suite
make test-all    # also round-trips every real param/*.param in the repo (slower, opt-in)
make lint        # ruff check
make fmt         # ruff format
```

Run the CLI directly with `./.venv/bin/geneid-train --help`, or `source
.venv/bin/activate` first and just call `geneid-train`.

### Installing into your own environment

Already have a virtualenv/conda env on Python ≥ 3.10? Skip `make venv` and install
into the active environment instead:

```bash
make install     # pip install --upgrade pip && pip install -e ".[dev]"
```

If you try to install under Python < 3.10, pip stops with a clear
`requires a different Python` error rather than failing mysteriously later.

## Try it

```bash
geneid-train param-info ../param/human1iso.param
geneid-train param-info ../param/drosophila.U12.070102.param   # U12-aware: yes
```

Build a validated training set from a GFF3 annotation + genomic FASTA:

```bash
# other flavors go through the converter first (GFF3 is the only ingested format)
geneid-train convert --from gff2 annotation.gff2 annotation.gff3
geneid-train prepare --gff annotation.gff3 --fastas genome.fa --min-aa 100 --results out
```

`convert` accepts `--from gff2` or `--from gtf` (GFF3 is the only ingested format).

Train a parameter file, then tune and evaluate it:

```bash
# GFF3 annotation + genome -> a complete geneid .param
# add --utr for a UTR-aware gene model (geneid -u with -S/-Y RNA-seq coverage)
geneid-train train --gff annotation.gff3 --fastas genome.fa \
    --species Genus_species --output Genus_species.param [--utr] [--u12] [--u2-branch]

# grid/search the exon weights against a held-out set
geneid-train optimize --param Genus_species.param \
    --eval-fastas eval.fa --eval-gff eval.gff3 --output Genus_species.optimized.param

# SN/SP of predictions vs annotation, and leave-group-out cross-validation
geneid-train evaluate predictions.gff annotation.gff
geneid-train jackknife --gff annotation.gff3 --fastas genome.fa \
    --species Genus_species --eval-fastas eval.fa --eval-gff eval.gff3
```

Retrofit an existing (CDS-only) param with newer capabilities, non-destructively:

```bash
# UTR-aware gene model + human soft intron-length model (weight 0 = inert) + U12
geneid-train retrofit Genus_species.param -o Genus_species.rnaseq.param \
    --utr --intron-length human --u12
```

`--utr` reuses the param's own trained intron range and sets the intergenic
minimum to 0 (neighbouring UTRs may abut/overlap). `--intron-length` takes
`human` or an explicit `mu,sigma`; the weight defaults to 0 (no prediction
change) — set `--intron-length-weight` or retrain per species to enable the
penalty. `--u12` injects the bundled pan-taxon U12 profile trio + acceptance
gates so `geneid -U` predicts U12 introns (the profiles carry their own geneid
geometry, so they drop into any param). Each injection is skipped if the param
already has it. (Multi-isochore params are a planned follow-up.)

`optimize`, `evaluate`, and `jackknife` need a compiled `geneid` binary on your
`PATH` (or pass `--geneid /path/to/geneid`); build it from the repo root with
`make`.
