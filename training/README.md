# geneid-train

Pure-Python training pipeline for the [geneid](../README.md) gene predictor,
with first-class U2/U12 intron support. Replaces the legacy Perl/AWK/C training
scripts. See [DESIGN.md](DESIGN.md) for the full plan and the frozen `.param`
format spec.

Status: **Phase 1** — package scaffold, CLI surface, and the `.param`
read/write core (byte-identical round-trip).

## Dev cycle

```bash
cd training
make install     # editable install + dev deps (pytest, ruff) into your env
make test        # fast unit suite
make test-all    # also round-trips every real param/*.param in the repo
make lint        # ruff check
make fmt         # ruff format
```

## Try it

```bash
geneid-train param-info ../param/human1iso.param
geneid-train param-info ../param/drosophila.U12.070102.param   # U12-aware: yes
```

The pipeline subcommands (`prepare`, `train`, `evaluate`, `jackknife`) are
stubs at this phase and report which phase implements them.
