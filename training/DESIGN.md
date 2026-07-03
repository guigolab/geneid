# geneid-train — design

A ground-up rewrite of geneid's training pipeline as a single pure-Python tool
that ships inside the geneid repository and is driven by one CLI (`geneid-train`).
It replaces the ~4,450-line Perl driver (`geneidTRAINer4docker.pl`), its ~12 AWK
scripts, the Perl modules (`Param.pm`, `Isocore.pm`, `GeneModel.pm`,
`geneidCEGMA.pm`, `SeqOp.pm`), and the compiled helpers `SSgff`, `pictogram`,
`cds2gff`, `gff2ps`, and `evaluation`.

The only compiled dependency that remains is the `geneid` predictor itself.

---

## 1. Goals & non-goals

**Goals**

- One installable Python package, one CLI, zero interactive STDIN prompts.
- Configuration through a single YAML file; every run fully reproducible (fixed
  seeds, recorded config).
- First-class **U12 (minor-spliceosome) intron** training. This is the original
  reason these scripts exist — U12 donor/acceptor/branch-point profiles and the
  U2-vs-U12 split are core, not an optional fungal add-on.
- Two training-set front-ends over a shared core:
  - **BUSCO** complete single-copy genes.
  - **RNA-seq** transcripts via PASA/TransDecoder, cross-checked against UniProt
    (DIAMOND, ≥90% coverage).
- Reimplement the SN/SP scorer (`evaluation`) in Python with richer, U2/U12-aware
  output.

**Non-goals**

- No reimplementation of the `geneid` predictor. It stays compiled and is invoked
  as the optimization/eval engine.
- **No sequence design of any kind.** Training is statistical parameter estimation
  from existing annotations: it reads gene models and emits a parameter file of
  matrices and numbers. It never edits, designs, or synthesizes sequence.

## 2. Decisions (Tyler, 2026-07-02)

| Decision | Choice |
|---|---|
| Compiled tools kept | `geneid` only; `evaluation` reimplemented in Python |
| Location / invocation | `training/` package inside the geneid repo, `geneid-train` CLI |
| Input front-ends | BUSCO **and** RNA-seq/TransDecoder, built in parallel on a shared core |
| Correctness stance | Redesign freely — legacy output is a *reference*, not a byte-repro gate |

**Defaults (Tyler, 2026-07-02):**

- **Single isochore by default.** Isochore modelling was essentially a human-only
  practice; on non-model organisms it needs more data and more training effort, so
  CNAG stopped doing it. Multi-isochore support stays as an optional/experimental
  path but is deprioritized.
- **Genome masking:** *prediction* must run on the masked genome; *training* does
  not require masking (unmasked is fine). Training candidates are already
  repeat-filtered upstream, so masking mainly affects flanks — either genome works
  for the training step.

**Annotation format: GFF3 only.** The tool ingests **GFF3 exclusively** as its
canonical internal format (the most robust flavor). CDS features are grouped into
transcripts by their `Parent` attribute. Other flavors (GFF2, GTF) are handled by
explicit converters in `core/convert.py` (`geneid-train convert --from gff2|gtf`),
which re-emit a well-formed gene→mRNA→CDS hierarchy; the pipeline itself never
parses them directly. Parsing geneid's own GFF prediction output is a separate
concern handled at the `engine.py` boundary.

**Two hard constraints survive "redesign freely":**

1. **The `.param` file format is frozen.** The compiled `geneid` binary parses it,
   so byte-level *format* compatibility is mandatory even though we are free to
   change how the *values* are computed. The format spec is section 4 below.
2. **geneid's GFF coordinate & exon-counting conventions** bound the site extractor
   and the Python SN/SP scorer. The scorer's metric definition is calibrated once
   against the C `evaluation` on a fixture set, then owned in Python.

The north-star metric is **held-out SN/SP equal-to-or-better than the current Perl
pipeline** on the test species, not identical intermediate numbers.

### Start-site scoring is weak by nature -- lean on coding potential, keep it swappable

Start-codon PWMs carry little information (confirmed empirically on xgXerMont:
order 0, the weakest profile class geneid supports) -- ATG context alone is a
poor discriminator. Two design consequences (Tyler, 2026-07-02):

- **Training data strategy:** rather than training the start profile to sharply
  separate real starts from decoys, treat every in-frame ATG in a training CDS as
  a start candidate and let the coding-potential (Markov) model carry the real
  discriminating weight; the start profile stays a weak prior on top, not the
  decision-maker. This matches how geneid's own gene assembly already combines
  site + coding scores -- it just means the *trainer* shouldn't over-invest in
  making the start PWM discriminative.
- **Architecture:** keep each site class (start/donor/acceptor/branch) behind a
  swappable scorer interface rather than hard-wiring the Markov/PWM math as the
  only option. The near-term implementation is the Markov/PWM pipeline in
  `stats/sites.py`; a future start-site scorer (e.g. a small NN) should be able to
  slot in without changing `stats/model.py`, `optimize.py`, or `evaluate.py` --
  those only need a profile object (or an external score), not the estimation
  method that produced it. Do not build this abstraction speculatively ahead of
  need; keep the interface narrow enough that swapping is possible when it
  actually comes up.

### Optional splice-class profiles (GC donors, U12)

geneid natively reads separate profiles for the full splice taxonomy —
`U2gta_Donor_profile` (bulk GT-AG), `U2gcag_Donor_profile` (GC-AG), the U12 donor/
acceptor/branch trio, `Poly_Pyrimidine_Tract_profile`, etc. (exact keywords in
§4). These are **optional**: the legacy trainer emitted only the generic
`Donor_profile`/`Acceptor_profile`, so geneid was categorically blind to GC-AG and
U12 introns.

The same optional treatment applies to the **U2 branch-point + poly-pyrimidine
tract** profiles (`Branch_point_profile`, `Poly_Pyrimidine_Tract_profile`) — these
are not a rare dinucleotide class (they're trainable from all U2 introns) but are
especially valuable for **fungal** genomes; reference/transplant source is
`param/human.070606.u2branch.ppt.param`. Unified rule for **every** optional
auxiliary profile (GC donor, U2 BP+PPT, the U12 trio(s), finer U2 gtg/gty donors):
**always train it, then include it only if it improves held-out F1 or SN/SP.**
Inclusion is decided at the optimize/evaluate stage (a train-and-test toggle), not
at `classify` — which only does the dinucleotide rare-class tally.

Design stance (Tyler, 2026-07-02): building these profiles is about **coverage,
not necessarily global accuracy**. Adding a rare-class profile sometimes *lowers*
overall SN/SP, but it lets geneid predict that class at all instead of leaving it
out categorically. Therefore the trainer:

- builds each optional class only when there are enough examples (else offers to
  **transplant** an existing cross-taxa profile — U12 params generalize well; the
  reference U12 source is `param/human3isoU12.param`, lifting the U12 profiles from
  isochore 1 as `convertParam2U12.pl` does. Note `human.070606.u2branch.ppt.param`
  is U2 branch/PPT, **not** U12);
- honours geneid's rule that a U12 subtype activates only if its **full
  donor+acceptor+branch trio** is present;
- treats each optional class as a toggle and **reports held-out SN/SP with vs.
  without it**, so the choice to include it is an informed, per-species decision
  rather than an unconditional default.

`prepare`/`classify` report donor-dinucleotide (GT/GC/AT) and U2/U12 candidate
tallies up front to drive the train-de-novo-vs-transplant-vs-omit decision.

## 3. Architecture

```
training/
  pyproject.toml
  Makefile                # install / test / lint / fmt
  DESIGN.md               # this file
  README.md
  src/geneid_train/
    cli.py                # geneid-train prepare | train | evaluate | jackknife | param-info
    config.py             # YAML run config + defaults, seed handling
    core/
      param.py            # Param model: read/write geneid .param (format authority)
      gff.py              # GFF3 read/write (the one canonical annotation format)
      convert.py          # GFF2 / GTF -> canonical GFF3 (the only on-ramp for other flavors)
      fasta.py            # FASTA + tbl I/O (replaces FastaToTbl/TblToFasta)
      seq.py              # translation, ORF/completeness checks, reverse-complement
    prepare/
      base.py             # shared filters: complete ORF, non-overlap, min-aa, flank extraction
      classify.py         # lightweight U2/U12 split for the training set (see note below)
      busco.py            # BUSCO full_table + genome -> validated gene-model IR
      transcripts.py      # PASA/TransDecoder + UniProt (DIAMOND) -> validated gene-model IR
    stats/
      sites.py            # site extraction (replaces SSgff) + PWM/order-1 matrices + boundaries
      coding.py           # order 4/5 Markov coding-vs-intron log-ratio (replaces geneidCEGMA)
      branch.py           # branch-point motif discovery (replaces MEME dependency)
      model.py            # intron/intergenic size stats, isocores, param assembly (U2+U12)
    optimize.py           # smart search over eWF/oWF/AccCtx driving a parallel geneid pool
    evaluate.py           # SN/SP scorer, U2/U12-aware (replaces C `evaluation`)
    engine.py             # subprocess wrapper around the one compiled dep: geneid
  tests/
    data/                 # small fixtures (synthetic mini.param, tiny gff/fasta)
    ...
```

Both front-ends emit the same **validated gene-model intermediate representation
(IR)**, so `prepare/base.py` holds all shared filtering and BUSCO vs. transcripts
are thin adapters.

**U2/U12 split for the trainer is deliberately lightweight** (Tyler, 2026-07-02):
the trainer does *not* need a full-fidelity intron classifier. `classify.py`
partitions a training set's introns by either (a) user-provided labels, or
(b) bootstrap-scoring each intron's donor/acceptor/branch against an existing,
cross-taxa-robust U12 parameter file (`param/drosophila.U12.070102.param` et al.),
with a coarse dinucleotide pre-filter (GT-AG / GC-AG vs AT-AC) narrowing candidates
first. The standalone high-fidelity classifier (`classifyIntrons.pl`, and a newer
version Tyler maintains) is a separate annotation-pipeline concern and is out of
scope here — pull it in only if bootstrap classification proves insufficient.

## 4. Parameter-file format (the frozen contract)

A `.param` file is a flat, section-oriented text file. Empirically:

- **Comment lines** start with `#`. **Blank lines** are allowed anywhere.
- **Section keywords** are always a bare single token on their own line matching
  `^[A-Za-z][A-Za-z0-9_]*$` with no leading whitespace. **Data lines never match
  this** (they contain spaces, digits, `+`, `:`, `.`), which is the parser's rule
  for telling headers from data.

Top-level / per-isocore sections, in order:

```
NO_SCORE                         <scalar>           # non-homology penalty
number_of_isochores              <N>
  # repeated per isochore:
  boundaries_of_isochore         <lo> <hi>          # %GC
  Absolute_cutoff_exons          <4 values>         # First Internal Terminal Single
  Coding_cutoff_oligos           <4 values>
  Site_factor                    <4 values>
  Exon_factor                    <4 values>
  HSP_factor                     <4 values>
  # --- U12 / RSS block (present only in U12-aware params) ---
  RSS_Markov_Score               <scalar>
  RSS_Donor_Score_Cutoff         <scalar>
  RSS_Acceptor_Score_Cutoff      <scalar>
  U12_Splice_Score_Threshold     <scalar>
  U12_Exon_Score_Threshold       <scalar>
  U12_Exon_weight                <scalar>
  # ---------------------------------------------------------
  Exon_weights                   <4 values>
  # profiles (PWA / order-k Markov). Header line is:
  #   length offset cutoff order [a b acc_context min_dist opt_dist pen_scale]
  # followed by rows "<position> <oligo(len order+1)> <score>"
  Start_profile
  U12_Branch_point_profile       # U12 only
  U12gtag_Acceptor_profile       # U12 only
  U12atac_Acceptor_profile       # U12 only
  Acceptor_profile
  U12gtag_Donor_profile          # U12 only
  U12atac_Donor_profile          # U12 only
  Donor_profile
  Stop_profile
  # coding model
  Markov_order                   <scalar>
  Markov_Initial_probability_matrix     # rows "<oligo> <i> <j> <logscore>"
  Markov_Transition_probability_matrix  # rows "<oligo> <i> <j> <logscore>"
maximum_number_of_donors_per_acceptor_site   <scalar>
General_Gene_Model               # GenAmic assembly rules; free-form lines, echoed verbatim
```

**U12 profile/attribute set** (from `convertParam2U12.pl`, the authority):

- attributes: `RSS_Markov_Score`, `RSS_Donor_Score_Cutoff`,
  `RSS_Acceptor_Score_Cutoff`, `U12_Splice_Score_Threshold`,
  `U12_Exon_Score_Threshold`, `U12_Exon_weight`
- profiles: `U12_Branch_point_profile`, `U12gtag_Acceptor_profile`,
  `U12atac_Acceptor_profile`, `U12gtag_Donor_profile`, `U12atac_Donor_profile`

`model.py` will emit these natively in U12 mode (`geneid-train train --u12`), so
U12 support is a training mode rather than a post-hoc conversion.
`convertParam2U12.pl` is retained only as a compatibility path for legacy params.

### Round-trip strategy (Phase 1 oracle)

`core/param.py` parses the file into ordered blocks (keyword + raw following
lines). Sections we do not yet model are echoed **verbatim**, guaranteeing a
byte-identical read→write round-trip on every real param in `param/`. Sections we
own (isocore scalars/vectors) are re-rendered from typed values; a modify→write→
re-read test proves a changed field lands while everything else is untouched. When
a `geneid` binary is available, an opt-in integration test asserts geneid accepts
the round-tripped file.

## 5. Legacy → new mapping

| Legacy | Replacement |
|---|---|
| `geneidTRAINer4docker.pl` | `cli.py` + module orchestration |
| `Param.pm`, `Isocore.pm` | `core/param.py` |
| `GeneModel.pm` | `stats/model.py` |
| `geneidCEGMA.pm` | `stats/coding.py` |
| `SeqOp.pm`, `FastaToTbl`, `TblToFasta` | `core/fasta.py`, `core/seq.py` |
| `SSgff` (C) | `stats/sites.py` |
| `Getkmatrix.awk`, `information.awk`, `submatrix*.awk`, `prepare*matrix*.awk` | `stats/sites.py` (harvest `GeneidTrainer_Darek/py_code`) |
| MEME + `parseMEME`/`extractRealBranches` | `stats/branch.py` |
| `evaluation` (C) | `evaluate.py` |
| `pictogram` (C) | `stats/*` + `logomaker` |
| `gff2ps` (C) | optional plotting, deferred |
| `classifyIntrons.pl` | `prepare/classify.py` (lightweight bootstrap only; not a full port) |
| grid search in Perl | `optimize.py` (parallel geneid pool, coordinate-descent / Bayesian) |

## 6. Phasing

1. **Scaffold + I/O + param format** — package, CLI, `core/param.py` with
   byte-identical round-trip on real params. *(current phase)*
2. **Shared prepare core + both adapters** — `base.py` filters, `classify.py`,
   BUSCO and TransDecoder/UniProt adapters → IR.
3. **Site + coding models** — `sites.py`, `coding.py`, `branch.py` (U2 & U12
   profile sets built in parallel).
4. **Param assembly** — full U2/U12 `.param` that `geneid` runs without error.
5. **Python `evaluate.py`** — calibrated once against C `evaluation`; U2/U12 split.
6. **Optimization** — parallel geneid pool, smarter-than-grid search; held-out SN/SP.
7. **Jackknife, plots, docs, tutorials** — Docker optional, never required.

## 7. Test & dev cycle

- `pytest` unit tests per module; fast synthetic fixtures under `tests/data/`.
- Param round-trip runs (opt-in) against the real `param/*.param` files in the repo.
- Bring-up fixtures may need a set with **labelled U12 introns**; alternatively the
  existing, cross-taxa-robust U12 parameters (`param/drosophila.U12.070102.param`,
  `param/human.070606.u2branch.ppt.param`) can **bootstrap** classification and
  serve as regression references.
- `make install` (editable + dev deps) / `make test` / `make lint` / `make fmt`
  (ruff). No network required to run the unit suite.

## 8. U2 branch-point discovery (PLANNED — not yet implemented)

**Status.** Design only. The U12 branch/donor/acceptor trio already ships (bundled,
IAOD-derived; `param/u12.py`) and the U12-vs-U2 donor screen is in `classify.py`.
This section is the plan for the *U2* branch-point (+ PPT) profiles, which are
trained **per genome** from that genome's own introns.

**Motivation.** The U2 `Branch_point_profile` (+ `Poly_Pyrimidine_Tract_profile`)
was historically discovered by hand with MEME — find the motif, search it, train a
PWM. The goal is to automate that as a self-supervised step that runs by default
and is **included only when it improves held-out accuracy** (§2 "Optional
splice-class profiles" toggle). Strongly-conserved-branch genomes (fungi) are the
clearest win; but since our genome portfolio is mostly non-fungal, the design
targets the **generic, degenerate-branch case** first (see engine choice below).

**What NOT to do.** Do *not* port BPP (Zhang et al. 2017, `github.com/zhqingit/BPP`)
as a tool. BPP is tuned to the *hard, degenerate human* branch problem and is
human-specific in ways that would actively mislead elsewhere: fixed offset windows
(BPS 21–34, background 187–200, PPT 3–16 nt upstream of the 3′SS) that **do not fit
short fungal introns**; a motif seeded from human U2 snRNP and trained on ~223k
human introns; and a PPT model (half the BPP score) explicitly co-evolved with
human U2AF65. Take only the *idea* — an EM mixture to discover a branch motif
without labels — and train everything on the target genome.

**Pipeline (`stats/branch.py`, new).**

1. **Adaptive search window.** Per U2 intron, scan the acceptor-upstream region
   defined *relative to intron length* with AG-exclusion (from the 3′SS back to the
   first upstream AG beyond a small gap), capped by the intron's own length so short
   introns still fit. No hardcoded human offsets.
2. **EM mixture, self-trained.** Two components — a branch-motif PWM (~7 nt,
   branch-A anchored) vs an intronic background — fit by EM over the candidate
   windows. **Seed** the motif from the *universal* U2-snRNA-complementary consensus
   (branch-A pairing → `yUNAy`/TACTAAC-like), **not** human weights; **train** only
   on the target genome so the result is species-appropriate. EM converges sharply
   for conserved (fungal) branches and still yields a usable model when degenerate.
   Chosen over a consensus/Hamming-to-U2-complement shortcut (Kupfer 2004) because
   the generic degenerate case covers most of our genomes; the conserved case is the
   easy sub-case EM also handles.
3. **Branch location + distances from data.** EM assigns each intron's branch-A;
   the branch header's `acc_context` / `min_dist` / `opt_dist` are set from
   *percentiles of the discovered branch-A→3′SS distances* (short for fungi, longer
   elsewhere) — not searched. This supersedes tuning branch *distance* in the
   optimiser; `optimize --tune-branch` is demoted to optional `pen_scale`-only
   refinement (`pen_scale` is not directly observable).
4. **PPT — deferred.** Because the PPT is species-specific and weak/inessential in
   fungi, model it later as a separate, independently-gated stage
   (`Poly_Pyrimidine_Tract_profile`); v1 is branch-PWM only.
5. **Emit + gate.** Reuse `param.u12`'s acceptor-side insertion for
   `Branch_point_profile`; assemble with vs without it and keep it only if held-out
   exon SNSP improves (`evaluate.py` / `jackknife.py`). Discovery runs by default;
   inclusion is earned and self-correcting across genome types.

**Validation.** Primary working test = **xgXerMont** (already wired end-to-end; a
metazoan with a degenerate U2 branch — exercises the generic EM path): check that
EM recovers a sensible branch consensus and the toggle behaves. This is a *weak*
positive control (a snail's branches aren't independently characterised); a genome
with **known/characterised U2 branch points** (e.g. a small fungal genome with the
sharp TACTAAC motif) is wanted for a stronger sanity check — TBD which.
