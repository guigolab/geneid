"""End-to-end training orchestration: gene models + genome -> geneid .param.

Ties the validated stages together — background model, site profiles (start /
donor / acceptor), coding-potential Markov matrices, gene-model length stats —
and assembles them into a complete parameter file. Each stage lives in its own
module and is validated independently against the reference run; this module is
the wiring.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from .param.assemble import assemble_param
from .prepare.base import GeneModel
from .prepare.sites import collect_sites
from .stats.background import from_genome
from .stats.coding import choose_orders, derive_coding_potential, format_markov_matrix
from .stats.genemodel import format_range, intron_range
from .stats.sites import (
    Matrix,
    format_profile,
    frequency,
    info_content,
    log_ratio,
    log_ratio_zero_order,
    mask_invariant_dinuc,
    position_matrix,
    profile_header,
    select_window,
    submatrix,
)

# order-selection thresholds for site profiles (from the legacy driver)
_MIN_SITES_ORDER1 = 1400  # donor/acceptor: order 1 above this, else order 0
_MIN_STARTS_ORDER2 = 5500  # start: order 2 above this, else order 0


def build_site_profile(
    seqs: Sequence[str],
    *,
    site: str,
    anchor: str,
    order: int,
    bg_freq: Matrix,
    bg_dimatrix: Matrix,
) -> list[str]:
    """Estimate one site profile and return its geneid param data lines.

    Order>=1 (donor/acceptor): dinucleotide log-ratio over the background matrix,
    windowed by information content, with the invariant ``anchor`` masked.
    Order 0 (start): single-nucleotide log-ratio (log_ratio_zero_order), windowed,
    no explicit mask (the invariant base masks itself via its 0/1 frequency).
    """
    freq = frequency(seqs)
    info = info_content(freq, bg_freq)
    window = select_window(info, site=site, order=order)
    if order >= 1:
        logm = log_ratio(position_matrix(seqs, order=order), bg_dimatrix)
        prof = mask_invariant_dinuc(
            submatrix(logm, window.start, window.end),
            window.st,
            window.nd,
            window.rd,
            anchor,
        )
    else:
        prof = submatrix(log_ratio_zero_order(freq, bg_freq), window.start, window.end)
    return format_profile(prof, profile_header(window, order))


def _site_order(n_donor: int, n_acceptor: int, n_start: int) -> tuple[int, int, int]:
    donor = 1 if n_donor >= _MIN_SITES_ORDER1 else 0
    acceptor = 1 if n_acceptor >= _MIN_SITES_ORDER1 else 0
    start = 2 if n_start >= _MIN_STARTS_ORDER2 else 0
    if donor == 0 or acceptor == 0:
        raise NotImplementedError(
            "order-0 donor/acceptor profiles are not yet supported "
            f"(donor={n_donor}, acceptor={n_acceptor} sites; need >= {_MIN_SITES_ORDER1}); "
            "provide a larger training set"
        )
    if start == 2:
        raise NotImplementedError(
            f"order-2 start profiles are not yet supported ({n_start} starts); "
            "order-2 masking is unvalidated against the reference"
        )
    return donor, acceptor, start


def train(
    models: Sequence[GeneModel],
    genome: Mapping[str, str],
    species: str,
    *,
    background: tuple[Matrix, Matrix] | None = None,
    seed: int = 0,
) -> str:
    """Train a geneid parameter file from complete, filtered gene ``models``.

    ``background`` may be supplied as ``(freq, dimatrix)`` to make a run fully
    deterministic (and to validate against a reference background); otherwise it
    is sampled from the genome.
    """
    if not models:
        raise ValueError("no gene models to train on")

    sites = collect_sites(models, genome)
    n_donor, n_acceptor, n_start = len(sites.donor), len(sites.acceptor), len(sites.start)
    donor_order, acceptor_order, start_order = _site_order(n_donor, n_acceptor, n_start)

    bg_freq, bg_dimatrix = background if background is not None else from_genome(
        genome.values(), seed=seed
    )

    start_lines = build_site_profile(
        sites.start, site="start", anchor="ATG", order=start_order,
        bg_freq=bg_freq, bg_dimatrix=bg_dimatrix,
    )
    donor_lines = build_site_profile(
        sites.donor, site="donor", anchor="GT", order=donor_order,
        bg_freq=bg_freq, bg_dimatrix=bg_dimatrix,
    )
    acceptor_lines = build_site_profile(
        sites.acceptor, site="acceptor", anchor="AG", order=acceptor_order,
        bg_freq=bg_freq, bg_dimatrix=bg_dimatrix,
    )

    cds_seqs = [m.cds(genome) for m in models]
    intron_seqs = [s for m in models for s in m.intron_seqs(genome)]
    init_logs, trans_logs, cbases, nbases = derive_coding_potential(cds_seqs, intron_seqs)
    _, transition_order = choose_orders(cbases, nbases)

    lo, hi = intron_range([len(s) for s in intron_seqs])

    return assemble_param(
        species=species,
        start_profile=start_lines,
        acceptor_profile=acceptor_lines,
        donor_profile=donor_lines,
        markov_order=transition_order,
        markov_initial=format_markov_matrix(init_logs),
        markov_transition=format_markov_matrix(trans_logs),
        intron_range=format_range(lo, hi),
        intergenic_range="200:Infinity",
    )
