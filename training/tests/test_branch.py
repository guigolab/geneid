import random

from geneid_train.prepare.base import Exon, GeneModel
from geneid_train.stats.branch import (
    BRANCH_A,
    branch_distances,
    branch_windows,
    distance_knobs,
    fit_branch_em,
    locate_branch,
)


def _planted(n=400, motif="TACTAAC", dist_from_end=20, seed=0):
    """n random windows each with `motif` planted at a fixed distance from the 3' end."""
    rng = random.Random(seed)
    wins = []
    for _ in range(n):
        length = rng.randint(35, 45)
        bases = [rng.choice("ACGT") for _ in range(length)]
        start = length - dist_from_end
        bases[start : start + len(motif)] = list(motif)
        wins.append("".join(bases))
    return wins


def test_em_recovers_planted_motif():
    wins = _planted(motif="TACTAAC")
    model = fit_branch_em(wins, width=7)
    # the per-position argmax base should spell the planted motif
    consensus = "".join(max("ACGT", key=lambda b: pos[b]) for pos in model.pwm)
    assert consensus == "TACTAAC"
    # the branch adenosine column is (nearly) invariant A
    assert model.pwm[BRANCH_A]["A"] > 0.9


def test_locate_branch_finds_the_plant():
    wins = _planted(n=50, motif="TACTAAC", dist_from_end=20)
    model = fit_branch_em(wins, width=7)
    j, score = locate_branch(wins[0], model)
    # motif planted at length-20; its start index is where locate should land
    assert wins[0][j : j + 7] == "TACTAAC"
    assert score > 0


def test_distances_cluster_at_the_plant():
    # motif at a fixed distance -> branch-A distance clusters at that value + anchor offset
    wins = _planted(n=200, motif="TACTAAC", dist_from_end=20)
    model = fit_branch_em(wins, width=7)
    dists = branch_distances(wins, model)
    # branch A is BRANCH_A into the motif; motif start is 20 from window end, plus _AG(2)
    # to the 3'SS -> branch-A distance = (20 - BRANCH_A) + 2
    expected = (20 - BRANCH_A) + 2
    mode = max(set(dists), key=dists.count)
    assert mode == expected
    knobs = distance_knobs(dists)
    assert knobs.min_dist <= knobs.opt_dist <= knobs.acc_context
    assert knobs.opt_dist == expected


def test_branch_profile_emission_format():
    from geneid_train.stats.branch import (
        BranchDistances,
        branch_profile_lines,
        branch_profile_section,
        branch_weight_scalar,
    )

    wins = _planted(n=100, motif="TACTAAC")
    model = fit_branch_em(wins, width=7)
    knobs = BranchDistances(acc_context=46, min_dist=5, opt_dist=27)
    lines = branch_profile_lines(model, knobs)
    # header: len offset cutoff order a b acc_context min_dist opt_dist pen_scale
    assert lines[0].split() == ["7", "5", "-20", "0", "0", "1", "46", "5", "27", "6"]
    assert lines[1].startswith("#")
    rows = {(int(p), b): float(v) for p, b, v in (ln.split() for ln in lines[2:])}
    assert len(rows) == 7 * 4  # every position x base
    # the invariant branch-A column is strongly positive (log-odds of A >> others)
    assert rows[(6, "A")] > rows[(6, "C")] and rows[(6, "A")] > rows[(6, "G")]
    # section + scalar wrappers
    assert branch_profile_section(model, knobs).startswith("Branch_point_profile\n")
    assert branch_weight_scalar(0) == "Branch_point_score_weight\n0\n"


def test_branch_windows_excludes_ag_and_is_adaptive():
    # intron: donor GT ... branch region ... acceptor AG
    intron = "GT" + "C" * 30 + "TACTAAC" + "CCCC" + "AG"
    exon = "AAAAAA"
    genome = {"c": exon + intron + exon}
    a2 = 6 + len(intron) + 1
    m = GeneModel("g", "c", "+", [Exon(1, 6), Exon(a2, a2 + 5)])
    wins = branch_windows([m], genome, max_window=45)
    assert len(wins) == 1
    # the terminal AG is excluded; the branch motif is inside the window
    assert not wins[0].endswith("AG")
    assert "TACTAAC" in wins[0]
