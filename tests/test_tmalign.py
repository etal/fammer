"""Properties of fammer.tmalign that need no external program."""

import itertools

import networkx
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fammer import tmalign

# Four structures with tied scores, so that several spanning trees share the
# maximum total and only the tie-breaking decides between them.
SCORES = {
    ("A", "B"): 0.8,
    ("A", "C"): 0.8,
    ("B", "C"): 0.8,
    ("C", "D"): 0.5,
    ("B", "D"): 0.5,
    ("A", "D"): 0.2,
}


def scored_pairs(order, flipped):
    """A fresh record pair per TM-align run, as align_structs makes them."""
    pairs = []
    for left, right in order:
        if (left, right) in flipped:
            left, right = right, left
        score = SCORES[tuple(sorted((left, right)))]
        pairs.append(
            (SeqRecord(Seq("ACDE"), id=left), SeqRecord(Seq("ACDE"), id=right), score)
        )
    return pairs


def edge_set(record_pairs):
    return frozenset(frozenset((left.id, right.id)) for left, right in record_pairs)


def max_total_score():
    """The largest total score of any spanning tree, by enumeration."""
    best = 0.0
    for edges in itertools.combinations(SCORES, 3):
        graph = networkx.Graph(list(edges))
        if graph.number_of_nodes() == 4 and networkx.is_tree(graph):
            best = max(best, sum(SCORES[edge] for edge in edges))
    return best


def test_mst_pairs_is_a_maximum_spanning_tree_independent_of_input_order():
    results = set()
    for order in itertools.permutations(SCORES):
        for flipped in ((), order[::2]):
            pairs = scored_pairs(order, set(flipped))
            selected = tmalign.mst_pairs(pairs)
            # The records returned are those of the input pairs, not copies.
            given = {id(rec) for left, right, _score in pairs for rec in (left, right)}
            assert all(id(rec) in given for pair in selected for rec in pair)
            results.add(edge_set(selected))
    assert len(results) == 1
    (edges,) = results
    assert len(edges) == 3
    assert networkx.is_tree(networkx.Graph([tuple(edge) for edge in edges]))
    total = sum(SCORES[tuple(sorted(edge))] for edge in edges)
    assert total == pytest.approx(max_total_score())
