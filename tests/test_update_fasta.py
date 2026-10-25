"""Tests of update_fasta.ungap_and_unique, the redundancy filter of update-fasta."""

import logging
import os

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fammer.update_fasta import ungap_and_unique

CDK_ALN = os.path.join(os.path.dirname(__file__), 'data', 'py27', 'e2e',
                       'build', 'SomeKinases', 'CMGC', 'CDK.aln')


def make_records(pairs):
    return [SeqRecord(Seq(seq), id=seq_id) for seq_id, seq in pairs]


def kept_ids(pairs):
    return [rec.id for rec in ungap_and_unique(make_records(pairs), 't.fasta')]


def purge_message(caplog, pairs):
    with caplog.at_level(logging.INFO):
        list(ungap_and_unique(make_records(pairs), 't.fasta'))
    return caplog.records[-1].getMessage()


def test_records_after_a_redundant_x_record_are_kept(caplog):
    """A redundant record containing X is skipped, not the end of the input."""
    pairs = [('long', 'ACDEFGHIKL'), ('withX', 'CDE-XGH'),
             ('short1', 'MNP'), ('short2', 'QRS')]
    assert kept_ids(pairs) == ['long', 'short1', 'short2']
    assert purge_message(caplog, pairs) == \
        'Updated t.fasta (purged 1 sequence)'


def test_purge_count_includes_both_kinds_of_redundant_record(caplog):
    """Redundant records with and without X are counted; kept ones are not."""
    pairs = [('long', 'ACDEFGHIKL'), ('novelX', 'WWXYYM'),
             ('redundantX', 'CDEXGH'), ('redundant', 'FGH'), ('novel', 'MNP')]
    assert kept_ids(pairs) == ['long', 'novelX', 'novel']
    assert purge_message(caplog, pairs) == \
        'Updated t.fasta (purged 2 sequences)'


def test_x_segments_found_in_different_records_are_not_redundant():
    pairs = [('first', 'AAAAACCCCC'), ('second', 'DDDDDEEEEE'),
             ('split', 'CCCXDDD')]
    assert kept_ids(pairs) == ['first', 'second', 'split']


@pytest.mark.parametrize('query, redundant', [
    ('CDEXGH', True),    # X stands for F
    ('XCDE', True),      # X stands for A
    ('GHIKX', True),     # X stands for L
    ('CXFG', False),     # X would have to stand for two residues, D and E
    ('XACD', False),     # no residue precedes A
    ('HIKLX', False),    # no residue follows L
    ('HIXCDE', False),   # both parts occur, but not in this order
    ('XXXXXXXXXX', True),  # as long as the kept sequence
])
def test_x_matches_exactly_one_residue(query, redundant):
    kept = kept_ids([('long', 'ACDEFGHIKL'), ('query', query)])
    assert kept == (['long'] if redundant else ['long', 'query'])


def test_x_in_a_kept_sequence_matches_only_x():
    assert kept_ids([('long', 'ACXEFGH'), ('plain', 'CDE')]) == \
        ['long', 'plain']
    assert kept_ids([('long', 'ACXEFGH'), ('withX', 'CXEF')]) == ['long']
    assert kept_ids([('a', 'ACXEF'), ('b', 'ACXEF')]) == ['a']


def test_other_characters_are_matched_literally():
    """Characters that are special in a regular expression match themselves."""
    assert kept_ids([('long', 'ACDEFGHIK'), ('stop', 'W*X')]) == \
        ['long', 'stop']
    assert kept_ids([('long', 'ACW*DEF'), ('stop', 'W*X')]) == ['long']


def test_alignment_after_an_edit_keeps_every_novel_record(caplog):
    """The built CDK alignment with its last record removed, as a curator
    might leave it: of 144 records, only the two X-containing structure
    sequences of kinases already in the set are redundant."""
    records = list(SeqIO.parse(CDK_ALN, 'clustal'))[:-1]
    assert len(records) == 144
    with caplog.at_level(logging.INFO):
        kept = list(ungap_and_unique(records, 'CDK.fasta'))
    assert len(kept) == 142
    assert {rec.id for rec in records} - {rec.id for rec in kept} == \
        {'gi|193506638|pdb|3BLH|A', 'gi|58176911|pdb|1UA2|A'}
    assert caplog.records[-1].getMessage() == \
        'Updated CDK.fasta (purged 2 sequences)'


@pytest.mark.parametrize('pairs', [
    [('long', 'ACDEFGHIKL'), ('withX', 'CDE-XGH'), ('a', 'MNP'), ('b', 'QRS')],
    [('a', 'AAAAACCCCC'), ('b', 'DDDDDEEEEE'), ('c', 'CCCXDDD'), ('d', 'MNP')],
    [('a', 'ACDEF'), ('b', 'ACDEF'), ('c', 'gh-ik.'), ('d', 'GHIK')],
    [('a', 'XXXX'), ('b', 'AC'), ('c', ''), ('d', '----')],
])
def test_filtering_twice_changes_nothing(pairs):
    once = [(rec.id, str(rec.seq)) for rec in
            ungap_and_unique(make_records(pairs), 't.fasta')]
    twice = [(rec.id, str(rec.seq)) for rec in
             ungap_and_unique(make_records(once), 't.fasta')]
    assert twice == once


def test_filtering_a_real_alignment_twice_changes_nothing():
    once = list(ungap_and_unique(SeqIO.parse(CDK_ALN, 'clustal'), 't.fasta'))
    twice = list(ungap_and_unique(
        [SeqRecord(Seq(str(rec.seq)), id=rec.id) for rec in once], 't.fasta'))
    assert [(r.id, str(r.seq)) for r in twice] == \
        [(r.id, str(r.seq)) for r in once]
