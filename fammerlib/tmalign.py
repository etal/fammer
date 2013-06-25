"""Align multiple structures with TMalign."""

import itertools
import logging
import os
import subprocess
import tempfile
from cStringIO import StringIO

import networkx

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biofrills import alnutils


def align_structs(pdb_fnames, seed_fnames=None):
    """Align multiple PDB structures using TM-align.

    Returns a list of aligned SeqRecords.
    """
    # 1. Align all-v-all structure pairs with TM-align
    allpairs = []
    for idx, ref_pdbfn in enumerate(pdb_fnames):
        assert ' ' not in ref_pdbfn
        for eqv_pdbfn in pdb_fnames[idx+1:]:
            assert eqv_pdbfn != ref_pdbfn
            logging.info("Aligning %s to %s", eqv_pdbfn, ref_pdbfn)
            try:
                tm_output = subprocess.check_output(['TMalign',
                                                     ref_pdbfn, eqv_pdbfn])
            except OSError:
                logging.warning("Failed command: TMalign %s %s",
                                ref_pdbfn, eqv_pdbfn)
                for fname in (ref_pdbfn, eqv_pdbfn):
                    if not os.path.isfile(fname):
                        logging.warning("Missing file: %s", fname)
                raise
            except subprocess.CalledProcessError, exc:
                raise RuntimeError("TMalign failed (returned %s):\n%s"
                                   % (exc.returncode, exc.output))
            tm_seqpair = read_tmalign_as_seqrec_pair(tm_output,
                                                     ref_pdbfn, eqv_pdbfn)
            allpairs.append(tm_seqpair)

    # In case of 2 structs, no need to combine alignments -- we're done
    if len(allpairs) == 1:
        recs = allpairs[0][:2]
        return alnutils.remove_empty_cols(recs)

    # 2. Resolve MST pairs & write seed tempfiles
    tmp_seed_fnames = []
    for seedpair in mst_pairs(allpairs):
        # fd, seedfn = tempfile.mkstemp(text=True)
        # SeqIO.write(seedpair, seedfn, 'fasta')
        # SeqIO.write(seedpair, os.fdopen(fd), 'fasta')
        with tempfile.NamedTemporaryFile('w+', delete=False) as handle:
            SeqIO.write(seedpair, handle, 'fasta')
            tmp_seed_fnames.append(handle.name)

    # 3. Use MAFFT to combine TMalign'd pairs into a multiple alignment;
    seq_fd, seq_fname = tempfile.mkstemp(text=True)
    # Create a blank file to appease MAFFT
    os.write(seq_fd, '')
    os.close(seq_fd)
    mafft_output = subprocess.check_output(['mafft',
        '--quiet', '--amino', '--localpair',
        '--maxiterate', '1000']
        + list(itertools.chain(*[('--seed', sfn)
                                 for sfn in (seed_fnames or []) + tmp_seed_fnames]))
        + [seq_fname])
    # Clean up
    os.remove(seq_fname)
    for sfn in tmp_seed_fnames:
        os.remove(sfn)

    # 4. Emit the aligned sequences
    recs = SeqIO.parse(StringIO(mafft_output), 'fasta')
    recs = clean_and_dedupe_seqs(recs)
    recs = alnutils.remove_empty_cols(recs)
    recs = purge_seqs(recs)
    return list(recs)


def read_tmalign_as_seqrec_pair(tm_output, ref_id, eqv_id):
    """Create a pair of SeqRecords from TMalign output."""
    lines = tm_output.splitlines()
    # Extract the TM-score (measure of structure similarity)
    # Take the mean of the (two) given TM-scores -- not sure which is reference
    tmscores = []
    for line in lines:
        if line.startswith('TM-score'):
            # TMalign v. 2012/05/07 or earlier
            tmscores.append(float(line.split(None, 2)[1]))
        elif 'TM-score=' in line:
            # TMalign v. 2013/05/11 or so
            tokens = line.split()
            for token in tokens:
                if token.startswith('TM-score='):
                    _key, _val = token.split('=')
                    tmscores.append(float(_val.rstrip(',')))
                    break
    tmscore = sum(tmscores) / len(tmscores)
    # Extract the sequence alignment
    lastlines = lines[-5:]
    assert lastlines[0].startswith('(":"') # (":" denotes the residues pairs
    assert not lastlines[-1].strip()
    refseq, eqvseq = lastlines[1].strip(), lastlines[3].strip()
    return (SeqRecord(Seq(refseq), id=ref_id,
                      description="TMalign TM-score=%f" % tmscore),
            SeqRecord(Seq(eqvseq), id=eqv_id,
                      description="TMalign TM-score=%f" % tmscore),
            tmscore)


def mst_pairs(pairs):
    """Given all pairwise distances, determine the minimal spanning subset.

    Convert pairwise distances to an undirected graph, determine the
    minumum spanning tree, and emit the minimal list of edges to connect all
    nodes.

    Input: iterable of (SeqRecord, SeqRecord, distance)
    Output: iterable of (SeqRecord, SeqRecord)
    """
    G = networkx.Graph()
    for left, right, score in pairs:
        G.add_edge(left, right, weight=1.0/score)
    mst = networkx.minimum_spanning_edges(G, data=False)
    return list(mst)


def tmscore_from_description(text):
    for token in text.split():
        if token.startswith('TM-score'):
            return float(token.split('=', 1)[1])


def clean_and_dedupe_seqs(records, best_score=False):
    """Remove the _seed_ prefix and omit duplicated records (by ID)."""
    if best_score:
        seen = {}
    else:
        seen = set()
    for record in records:
        # Remove the _seed_ prefix from each sequence ID
        if record.id.startswith('_seed_'):
            record.id = record.id[len('_seed_'):]
            record.name = record.id
        # Check for duplicates.
        if best_score:
            # If a previously seen PDB was aligned better (per the TM-score),
            # defer to that one
            tmscore = tmscore_from_description(record.description)
            if record.id in seen and seen[record.id] >= tmscore:
                # This PDB was aligned better previously; skip
                continue
            seen[record.id] = tmscore
        else:
            # Keep a duplicate sequence if it was aligned differently
            ident = (record.id, str(record.seq))
            if ident in seen:
                continue
            seen.add(ident)
        yield record


def purge_seqs(records):
    """Drop duplicated records by identical sequence."""
    seen = set()
    for rec in records:
        seq = str(rec.seq)
        if seq not in seen:
            yield rec
            seen.add(seq)

