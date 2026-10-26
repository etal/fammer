"""Use TMalign for multiple structure alignment.

Procedure:

1. Run all-vs-all pairwise alignments with TMalign; keep the alignments and
   also the TM-score (a metric of structural similarity/alignability).
2. Build a graph with one node per structure and one edge per pair, weighted
   by TM-score, and determine its maximum spanning tree.
3. Take the pairwise alignments corresponding to the edges of the tree, N-1 of
   them for N structures, and use them as seeds in MAFFT (--localpair).
4. Remove duplicated sequences and all-gap columns from the resulting alignment.

Each structure must be listed only once, because a structure listed twice would
be aligned to itself. No PDB file name may contain whitespace, because the name
becomes a sequence identifier in FASTA files, and a FASTA identifier ends at
the first whitespace. No seed file name may contain whitespace, because MAFFT
splits seed file paths at whitespace.

"""

import collections
import itertools
import logging
import math
import os
import subprocess
import tempfile
from io import StringIO

import networkx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biofrills import alnutils


def _check_file_names(pdb_fnames, seed_fnames):
    """Reject file names that the alignment procedure cannot carry.

    Each PDB file name becomes the identifier of a sequence in the alignment,
    and build reads it back from the .pdb.seq file as a path. A FASTA
    identifier ends at the first whitespace, so a name containing whitespace
    would be truncated; MAFFT likewise splits a seed file path at whitespace.
    A structure listed twice would be aligned to itself.
    """
    split = [fname for fname in [*pdb_fnames, *seed_fnames] if fname.split() != [fname]]
    if split:
        raise ValueError(
            "File names must not contain whitespace: " + ", ".join(map(repr, split))
        )
    counts = collections.Counter(pdb_fnames)
    repeated = sorted(fname for fname, count in counts.items() if count > 1)
    if repeated:
        raise ValueError(
            "Each structure must be given once: " + ", ".join(map(repr, repeated))
        )


def align_structs(pdb_fnames, seed_fnames=None):
    """Align multiple PDB structures using TM-align.

    Returns a list of aligned SeqRecords.
    """
    seed_fnames = list(seed_fnames or [])
    _check_file_names(pdb_fnames, seed_fnames)
    # 1. Align all-v-all structure pairs with TM-align
    allpairs = []
    for idx, ref_pdbfn in enumerate(pdb_fnames):
        for eqv_pdbfn in pdb_fnames[idx + 1 :]:
            logging.info("Aligning %s to %s", eqv_pdbfn, ref_pdbfn)
            try:
                tm_output = subprocess.check_output(
                    ["TMalign", ref_pdbfn, eqv_pdbfn], text=True
                )
            except OSError:
                logging.warning("Failed command: TMalign %s %s", ref_pdbfn, eqv_pdbfn)
                for fname in (ref_pdbfn, eqv_pdbfn):
                    if not os.path.isfile(fname):
                        logging.warning("Missing file: %s", fname)
                raise
            except subprocess.CalledProcessError as exc:
                raise RuntimeError(
                    f"TMalign failed (returned {exc.returncode}):\n{exc.output}"
                ) from exc
            tm_seqpair = read_tmalign_as_seqrec_pair(tm_output, ref_pdbfn, eqv_pdbfn)
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
        with tempfile.NamedTemporaryFile("w+", delete=False) as handle:
            SeqIO.write(seedpair, handle, "fasta")
            tmp_seed_fnames.append(handle.name)

    # 3. Use MAFFT to combine TMalign'd pairs into a multiple alignment;
    seq_fd, seq_fname = tempfile.mkstemp(text=True)
    # Create a blank file to appease MAFFT
    os.write(seq_fd, b"")
    os.close(seq_fd)
    mafft_output = subprocess.check_output(
        ["mafft", "--quiet", "--amino", "--localpair", "--maxiterate", "1000"]
        + list(
            itertools.chain(*[("--seed", sfn) for sfn in seed_fnames + tmp_seed_fnames])
        )
        + [seq_fname],
        text=True,
    )
    # Clean up
    os.remove(seq_fname)
    for sfn in tmp_seed_fnames:
        os.remove(sfn)

    # 4. Emit the aligned sequences
    recs = SeqIO.parse(StringIO(mafft_output), "fasta")
    recs = clean_and_dedupe_seqs(recs)
    recs = alnutils.remove_empty_cols(recs)
    recs = purge_seqs(recs)
    return list(recs)


def read_tmalign_as_seqrec_pair(tm_output, ref_id, eqv_id):
    """Parse TMalign output into a pair of SeqRecords and their mean TM-score."""
    lines = tm_output.splitlines()
    # Extract the TM-score (measure of structure similarity)
    # Take the mean of the (two) given TM-scores -- not sure which is reference
    tmscores = []
    for line in lines:
        if line.startswith("TM-score"):
            # TMalign v. 2012/05/07 or earlier
            tmscores.append(float(line.split(None, 2)[1]))
        elif "TM-score=" in line:
            # TMalign v. 2013/05/11 or so
            tokens = line.split()
            for token in tokens:
                if token.startswith("TM-score="):
                    _key, _val = token.split("=")
                    tmscores.append(float(_val.rstrip(",")))
                    break
    if not tmscores:
        raise ValueError(f"No TM-score in TM-align output for {ref_id}, {eqv_id}")
    tmscore = math.fsum(tmscores) / len(tmscores)
    # Extract the sequence alignment, which follows the line explaining the
    # ":" and "." markers; newer versions print more lines after it.
    header = next(
        (idx for idx, line in enumerate(lines) if line.startswith('(":"')), None
    )
    if header is None or header + 3 >= len(lines):
        raise ValueError(f"No alignment in TM-align output for {ref_id}, {eqv_id}")
    refseq, eqvseq = lines[header + 1].strip(), lines[header + 3].strip()
    return (
        SeqRecord(Seq(refseq), id=ref_id, description=f"TMalign TM-score={tmscore:f}"),
        SeqRecord(Seq(eqvseq), id=eqv_id, description=f"TMalign TM-score={tmscore:f}"),
        tmscore,
    )


def mst_pairs(pairs):
    """Select the pairwise alignments that connect all structures best.

    Build an undirected graph with one node per structure, keyed by record ID,
    and one edge per pair, weighted by TM-score. Return the record pairs of the
    edges of its maximum spanning tree: N-1 pairs for N structures.

    Pairs are added in order of decreasing score and then of ID, each with its
    two IDs in sorted order, so the graph, and with it the selected tree, does
    not depend on the order of the pairs or of the two records in a pair.

    Input: iterable of (SeqRecord, SeqRecord, TM-score)
    Output: list of (SeqRecord, SeqRecord)
    """
    edges = [
        (score, *sorted((left.id, right.id)), (left, right))
        for left, right, score in pairs
    ]
    G = networkx.Graph()
    for score, node1, node2, records in sorted(
        edges, key=lambda edge: (-edge[0], edge[1], edge[2])
    ):
        G.add_edge(node1, node2, weight=score, records=records)
    return [
        data["records"]
        for _left, _right, data in networkx.maximum_spanning_edges(G, data=True)
    ]


def tmscore_from_description(text):
    for token in text.split():
        if token.startswith("TM-score"):
            return float(token.split("=", 1)[1])


def clean_and_dedupe_seqs(records, best_score=False):
    """Remove the _seed_ prefix and omit duplicated records (by ID)."""
    if best_score:
        seen = {}
    else:
        seen = set()
    for record in records:
        # Remove the _seed_ prefix from each sequence ID
        if record.id.startswith("_seed_"):
            record.id = record.id[len("_seed_") :]
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
