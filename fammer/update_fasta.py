"Update FASTA: Overwrite the original FASTA files with ungapped .aln sequences."

import logging
import os
import re
from os.path import isdir, isfile, join

from Bio import SeqIO
from Bio.Seq import Seq

from ._share import write_fasta
from .tasks import Task


def cmd_update_fasta(args):
    if isdir(args.basedir):
        for task in find_update_tasks(args.basedir):
            task.build()
    elif isfile(args.basedir) and args.basedir.endswith(".aln"):
        task = one_update_task(args.basedir)
        if task is not None:
            task.build()
    else:
        raise RuntimeError(f"Target must be a directory or .aln file: {args.basedir}")


def find_update_tasks(topdir):
    """Locate and update all .fasta sequences."""
    # if <base>.aln and <base>.fasta exist (same base)
    # and .aln is newer than .fasta (handled by Task)
    # [ENH: and same # seqs?]
    # then update_fasta()
    for dirpath, _subdirs, fnames in os.walk(topdir):
        todo_bases = [
            join(dirpath, fn[:-4])
            for fn in fnames
            if fn.endswith(".aln") and fn[:-4] + ".fasta" in fnames
        ]
        for base in todo_bases:
            if not is_same_aln_fasta(base + ".aln", base + ".fasta"):
                yield Task(base + ".fasta", action=update_fasta, depends=base + ".aln")


def one_update_task(aln_fname):
    fasta_fname = aln_fname[:-4] + ".fasta"
    if not isfile(fasta_fname):
        raise RuntimeError(
            f"No corresponding .fasta file for {aln_fname}: need {fasta_fname}"
        )
    if is_same_aln_fasta(aln_fname, fasta_fname):
        logging.info("%s already matches %s.", fasta_fname, aln_fname)
    else:
        return Task(fasta_fname, action=update_fasta, depends=aln_fname)


def update_fasta(task):
    """Overwrite the .fasta with ungapped unique .aln sequences."""
    # ungap_and_unique reads the alignment while write_fasta writes, so the
    # alignment stays open until the FASTA file is written.
    with SeqIO.parse(str(task.depends[0]), "clustal") as parsed:
        write_fasta(ungap_and_unique(parsed, task.target), task.target)
    # SeqIO.write(records, task.target, 'fasta')


def ungap_and_unique(records, target_name):
    """Remove gap characters, drop redundant SeqRecords, and sort by length.

    Records are visited in order of decreasing ungapped length, and a record
    is dropped if its ungapped sequence occurs within a sequence already
    kept.  An 'X' in the record's sequence matches any single residue; an 'X'
    in a kept sequence matches only an 'X'.
    """
    # ENH: warn about:
    # - Sequences containing 'X' (log the number of X's)
    # - Short sequences (e.g. if seq_length < .5 * alignment_width)

    seen = set()
    counter = 0
    for rec in sorted(
        records,
        key=lambda r: len(r.seq) - r.seq.count("-") - r.seq.count("."),
        reverse=True,
    ):
        # Remove gaps ('-')
        thisseq = str(rec.seq).replace("-", "").replace(".", "").upper()
        # Drop if this sequence is completely contained by another
        if "X" in thisseq:
            # Each ambiguous residue matches any single residue
            pattern = re.compile(re.escape(thisseq).replace("X", "."))
            redundant = any(pattern.search(longseq) for longseq in seen)
        else:
            redundant = any(thisseq in longseq for longseq in seen)
        if redundant:
            counter += 1
            continue
        seen.add(thisseq)
        # Apply the ungapping from above
        rec.seq = Seq(thisseq)
        yield rec

    # I/O
    msg = f"Updated {target_name}"
    if counter:
        msg += f" (purged {counter} sequence{'' if counter == 1 else 's'})"
    logging.info(msg)


def is_same_aln_fasta(aln_fn, fasta_fn):
    """Compare ungapped seqs in aln and fasta; are they all the same?"""

    def rec2tuple(rec, do_ungap):
        if do_ungap:
            strseq = str(rec.seq).replace("-", "").replace(".", "").upper()
        else:
            strseq = str(rec.seq).upper()
        return (rec.id, strseq)

    with SeqIO.parse(aln_fn, "clustal") as parsed:
        aln_hash = {rec2tuple(rec, True) for rec in parsed}
    with SeqIO.parse(fasta_fn, "fasta") as parsed:
        fa_hash = {rec2tuple(rec, False) for rec in parsed}
    return aln_hash == fa_hash
