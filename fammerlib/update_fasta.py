"Update FASTA: Overwrite the original FASTA files with ungapped .aln sequences."

import logging
import os
from os.path import isdir, isfile, join

from Bio import SeqIO

from ._share import write_fasta
from .tasks import Task


def cmd_update_fasta(args):
    if isdir(args.basedir):
        for task in find_update_tasks(args.basedir):
            task.build()
    elif isfile(args.basedir) and args.basedir.endswith('.aln'):
        task = one_update_task(args.basedir)
        if task is not None:
            task.build()
    else:
        raise RuntimeError("Target must be a directory or .aln file: %s"
                         % args.basedir)


def find_update_tasks(topdir):
    """Locate and update all .fasta sequences."""
    # if <base>.aln and <base>.fasta exist (same base)
    # and .aln is newer than .fasta (handled by Task)
    # [ENH: and same # seqs?]
    # then update_fasta()
    for dirpath, _subdirs, fnames in os.walk(topdir):
        todo_bases = [join(dirpath, fn[:-4])
                      for fn in fnames
                      if fn.endswith('.aln') and fn[:-4] + '.fasta' in fnames]
        for base in todo_bases:
            if not is_same_aln_fasta(base + '.aln', base + '.fasta'):
                yield Task(base + '.fasta',
                    action=update_fasta,
                    depends=base + '.aln')


def one_update_task(aln_fname):
    fasta_fname = aln_fname[:-4] + '.fasta'
    if not isfile(fasta_fname):
        raise RuntimeError("No corresponding .fasta file for %s: need %s"
                           % (aln_fname, fasta_fname))
    if is_same_aln_fasta(aln_fname, fasta_fname):
        logging.info("%s already matches %s.", fasta_fname, aln_fname)
    else:
        return Task(fasta_fname,
                    action=update_fasta,
                    depends=aln_fname)


def update_fasta(task):
    """Overwrite the .fasta with ungapped unique .aln sequences."""
    records = ungap_and_unique(SeqIO.parse(str(task.depends[0]), 'clustal'),
                               task.target)
    write_fasta(records, task.target)
    # SeqIO.write(records, task.target, 'fasta')


def ungap_and_unique(records, target_name):
    """Remove gap characters, drop duplicate SeqRecords, and sort by length.

    Equivalent to:
    seqioconvert clustal fasta < t.aln | purge.py -c | press.py -u > t.fasta
    """
    # ENH: warn about:
    # - Sequences containing 'X' (log the number of X's)
    # - Short sequences (e.g. if seq_length < .5 * alignment_width)

    seen = set()
    counter = 0
    for rec in sorted(records,
            key=lambda r: len(r.seq) - r.seq.count('-') - r.seq.count('.'),
            reverse=True):
        # Remove gaps ('-')
        thisseq = str(rec.seq).replace('-', '').replace('.', '').upper()
        # Drop if this sequence is completely contained by another
        if 'X' in thisseq:
            # Ignore ambiguous chars by checking each unambigous segment
            for segment in thisseq.split('X'):
                for longseq in seen:
                    if segment in longseq:
                        # Check the next segment
                        break
                else:
                    # No match -- thisseq is novel
                    seen.add(thisseq)
                    rec.seq._data = thisseq
                    yield rec
                    break
            else:
                # All fragments have been "seen"
                counter += 1
                break
        else:
            for longseq in seen:
                if thisseq in longseq:
                    counter += 1
                    break
            else:
                seen.add(thisseq)
                # Apply the ungapping from above
                rec.seq._data = thisseq
                yield rec

    # I/O
    msg = "Updated %s" % target_name
    if counter:
        msg += (" (purged %d sequence%s)"
                % (counter, '' if counter == 1 else 's'))
    logging.info(msg)


def is_same_aln_fasta(aln_fn, fasta_fn):
    """Compare ungapped seqs in aln and fasta; are they all the same?"""
    def rec2tuple(rec, do_ungap):
        if do_ungap:
            strseq = str(rec.seq).replace('-', '').replace('.', '').upper()
        else:
            strseq = str(rec.seq).upper()
        return (rec.id, strseq)

    aln_hash = set([rec2tuple(rec, True)
                    for rec in SeqIO.parse(aln_fn, 'clustal')])
    fa_hash = set([rec2tuple(rec, False)
                   for rec in SeqIO.parse(fasta_fn, 'fasta')])
    return aln_hash == fa_hash



