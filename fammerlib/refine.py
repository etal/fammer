"Refine: Leave-one-out validation and refinement of each profile."

import logging
import os
import tempfile
from glob import glob
from os.path import abspath, basename, isdir, isfile, join

from Bio import SeqIO

from .tasks import ext, sh
from ._share import get_hmm_evalues


def cmd_refine(args):
    assert isdir(args.basedir)

    # Defaults: locate args.unclassified or set it
    if args.unclassified:
        unclass_fa = args.unclassified
    else:
        topname = basename(abspath(args.basedir).rstrip('/'))
        unclass_fa = join(args.basedir, topname + '-Unique.fasta')
    if not isfile(unclass_fa):
        raise ValueError("Can't find the unclassified sequence set "
                            "(%s does not exist)" % unclass_fa)

    all_fa = glob(join(args.basedir, '*.fasta'))
    all_fa.remove(unclass_fa)
    if args.exclude:
        if args.exclude in all_fa:
            all_fa.remove(args.exclude)
        elif args.exclude + '.fasta' in all_fa:
            all_fa.remove(args.exclude + '.fasta')
        else:
            logging.warn("Didn't find excluded profile %s in %s",
                         args.exclude, args.basedir)
    # Prerequisites: HMM profiles for each sequence set must be built
    for fa in all_fa:
        if not isfile(ext(fa, 'hmm')):
            raise ValueError("Unbuilt HMM for sequence set: %s" % fa)

    for this_fa in all_fa:
        this_hmm = ext(this_fa, 'hmm')
        logging.info("Refining %s...", this_fa)

        # Scan the seed sequences in `this_fa` w/ `this_hmm`
        # -> {seq_id: top_dom_eval}
        seq_evals = get_hmm_evalues(this_hmm, this_fa)

        # Starting with the lowest-scoring seed, do leave-one-out.
        # Leave at least 2 sequences in the profile, otherwise MAFFT will fail
        # (and a set of only 2 seqs should be merged into Unclassified anyway)
        for seq_id, evalue in sorted(seq_evals.items()[:-2],
                                    key=lambda kv: kv[1]):
            # Copy all seqs in this_fa, except seq_id, to a temp file
            with tempfile.NamedTemporaryFile(suffix='.fa') as tempseqfile:
                for rec in SeqIO.parse(this_fa, 'fasta'):
                    if rec.id == seq_id:
                        # Put this seq in a separate file
                        SeqIO.write(rec, tempseqfile.name + '.seed', 'fasta')
                    else:
                        SeqIO.write(rec, tempseqfile, 'fasta')
                tempseqfile.flush()
                # Align the "lean" profile w/ MAFFT (as usual)
                sh("mafft --quiet --amino --reorder --maxiterate 1000 "
                   "--genafpair --ep 0.123 %s > %s.seq"
                   % (tempseqfile.name, tempseqfile.name))
                # Build the "lean" HMM w/ hmmbuild
                SeqIO.convert(tempseqfile.name + '.seq', 'fasta',
                              tempseqfile.name + '.stk', 'stockholm')
                sh('hmmbuild %s %s > /dev/null'
                   % (tempseqfile.name + '.hmm', tempseqfile.name + '.stk'))
                os.remove(tempseqfile.name + '.seq')
                os.remove(tempseqfile.name + '.stk')

                # Scan the dropped seed seq w/ the "lean" HMM
                seed_eval = get_hmm_evalues(tempseqfile.name + '.hmm',
                                            tempseqfile.name + '.seed'
                                           )[seq_id]
                # Scan the unclassified seqs w/ the "lean" HMM
                # Keep the top evalue only; that's our threshold
                db_evals = get_hmm_evalues(tempseqfile.name + '.hmm',
                                           unclass_fa)
                try:
                    out_eval = max(db_evals.values())
                except ValueError:
                    # No HMM hits!
                    logging.warn("No unclassified hits for %s!", this_fa)
                    break

                if out_eval <= seed_eval:
                    if args.dry_run:
                        logging.info("* drop %s (evalue %s vs. %s)",
                                     seq_id, seed_eval, out_eval)
                    else:
                        # TODO:
                        # - append the dropped seq to unclass_fa
                        # - replace this_fa w/ the "lean" profile
                        # - Logging: print the name of the dropped seq
                        logging.info(
                            "Dropped %s from %s into %s (evalue %s vs. %s)",
                            seq_id, this_fa, unclass_fa, seed_eval, out_eval)
                else:
                    # Since the seeds are sorted by evalue, if this sequence is a
                    # keeper, then so are the rest. On to the next profile!
                    if args.dry_run:
                        logging.info("* keep %s (evalue %s vs. %s)",
                                     seq_id, seed_eval, out_eval)
                    break

