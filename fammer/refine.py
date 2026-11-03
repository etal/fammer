"Refine: Leave-one-out validation and refinement of each profile."

import logging
import os
import re
import shutil
import tempfile
from glob import glob
from os.path import abspath, basename, isdir, isfile, join, samefile

from Bio import SeqIO

from ._share import (
    MEMBER_EVALUE,
    UsageError,
    as_target,
    get_hmm_scores,
    is_catch_all,
    leave_one_out_hmms,
    member_names,
    write_fasta,
)
from .tasks import ext, noext

logger = logging.getLogger(__name__)

# refine leaves at least this many sequences in a profile.
MIN_MEMBERS = 2


def cmd_refine(args):
    if not isdir(args.basedir):
        raise UsageError(f"Not a directory: {args.basedir}")

    # Defaults: locate args.unclassified or set it
    if args.unclassified:
        unclass_fa = args.unclassified
    else:
        topname = basename(abspath(args.basedir).rstrip("/"))
        unclass_fa = join(args.basedir, topname + "-Unique.fasta")
    if not isfile(unclass_fa):
        raise UsageError(
            f"Can't find the unclassified sequence set ({unclass_fa} does not exist)"
        )

    # The subfamily profiles directly in basedir, sorted so that the result
    # does not depend on the order the file system lists them in
    profiles = []
    excluded = False
    for fasta in sorted(glob(join(args.basedir, "*.fasta"))):
        if samefile(fasta, unclass_fa) or is_catch_all(basename(noext(fasta))):
            continue
        if args.exclude in (fasta, noext(fasta), basename(noext(fasta))):
            excluded = True
            continue
        profiles.append(fasta)
    if args.exclude and not excluded:
        logger.warning(
            "Didn't find excluded profile %s in %s", args.exclude, args.basedir
        )
    # Prerequisites, checked before any file changes: each profile's alignment,
    # from the last build, holds exactly the sequences of its FASTA file, and
    # none of them is also an unclassified sequence
    with open(unclass_fa) as handle:
        unclassified = [as_target(rec) for rec in SeqIO.parse(handle, "fasta")]
    unclassified_names = {rec.id for rec in unclassified}
    for fasta in profiles:
        aln = ext(fasta, "aln")
        names = member_names(fasta)
        if not isfile(aln) or aligned_names(aln) != names:
            raise UsageError(
                f"{aln} does not hold exactly the sequences of {fasta}; "
                "run fammer build first"
            )
        both = names & unclassified_names
        if both:
            raise UsageError(
                f"Sequences in both {fasta} and {unclass_fa}: {', '.join(sorted(both))}"
            )

    for fasta in profiles:
        logger.info("Refining %s ...", fasta)
        dropped = refine_profile(fasta, unclassified)
        if not dropped:
            continue
        # Later profiles are compared with the unclassified set as it is after
        # the move, whether or not the files are changed
        unclassified.extend(dropped)
        if args.dry_run:
            continue
        move_records(fasta, unclass_fa, {rec.id for rec in dropped})
        logger.info(
            "Moved %d sequences from %s into %s; run fammer build to realign them",
            len(dropped),
            fasta,
            unclass_fa,
        )


def refine_profile(fasta, unclassified):
    """Select the sequences of a profile that do not score above unclassified ones.

    Each sequence of the profile is held out in turn: an HMM is built from the
    profile's alignment without its row, and the held-out sequence and the
    unclassified sequences are scored against it, in bits, counting the hits
    at MEMBER_EVALUE or better. A sequence fails if its own score is no higher
    than the best unclassified score, or if the HMM does not hit it while it
    hits an unclassified sequence. The lowest-scoring sequence that fails is
    dropped and joins the unclassified sequences, and the remaining sequences
    are tested again, until none fails or MIN_MEMBERS remain.

    Returns the dropped sequences, as ungapped records, in the order they were
    dropped.
    """
    with SeqIO.parse(ext(fasta, "aln"), "clustal") as parsed:
        rows = list(parsed)
    if len(rows) <= MIN_MEMBERS:
        logger.warning(
            "Skipping %s: it has %d sequences, and refine leaves at least %d",
            fasta,
            len(rows),
            MIN_MEMBERS,
        )
        return []
    targets = list(unclassified)
    dropped = []
    with tempfile.TemporaryDirectory(prefix="fammer-") as workdir:
        seqfname = join(workdir, "targets.fasta")
        while True:
            scores = {}
            by_name = {row.id: row for row in rows}
            for name, hmm in leave_one_out_hmms(rows, by_name.keys(), workdir):
                write_fasta([as_target(by_name[name]), *targets], seqfname)
                hits = get_hmm_scores(hmm, seqfname)
                own = hits.pop(name, None)
                scores[name] = (own, max(hits.values(), default=None))
            failing = [
                name
                for name, (own, best_out) in scores.items()
                if best_out is not None and (own is None or own <= best_out)
            ]
            if not failing:
                log_closest_keeper(fasta, scores)
                break
            worst = min(failing, key=lambda name: own_score(scores[name][0]))
            own, best_out = scores[worst]
            logger.info("* drop %s (%s vs. %.1f bits)", worst, describe(own), best_out)
            record = as_target(by_name[worst])
            rows = [row for row in rows if row.id != worst]
            dropped.append(record)
            targets.append(record)
            if len(rows) <= MIN_MEMBERS:
                logger.warning(
                    "Stopped refining %s at %d sequences", fasta, MIN_MEMBERS
                )
                break
    return dropped


def log_closest_keeper(fasta, scores):
    """Log the kept sequence of lowest score, or that nothing unclassified scored."""
    compared = {name: pair for name, pair in scores.items() if pair[1] is not None}
    if not compared:
        logger.warning(
            "No unclassified sequence scores at e-value %s or better against %s",
            MEMBER_EVALUE,
            fasta,
        )
        return
    name = min(compared, key=lambda name: own_score(compared[name][0]))
    own, best_out = compared[name]
    logger.info("* keep %s (%s vs. %.1f bits)", name, describe(own), best_out)


def own_score(score):
    """A held-out score for ordering, where no hit is the lowest."""
    return float("-inf") if score is None else score


def describe(score):
    """A held-out score for a log line."""
    return "no hit" if score is None else f"{score:.1f}"


def aligned_names(aln):
    """The identifiers of the rows of a Clustal alignment."""
    with SeqIO.parse(aln, "clustal") as parsed:
        return {row.id for row in parsed}


def move_records(fasta, dest, names):
    """Move the FASTA records with the given identifiers from one file to another.

    The records are appended to dest in the order and with the text they have
    in fasta, and then fasta is rewritten without them, keeping its
    permissions; every other line of either file is left as it was, with its
    line ending. A line ending that has to be added, after the last line of
    dest or of a moved record, is the one dest uses, or the one fasta uses
    where dest has none; a file with any CRLF line ending counts as CRLF.
    """
    with open(fasta, newline="") as handle:
        text = handle.read()
    with open(dest, newline="") as handle:
        dest_text = handle.read()
    newline = line_ending(dest_text, line_ending(text, "\n"))
    head, *records = re.split(r"(?m)^(?=>)", text)
    moved = [rec for rec in records if record_name(rec) in names]
    kept = [rec for rec in records if record_name(rec) not in names]
    # Write the rewritten profile to a temporary file first, so that a failure
    # while writing it changes neither file, and replace the profile with it
    # only after the records are appended; on any failure the temporary file
    # is removed
    rewritten = tempfile.NamedTemporaryFile(
        "w", dir=os.path.dirname(abspath(fasta)), delete=False, newline=""
    )
    try:
        with rewritten:
            rewritten.write(head + "".join(kept))
        shutil.copymode(fasta, rewritten.name)
        with open(dest, "a", newline="") as handle:
            if dest_text and not dest_text.endswith("\n"):
                handle.write(newline)
            handle.writelines(
                rec if rec.endswith("\n") else rec + newline for rec in moved
            )
        os.replace(rewritten.name, fasta)
    except BaseException:
        os.unlink(rewritten.name)
        raise


def line_ending(text, default):
    """The line ending that text uses, CRLF or LF, or default if it has none."""
    if "\r\n" in text:
        return "\r\n"
    return "\n" if "\n" in text else default


def record_name(text):
    """The identifier of a FASTA record's text, as Biopython reads it."""
    words = text[1:].split("\n", 1)[0].split(None, 1)
    return words[0] if words else ""
