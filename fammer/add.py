"Add sequences from a scanned database to the profile tree."

import functools
import logging
import os
import tempfile
from collections import defaultdict
from contextlib import closing
from os.path import basename, commonpath, dirname, isdir, isfile, join

from Bio import SeqIO
from Bio.Seq import Seq

from ._share import (
    BORDERLINE_MARGIN,
    MEMBER_EVALUE,
    UsageError,
    as_target,
    best_hits,
    get_hmm_rows,
    get_hmm_scores,
    is_catch_all,
    is_consensus,
    leave_one_out_hmms,
    member_names,
    ungapped,
    write_fasta,
)
from .tasks import ends_with_newline, ext, noext, sh

logger = logging.getLogger(__name__)

# A profile of at most this many sequences scores its weakest member held out;
# see weakest_member_score.
LEAVE_ONE_OUT_MAX = 10


def cmd_add(args):
    if not isdir(args.basedir):
        raise UsageError(f"Not a directory: {args.basedir}")
    hmmdb = args.basedir.rstrip("/") + "_all.hmm"
    if not isfile(hmmdb):
        raise UsageError(f"No HMM profile database at {hmmdb}")
    if not isfile(args.sequence_db):
        raise UsageError(f"No sequence file at {args.sequence_db}")
    if not 0 < args.min_length <= 1:
        raise UsageError(
            "Minimum hit length (-L, --min-length) must be between 0 and 1."
        )

    # Traverse the directory tree -- if there's an unbuilt HMM that we'll need,
    # best to identify it now and raise an error
    profiles = {}
    for fasta, hmm in walk_profiles(args.basedir):
        name = basename(noext(fasta))
        if name in profiles:
            raise UsageError(
                f"Two profiles are named {name}: {profiles[name][0]} and {fasta}"
            )
        profiles[name] = (fasta, hmm)
    names_of = {name: member_names(fasta) for name, (fasta, _hmm) in profiles.items()}
    # A sequence already in any profile is not added again, to that profile or
    # to another, so a repeated run adds each sequence once
    present = set().union(*names_of.values())
    # Scan the target database with every profile of the tree, as cmd_scan
    # does, but assign each sequence only among the profiles that have a FASTA
    # file to receive it
    rows = get_hmm_rows(hmmdb, args.sequence_db, args.evalue)
    if not rows:
        raise UsageError(f"No profile hits in {args.sequence_db}")
    rows = [row for row in rows if row[0] not in present and row[1] in profiles]
    with tempfile.TemporaryDirectory(prefix="fammer-") as workdir:
        threshold = acceptance_thresholds(profiles, workdir)
        dbhits = best_hits(rows)
        # A hit that its profile does not clearly separate is reported for the
        # curator to examine, and not written
        borderline = borderline_hits(dbhits, rows, threshold)
        report_borderline(borderline)
        dbhits = {s: hit for s, hit in dbhits.items() if s not in borderline}
        appended = set()
        with closing(SeqIO.index(args.sequence_db, "fasta")) as index:
            for name, (fasta, hmm) in profiles.items():
                logger.info("Processing %s ...", noext(fasta))
                hits = list(filter_profile(dbhits, name))
                if not hits:
                    continue
                min_score = threshold(name)
                if min_score is None:
                    warn_without_threshold(fasta)
                    continue
                appended |= add_hits(
                    fasta, hmm, hits, min_score, index, args.min_length, args.dry_run
                )
        if args.report_promotions:
            # A record that this run appended is left out, since in a dry run
            # it is not in the files, and so is a record that a subfamily also
            # holds
            in_subfamilies = set().union(
                *(names for name, names in names_of.items() if not is_catch_all(name))
            )
            report_promotions(
                profiles,
                hmmdb,
                threshold,
                args.evalue,
                args.min_length,
                skip=appended | in_subfamilies,
                workdir=workdir,
            )


def profile_directories(profiles):
    """The directory of each profile's FASTA file, by profile name."""
    return {name: dirname(fasta) for name, (fasta, _hmm) in profiles.items()}


def warn_without_threshold(fasta):
    """Log that a profile is skipped for want of an acceptance threshold."""
    logger.warning(
        "Skipping %s: its HMM profile hits none of its own sequences "
        "at e-value %s or better, so there is no threshold to accept "
        "a hit; run fammer build if its alignment is out of date",
        fasta,
        MEMBER_EVALUE,
    )


def inside(directory, top):
    """True if a directory is top or lies below it, however either is spelled."""
    return commonpath([directory, top]) == commonpath([top])


def accepting_subfamily(rows, top, threshold, directories):
    """The subfamily profile at or below a directory that accepts a sequence.

    The rows are the hmmsearch table rows of one sequence, tuples of (sequence
    name, profile name, e-value, bit score); threshold maps a profile name to
    its acceptance threshold in bits, or to None if it has none, and
    directories maps a profile name to the directory of its FASTA file. Of the
    profiles that are not catch-alls, lie in the directory top or below it,
    and have a threshold that the sequence's bit score reaches, the one whose
    best domain has the smallest e-value is selected, as best_hits selects it.

    Returns (profile name, e-value, bit score), or None if no profile accepts
    the sequence.
    """
    accepted = [
        row
        for row in rows
        if not is_catch_all(row[1])
        and inside(directories[row[1]], top)
        and (min_score := threshold(row[1])) is not None
        and row[3] >= min_score
    ]
    return next(iter(best_hits(accepted).values()), None)


def rows_by_sequence(rows):
    """Group hmmsearch table rows by their sequence name, the first field."""
    rows_of = defaultdict(list)
    for row in rows:
        rows_of[row[0]].append(row)
    return rows_of


def borderline_hits(hits, rows, threshold):
    """The assigned hits that their profile does not clearly separate.

    The hits are {sequence name: (profile name, e-value, bit score)}, as
    best_hits assigns them from the rows, and threshold maps a profile name to
    its acceptance threshold in bits, or to None. A hit that reaches its
    profile's threshold is borderline when it reaches it by less than
    BORDERLINE_MARGIN bits, or when it does not lead the closest other profile
    among the rows by BORDERLINE_MARGIN bits: a weak best match, or two
    profiles that do not separate it. The other profile may score the hit
    higher, since best_hits assigns by e-value. A catch-all and a subfamily of
    its group are rivals like any two profiles. A hit below its threshold is
    not borderline; its profile rejects it.

    Returns {sequence name: description for the report}.
    """
    rows_of = rows_by_sequence(rows)
    borderline = {}
    for seqname, (profname, _evalue, score) in hits.items():
        min_score = threshold(profname)
        if min_score is None or score < min_score:
            continue
        rival = max(
            (row for row in rows_of[seqname] if row[1] != profname),
            key=lambda row: row[3],
            default=None,
        )
        if score - min_score < BORDERLINE_MARGIN or (
            rival is not None and score - rival[3] < BORDERLINE_MARGIN
        ):
            description = (
                f"{seqname} for {profname} ({score:.1f} vs. {min_score:.1f} bits"
            )
            if rival is not None:
                description += f"; {rival[1]} {rival[3]:.1f}"
            borderline[seqname] = description + ")"
    return borderline


def report_borderline(borderline):
    """Log each borderline hit as a warning, which --quiet does not hide."""
    for seqname in sorted(borderline):
        logger.warning("* borderline %s", borderline[seqname])


def add_hits(fasta, hmm, hits, min_score, index, min_length, dry_run):
    """Append the hits that a profile accepts to its FASTA file.

    The hits are tuples of (sequence name, bit score), and min_score is the
    profile's acceptance threshold in bits (see acceptance_thresholds). If
    dry_run is set, nothing is written.

    Returns the names of the sequences appended, or that would be.
    """
    hit_seqs = accepted_hits(hmm, hits, min_score, index, min_length)
    added = {record.id for record in hit_seqs}
    if not hit_seqs:
        return added

    # ENH: extend hits that are close to the edge of the profile (e.g. 5aa) to
    # fill sequence to the edges (no gaps)
    #   e.g. if hmmstart = 4, add another 3aa sliced from the full sequence
    #   to the aligned bit -- if only 2aa available, then take that much
    # - can I tell from the alignment? (# terminal gaps)

    if dry_run:
        logger.info("Would add %d hits to profile %s", len(hit_seqs), fasta)
        return added
    # Add the accepted hits to the FASTA profile. A profile whose last line
    # has no line ending gets one first, so the first header starts a line.
    logger.info("Adding %d hits to profile %s", len(hit_seqs), fasta)
    needs_newline = not ends_with_newline(fasta)
    with open(fasta, "a") as fafile:
        if needs_newline:
            fafile.write("\n")
        write_fasta(hit_seqs, fafile)
    return added


def accepted_hits(hmm, hits, min_score, index, min_length):
    """The hits that a profile accepts, each as its longest aligned domain.

    The hits are tuples of (sequence name, bit score), and index maps a
    sequence name to its record. A hit is accepted when its bit score is at
    least min_score and the longest domain that the HMM aligns covers at least
    the fraction min_length of the HMM's length.

    Returns a list of records, each the ungapped domain under the sequence
    name.
    """
    hits = list(filter_score(hits, min_score))
    if not hits:
        return []

    # Align the hit sequences (retrieved from the target DB)
    aln = get_hmm_alignment(hmm, hits, index)

    # Calculate minimum sequence length based on HMM profile length
    hmm_len = get_hmm_len(hmm)
    min_len_aa = int(hmm_len * min_length)

    accepted = []
    for acc, record in longest_domains(aln, min_len_aa):
        # Remove the "/123-456" suffix from aligned record keys
        record.id = acc
        # Remove gaps from the hit sequence. The per-residue posterior
        # probabilities from hmmsearch -A do not fit the ungapped
        # sequence, so drop them first.
        record.letter_annotations = {}
        record.seq = Seq(ungapped(record.seq))
        accepted.append(record)
    return accepted


def report_promotions(profiles, hmmdb, threshold, evalue, min_length, skip, workdir):
    """Report each catch-all member that a subfamily of its group would accept.

    The members of every catch-all, except the consensus records and the
    names in skip, are scanned with the profile database hmmdb, and each is
    assigned to the profile whose best domain has the smallest e-value, as
    best_hits selects it; a member whose best profile is a catch-all, as a
    member's own catch-all usually is, is assigned instead to the subfamily
    that accepting_subfamily selects at or below that catch-all's directory,
    if there is one. A member assigned to a subfamily in its own catch-all's
    directory or below it, which accepts it by threshold and by min_length as
    add_hits does, is reported with its score against that subfamily, the
    subfamily's threshold, and its score against the catch-all's HMM built
    without it, since against the catch-all's own HMM, built partly from it,
    it scores higher than it would as an unseen sequence. Nothing is moved:
    whether a member belongs to the subfamily is the curator's decision. A
    member that two catch-alls hold is considered only as a member of the
    first of them in the order of their names.
    """
    directories = profile_directories(profiles)
    source_of, records = {}, {}
    for name in sorted(profiles):
        if not is_catch_all(name):
            continue
        with SeqIO.parse(profiles[name][0], "fasta") as parsed:
            for record in parsed:
                if is_consensus(record.id) or record.id in skip or record.id in records:
                    continue
                source_of[record.id] = name
                records[record.id] = as_target(record)
    if not records:
        return
    seqfname = join(workdir, "catch-all-members.fasta")
    write_fasta(records.values(), seqfname)
    rows = [row for row in get_hmm_rows(hmmdb, seqfname, evalue) if row[1] in profiles]
    rows_of = rows_by_sequence(rows)
    hits_of = defaultdict(dict)
    for seqname, (profname, _evalue, score) in best_hits(rows).items():
        if is_catch_all(profname):
            selected = accepting_subfamily(
                rows_of[seqname], directories[profname], threshold, directories
            )
            if selected is None:
                continue
            profname, _evalue, score = selected
        if inside(directories[profname], directories[source_of[seqname]]):
            hits_of[profname][seqname] = score
    candidates = defaultdict(list)
    for dest in sorted(hits_of):
        fasta, hmm = profiles[dest]
        min_score = threshold(dest)
        if min_score is None:
            warn_without_threshold(fasta)
            continue
        hits = hits_of[dest]
        for record in accepted_hits(hmm, hits.items(), min_score, records, min_length):
            candidates[source_of[record.id]].append(
                (record.id, dest, hits[record.id], min_score)
            )
    for source in sorted(candidates):
        held_out = held_out_scores(
            *profiles[source],
            [records[seqname] for seqname, *_rest in candidates[source]],
            workdir,
        )
        for seqname, dest, score, min_score in candidates[source]:
            if seqname not in held_out:
                own = "not scored"
            elif held_out[seqname] is None:
                own = "no hit"
            else:
                own = f"{held_out[seqname]:.1f}"
            logger.info(
                "* candidate %s in %s for %s (%.1f vs. %.1f bits; %s without it %s)",
                seqname,
                source,
                dest,
                score,
                min_score,
                source,
                own,
            )


def held_out_scores(fasta, hmm, records, workdir):
    """Score each record against its profile's HMM built without it.

    The records are sequences of the profile; each that is a row of the
    profile's alignment is scored against an HMM built from the other rows,
    and each that is not, having been added since the last build, against the
    profile's own HMM, which was not built from it. A record is not scored
    when there is no such HMM: the profile has no alignment, or the record is
    its only row. Unlike weakest_member_score, every row is held out, however
    many the profile has. The scores are best-domain bit scores at
    MEMBER_EVALUE or better.

    Returns {record identifier: bit score, or None for no hit}, without the
    records that are not scored.
    """
    rows, _members = aligned_members(fasta)
    if not rows:
        return {}
    row_ids = {row.id for row in rows}
    seqfname = join(workdir, "held-out.fasta")
    write_fasta(records, seqfname)
    scores = {}
    if len(rows) > 1:
        held = {record.id for record in records} & row_ids
        for name, lean in leave_one_out_hmms(rows, held, workdir):
            scores[name] = get_hmm_scores(lean, seqfname).get(name)
    added = [record.id for record in records if record.id not in row_ids]
    if added:
        own = get_hmm_scores(hmm, seqfname)
        scores.update((name, own.get(name)) for name in added)
    return scores


def walk_profiles(topdir):
    """Iterate through paired FASTA and HMM profiles in the directory tree."""
    for dirpath, _subdirs, fnames in os.walk(topdir):
        fasta_fnames = [join(dirpath, fn) for fn in fnames if fn.endswith(".fasta")]
        for fasta in fasta_fnames:
            hmm_fname = ext(fasta, "hmm")
            if not isfile(hmm_fname):
                raise UsageError(
                    f"No .hmm profile for {fasta}; run fammer build --hmmer first"
                )
            yield fasta, hmm_fname


def acceptance_thresholds(profiles, workdir):
    """The bit score that a hit must reach to join each profile.

    The profiles map a profile name to its FASTA and HMM files. A profile's
    threshold is the higher of two scores: the score of its weakest member
    (see weakest_member_score), and the highest best-domain bit score, at
    MEMBER_EVALUE or better, that its HMM gives to a member of a subfamily in
    another directory. A sequence of a family that the tree does not hold can
    outscore the weakest member of a profile, and so can a member of another
    group; a hit must score at least as high as such members do, so that a
    family missing from the tree is not taken in by the profile whose weakest
    member scores lowest. The subfamilies of the profile's own directory are
    left out, since related subfamilies are not separated by one score: their
    members commonly outscore each other's weakest members. For the second
    score, the members of a subfamily are those that weakest_member_score
    scores, without its consensus records, and the members of catch-alls do
    not count. A catch-all's own consensus records do count for the first
    score, its weakest member. A profile whose weakest member has no score has
    no threshold, and its value is None.

    Returns a function of a profile name, which computes each threshold once;
    the members of all subfamilies are written once to a file in workdir,
    which each profile's HMM searches.
    """
    directories = profile_directories(profiles)
    directories_of = defaultdict(set)
    subfamily_members = {}
    for name, (fasta, _hmm) in sorted(profiles.items()):
        if is_catch_all(name):
            continue
        _rows, members = aligned_members(fasta)
        for record in members:
            if not is_consensus(record.id):
                directories_of[record.id].add(directories[name])
                subfamily_members.setdefault(record.id, record)
    seqfname = join(workdir, "subfamily-members.fasta")
    write_fasta(subfamily_members.values(), seqfname)

    @functools.cache
    def threshold(name):
        fasta, hmm = profiles[name]
        weakest = weakest_member_score(fasta, hmm)
        if weakest is None:
            return None
        others = [
            score
            for seqname, score in get_hmm_scores(hmm, seqfname).items()
            if directories[name] not in directories_of[seqname]
        ]
        return max([weakest, *others])

    return threshold


def aligned_members(fasta):
    """The alignment rows of a profile and its members, as sequences to search.

    The members are the rows of the profile's alignment (.aln), from which
    fammer build made its HMM, that are also in its FASTA file, without gaps;
    a profile without an alignment has no rows, and its members are the
    records of its FASTA file.

    Returns a tuple of the rows and the members.
    """
    alnfname = ext(fasta, "aln")
    if not isfile(alnfname):
        with SeqIO.parse(fasta, "fasta") as records:
            return [], [as_target(record) for record in records]
    names = member_names(fasta)
    with SeqIO.parse(alnfname, "clustal") as parsed:
        rows = list(parsed)
    return rows, [as_target(row) for row in rows if row.id in names]


def weakest_member_score(fasta, hmm):
    """The lowest best-domain bit score among a profile's own sequences.

    The sequences are its members as aligned_members finds them. A record
    that add appended since the last build is not a row of the alignment, so
    it cannot lower the score before the HMM is rebuilt with it.

    A profile of at most LEAVE_ONE_OUT_MAX sequences scores each against an HMM
    built from the alignment without its row, since against an HMM built
    partly from itself a sequence of a small profile scores far higher than an
    unseen member of the family does. A larger profile, or one without an
    alignment, scores them against its own HMM. Only the sequences hit at
    MEMBER_EVALUE or better count; without any, the result is None. The scores
    are bit scores, since a bit score, unlike an e-value, does not depend on
    the number of sequences searched, so the profile's own scores and those of
    a target are comparable.
    """
    rows, members = aligned_members(fasta)
    with tempfile.TemporaryDirectory(prefix="fammer-") as workdir:
        seqfname = join(workdir, "members.fasta")
        write_fasta(members, seqfname)
        if len(rows) < 2 or len(members) > LEAVE_ONE_OUT_MAX:
            scores = get_hmm_scores(hmm, seqfname)
        else:
            scores = {}
            held_out = {member.id for member in members}
            for name, lean in leave_one_out_hmms(rows, held_out, workdir):
                score = get_hmm_scores(lean, seqfname).get(name)
                if score is not None:
                    scores[name] = score
    return min(scores.values(), default=None)


def get_hmm_alignment(hmmfname, hits, index):
    """Align the hit sequences to the HMM profile with hmmsearch -A.

    Returns the Stockholm records, one per aligned domain, each named
    ``<sequence name>/<start>-<end>``.
    """
    fullseqs = [index[hit[0]] for hit in hits]
    with (
        tempfile.NamedTemporaryFile(mode="w+", suffix=".fa") as seqfile,
        tempfile.NamedTemporaryFile(mode="w+", suffix=".stk") as alnfile,
    ):
        write_fasta(fullseqs, seqfile)
        seqfile.flush()
        sh(f"hmmsearch --noali -A {alnfile.name} {hmmfname} {seqfile.name} > /dev/null")
        alnfile.seek(0)
        return list(SeqIO.parse(alnfile, "stockholm"))


def get_hmm_len(hmmfname):
    """Read the HMM profile length from an HMM file."""
    with open(hmmfname) as hmmfile:
        for line in hmmfile:
            if line.startswith("LENG"):
                key, length = line.split()
                assert key == "LENG"
                return int(length)


def filter_profile(hits, profile_name):
    """Select the hits matching the given profile name.

    Hits are a dict of: {sequence name: (best profile, e-value, bit score)}

    Generates tuples of: (sequence name, bit score)
    """
    for seqname, (profname, _evalue, score) in hits.items():
        if profname == profile_name:
            yield (seqname, score)


def longest_domains(records, min_len_aa):
    """Select the longest domain of each sequence that hmmsearch -A aligned.

    The records are named ``<sequence name>/<start>-<end>``, one per aligned
    domain, so a sequence with several domains appears several times. Domains
    shorter than min_len_aa are skipped; of two equally long domains, the
    first is kept.

    Generates tuples of: (sequence name, record)
    """
    longest = {}
    for record in records:
        # Sequence names may themselves contain "/"
        acc, offset = record.id.rsplit("/", 1)
        start, end = map(int, offset.split("-"))
        hit_len = end - start + 1  # HMMer uses 1-based numbering
        if hit_len < min_len_aa:
            continue
        if acc in longest and longest[acc][0] >= hit_len:
            continue
        longest[acc] = (hit_len, record)
    for acc, (_hit_len, record) in longest.items():
        yield acc, record


def filter_score(hits, min_score):
    """Skip hits whose bit score is below the given one."""
    for hit in hits:
        if hit[1] >= min_score:
            yield hit
