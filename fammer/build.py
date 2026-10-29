"Build alignments and profiles in a directory tree of sequence sets."

import logging
import os
import shutil
import subprocess
from glob import glob
from os.path import basename, isdir, isfile, join

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biofrills import alnutils, consensus

from . import tmalign
from ._share import UsageError, write_fasta
from .tasks import Task, ext, is_empty, noext, sh, which

logger = logging.getLogger(__name__)


class Result:
    """Output of taskify_subdirs. Collection of tasks.

    Contains references to the tasks for the main alignment, the PDB seed
    alignment, and, optionally, the corresponding HMMer .hmm profile and/or
    MAPGAPS .cma and .tpl.
    """

    def __init__(self, aln, pdbseq=None, hmm=None, cma=None, tpl=None):
        self.aln = aln
        self.pdbseq = pdbseq
        self.hmm = hmm
        self.cma = cma
        self.tpl = tpl

    def __str__(self):
        return str(self.aln)


def cmd_build(args):
    if not isdir(args.basedir):
        raise UsageError(f"Not a directory: {args.basedir}")
    base = args.basedir.rstrip("/")

    # Write the directory layout as a Newick tree
    if args.tree:
        tree_str = treeify_subdirs(base) + ";\n"
        with open(base + ".nwk", "w+") as outfile:
            outfile.write(tree_str)
        logger.info("Wrote %s.nwk", base)

    # Build all alignments, plus HMM .hmm or MAPGAPS .cma/.tpl
    base_result = taskify_subdirs(base, args.hmmer, args.mapgaps, args.pdb, 0)

    # HMM profiles
    if args.hmmer:
        T_base_all_hmm = Task(
            base + "_all.hmm",
            action=all_hmm,
            depends=base_result.hmm,
        )
        T_base_all_hmm.build()
        # Clean everything
        if args.clean:
            T_base_all_hmm.clean()

    # MAPGAPS profiles
    if args.mapgaps:
        T_mg_run_map = Task(
            base + ".mpa", action=mg_run_map, depends=[base_result.cma, base_result.tpl]
        )
        T_mg_run_map.build()
        if args.clean:
            T_mg_run_map.clean()

    if not args.hmmer and not args.mapgaps:
        # Just build the alignments
        base_result.aln.build()
        if args.clean:
            base_result.aln.clean()


# To build the HMM (and CMA) profiles


def taskify_subdirs(topdir, hmmer, mapgaps, use_pdb, level):
    """Recursively define tasks for a directory tree.

    Return a task to build the HMM profile for topdir.

    this.hmm
        this.aln
            this.fasta -- consenses from each group/family
                this/subdir1.hmm (recurse)
                this/subdir2.hmm (recurse)
                this/...
                this/[others.hmm]
                --> task: hmmemit from each of deez

    How do we know which sub-hmms to build?
        subdirs of this/ --> group names --> [this/groupname.hmm]
        this/*.fasta --> [this/family.hmm]

    If mapgaps is True, build CMA profiles and templates too.
    """
    this = topdir.rstrip("/")
    subdirs = list(
        filter(isdir, [join(topdir, sd) for sd in sorted(os.listdir(topdir))])
    )
    # Each subtask pair: (HMM, CMA)-making tasks for the subdir
    # Groups with their own families within -- recurse
    subtask_group_results = [
        taskify_subdirs(sd, hmmer, mapgaps, use_pdb, level + 1)
        for sd in sorted(subdirs)
    ]
    # Families / tips of the profile tree -- build from scratch
    subtask_family_results = []
    for subfa in sorted(glob(join(topdir, "*.fasta"))):
        # Skip the group sequence sets, they're already covered
        if noext(subfa) in subdirs:
            continue
        subresult = Result(
            Task(
                ext(subfa, "aln"),
                action=align_fasta_mafft,
                depends=subfa,
                cleans=ext(subfa, "seq"),
            )
        )
        if hmmer:
            subresult.hmm = Task(
                ext(subfa, "hmm"),
                action=aln2hmm,
                depends=subresult.aln,
                cleans=ext(subfa, "stk"),
            )
        if mapgaps:
            subresult.cma = Task(
                ext(subfa, "cma"),
                action=mg_aln2cma,
                kwargs={"level": level + 1},
                depends=subresult.aln,
                cleans=[ext(subfa, "cons.cma"), ext(subfa, "cons_iron.cma")],
            )
        subtask_family_results.append(subresult)

    subtask_family_results.sort(key=str)  # Needed? We sort FASTAs above

    # Structural alignment of PDBs in this dir; reuse subfamily PDB alignments
    these_pdbs = sorted(glob(join(topdir, "*.pdb")))
    sub_pdb_seqs = [sgr.pdbseq for sgr in subtask_group_results]
    task_pdbseq = Task(
        this + ".pdb.seq",
        action=align_pdbs,
        kwargs={"use_pdb": use_pdb},
        depends=these_pdbs + sub_pdb_seqs,
    )

    # Aggregate those profile consensus sequences & make a meta-profile
    result = Result(
        Task(
            this + ".aln",
            action=align_profiles,
            kwargs={"use_pdb": use_pdb},
            depends=[r.aln for r in (subtask_group_results + subtask_family_results)]
            + [task_pdbseq],
            cleans=ext(
                map(str, subtask_group_results + subtask_family_results), "cons.seq"
            )
            + [this + ".families.fa", this + ".families.seq", this + ".seq"],
        ),
        pdbseq=task_pdbseq,
    )
    if hmmer:
        result.hmm = Task(
            this + ".hmm",
            action=aln2hmm,
            depends=(
                [result.aln]
                + [sgr.hmm for sgr in subtask_group_results]
                + [sfr.hmm for sfr in subtask_family_results]
            ),
            cleans=this + ".stk",
        )
    if mapgaps:
        result.tpl = Task(
            this + ".tpl",
            action=mg_aln2cma,
            kwargs={"level": level},
            depends=result.aln,
            cleans=[this + ".cons.cma", this + ".cons_iron.cma"],
        )
        result.cma = Task(
            this + ".cma",
            action=mg_cat_cma,
            depends=[
                task for sgr in subtask_group_results for task in (sgr.tpl, sgr.cma)
            ]
            + [sfr.cma for sfr in subtask_family_results],
        )

    return result


# Actions


def align_fasta_mafft(task):
    """Align a FASTA file with MAFFT. Clustal output."""
    seq = ext(task.depends[0], "seq")
    sh(
        "mafft --quiet --amino --reorder --maxiterate 1000 "
        f"--genafpair --ep 0.123 {task.depends[0]} > {seq}"
    )
    # Convert FASTA to "pressed" (single-row) Clustal
    with SeqIO.parse(seq, "fasta") as parsed:
        records = list(parsed)
    # Check for 'X' characters in the sequences -- these cause problems
    for rec in records:
        if "X" in str(rec.seq):
            logger.warning("Sequence %r contains unknown residue X", rec.id)
    max_id_len = max(len(r.id) for r in records)
    with open(task.target, "w+") as outfile:
        outfile.write("CLUSTAL X (-like) multiple sequence alignment\n\n")
        outfile.writelines(
            [f"{rec.id.ljust(max_id_len)} {rec.seq}\n" for rec in records]
        )


def align_pdbs(task, use_pdb=None):
    """Create a structure-based sequence alignment from PDB files.

    Inputs are PDB files and FASTA alignments (of previously aligned PDBs).
    """
    if not use_pdb or not task.depends:
        # Just touch the '.pdb.seq' file; don't use TM-align
        with open(task.target, "a"):
            return

    pdbs = []
    sub_pdb_seqs = []
    for elem in task.depends:
        if str(elem).endswith(".pdb"):
            pdbs.append(elem)
        else:
            sub_pdb_seqs.append(str(elem))

    # Scan existing PDB alignments to choose a reference PDB from each, and
    # pass each non-empty one on to MAFFT as a seed.
    pdbseedfnames = []
    for sub_pdb_fname in sub_pdb_seqs:
        if is_empty(sub_pdb_fname):
            # Empty: that directory has fewer than two structures to align, so
            # take its structures here
            #   ENH - recursively
            sub_dir = sub_pdb_fname.removesuffix(".pdb.seq")
            sub_pdbs = sorted(glob(join(sub_dir, "*.pdb")))
            if sub_pdbs:
                logger.info(
                    "Picked up %d neglected PDBs: %s", len(sub_pdbs), " ".join(sub_pdbs)
                )
                pdbs.extend(sub_pdbs)
            continue
        pdbseedfnames.append(sub_pdb_fname)
        best_tmscore = -1
        best_pdb = None
        with SeqIO.parse(sub_pdb_fname, "fasta") as parsed:
            for rec in parsed:
                # The rows of a MUSTANG seed alignment carry no TM-score.
                tmscore = tmalign.tmscore_from_description(rec.description)
                if tmscore is not None and tmscore > best_tmscore:
                    best_tmscore = tmscore
                    best_pdb = rec.id
        if best_pdb is None:
            logger.warning("Empty PDB alignment: %s", sub_pdb_fname)
        else:
            logger.info("Best PDB of %s: %s", sub_pdb_fname, best_pdb)
            pdbs.append(best_pdb)

    try:
        mustang_tmpfname = "_tmp_mustang.afasta"
        if len(pdbs) > 1 and which(["mustang"]):
            # Align PDBs with MUSTANG.
            subprocess.check_call(
                ["mustang", "-o", "_tmp_mustang", "-F", "fasta", "-s", "OFF", "-i"]
                + pdbs
            )
            pdbseedfnames.append(mustang_tmpfname)

        # This is where the magic happens.
        records = tmalign.align_structs(pdbs, pdbseedfnames)

    finally:
        if isfile(mustang_tmpfname):
            os.remove(mustang_tmpfname)

    write_fasta(records, task.target)


def align_profiles(task, use_pdb=None):
    """Align several FASTA files with MAFFT. Clustal output.

    Cleans: [depends].cons.seq, [target].families.fa, [target].families.seq,
        [target].seq
    """
    seeds, singles = [], []
    # PDB alignment -- include as a seed profile if requested
    subalignments, pdb_seed = task.depends[:-1], str(task.depends[-1])
    if use_pdb and not is_empty(pdb_seed):
        seeds.append(pdb_seed)
    else:
        logger.info("Empty PDB alignment: %s", pdb_seed)

    # Get subfamily and subgroup consensus sequences/profiles
    for subaln in subalignments:
        aln = AlignIO.read(str(subaln), "clustal")
        with open(ext(subaln, "cons.seq"), "w+") as outfile:
            outfile.write(f">{basename(noext(subaln))} consensus\n")
            cons_seq = consensus.consensus(aln, trim_ends=False, gap_threshold=0.6)
            if isdir(noext(subaln)):
                # Group profiles: include the subfamily consenses, too
                outfile.write(cons_seq + "\n")
                for record in aln:
                    outfile.write(f">{record.id}\n")
                    outfile.write(f"{record.seq}\n")
            else:
                # Ungapped family consensus sequences
                outfile.write(cons_seq.replace("-", "") + "\n")

    # Merge the sequences and profiles
    for subconsseq in ext(subalignments, "cons.seq"):
        if isdir(subconsseq[:-9]):
            # Group
            seeds.append(subconsseq)
        else:
            singles.append(subconsseq)
    famfa = ext(task.target, "families.fa")
    allseq = ext(task.target, "seq")
    if not singles and not seeds:
        raise UsageError(f"No .fasta files found to build {task.target}")
    # Written even when empty, so that MAFFT never reads a file left by an
    # earlier build; with seeds, MAFFT accepts an empty input.
    with open(famfa, "wb") as outfile:
        for single in singles:
            with open(single, "rb") as infile:
                shutil.copyfileobj(infile, outfile)
    # Align the families, with the group profiles and the requested PDB
    # alignment, if any, as seeds.
    #
    # G-INS-i (--globalpair) aligns the trimmed consensus rows end to end.
    # The 2012-2013 sources passed --globalgenafpair, which MAFFT 6.903 to
    # 7.123 computed as --globalpair and later releases reject.
    seed_args = "".join(f"--seed {s} " for s in seeds)
    sh(
        "mafft --quiet --amino --globalpair --maxiterate 1000 "
        f"{seed_args}{famfa} > {allseq}"
    )
    # Convert FASTA to "pressed" (single-row) Clustal
    with SeqIO.parse(allseq, "fasta") as parsed:
        records = [
            rec
            for rec in parsed
            # Drop PDB-derived sequences
            if "TMalign" not in rec.description
            and "TM-score" not in rec.description
            and not rec.id.endswith(".pdb")
        ]
    records = list(alnutils.remove_empty_cols(records))
    if seeds:
        # MAFFT prefixes seed alignments with '_seed_' -- get rid of that
        for rec in records:
            if rec.id.startswith("_seed_"):
                rec.id = rec.id[6:]
    try:
        max_id_len = max(len(r.id) for r in records)
    except ValueError as exc:
        # No records remain: MAFFT wrote none, or all were PDB-derived.
        inputs = " ".join(map(str, task.depends))
        raise ValueError(
            f"Profile alignment failed for {task.target}.\nInputs: {inputs}"
        ) from exc

    with open(task.target, "w+") as outfile:
        outfile.write("CLUSTAL X (-like) multiple sequence alignment\n\n")
        outfile.writelines(
            [f"{rec.id.ljust(max_id_len)} {rec.seq}\n" for rec in records]
        )


def aln2hmm(task):
    """Convert a Clustal alignment to an HMM profile.

    Cleans: .stk
    """
    stk = ext(task.depends[0], "stk")
    with open(str(task.depends[0])) as aln:
        SeqIO.convert(aln, "clustal", stk, "stockholm")
    sh(f"hmmbuild {task.target} {stk}")


def all_hmm(task):
    """Concatenate all HMM profiles into a database & press."""
    base = noext(task.depends[-1])
    sh(
        f"cat {base}.hmm `find {base}/ -name '*.hmm' | grep -v '_all.hmm'` "
        f"> {task.target}"
    )
    sh(f"hmmpress -f {task.target}")


# MAPGAPS actions


def mg_aln2cma(task, level=None):
    """Convert an alignment profile to CMA (or .tpl).

    Depends: .aln
    Cleans: .cons.cma, .cons_iron.cma
    """
    try:
        from biocma import biocma
    except ImportError as exc:
        raise ImportError(
            "build --mapgaps requires the Python package biocma 0.3, which is "
            "not yet on PyPI: install it from https://github.com/etal/biocma, "
            f"then install fammer with its mapgaps extra ({exc})"
        ) from exc
    base = noext(task.target)
    name = basename(base)
    # Add consensus back to the subfamily-consensus seq set (.aln)
    # to produce a CMA (.cons.cma)
    aln = AlignIO.read(str(task.depends[0]), "clustal")
    cons_rec = SeqRecord(
        Seq(consensus.consensus(aln, trim_ends=False, gap_threshold=0.6)),
        id=name,
        description=name + " consensus",
    )
    aln = MultipleSeqAlignment([cons_rec] + list(aln))
    # Tidy up the CMA
    cmaln = biocma.ChainMultiAlignment(aln, level=level)
    biocma.write([cmaln], task.target, do_iron=True)


def mg_cat_cma(task):
    """Concatenate subfamily MAPGAPS profiles.

    Depends: .cma of each subfamily, .tpl of this & sub-groups.
    """
    assert task.depends, "No CMA files were given"
    sh(f"cat {' '.join(map(str, task.depends))} > {task.target}")


def mg_run_map(task):
    """Build/compile the complete set of MAPGAPS profiles."""
    sh(f"run_map {noext(task.target)}")


# Tree


def treeify_subdirs(topdir):
    """Build a Newick string from the directory tree structure.

    Internal nodes = names of non-empty dirs
    External nodes = names of .fasta files
    """
    # Full paths of directories under topdir
    subdirs = list(
        filter(isdir, [join(topdir, sd) for sd in sorted(os.listdir(topdir))])
    )
    # Do internal nodes first, then the tips
    # Internal nodes: subtree Newick strings
    subtree_strs = [treeify_subdirs(sd) for sd in sorted(subdirs)]
    # Tips: string names (basename minus trailing .fasta)
    tip_names = [
        basename(fafname)[: -len(".fasta")]
        for fafname in glob(join(topdir, "*.fasta"))
        if noext(fafname) not in subdirs
    ]
    tip_names.sort(key=str)
    return f"({','.join(subtree_strs + tip_names)}){basename(topdir.rstrip('/'))}"
