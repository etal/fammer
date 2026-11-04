======
Fammer
======

Utilities for curating a hierarchical set of sequence profiles representing a
protein superfamily.

If you use this software in a publication, please cite the paper that describes it:

    Talevich, E. & Kannan, N. (2013) Structural and evolutionary adaptation of
    rhoptry kinases and pseudokinases, a family of coccidian virulence factors.
    *BMC Evolutionary Biology* 13:117.
    https://doi.org/10.1186/1471-2148-13-117

Freely distributed under the permissive BSD 2-clause license (see LICENSE).

Purpose and design
------------------

Fammer builds and curates a hierarchical set of sequence profiles for one protein superfamily: a profile for each family or subfamily, and above them a profile for each group and one for the superfamily, arranged as a directory tree.
The profile set has three uses.
It detects members of the superfamily in large sequence databases, such as the proteome of a newly sequenced genome.
It classifies each member into the most specific level of the hierarchy that fits it.
It aligns all members to one another through a template alignment of the consensus sequences of the profiles, as MAPGAPS does [Neuwald2009]_.
Such alignments are the input to comparative analyses of the sequence patterns that distinguish one subfamily from its relatives, such as CHAIN analysis and CladeCompare.

Fammer was written to partially automate the construction and curation of such profile sets, and was used to classify the rhoptry kinases of coccidian parasites into 42 subfamilies [Talevich2013]_.
Its design follows the iterative, partly manual procedure of [Kannan2007]_, and Fammer automates the repetitive steps of that procedure.
The procedures of these papers, and the rules below that follow from them, define what Fammer is meant to do.

- **A profile models closely related sequences.**
  Each family or subfamily set holds sequences that standard methods align accurately; the relationships between distant subfamilies are represented by the hierarchy and the alignments of consensus sequences above it, not by any one alignment.
- **The profile set covers the widespread subfamilies.**
  The method assumes that most of the phylogenetically widespread families and subfamilies of the superfamily have a profile, so that a sequence that fits none of them is unusual: a member of a sparsely populated, lineage-specific subgroup, or a divergent outlier.
  Where widespread subfamilies are missing, their members fall to the nearest profile or catch-all that the tree does hold, and classifications and acceptance rules can be judged only weakly on such a tree.
  The example tree ``profiles/SomeKinases`` is such a tree: it is a deliberately incomplete subset of a few groups of non-tyrosine kinases, which omits many of the families of the KinBase classification and possibly some subfamilies of the families it holds, so it serves as test data rather than as a benchmark.
- **A sequence belongs to the profile that scores it highest.**
  Related subfamilies are not separated by an absolute score; they are told apart by which profile scores a sequence highest.
  MAPGAPS classifies a sequence whose best match is a group profile, or the superfamily profile, at that level, as a member of a subgroup that no profile models yet, and a set of such sequences is a candidate for a new profile [Neuwald2009]_.
  ``add`` does not: it considers such a sequence for the best-matching profile among those that have a FASTA file (see `add`_).
- **A catch-all set is a profile of its group.**
  A file named ``*-Unclassified``, ``*-Unclass`` or ``*-Unique`` holds the members of a group that fit none of its subfamilies, together with any consensus records that make up its generic profile (see `Directory tree is the superfamily hierarchy`_).
  These members include orphan sequences and members of lineage-specific subfamilies too sparsely populated to have their own profile, which the curator collapsed into the catch-all.
  A sequence in a catch-all has a classification at the level of its group; it is not necessarily waiting for a subfamily.
  A catch-all corresponds to the "unique" set of [Talevich2013]_.
- **Related but irrelevant sequences are excluded by competing profiles.**
  A profile of the sequences to be excluded, such as the profile of other protein kinases that [Talevich2013]_ included beside the rhoptry kinase profiles, or the excluded profiles of MAPGAPS, takes those sequences as its best match.
  A score threshold alone cannot tell a member of a missing family from a divergent member of a present one.
- **The curator decides; Fammer proposes.**
  In [Talevich2013]_, subfamilies were defined from clusters of orthologs and trees of consensus sequences, outliers were moved to the "unique" set by inspection of conserved motifs, the "unique" set was reclassified with ``scan``, and each subfamily was validated against the "unique" set by ``refine``.
  The leave-one-out validation of ``refine`` follows [Hedlund2010]_.
  [Kannan2007]_ examined by hand every sequence with a borderline score, that is, a weak best match or a small difference between its best two models.
  A missed member costs one more round of curation; a wrong member changes the profile, and with it every later search.
  A command that writes to curated sequence sets should therefore report a borderline case for the curator to examine instead of writing it, and should be able to show its changes without making them.
  ``add`` reports borderline hits instead of adding them, ``refine`` marks its borderline drops in the log but still makes them, and both show their changes without making them with ``-n``.
- **Iteration converges.**
  Alternating ``build`` and ``add`` repeats the build, search and refine cycle of [Kannan2007]_; as a requirement of Fammer's own, the rounds should reach one that adds nothing, rather than grow without end.

.. [Talevich2013] Talevich E, Kannan N (2013) Structural and evolutionary adaptation of rhoptry kinases and pseudokinases, a family of coccidian virulence factors. *BMC Evolutionary Biology* 13:117.
   The procedure is in Methods, "Subfamily classification" (https://link.springer.com/article/10.1186/1471-2148-13-117#Sec17).
   The resulting HMM profile database, with the profile of other protein kinases, is Additional file 3 (https://link.springer.com/article/10.1186/1471-2148-13-117#Sec23).
.. [Neuwald2009] Neuwald AF (2009) Rapid detection, classification and accurate alignment of up to a million or more related protein sequences. *Bioinformatics* 25(15):1869–1875. https://doi.org/10.1093/bioinformatics/btp342
.. [Kannan2007] Kannan N, Taylor SS, Zhai Y, Venter JC, Manning G (2007) Structural and functional diversity of the microbial kinome. *PLoS Biology* 5(3):e17. https://doi.org/10.1371/journal.pbio.0050017
.. [Hedlund2010] Hedlund J, Jörnvall H, Persson B (2010) Subdivision of the MDR superfamily of medium-chain dehydrogenases/reductases through iterative hidden Markov model refinement. *BMC Bioinformatics* 11:534. https://doi.org/10.1186/1471-2105-11-534


Installation
------------

The installation provides a Python library, ``fammer``, and two commands, ``fammer`` and ``fammer-tmalign``.
Fammer requires Python 3.12 or later.
This Python 3 version of Fammer is not yet released on PyPI, where the package named ``fammer`` is the earlier Python 2 release, 0.2, so install it from the source code::

    git clone https://github.com/etal/fammer.git
    cd fammer

Fammer runs the following programs, which must be on ``PATH`` under the names given:

- MAFFT_ (``mafft``), for ``build``, ``fammer-tmalign``, and ``cluster`` with an unaligned sequence file.
- HMMER_ 3 (``hmmbuild``, ``hmmpress`` and ``hmmsearch``), for ``build --hmmer``, ``scan``, ``add`` and ``refine``.
- TM-align_ (``TMalign``), for ``fammer-tmalign``, and for ``build`` of a tree that holds PDB files unless ``build`` is given ``--no-pdb``.
- MUSTANG_ (``mustang``), optional: where it is installed, ``build`` also aligns the structures of a tree with it (see `build`_).
- FastTree_ (``fasttree`` or ``FastTree``), for ``cluster``, unless ``cluster`` is given a tree with ``--tree``.
- MAPGAPS_ (``run_map``), optional, for ``build --mapgaps``.

With pixi
`````````

pixi_ installs Python, Fammer, its Python dependencies and all of these programs except MAPGAPS into an environment inside the source directory, at the versions recorded in ``pixi.lock``.
The lock file covers Linux (x86-64) and macOS (Apple silicon); on other systems, install with pip.
From the source directory::

    pixi install
    pixi run fammer --help

``pixi run`` runs one command in the environment, and ``pixi shell`` starts a shell in which ``fammer`` and the programs it runs are on ``PATH``.
The environment also holds the development tools (see `Development`_).

With pip
````````

Install Fammer and its Python dependencies, Biopython_, biofrills_, networkx_ and NumPy_, with pip into a virtual environment or a conda environment of Python 3.12 or later; the system Python of Debian and Ubuntu refuses ``pip install`` outside a virtual environment.
Fammer requires biofrills 0.4, which is not yet released on PyPI, so first install biofrills from its Git repository, at the commit that ``[tool.uv.sources]`` in ``pyproject.toml`` names; this needs a ``git`` executable.
Then, from the source directory::

    pip install git+https://github.com/etal/biofrills@0a3f0efeaee20caff07c955ea71ca5ea6173010f
    pip install .

Then install the external programs.
conda-forge and bioconda provide all of them except MAPGAPS; in a conda environment::

    conda install -c conda-forge -c bioconda mafft hmmer tmalign mustang fasttree

On Debian or Ubuntu, the system packages provide them::

    sudo apt-get install mafft hmmer tm-align mustang fasttree

MAPGAPS support
```````````````

The ``--mapgaps`` option of ``build`` also needs biocma_ 0.3, which is likewise not yet on PyPI and is not installed by default.
The pixi environment includes it; with pip, install it from its repository at the pinned commit, and then install Fammer with its ``mapgaps`` extra::

    pip install git+https://github.com/etal/biocma@307d607b1739f647610918113176bc7dcd71dbbc
    pip install '.[mapgaps]'

No package manager provides MAPGAPS itself.
Version 1.0.1 is distributed from the MAPGAPS_ site as a tarball of programs compiled for Linux x86-64; unpack it and put its ``bin`` directory on ``PATH``.

.. _MAFFT: https://mafft.cbrc.jp/alignment/software/
.. _HMMER: https://hmmer.org/
.. _TM-align: https://zhanggroup.org/TM-align/
.. _MUSTANG: https://lcb.infotech.monash.edu/mustang/
.. _FastTree: https://morgannprice.github.io/fasttree/
.. _MAPGAPS: https://mapgaps.igs.umaryland.edu/
.. _Biopython: https://biopython.org/
.. _biofrills: https://github.com/etal/biofrills
.. _biocma: https://github.com/etal/biocma
.. _networkx: https://networkx.org/
.. _NumPy: https://numpy.org/


Basic usage
-----------

Global options of ``fammer``:

  ``-h``, ``--help``
      Show a help message and basic usage.
  ``-q``, ``--quiet``
      Don't print status messages, only warnings and errors.

Sub-commands:

    `build`_
        Build profiles from a given directory tree.
    `update-fasta`_
        Replace original FASTA sequence sets with the ungapped sequences from
        the corresponding alignment (.aln) files, sorted by decreasing length.
    `scan`_
        Scan and classify input sequences using a set of profiles.
    `add`_
        Scan a target database with the given HMM profile set.  Add hits that
        meet acceptance thresholds to the profile FASTA files, and report
        borderline hits instead of adding them.
    `refine`_
        Leave-one-out validation of the subfamily profiles in one directory.
        Move the sequences that score no higher than an unclassified sequence
        into the unclassified set.
    `cluster`_
        Split a sequence set into clusters (based on phylogeny).

Each sub-command takes ``-h`` for its own options.
The separate command `fammer-tmalign`_ aligns protein structures.

To try Fammer on the example data, build the alignments and HMM profiles of a copy of the example tree ``profiles/SomeKinases``, then classify the predicted proteome of *Toxoplasma gondii* ME49 in ``test/TgME49.pep`` with them::

    cp -R profiles/SomeKinases /tmp/
    fammer build --hmmer /tmp/SomeKinases
    fammer scan /tmp/SomeKinases test/TgME49.pep

In the pixi environment, run each command as ``pixi run fammer ...``, or in ``pixi shell``.


Commands
--------

build
`````

Construct a profile database from a directory tree of family profile alignments.

Given a directory tree under ``Superfamily/`` laid out as in `Directory tree is the superfamily hierarchy`_, run ``fammer build Superfamily`` to align each sequence file with MAFFT and then, recursively upward, the consensus sequences of the subfamilies in each directory::

    Superfamily/
        Group1/
            subfam1_Group1.fasta
            subfam1_Group1.aln
            subfam2_Group1.fasta
            subfam2_Group1.aln
            subfam3_Group1.fasta
            subfam3_Group1.aln
            ...
            Group1-Unclassified.fasta
            Group1-Unclassified.aln
        Group1.aln
        ...
        Superfamily-Unclassified.fasta
        Superfamily-Unclassified.aln
    Superfamily.aln

The alignments are in un-wrapped Clustal format.

You can manually adjust the alignments and rebuild, if desired, perhaps
iteratively. Only the "parent" family alignments will be rebuilt as needed, e.g.
if ``subfam1_Group1.aln`` is edited, then only ``Group1.aln`` and
``Superfamily.aln`` will be rebuilt the next time ``fammer build
Superfamily`` is called because the consensus sequences that constitute those
alignments may have changed. (It's like Make.)

Finally, use the option ``--hmmer`` to build profiles::

    Superfamily/
        Group1/
            subfam1_Group1.fasta
            subfam1_Group1.aln
            subfam2_Group1.hmm
            ...
        Group1.aln
        Group1.hmm
        ...
    Superfamily.aln
    Superfamily.hmm
    Superfamily_all.hmm     # concatenated profiles
    Superfamily_all.hmm.{h3f,h3i,h3m,h3p}   # indexes from hmmpress

The ``--mapgaps`` option works similarly, if MAPGAPS and the Python library biocma_ are installed (see `MAPGAPS support`_).

The ``--clean`` option can be included with any of the above commands to remove
intermediate files.

If the directory tree holds PDB files, ``fammer build`` aligns them with TM-align_, which must then be installed, unless given the ``--no-pdb`` option.
It first creates a structural alignment of the PDBs in each directory, then uses that alignment as the seed for higher-up alignments::

    Superfamily/
        Group1/
            subfam1_Group1.fasta
            subfam1_Group1.aln
            1ATP.pdb
            1O6K.pdb
            3C4X.pdb
            ...
        Group1.pdb.seq  # Alignment of 1ATP, 1O6K, 3C4X
        Group1.aln
        ...
    Superfamily.aln

In this example, the alignment generated by aligning the structures 1ATP, 1O6K
and 3C4X is passed to MAFFT as a seed for ``Group1.aln``, along with the
unaligned consensus sequences of each subfamily of Group1 (subfam1, subfam2,
...). The seed sequences are removed from Group1.aln after the alignment of
consensus sequences is completed. This can help correctly align the more
divergent families and groups to each other.

The structures of each directory are aligned as `fammer-tmalign`_ aligns them, from the pairwise alignments of TM-align.
Where MUSTANG_ is installed, ``build`` also aligns the structures of each directory that has at least two with MUSTANG, and passes that multiple alignment to MAFFT as one more seed beside the pairwise alignments.
A directory with exactly two structures to align uses no seed, so its structural alignment is the TM-align alignment of the pair (see `fammer-tmalign`_).

A directory with fewer than two structures to align leaves its ``.pdb.seq`` file empty, and the parent directory then aligns the structure files found directly in that directory with its own.
Structures further down are not passed on: if the parent directory has no other structure either, a single structure is not used at any level.

For nested directory trees, the option ``--tree`` generates a Newick file
representing the structure of the directory tree. A tree based on the above
examples would look something like this (ignoring whitespace), created as
``Superfamily.nwk``::

    ((subfam1_Group1, subfam2_Group1, subfam3_Group1,
      Group1-Unclassified)Group1,
     (subfam1_Group2, subfam2_Group2, subfam3_Group2,
      Group2-Unclassified)Group2,
     (subfam1_Group3, subfam2_Group3, subfam3_Group3,
      Group3-Unclassified)Group3,
     Superfamily-Unclassified)Superfamily;

This tree could be passed to RAxML as a constraint tree in an effort to identify
deeper subfamilies, for example.


update-fasta
````````````

Convert the contents of the ``.aln`` sequence alignment files back to unaligned
FASTA format, overwriting the corresponding ``.fasta`` files.

After initially building a tree of sequence alignments, you might edit the
Clustal alignments, deleting spurious sequences or trimming the alignment to the
edges of a conserved domain. With ``update-fasta``, you update the contents of
the unaligned sequence files to match the ``.aln`` files.

The next step is usually to either (a) do some sequence processing unrelated to
fammer, e.g.  clustering, or (b) realign everything. Since you've presumably
removed some junk from the input sequences, the resulting alignments may be
better.


scan
````

Scan/search a set of sequences (FASTA) with the HMM profile database and assign a
classification to each hit.

The profile database is either the ``_all.hmm`` file that ``build --hmmer`` writes, or the top directory of the tree, in which case ``scan`` reads the ``_all.hmm`` file beside it: ``fammer scan Superfamily seqs.fasta`` reads ``Superfamily_all.hmm``.
A sequence is a hit when ``hmmsearch`` reports it at an E-value of at most ``-E``, by default 0.0001, and it is classified by the profile whose best domain gives it the smallest E-value.

This is essentially a set of wrappers to process the output of ``hmmsearch``,
simplifying the results for common use cases. The three output forms are:

    **summary** (default):
        Print two formated columns for each profile in the given HMM profile
        database that matched at least one hit: the name of the profile and the
        number of hits for which it was the best match.
    **table** (``--table``):
        For each sequence in the target sequence set that matched a profile in
        the HMM profile database, print the sequence ID/accession and the name
        of the best-matching profile, separated by a tab character.
    **sequence sets** (``--seqsets``):
        For each profile and matching sequence set (as they'd appear in summary
        output), write a file containing the matching sequences. The output
        files, written to the current directory, are named after the sequence
        file and the profile: ``twenty.fasta.CDPK.fasta`` for the sequences of
        ``twenty.fasta`` that match ``CDPK`` best.

Note that ``--table`` and ``--seqsets`` can be combined.
``-o FILE`` (``--seqout``) also writes every matching sequence to one FASTA file.
``-n`` (``--include``) keeps only the sequences whose best-matching profile has a name containing one of the given strings, separated by commas.
``-x`` (``--exclude``) leaves out the sequences whose best-matching profile has a name containing one of them.
A sequence left out is not reassigned to its next best profile.

``scan`` reads the target as FASTA, whatever its file name, and logs a warning and succeeds when no sequence matches any profile, including when the target is empty or holds only whitespace.


add
```

Scan a target database with the given HMM profile set and add hits that meet a series of acceptance thresholds to the profile FASTA files.
``add`` reads the HMM profile database ``Superfamily_all.hmm`` beside the tree ``Superfamily``, which ``build --hmmer`` writes.

Once you've constructed profiles from a collection of carefully selected sequences representing each subfamily, you can use this command to scan another sequence set and automatically add strong hits to the corresponding profile sequence sets.
The target database could be a proteome, or a larger sequence database like **refseq_proteins**, if you're confident in your coverage of the superfamily and want to improve the sensitivity of your profiles.
A sequence that is already in the tree is never added again, so the tree's own ``*-Unclassified.fasta`` files are of no use as a target; ``--report-promotions``, described below, lists the catch-all members that a subfamily would accept instead.

Each hit is considered for one of the profiles that have a FASTA file of their own: the one whose best-matching domain gives it the smallest E-value in the target database.
The group-level profiles do not take part, so a hit that matches a group profile best is considered for the best-matching profile among those with a FASTA file.
Hits are assigned to profiles by name, so two profiles of the tree cannot share a name.
A catch-all profile (see `Directory tree is the superfamily hierarchy`_) does not compete with the subfamilies of its group: a hit that matches a catch-all best is considered instead for the subfamily profile with the smallest E-value among those in the catch-all's directory or below it whose threshold its bit score reaches, and stays with the catch-all if there is none, or if that subfamily rejects it for its length.
A hit is accepted when it passes the significance cutoff ``-E`` (by default 1e-20) in the target database, when the bit score of its best domain is at least the threshold of the profile, and when the aligned region covers at least the fraction ``-L`` (by default 0.8) of the profile length.

The threshold of a profile is the higher of two scores, both counting only the sequences that its HMM hits at an E-value of 0.001 or better.
The first is the lowest best-domain bit score among the profile's own sequences.
These are the sequences of the alignment from which ``fammer build`` made the HMM, so a sequence that ``add`` appended since the last build does not count until the HMM is rebuilt with it.
A profile of at most ten sequences scores each of them against an HMM built from its alignment without that sequence, since against an HMM built partly from itself a sequence of a small profile scores far higher than an unseen member of the family does; a larger profile scores its sequences against its own HMM.
The second is the highest best-domain bit score of any aligned sequence of a subfamily profile in another directory, scored against the profile's own HMM; consensus records and the sequences of catch-alls do not count.
Without it, a profile whose weakest sequence scores low accepts kinases of families that the tree has no profile for, since such a kinase can outscore that weakest sequence as easily as a member of another group can.
The subfamilies of the profile's own directory are left out, because related subfamilies are not separated by a single score: in the built ``profiles/SomeKinases``, a member of a sister subfamily outscores the weakest sequence of 20 of the 32 subfamilies.
For example, the weakest sequence of ``CAMK1`` scores 101.2 bits against it, a member of a sister subfamily 310.1, and a member of ``Aurora``, in another directory, 242.9, which is ``CAMK1``'s threshold.
The threshold uses bit scores because an E-value grows with the number of sequences searched, while a bit score does not, so the threshold is the same whatever the size of the target.
The second score takes one more search of the tree's subfamily sequences for each profile whose threshold is needed: one that a hit is assigned to, or one that competes for a hit that matches a catch-all best.

A hit is borderline when it clears its profile's threshold by less than 10 bits, or when it does not lead the closest other profile with a FASTA file by 10 bits, not counting a catch-all that passed it on; the other profile may even score it higher, since a hit is assigned by its E-value.
A near tie between a subfamily and the catch-all of its group is therefore borderline when the subfamily matches the hit best, and not when the catch-all matches it best and passes it to the subfamily.
A hit that a catch-all passed on and receives back, because the subfamily rejected it for its length, is tested again as a hit of the catch-all.
The borderline rule follows the manual examination of borderline scores in [Kannan2007]_.
``add`` does not write a borderline hit; it reports it as a warning, which ``--quiet`` does not hide, with its score, the threshold, and the closest other profile and its score, for example (without the log prefix)::

    * borderline TGME49_291050 for CDPK (262.8 vs. 242.8 bits; TSSK 253.2)

A borderline hit is reported again by every run until a rebuild separates it or the curator places it by hand.

Only the aligned region is added, not the full-length sequence.
A sequence with several domains that pass these thresholds contributes its longest domain once, and a sequence whose identifier is already in any profile of the tree is not added again.
A second run with the same target and the same HMMs therefore adds nothing.
After ``fammer build`` has rebuilt the HMMs from the enlarged sequence sets, another run may add a few more sequences.
On a copy of ``profiles/SomeKinases`` with ``test/TgME49.pep`` as the target, the first run adds 21 sequences to nine profiles and reports 7 as borderline, the run after the first rebuild adds 4 more, and the run after the second rebuild adds none and still reports 4; on a copy of ``profiles/ROPK``, the first run adds 15 sequences to ten profiles and reports none.
A family with no profile in the tree and no member in another group can still be taken in, so review what each round adds before rebuilding.
The second score also rejects some true members that score low: in that run, the Aurora kinase TGME49_203010 scores 254.0 bits against ``Aurora``, above its weakest sequence (250.6) but below a member of ``CAMK1`` (259.9).
Holding each sequence out does not help a small profile whose sequences are nearly identical: each still scores far higher against the others than an unseen member of the family does, so the profile may accept nothing at all.
If a profile has hits that it does not yet hold, but its HMM finds none of its own sequences at 0.001 or better, there is no threshold to compare against, and the profile is skipped with a warning.
A target that is empty or holds only whitespace has no hits, which ``add`` reports as an error.
When a profile file does not end with a line ending, as a hand edit can leave it, ``add`` writes one before the first added record.

With ``--report-promotions``, ``add`` also reports the members of each catch-all that a subfamily of its group would accept, and moves nothing; whether such a member belongs in the subfamily is the curator's decision.
The members of every catch-all are scanned with the tree's profiles and assigned as a hit of the target is, by the same thresholds and the same ``-L``; a member assigned to a subfamily in its catch-all's directory or below it is reported.
Consensus records are never reported, nor are the sequences that this run appended or that a subfamily also holds.
Each report gives the member's score against the subfamily, the subfamily's threshold, and the member's score against an HMM of the catch-all built without it, since against the catch-all's own HMM, built partly from it, a member scores higher than an unseen sequence would; a member that is the catch-all's only aligned row, or that belongs to a catch-all without an alignment, is not scored that way.
The report is logged at the same level as the progress messages, so ``--quiet`` hides it::

    * candidate tr|B9N382|B9N382_POPTR in STE-Unclassified for STE7 (256.3 vs. 208.1 bits; STE-Unclassified without it 352.0)

On the same copy of ``profiles/SomeKinases``, 20 members of ``STE-Unclassified`` are reported: 15 plant sequences for ``STE7``, each of which still scores 41 to 101 bits higher against ``STE-Unclassified`` held out, and 5 *Dictyostelium* sequences for ``STE11``, which score 5 to 29 bits higher against ``STE11``.
The curator keeps both in the catch-all, as lineages distinct from those subfamilies.
On ``profiles/ROPK``, 7 members of ``ROPK-Unique`` are reported, 5 of them for ``PKL``, the profile of protein kinases other than rhoptry kinases.
With ``-n`` (``--dry-run``) ``add`` writes nothing and logs what it would add.


refine
``````

Leave-one-out validation of the subfamily profiles in one directory, so that each subfamily holds only sequences that its profile separates from the unclassified ones.
Unlike the other commands, this is non-recursive.

Given a directory and its set of unclassified sequences (``-u``; by default *dirname*/*dirname*-Unique.fasta), ``refine`` tests each subfamily profile in the directory except the unclassified set, any catch-all profile, and the profile named by ``-x``.
Each sequence of a profile is held out in turn: an HMM is built from the profile's alignment without it, and the held-out sequence and the unclassified sequences are scored against that HMM, in bits, counting the hits at an E-value of 0.001 or better.
A sequence fails if it scores no higher than the highest-scoring unclassified sequence, or if the HMM misses it while it hits an unclassified sequence.
The lowest-scoring sequence that fails is dropped and joins the unclassified sequences, and the remaining ones are tested again, until none fails or two remain.
Each dropped sequence is logged with its score and the best unclassified score, and so is the kept sequence of lowest score; a drop whose score is within 10 bits of the best unclassified score is marked ``-- borderline``, but still made.
With ``-n`` (``--dry-run``) nothing is written; otherwise the dropped records are moved, with their text unchanged, from the profile's FASTA file to the unclassified file, and ``fammer build`` then realigns both.
``refine`` reads the alignments of the last build, so a profile whose ``.aln`` file does not hold exactly the sequences of its FASTA file is an error until the tree is rebuilt, as is a sequence that is both in a profile and in the unclassified set; both are checked before any file changes.
Each round of testing builds one HMM per remaining sequence, so a profile that loses many sequences takes correspondingly longer.

Note that if a member of a known subfamily was mistakenly placed in the unclassified set (i.e. was missed by the initial classification), then many of the legitimate members of the subfamily profile could score worse than this high-scoring "unclassified" sequence and be erroneously marked for removal from the profile.
This is easy to spot in the logged output, where nearly the same unclassified score recurs against many dropped sequences.
Run with ``-n`` first, and move such a sequence into its subfamily before refining; ``fammer add --report-promotions`` lists the members of a catch-all that a subfamily of its group would accept.
A consensus record in the unclassified set, one whose identifier ends with ``consensus``, is not an unclassified sequence and is left out of the comparison, since it would outscore many members of the subfamily it summarizes.

cluster
```````

Extract clusters from a sequence set based on phylogenetic relationships.

Uses FastTree to quickly build a tree with branch support values, or reads a given Newick tree, then extracts well-supported clades from the tree and writes the corresponding sequence sets to FASTA files.
Unclustered sequences are written to another "Unique" file.
Also writes a phyloXML tree file (.xml) showing clusters as colorized clades.

The output files are written to the current directory, not beside the input, and are named after the base name of the input file: ``cluster path/to/fam.aln`` writes ``fam.aln.0``, ``fam.aln.1`` and so on, ``fam.aln.Unique`` and ``fam.aln.xml``.
Clusters are numbered by decreasing size, and clusters of equal size by their sorted sequence identifiers.
Sequences read from an alignment (``.aln``) are written without gaps; an unaligned input (``.fasta``) is first aligned with MAFFT.
``cluster`` gives FastTree only the columns in which the sequences that hold a residue carry at least 40% of the total sequence weight, where the weights count similar sequences less, and ``-`` and ``.`` both count as gaps; insert columns are left out (see `Letter case in alignments`_).
FastTree is run as ``fasttree`` or as ``FastTree``, whichever is found on ``PATH``.
Each sequence identifier must occur once in the input.

Three options set the tree and the thresholds that choose the clusters:

``--size N``
    The smallest number of sequences in a cluster; the sequences of a smaller clade go to the Unique file.
    By default, the larger of 4 and the number of sequences to the power 0.4, rounded to an integer.

``--confidence X``
    Collapse the splits whose branch support is below X before clustering, on the scale of the tree's own support values: FastTree writes values between 0 and 1, while bootstrap percentages run from 0 to 100.
    The splits are the internal clades other than the root, which is never collapsed.
    A split without a support value counts as weakly supported.
    By default, X is the mean of the splits' support values, and a tree whose splits have no support values is not collapsed.
    An explicit X is reported as an error if no split has a support value, if X is above every split's support, or if X is above 0 and up to 1 while the tree's supports run above 1, as bootstrap percentages do, since X would then collapse almost no split.
    X of 0, which collapses only the splits without a support value, is allowed on either scale.

``--tree FILE``
    Cluster this Newick tree instead of building one, so that neither MAFFT nor FastTree is run.
    The tips must be named by exactly the identifiers of the sequence file, each once.
    An identifier containing whitespace or any of the characters ``:,;()[]'`` must be quoted, as in ``'TAIR|locus:2011201'``; unquoted, Newick reads it as a shorter name, which is then reported as naming no input sequence.
    The sequence file is still required, for the sequences to write.

After writing the files, ``cluster`` logs the silhouette of the clusters, in lines such as::

    Mean silhouette 0.630 over the 53 sequences of 6 clusters, leaving out 11 Unique sequences
    Silhouette of DYRK.fasta.0: 0.592

The silhouette of a sequence compares its mean distance to the other members of its cluster with its mean distance to the members of the nearest other cluster.
It ranges from -1 to 1, and higher values mean tighter and better separated clusters.
The first line gives the mean over all clustered sequences, as ``sklearn.metrics.silhouette_score`` computes it, and each following line the mean over the members of one cluster.
The distances are patristic distances, the summed branch lengths between two tips, in the tree as FastTree built it or as given with ``--tree``.
They are measured before weakly supported splits are collapsed, because collapsing a split lengthens the path between the tips below it, which would make the distances depend on ``--confidence``.
The Unique sequences are not a cluster, so they are left out of the score.
With fewer than two clusters the silhouette is undefined, and a line saying so is logged instead.
A larger ``--size`` or ``--confidence`` moves sequences to the Unique file, and leaving out the sequences that fit no cluster tends to raise the score, so compare settings by the score together with the number of sequences it covers.


fammer-tmalign
``````````````

Align protein structures, given as PDB files, into a multiple sequence alignment in FASTA format::

    fammer-tmalign [-s SEED.fasta]... [-o OUTPUT.fasta] STRUCTURE.pdb...

``fammer build`` aligns the structures of each directory of a tree in the same way, and also with MUSTANG where it is installed (see `build`_).
Every pair of structures is aligned with TM-align_.
The pairwise alignments that form the spanning tree of maximum total TM-score are passed to MAFFT as seed alignments, which MAFFT combines into a multiple sequence alignment.
A spanning tree over N structures has N-1 edges, so this passes N-1 seeds rather than one for every pair of structures.
``-s`` (``--seed``) passes an existing FASTA alignment to MAFFT as one more seed and includes its rows in the output; it can be given more than once.
With exactly two structures, the output is the TM-align alignment of the pair: MAFFT is not run, and the rows of a seed alignment are not included.
A structure that occurs in two of the seed alignments can occur in the output more than once, as differently gapped rows, because MAFFT does not keep a sequence that recurs across seed alignments fixed.
The alignment is written to standard output, or with ``-o`` (``--output``) to a file.

Each structure's file name becomes the identifier of its sequence in the alignment, and a FASTA identifier ends at the first whitespace, so structure file names must not contain whitespace.
Seed file names must not contain whitespace either, because MAFFT splits seed file paths at whitespace.
Each structure must be listed only once, because a structure listed twice would be aligned to itself.
TM-align reports no error for a missing, empty or malformed structure file, so ``fammer-tmalign`` and ``fammer build`` check that the structure and seed files exist before aligning any pair, and report TM-align output without a TM-score or alignment as an input error.


Bundled modules
---------------

tasks
`````

Serves the basic purpose of a build tool like Make or Rake: compare the time
stamps of input and output files at each step of the `build`_ process, and only
rebuild the outputs that are out of date. Also track intermediate files to clean
up after the process successfully completes. See `this blog post about it
<http://etalog.blogspot.com/2012/01/building-analysis-how-to-avoid.html>`_.


How to curate profiles
----------------------

Directory tree is the superfamily hierarchy
```````````````````````````````````````````

Begin by creating a directory tree with each subfamily's representative
sequences in an unaligned FASTA file.  The FASTA file names must end in
``.fasta``.

A simple un-nested layout looks like::

    Superfamily/
        subfam1.fasta
        subfam2.fasta
        subfam3.fasta
        ...
        Superfamily-Unclassified.fasta

Typically, a protein superfamily will have some members that cannot be cleanly classified into subfamilies.
Put them in a catch-all file whose name ends with ``-Unclassified``, ``-Unclass`` or ``-Unique``; ``add`` and ``refine`` treat such a file as the members of a group that fit none of its subfamilies, rather than as a subfamily; it may also hold consensus records (below).

A catch-all's profile matches its group as a whole, without matching any one subfamily specifically.
Unclassified sequences from annotated genomes are scarce, sometimes unreliable, and sometimes incomplete, so the catch-all sets of the ``CAMK``, ``CMGC`` and ``STE`` groups in ``profiles/SomeKinases`` are supplemented with a consensus sequence of each known subfamily of the group, which together make a generic profile of the group: the consensus sequences of different subfamilies are intended to dilute each other, so that a sequence of a real subfamily still matches that subfamily more strongly.
A record whose identifier ends with ``consensus`` is such a consensus sequence.
``add --report-promotions`` never reports one, ``refine`` leaves it out of the unclassified set, and the threshold of ``add`` does not count it as a sequence of another group.
It does count toward the weakest-sequence score of the catch-all that holds it: ``CAMK-Unclass`` and ``CMGC-Unclass`` hold nothing else, so that score comes from their consensus records alone.

A catch-all can also hold the members of lineage-specific subfamilies that are too sparsely populated to support a profile of their own.
In ``profiles/SomeKinases``, ``STE-Unclassified.fasta`` holds plant and *Dictyostelium* STE kinases, represented by the records ``STE-plant1-consensus``, ``Dicty2-consensus`` and ``Dicty3-consensus``; these lineages are distinct from the ``STE7`` and ``STE11`` subfamilies.
``CAMK-Unclass.fasta`` also holds consensus records of CAMK families that have no profile in that tree, such as CAMKL, the family of the SNF1/AMPK kinases.
Whether a catch-all member belongs to one of its group's subfamilies is a question for phylogenetic analysis, not for a score threshold.

Recursive nesting is also allowed (in fact, that's the point of this project).
It looks like this::

    Superfamily/
        Group1/
            subfam1_Group1.fasta
            subfam2_Group1.fasta
            subfam3_Group1.fasta
            ...
            Group1-Unclassified.fasta
        Group2/
            subfam1_Group2.fasta
            subfam2_Group2.fasta
            subfam3_Group2.fasta
            ...
            Group2-Unclassified.fasta
        Group3/
            subfam1_Group3.fasta
            subfam2_Group3.fasta
            subfam3_Group3.fasta
            ...
            Group3-Unclassified.fasta
        ...
        Superfamily-Unclassified.fasta


At least 2 sequences are needed to define each subfamily.  The sequences are
initially retrieved and organized manually. External databases can help; for
example, KinBase_ provides a classification scheme and representative sequences
for the protein kinase superfamily.  The directory hierarchy defines the higher
levels of organization of the superfamily.

Now build the initial alignments::

    fammer build --clean Superfamily/


Trim the tails
``````````````

For best results, sequences should all cover the same conserved domain region
and align exactly at N- and C-termini. For example, the eukaryotic protein
kinase domain begins 7 residues before the glycine-rich loop (GxGxxG motif) and
ends 12 residues after a conserved arginine in subdomain XI (RP[TS] motif).
The consensus function used in Fammer includes the flanking sequences before
the first conserved block or after the last conserved block, regardless of
"gappiness", so even one sequence in a set can be used to define domain
boundaries for the whole set.

Edit the subfamily-level alignments (.aln) to trim the sequences to the domain
boundaries. I use Vim's block-selection mode (Ctrl-v); you might prefer
Jalview_. Do not edit higher-level .aln files; just look at them to see which
subfamily-level alignments have extended tails, then edit those .aln files.

When you think you've trimmed all the .aln files to the conserved domain
region, update the corresponding unaligned sequence files (.fasta) to match::

    fammer update-fasta Superfamily/

Rebuild the alignments whose source sequence sets have changed::

    fammer build --clean Superfamily/

MAFFT may create better alignments for sets where the unalignable regions have
now been removed. When this step completes, examine the top-level alignment
(e.g. Superfamily.aln) -- are the domain boundaries aligned cleanly across all
sequence sets? Repeat the process if any subfamilies need to be trimmed further.

If you realize you've trimmed a profile too far, use your version control
system (you are using Git or Mercurial to take snapshots of the tree, right?)
to revert the .fasta file to an earlier version, then rebuild and try trimming
again.

If you never had the full-length sequences for a subfamily to begin with, try
to find one representative full-length sequence from a database like UniProt,
add it to the .fasta file, realign, and trim that sequence to the region it
should cover. This won't improve the HMM or MAPGAPS profile for that subfamily
much, but it will help the higher-level alignments that include the consensus
sequence of that subfamily.


Letter case in alignments
`````````````````````````

HMMER and CHAIN (Neuwald et al. 2003, Genome Research 13:673) mark the state of each residue by its case: an uppercase residue is aligned to a match column of a profile and a lowercase residue is an insertion relative to it.
In HMMER alignments, ``-`` is a deletion in a match column and ``.`` a gap in an insert column.
Alignments in this convention come from ``hmmalign`` and ``hmmsearch -A``, and from the CMA files that CHAIN and MAPGAPS read and write.
BLAST and PSI-BLAST use lowercase for masking instead: with ``-lcase_masking``, lowercase residues in a query or subject sequence are masked, so convert such residues to uppercase before placing them in an alignment.

MAFFT writes protein residues in uppercase, so every ``.aln`` file that ``fammer build`` creates is in uppercase, whatever the case of the input sequences.
A family ``.aln`` file that you edit by hand is used unchanged by ``fammer build`` as long as it is newer than its ``.fasta`` file.
In such a file, Fammer reads letter case as HMMER does: a column whose residues are all lowercase is an insert column, and the family consensus has a gap there; a column that mixes lowercase and uppercase residues is read as uppercase, with a warning.
An alignment that holds no uppercase residue at all, such as the nucleotide alignments that MAFFT writes, marks nothing by case and is read without regard to it.
``fammer cluster`` builds its tree only from columns that hold an aligned residue, so it leaves out the insert columns of an alignment given to it, whatever their proportion of residues.
``hmmbuild`` does not read letter case in the alignments that Fammer gives it, so the match states of an HMM profile are chosen by the proportion of residues in each column, whatever their case.
When ``fammer update-fasta`` rewrites a ``.fasta`` file from its edited ``.aln`` file, it writes the sequences in uppercase, and ``fammer build`` then realigns the family with MAFFT, so lowercase in that ``.aln`` file does not survive the rebuild.

``fammer build --mapgaps`` writes a ``.cma`` profile for each family and a ``.tpl`` template for each group, in the CMA format of CHAIN and MAPGAPS, and these files use the convention themselves.
A column of an alignment becomes a match column of its profile where the consensus of the alignment has a residue, and the residues of the other columns are written in lowercase, except that an insertion next to a deletion is written in uppercase as a match.
A profile can therefore have match columns that its parent profile lacks, and in the parent's template the row of its consensus sequence holds lowercase residues at most of those columns.
For example, built from ``profiles/SomeKinases``, the CAMK row of ``SomeKinases.tpl`` reads ``KDKlDREQ`` and ``LLDskDPD`` where ``CAMK.tpl`` reads ``KDKLDREQ`` and ``LLDSKDPD``.


PDB files for structural alignment
``````````````````````````````````

While MAFFT typically creates good alignments within a subfamily, for
high-level profiles it may struggle to align the consensus sequences of very
divergent families or groups to each other.  In these cases, seeding with a
structural alignment can help line up homologous regions.

Manually identify the high-quality solved crystal structures that correspond to
families in your tree, and place those PDB files in the directory tree at the
same level as the subfamily they represent.

Open the PDB file (.pdb) in a text editor and:

- Remove the ATOM records corresponding to residues outside the conserved
  domain.
- If multiple chains are present, choose the best, most complete chain and
  delete the others. (Otherwise, Fammer will take the first chain by default.)

To determine the domain boundaries, load a "reference" PDB (e.g. 1ATP for
kinases) and the other PDB together in PyMOL and align using the command
"cealign" or "fit". Visually find which residues correspond to the start and
end of the conserved domain.  If your PDB structure diverges from the reference
structure drastically before or after a certain point (i.e. N- or C-terminal
region of the domain is non-homologous -- excluding inserts), it may be best to
truncate the PDB to remove the non-homologous portion as it cannot be aligned
accurately.

Save the edited .pdb file, and rebuild the higher-level profiles::

    fammer build --clean Superfamily/

Now that PDB files are present in the directory tree, Fammer writes
structural alignments in FASTA format to the .pdb.seq files.

Once complete, examine the top-level alignment (e.g. Superfamily.aln) for
misalignments. To fix these, first edit the corresponding .pdb.seq file.
Usually there is just one structure that was misaligned.  Rebuild and re-edit
as necessary. If a particular PDB is very poorly aligned to the others, it may
be best to just remove it altogether -- it may have a different conformation
from the others.

Development
-----------

Two development workflows are supported from the same ``pyproject.toml``, and both lock files are committed.
Both install Fammer in editable mode with its ``mapgaps`` extra and the tools in the ``dev`` dependency group: pytest, ruff, and scikit-learn, which the tests use as an independent implementation of the silhouette.
Both also need a ``git`` executable, because biofrills and biocma are installed from their repositories until they are released on PyPI.

uv_ takes the other Python packages from PyPI, at the versions recorded in ``uv.lock``::

    uv sync --all-extras
    uv run pytest
    uv run ruff check
    uv run ruff format

pixi_ takes Python and the other Python packages from conda-forge, and also installs the programs that Fammer runs: MAFFT from conda-forge, and HMMER, TM-align, MUSTANG and FastTree from bioconda, at the versions recorded in ``pixi.lock``.
The ``default`` environment uses the latest Python that conda-forge provides, and ``py312`` the oldest supported one::

    pixi run test
    pixi run -e py312 test
    pixi run test-integration
    pixi run lint            # ruff check
    pixi run fmt             # ruff format
    pixi run fmt-check       # ruff format --check

The unit tests compare the library with reference outputs recorded under Python 2.7 in ``tests/data``, replay the recorded ``cluster`` runs with the recorded output of MAFFT and FastTree, and check properties that need no reference.
The tests marked ``integration`` run ``fammer build`` on a copy of ``profiles/SomeKinases`` and compare the result with the reference outputs, then run ``scan``, ``update-fasta`` and ``cluster`` with the external programs.
``pixi run test`` leaves them out, and elsewhere they are skipped unless the programs are on ``PATH`` and the tests run from a Git checkout, so ``uv run pytest`` runs the unit tests alone unless those programs are already installed; ``pixi run test-integration`` fails instead of skipping.
MAPGAPS is not packaged for either workflow, so ``--mapgaps`` cannot be run end to end.

The GitHub Actions workflow in ``.github/workflows/ci.yml`` runs on every pull request and every push to ``master``, and can also be started by hand.
It checks the code with ``ruff check`` and ``ruff format --check``, runs the unit tests with uv on every supported Python version on Linux and macOS, and runs ``pixi run test-integration`` on both systems.
Dependabot, configured in ``.github/dependabot.yml``, proposes new releases of the actions that the workflow uses in one pull request a month.
It does not update the pixi release that the integration job installs, and it is not configured for the Python dependencies, because ``uv.lock`` and ``pixi.lock`` must change together.

Code is formatted with ``ruff format``, and the formatting-only commits are listed in ``.git-blame-ignore-revs``.
To skip them in ``git blame``::

    git config blame.ignoreRevsFile .git-blame-ignore-revs

.. _uv: https://docs.astral.sh/uv/
.. _pixi: https://pixi.sh/
.. _KinBase: http://kinase.com/web/current/kinbase/
.. _Jalview: https://www.jalview.org/

