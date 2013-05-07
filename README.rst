======
Fammer
======

Utilities for curating a hierarchical set of sequence profiles representing a
protein superfamily.


Installation
------------

Dependencies:

- Python_ 2.7
- Python libraries Biopython_, biofrills_, biocma_, networkx_
- MAFFT_
- PRANK_
- Exonerate_ (optional, used by PRANK for speed)
- HMMer_ 3.0
- MAPGAPS_ (optional)
- TMalign_ (optional) -- for structural alignments
- FastTree_ (optional) -- for clustering

.. _Python: http://www.python.org/download/
.. _Biopython: http://biopython.org/wiki/Download
.. _biofrills: https://github.com/etal/biofrills
.. _biocma: https://github.com/etal/biocma
.. _networkx: http://networkx.lanl.gov/
.. _MAFFT: http://mafft.cbrc.jp/alignment/software/
.. _PRANK: http://code.google.com/p/prank-msa/
.. _Exonerate: http://www.ebi.ac.uk/~guy/exonerate/
.. _HMMer: http://hmmer.janelia.org/
.. _MAPGAPS: http://mapgaps.igs.umaryland.edu/
.. _TMalign: http://cssb.biology.gatech.edu/skolnick/webservice/TM-align/index.shtml
.. _FastTree: http://www.microbesonline.org/fasttree/


I recommend creating a symbolic link to ``fammer.py`` in your ``$PATH``, e.g.::

    ln -s `pwd`/fammer.py ~/bin/fammer


Basic usage
-----------

Global options:

  ``-h``, ``--help``
      Show a help message and basic usage.
  ``--quiet``
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
        meet acceptance thresholds to the profile FASTA files.
    `refine`_
        Leave-one-out validation of HMM profiles.
    `cluster`_
        Split a sequence set into clusters (based on phylogeny).


Directory tree is the superfamily hierarchy
-------------------------------------------

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

Typically, a protein superfamily will have some members that cannot be cleanly
classified into subfamilies.

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


Commands
--------

build
`````

Construct a profile database from a directory tree of family profile alignments.

Assume we have a directory tree set up under ``Superfamily/`` as above.
Next, run ``fammer build Superfamily`` to align all sequence files with MAFFT,
and (recursively up) align the consensus sequences of each subfamily together::

    Superfamily/
        Group1/
            subfam1_Group1.fasta
            subfam1_Group1.aln
            subfam2_Group1.fasta
            subfam2_Group1.aln
            subfam3_Group1.fasta
            subfam2_Group1.aln
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
``Superfamily.aln`` will be rebuilt the next time ``fammer build Superfamily``
is called because the consensus sequences that constitute those alignments may
have changed. (It's like Make.)

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

The ``--mapgaps`` option works similarly, if you have the necessary programs
installed.

The ``--clean`` option can be included with any of the above commands to remove
intermediate files.

If you have included PDB structures in your directory tree and have a structure
alignment program installed, the ``--pdb`` option will first create a structural
alignment of the PDBs in the directory, then use that alignment as the seed for
higher-up alignments::

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
        filenames indicate the name of the source sequence file name and the
        matching HMM profile names.

Note that ``--table`` and ``--seqsets`` can be combined.


add
```

Scan a target database with the given HMM profile set and add hits that meet
a series of acceptance thresholds to the profile FASTA files.

Once you've constructed profiles from a collection of carefully selected
sequences representing each subfamily, you can use this command to scan another
sequence set and automatically add strong hits to the corresponding profile
sequence sets. The target database could be the ``*-Unclassified.fasta``
sequence sets, to catch any classifiable members that were not noticed
initially, or a larger sequence database like **refseq_proteins**, if you're
confident in your coverage of the superfamily and want to improve the
sensitivity of your profiles.


refine
``````

Leave-one-out validation of sequence profiles.
Unlike the other commands, this is non-recursive.

Given a target subdirectory and the name of the subdirectory's
``*-Unclassified.fasta`` file (if not specified, it looks for
*dirname*-Unclassified.fasta), scan each subfamily's sequence set (``.fasta``)
with the corresponding HMM profile (``.hmm``), and also scan the
``-Unclassified.fasta`` file with all the HMMs to obtain scores for each
sequence and each profile. Then, compare the scores of sequences in a subfamily,
starting with the worst-scoring sequence, to the highest-scoring "unclassified"
sequence by the same profile. If, for a given profile, a classified sequence
scores worst than an unclassified one, mark the classified one for removal from
the sequence profile.

Note that if a member of a known subfamily was mistakenly placed in
``-Unclassified.fasta`` (i.e. was missed by the initial classification), then
many of the legitimate members of the subfamily profile could score worse than
this high-scoring "unclassified" sequence and be erroneously marked for removal
from the profile. This is easy enough to spot in the logged output. One way to
avoid it is to first use the ``add`` command with  ``-Unclassified.fasta`` as
the target, to catch and classify such sequences beforehand.

cluster
```````

Extract clusters from a sequence set based on phylogenetic relationships.

Uses FastTree to quickly build a tree with branch support values, then extracts
well-supported clades from the tree and writes the corresponding sequence sets
to FASTA files. Unclustered sequences are written to another "Unique" file.
Also writes a phyloXML tree file (.xml) showing clusters as colorized clades.


Bundled modules
---------------

tasks
`````

Serves the basic purpose of a build tool like Make or Rake: compare the time
stamps of input and output files at each step of the `build`_ process, and only
rebuild the outputs that are out of date. Also track intermediate files to clean
up after the process successfully completes. See `this blog post about it
<http://etalog.blogspot.com/2012/01/building-analysis-how-to-avoid.html>`_.

tmalign
```````

Align multiple structures using TMalign_ for pairwise alignments and a minimum
spanning tree constructed from the pairwise TM-scores to assemble the pairwise
alignments into a multiple sequence alignment. This module can also be used as a
command-line script.

