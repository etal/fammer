======
Fammer
======

Utilities for curating a hierarchical set of sequence profiles representing a
protein superfamily.

If you use this software in a publication, please cite our paper that describes
it:

    Talevich, E. & Kannan, N. (2013) Structural and evolutionary adaptation of
    rhoptry kinases and pseudokinases, a family of coccidian virulence factors.
    *BMC Evolutionary Biology* 13:117
    doi:10.1186/1471-2148-13-117

    Available at: http://www.biomedcentral.com/1471-2148/13/117

Freely distributed under the permissive BSD 2-clause license (see LICENSE).

Installation
------------

The installation consists of a Python library, ``fammerlib``, and two scripts,
``fammer.py`` and ``tmalign.py``.

Download the .zip file and unpack it, or clone this Git repository, to get the
source code.

To use all the features of Fammer, you'll need the following third-party
programs installed:

- Python_ 2.7
- MAFFT_
- HMMer_ 3.0
- MAPGAPS_ (optional)
- TMalign_ (optional) -- for structural alignments
- FastTree_ (optional) -- for clustering

.. _Python: http://www.python.org/download/
.. _MAFFT: http://mafft.cbrc.jp/alignment/software/
.. _HMMer: http://hmmer.janelia.org/
.. _MAPGAPS: http://mapgaps.igs.umaryland.edu/
.. _TMalign: http://cssb.biology.gatech.edu/skolnick/webservice/TM-align/index.shtml
.. _FastTree: http://www.microbesonline.org/fasttree/

.. For hackers, also PRANK: http://code.google.com/p/prank-msa/

If you're on a Debian-based Linux system, check your package manager for these
first to save yourself some time::

    sudo apt-get install mafft hmmer tm-align fasttree python-pip

Then, install the Python library dependencies and Fammer itself as follows.

Recommended:
````````````

Install the Python packaging system pip or setuptools. Then run the setup
script, and all Python dependencies will be pulled in::

    python setup.py build
    python setup.py install

(You might need root privileges for the last step.)

Manual:
```````

Install the Python libraries Biopython_, biofrills_, biocma_ and networkx_.
Then run the setup script as above.

.. _Biopython: http://biopython.org/wiki/Download
.. _biofrills: https://github.com/etal/biofrills
.. _biocma: https://github.com/etal/biocma
.. _networkx: http://networkx.lanl.gov/



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


Commands
--------

build
`````

Construct a profile database from a directory tree of family profile alignments.

Assume we have a directory tree set up under ``Superfamily/`` as above.
Next, run ``fammer.py build Superfamily`` to align all sequence files with
MAFFT, and (recursively up) align the consensus sequences of each subfamily
together::

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
``Superfamily.aln`` will be rebuilt the next time ``fammer.py build
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


At least 2 sequences are needed to define each subfamily.  The sequences are
initially retrieved and organized manually. External databases can help; for
example, KinBase provides a classification scheme and representative sequences
for the protein kinase superfamily.  The directory hierarchy defines the higher
levels of organization of the superfamily.

Now build the initial alignments::

    fammer.py build --clean Superfamily/


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
JalView. Do not edit higher-level .aln files; just look at them to see which
subfamily-level alignments have extended tails, then edit those .aln files.

When you think you've trimmed all the .aln files to the conserved domain
region, update the corresponding unaligned sequence files (.fasta) to match::

    fammer.py update-fasta Superfamily/

Rebuild the alignments whose source sequence sets have changed::

    fammer.py build --clean Superfamily/

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

    fammer.py build --clean Superfamily/

Now that PDB files are present in the directory tree, Fammer will writes
structural alignments in FASTA format to the .pdb.seq files.

Once complete, examine the top-level alignment (e.g. Superfamily.aln) for
misalignments. To fix these, first edit the corresponding .pdb.seq file.
Usually there is just one structure that was misaligned.  Rebuild and re-edit
as necessary. If a particular PDB is very poorly aligned to the others, it may
be best to just remove it altogether -- it may have a different conformation
from the others.

