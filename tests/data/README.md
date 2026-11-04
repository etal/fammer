# Reference outputs from the Python 2 code base

The files under `py27/` record what fammer computed on Python 2.7, before the port to Python 3.
They are the reference against which the port is checked: a ported function or command should reproduce these results, or the difference should be explained by a documented, deliberate change.
The files under `tool-outputs/` are raw output of the external programs fammer runs, kept so that fammer's parsers can be tested without those programs installed.
`differences.json` lists the recorded cases that the current code deliberately computes differently, as described under "Deliberate differences".

## Producing versions

`py27/versions.json` is written by the capture script and is the authoritative record.
At capture time it read:

| Component | Version |
|---|---|
| Python | 2.7.15 (CPython, conda-forge build) |
| Platform | Darwin 25.6.0, x86_64 under Rosetta 2 on Apple silicon |
| Biopython | 1.76 |
| NetworkX | 2.2 |
| ReportLab | 3.5.42 (`fammerlib.cluster` imports `Bio.Graphics`, which requires it) |
| biofrills | commit `36684bb`, its last Python 2 commit |
| biocma | commit `eac0c57` |
| fammer | commit `5bcf8eb`, sources unmodified |
| MAFFT | 7.526 (conda-forge, osx-arm64) |
| HMMER | 3.4 (bioconda, osx-arm64) |
| FastTree | 2.2.0, double precision (bioconda, osx-arm64) |
| TM-align | 20240303 (bioconda, osx-arm64) |

biofrills and biocma are extracted with `git archive` from sibling checkouts at the commits above, so later work in those repositories does not change the capture.

## Regenerating

conda-forge publishes no Python 2.7 build for osx-arm64, but its osx-64 build runs under Rosetta.
The external programs run natively:

```sh
CONDA_SUBDIR=osx-64 conda create -p /tmp/fammer-py27 --override-channels -c conda-forge \
    python=2.7 biopython=1.76 numpy networkx=2.2 reportlab
conda create -p /tmp/fammer-tools --override-channels -c conda-forge -c bioconda \
    mafft hmmer tmalign fasttree
```

The capture must run the Python 2 sources of commit `5bcf8eb`, and the script refuses any other version of them.
Later commits change those sources, and even a change that preserves behavior moves the source line numbers recorded in the standard error of the end-to-end runs.
So run it from a Git worktree of `86ac15b`, the first commit that holds the capture script, with the current `tests/data/` copied in.
From the repository root, with `../biofrills` and `../biocma` checked out (or `BIOFRILLS_REPO` and `BIOCMA_REPO` naming them):

```sh
git worktree add --detach ../fammer-capture 86ac15b
cp -R tests/data/. ../fammer-capture/tests/data/
PATH=/tmp/fammer-tools/bin:$PATH /tmp/fammer-py27/bin/python -B ../fammer-capture/tests/data/capture_py27.py &&
  rm -r tests/data/py27 tests/data/tool-outputs &&
  cp -R ../fammer-capture/tests/data/py27 ../fammer-capture/tests/data/tool-outputs tests/data/
git worktree remove --force ../fammer-capture
```

The worktree is a sibling of the repository so that the script finds `../biofrills` and `../biocma` from it too.
The script deletes and rewrites `py27/` and `tool-outputs/` inside the worktree, and takes about two and a half minutes.
If it stops early, nothing in the repository itself has changed; remove the worktree with the last command before trying again.
It refuses to run if the fammer sources or the input directories `profiles/` and `test/` differ from the committed state of the worktree, and `versions.json` records the sources commit and the last commit that changed the inputs.
Two consecutive captures on the same machine produced identical files.
A capture on a different kind of file system produces the order-only differences described under "Interpreting a mismatch".
On Linux, create the first environment without `CONDA_SUBDIR`.
The recorded sources call the FastTree program `fasttree`, whereas bioconda installs it as `FastTree`; the case-insensitive default macOS file system resolves one name to the other, but on Linux the capture needs a symbolic link named `fasttree`.
The current code looks for either name.

### The MAFFT wrapper

The captured sources of commit `5bcf8eb` align the family consensus sequences of every directory in `build.align_profiles` with `mafft --globalgenafpair`.
MAFFT stopped implementing that option in 2012 and stopped accepting it at the turn of 2014, as its source releases show:

| MAFFT | `--globalgenafpair` |
|---|---|
| 6.864 (November 2011) and earlier | global pairwise alignment with generalized affine gap costs |
| 6.903 (May 2012) to 7.123 (October 2013) | accepted, but the pairwise step ignores it and runs the same global algorithm as `--globalpair` |
| 7.130 (December 2013) | prints `No longer supported!` and exits without an alignment |
| 7.141 (February 2014) and later, including every packaged build | exits with status 1 and `--globalgenaf is no longer supported!` |

Those sources therefore cannot build any profile tree with a MAFFT that can be installed today.
The capture puts `tool-shims/` first on `PATH`; its `mafft` script replaces `--globalgenafpair` with `--globalpair` and runs the real MAFFT, and every other MAFFT call passes through unchanged.
For the MAFFT releases contemporary with fammer's development in 2012 and 2013, this replacement computes what the original option computed.
The current code passes `--globalpair` itself, so it needs no wrapper, and with MAFFT 7.526, the version recorded in `py27/versions.json`, it writes files identical to those recorded in `e2e/build.json`.
The wrapper remains for the capture, which runs the unmodified sources.

## Conventions

Every function call is recorded as `{"value": ...}` or, where the original code raised an exception, `{"error": {"type": ..., "message": ...}}`.
Exceptions are part of the reference.
The type is the reference and the message is the interpreter's: Python 3 raises `FileNotFoundError`, a subclass of `OSError`, where these files record `OSError`, and words several messages differently (for example `could not convert string to float: 'abc'`), so compare types with `isinstance` and treat messages as informational.
Calls that log record a `log` list of `[level, module, message]` rows at level INFO and above.
Floats are written with Python 2.7's `repr`, the shortest decimal string that reads back as the identical double.
Tuples are written as JSON lists.

Machine-specific paths are replaced by labels; the current files contain `<repo>` (the repository), `<taskdir>` (the scratch directory of the task scenarios) and `<tmpfile>` (a randomly named temporary file).
The script defines further labels for its other scratch directories and the Python installation, which appear only if such a path reaches a recorded string.
The standard error of a command-line run uses fammer's log format, `module [@line]: message`, and may hold a traceback; both carry line numbers of the Python 2 sources, which differ from those of the current sources.

## Layout

### `py27/library/`

Pure-Python functions, called in-process.

- `tasks.json`: `ext` and `noext` on several path shapes (note `noext('dir.v2/fam') == 'dir'`), `is_empty` by file size, `Task.build`, `is_up_to_date` and `clean` over files with controlled modification times (including a missing dependency and a failing action), `sh` on succeeding and failing commands, and `which` on a missing program.
- `build-SomeKinases.json`, `build-ROPK.json`: `treeify_subdirs` and the task graph that `taskify_subdirs` defines for every combination of `hmmer`, `mapgaps` and `use_pdb`, run on a copy of the tracked files of `profiles/<tree>`.
  Each task is recorded by target, with its action's name, keyword arguments, dependencies (`{"task": target}` or `{"file": path}`) and files to clean.
- `build-SomeKinases-group-fasta.json`: the same for `profiles/SomeKinases` with a FASTA file added beside each group directory (`CAMK.fasta`, `CMGC.fasta`, `STE.fasta`).
  Both functions skip a FASTA file named after a sibling directory, and the bundled trees contain no such file, so only this case detects a port that turns their list of subdirectories into a one-shot iterator.
- `mapgaps.json` and `mapgaps/`: `build.mg_aln2cma`, which writes MAPGAPS profiles through biocma, on three alignments from `e2e/build/`, one per level of the profile tree: `SomeKinases/CAMK/CDPK.cma` (level 2), `SomeKinases/CAMK.tpl` (level 1) and `SomeKinases.tpl` (level 0).
  The written files are kept under `mapgaps/` at their paths in the build tree.
  Each case also records the consensus sequence that `mg_aln2cma` prepends to the alignment, so that the biocma writer can be compared on its own when the biofrills consensus differs.
- `share-scan-add.json`: `_share.parse_scanned_hits` (since replaced by `parse_scanned_rows` and `best_hits`; see "Assignment of `add` hits" under "Deliberate differences") on the example in its docstring, on a synthetic table with ties, and on `tool-outputs/hmmsearch/twenty.fasta.SomeKinases_all.tbl`; `scan.write_summary` and `scan.write_table` output; `add.filter_profile` and `add.filter_evalue` (since replaced; see "Acceptance threshold of `add`" under "Deliberate differences"); `_share.write_fasta`.
- `update_fasta.json`: `ungap_and_unique` on synthetic record sets and on `profiles/ROPK/PKL.aln`, and `is_same_aln_fasta` against several FASTA variants of one alignment.
- `tmalign.json`: `read_tmalign_as_seqrec_pair` on each file in `tool-outputs/tmalign/`, as written and with its final CPU-time line removed; `tmscore_from_description`, `clean_and_dedupe_seqs`, `purge_seqs` and `mst_pairs`.
- `cluster.json`: `geom_mean`, `cluster.write_fasta` (since removed; see "FASTA headers" under "Deliberate differences"), and `load_tree` on input that is neither `.aln` nor `.fasta`.
- `hmmer.json`: `_share.get_hmm_evalues` (since removed; see "`refine`" under "Deliberate differences"), `add.get_hmm_len` and `add.walk_profiles`, which need HMMER or built profiles and so run in the end-to-end workspace after `fammer build`.

### `py27/e2e/`

Command-line runs of `fammer.py` and `tmalign.py` in a scratch directory holding a copy of the tracked files of `profiles/SomeKinases` and the sequence files of `test/`.
Each run records its arguments, exit status, standard output, standard error, and the files it created, modified or deleted.

- `build.json`: `fammer build SomeKinases --hmmer --no-pdb --tree`, with a manifest of the 213 files written (size and SHA-256; `.hmm` files hashed without their `DATE` line; HMMER's binary `.h3?` files by size only, because they embed the date).
  The `.aln`, `.cons.seq`, `.families.fa` and `.nwk` files are kept in full under `e2e/build/`.
- `build-clean.json`: `fammer build ... --clean` on the built tree, which rebuilds nothing and deletes the intermediate files listed in the task graph.
- `scan.json`: `fammer scan` on `test/twenty.fasta` and `test/TgME49.pep`, with the summary and table outputs, the `-n`, `-x`, `-o` and `-O` options, and a stricter E-value.
  Each run also left the hmmsearch table `<target>.tbl` in the working directory; the current code writes the table to a temporary file and removes it, so `test_scan_reproduces_recorded_run` expects every recorded file except that one.
- `update-fasta.json`: `fammer update-fasta` on the built tree and on one alignment, each before and after the last record of an alignment (`CDK.aln`, then `DYRK.aln`) was deleted.
- `add.json`: `fammer add SomeKinases TgME49.pep` on the built tree, with the text appended to each FASTA file.
  `test_add_appends_the_recorded_hits_once` compares the appended records, header and sequence, as a multiset per file; see "Repeated runs of `add`", "Acceptance threshold of `add`" and "Assignment of `add` hits" under "Deliberate differences".
- `refine.json`: `fammer refine SomeKinases -n -u SomeKinases/Testfam-Unclassified.fasta`.
  `test_refine_dry_run_keeps_every_sequence_and_writes_nothing` compares the profiles it tests and the files it writes; see "`refine`" under "Deliberate differences".
- `pdb.json`: `fammer build` with structures enabled, and `tmalign.py` on two and on three structures.
  Both `tmalign.py` runs fail; the first lists `tm.out.fasta` as created because the recorded sources opened the `-o` file before aligning, and the second lists it as neither created nor modified because the first had left it empty.
  The current `fammer-tmalign` writes the file only after a successful alignment, so neither run would leave it.
- `cluster.json`: `cluster.cmd_cluster` on the built `CDK.aln`, on the unaligned `DYRK.fasta`, and on `profiles/ROPK/PKL.aln` (631 sequences), run in-process so that the MAFFT and FastTree output it read could be kept under `tool-outputs/cluster/`.

### `tool-outputs/`

- `tmalign/`: TM-align output for three structure pairs, with the reported CPU time set to zero so that regeneration is reproducible.
- `hmmsearch/`: the `--tblout` table of `test/twenty.fasta` against the built profile database, without its `# Date` and `# Current dir` lines.
- `cluster/`: the standard output of each MAFFT and FastTree call made by `cluster`, in call order.

## Recorded defects

The capture records what the code did, including results that are wrong.
A port that reproduces these cases is faithful; the fixes are separate, deliberate changes, and the cases they move are listed under "Deliberate differences".

- **`update-fasta` drops sequences.**
  In `ungap_and_unique`, an X-containing sequence whose X-separated segments all occur in longer kept sequences ends the loop, and every shorter sequence is lost.
  The synthetic case `x_redundant_then_unique` shows it on four records.
  The end-to-end case `single_aln_after_edit` shows it on real data: after one record is deleted from `CDK.aln`, the rewritten `CDK.fasta` holds 24 of the 144 remaining sequences, because the 25th, `gi|193506638|pdb|3BLH|A`, contains an X.
  The current code keeps those sequences.
- **Duplicated identifiers in FASTA headers.**
  `_share.write_fasta` and `cluster.write_fasta` write the identifier before a description that already begins with it (`>seq1 seq1 some description`), and `cluster.write_fasta` writes an empty line after a header without a description.
  The current code writes the identifier once; see "FASTA headers" under "Deliberate differences".
- **The structure-alignment path fails.**
  TM-align 20240303 ends its output with a CPU-time line, so `read_tmalign_as_seqrec_pair` fails its assertion on all three recorded pairs; with that line removed it parses them.
  `mst_pairs` raises `TypeError` because Biopython 1.67 and later cannot hash a `SeqRecord`.
  Its `identity_hash_emulation` case is not this environment's behavior: it gives `SeqRecord` the identity hash of Biopython before 1.67 and shows that the function then returns all six pairs of four structures, not a three-edge spanning tree.
  The current code reads the output as written and returns a spanning tree; see "Structure alignment" under "Deliberate differences".
- **`add` duplicates sequences.**
  The recorded `add` appends every accepted hit without checking the identifiers the profile already holds, so a second run appends the same 30 records again, and `build` then fails on the duplicate identifiers when it writes the Stockholm alignment.
  It also appends one record for each domain that `hmmsearch -A` aligns in a sequence, under the same identifier; at the default `-L 0.8` no recorded sequence has two domains that pass the length filter, so the recorded run is free of this defect.
  A profile whose HMM hits none of its own sequences raises `ValueError` after the profiles visited before it were already appended to, and so does a hit whose sequence name contains `/`, since the record names of `hmmsearch -A` are split at every `/`; the recorded run exercises neither, since every profile's HMM hits some of its own sequences and no recorded identifier contains a slash.
  The temporary `<profile>.hmm.fa` file, written into the profile tree, is left there when `hmmsearch` exits with an error.
- **`cluster` fails on identifiers containing `:`.**
  Newick reserves `:`, so tip names such as `TAIR|locus:2011201` are read back truncated, and the `CDK.aln` run raises `KeyError` after writing its first clusters.
  The current code gives FastTree placeholder names; see "Identifiers containing a colon" under "Deliberate differences".
- **`cluster.geom_mean` underflows.**
  The product of 249 or more branch lengths of 0.05 is zero, as is any product containing a zero-length branch.
  FastTree gives 9 of the 631 tips of `PKL.aln` a length of zero, but collapsing the weakly supported splits above them adds the length of each split to its tips, so no tip length is zero when the threshold is computed.
  The product of the 631 tip lengths underflows to zero nonetheless, so the threshold is zero and all 418 internal clades become candidates to cut.
  The current code computes the mean in log space; see "Branch-length threshold" under "Deliberate differences".

## Interpreting a mismatch

- **Consensus sequences.**
  `fammer build` writes biofrills consensus sequences into `.cons.seq` files and aligns them, so the alignments and profiles above every family depend on biofrills.
  The Python 3 port of biofrills replaces its dictionary-order tie-breaking with an explicit rule, which changes some consensus sequences.
  Compare the `.cons.seq` files under `e2e/build/` first: a difference there comes from biofrills, and differences in the alignments that consume them follow from it.
- **External programs.**
  End-to-end results and `tool-outputs/` depend on the program versions above; a different MAFFT, HMMER or FastTree version is expected to change them.
- **Cluster numbering.**
  The recorded `cluster` numbers clusters of equal size in the dictionary order of `Bio.Phylo` clade objects, which Python 2 hashes by identity, so the numbering and the colors in its phyloXML output varied between runs.
  `e2e/cluster.json` therefore records the clusters as membership lists ordered by size and content, not by output file name, and does not record the phyloXML file.
  The current code numbers clusters of equal size by their sorted member identifiers, which gives the order of the recording on all three runs, so the replay compares the clusters in the order of their numbers.
- **File-system enumeration order.**
  Several code paths visit files in the order the file system lists them: the recorded `taskify_subdirs` lists PDB files with an unsorted `glob`, which the current code sorts, `build.all_hmm` concatenates profiles in the order `find` returns them, and `add`, `update-fasta` and `refine` walk the tree with `os.walk` or `glob`.
  These fixtures follow the order of the case-insensitive APFS volume that ran the capture.
  On another file system, the recorded order of the `.pdb` dependencies in `build-*.json`, the SHA-256 of `SomeKinases_all.hmm`, the row order of the hmmsearch table and its checksum in `scan.json`, and the order of the log lines of `add`, `refine`, `update-fasta` and the structure build all differ while their content does not.
  Compare those lists as sets, the table rows and log lines as multisets, and treat the two checksums as specific to this file system.
- **Python 2 dictionary order.**
  `scan -o` and `scan -O` write sequences in the iteration order of a dictionary, which is neither the order of the hmmsearch table nor sorted; compare those files as sets of records.
  `refine` drops two entries of a dictionary before sorting the rest, so which seeds it evaluates depends on Python 2 string hashing; the capture runs with `PYTHONHASHSEED=0`.
- **Standard error.**
  Remove the `[@line]` part of each log line and reduce a traceback to its final exception line before comparing.
- **Not captured.**
  A successful structure alignment, which no current toolchain can produce from the recorded Python 2 sources; `test_build_hmmer_with_structures` and the "Structure alignment" cases under "Deliberate differences" cover the current code instead.
  `--mapgaps` end to end, because the MAPGAPS programs are not packaged; its task graph is in `build-*.json`, and `mapgaps/` holds three of the profile files it writes before running `run_map`.
  `scan -A` and `build --align-script`, which had no effect in the recorded sources and which the current code no longer accepts.

## Deliberate differences

`replay.py` runs the recorded cases again with the current sources of this repository and writes each case whose outcome differs from the recorded one to `differences.json`.
It replays every case under `library/` except those of `hmmer.json`, which need HMMER and built profiles, the runs of `e2e/update-fasta.json`, and the runs of `e2e/cluster.json`.
The cases, and the way an outcome is compared with a recorded one, are defined in `tests/reference_cases.py`, which `tests/test_reference.py` uses as well.
A library entry gives the fixture file relative to `py27/`, the section, the case name or position, and the recorded and current outcomes.
Values are compared exactly, and an exception by its type only; a recorded built-in exception type also accepts its subclasses, so the `FileNotFoundError` of Python 3 matches a recorded `OSError`.
The cases of `library/mapgaps.json` form three sections: `consensus` compares the recorded consensus sequence with the one biofrills now computes, `mg_aln2cma` runs the function, and `biocma.write` runs biocma on the alignment with the recorded consensus prepended; the last two compare the SHA-256 of the written profile with that of the file kept under `library/mapgaps/`.
An update-fasta run is replayed in-process on a copy of `profiles/SomeKinases` with the alignments under `e2e/build/`, so it needs no external program; its entry gives the fields that differ among the log lines (without line numbers, in sorted order) and the lists of created, modified and deleted files, and, for each file whose contents differ, the record counts and checksums of both versions with the identifiers found in only one of them.
A cluster run is replayed in-process with the recorded output of MAFFT and FastTree; its entry gives the fields that differ among the outcome, the log, the files written, the members of each cluster in the order of the cluster numbers, and the members of the Unique file.
Every entry is deliberate and is explained below; a test of the current code expects the fixture value except where `differences.json` records another, and `tests/test_reference.py` fails on an entry whose case no longer differs.
After a change that moves a recorded output, regenerate the file from the repository root and explain the new entries here:

```sh
pixi run python tests/data/replay.py
```

The script needs a Git checkout, for the copies of `profiles/`, and the optional dependency biocma, which the default pixi environment installs.
The sources of the Python 3 port, commit `a0a5233`, produce an empty list.

- **Redundant sequences in `update-fasta`** (3 cases).
  `ungap_and_unique` no longer stops at the first redundant sequence that contains an X, and it tests such a sequence differently.
  An X now matches exactly one residue, and the whole sequence must occur in a single kept sequence, as a sequence without X must; the recorded code instead required each X-separated segment to occur somewhere among the kept sequences, in any of them and in any order.
  An X in a kept sequence still matches only an X.
  `x_redundant_then_unique` keeps `short_novel1` and `short_novel2`, with the same purge count.
  `x_segments_in_different_records` keeps `split`, whose segments `CCC` and `DDD` occur only in two different records, and `tail`, which the recorded code lost; nothing is purged.
  `single_aln_after_edit` writes 142 of the 144 sequences instead of 24, and purges the two structure sequences `gi|193506638|pdb|3BLH|A` and `gi|58176911|pdb|1UA2|A`, which match the full-length sequences `CDK9_Hsap` and `CDK7_Hsap` with their single X in place of a threonine.
  On every alignment under `e2e/build/` and on `profiles/ROPK/PKL.aln`, the new rule and the recorded rule with only the stopping defect corrected keep the same records.
- **Structure alignment** (5 cases).
  `read_tmalign_as_seqrec_pair` finds the alignment by the line that explains the `:` and `.` markers, instead of taking the last five lines, so it reads the output of TM-align 20240303 as written.
  On each of the three recorded pairs it returns what the recorded code returned once the final CPU-time line was removed.
  `mst_pairs` keys the graph by record identifier, the PDB file name, instead of by `SeqRecord`, and returns the record pairs on the edges of the maximum spanning tree over TM-score: N-1 pairs for N structures.
  The recorded code asked for a minimum spanning tree over the reciprocal TM-score, which selects the same edges but divides by zero on a score of zero.
  Under Biopython before 1.67 each pair of records formed its own component, so every one of the N(N-1)/2 pairs became a MAFFT seed, and under later versions it raises `TypeError`.
  Both `mst_pairs` cases now return the three edges from `A.pdb`, which carry the three highest scores, whether or not the records can be hashed.
- **FASTA headers** (7 cases).
  `update-fasta`, `build` (for `.pdb.seq` files) and `cluster` now write unwrapped FASTA through the one helper `_share.write_fasta`, which uses Biopython's `fasta-2line` format.
  A header is the record's description when the first word of the description is the identifier, as it is for every record a FASTA or Clustal parser returns, and otherwise the identifier followed by the description; sequences stay on one line.
  `cluster.write_fasta` is gone: `cluster` now ungaps the records of an alignment and passes them to `_share.write_fasta`, and its three recorded cases are replayed that way.
  `fasta_parsed`, `clustal_parsed`, `with_description` and `with_description_ungap` lose the repeated identifier, and `empty_description` loses its empty line.
  The `constructed` case of `_share.write_fasta` now raises `IndexError`, because Biopython's writer fails on its `blankdesc` record, whose description is two spaces; no parser produces such a description, since both strip the title, and the recorded code wrote `>blankdesc`.
  On the end-to-end runs, `tree_after_edit` rewrites `DYRK.fasta` with the same 63 records and headers that name each identifier once, and the current checksum of `CDK.fasta` in `single_aln_after_edit` reflects the same change.
  A structure alignment is written with headers such as `>1JST.pdb TMalign TM-score=0.826490`: the identifier is written before a description of the form `TMalign TM-score=…`.
  For three or more structures, `clean_and_dedupe_seqs` also removes the `_seed_` prefix that MAFFT adds to the first word of a seed row's description, so that word is the identifier and the header names the structure once.
- **Branch-length threshold** (6 cases and one run).
  `geom_mean` computes the mean of the logarithms, so a long input no longer underflows: `n250_of_0.05` and `n500_of_0.05` give `0.05000000000000001` instead of `0.0`, and `three_values` and `n100_of_0.05` differ from the recorded values in the last bit.
  It accepts positive numbers only, so `with_zero` raises `ValueError`, and `empty` raises `ZeroDivisionError` instead of `TypeError`.
  `get_clades_to_cut` leaves tips of length zero, which FastTree gives to a sequence identical to another, out of the mean, and never cuts a branch of length zero.
  `statistics.geometric_mean` is not used, because Python 3.12 raises an error on a zero and Python 3.13 and later return zero.
  On `PKL.aln` the threshold is 0.0701 instead of zero, and 369 of the 418 internal clades are candidates to cut; the clusters and the Unique file are unchanged, and of the lines the recorded run logged, only the count differs.
  The trees of `CDK.aln` and `DYRK.fasta`, whose products did not underflow, keep 100 and 42 candidates.
- **Identifiers containing a colon** (one run).
  `cluster` gives FastTree the names `seq0`, `seq1` and so on, in alignment order, and renames the tips of the tree it returns to the original identifiers.
  The `CDK.aln` run now succeeds: after the same eight clusters it writes `CDK.aln.Unique` with 45 sequences, including `TAIR|locus:2011761` and `TAIR|Locus:2037410`, and `CDK.aln.xml`, and logs both files.
  The recorded FastTree output names each tip by the identifier FastTree read, truncated at the first `:` or other character that ends a name, so the replay renames those tips to the names the current code gives FastTree, matching the sequences by position.
  FastTree 2.2.0 on the `CDK.aln` blocks with the placeholder names returns the recorded topology.
- **Silhouette of the clusters** (three runs).
  `cluster` now logs, after the files it wrote, the mean silhouette of the clustered sequences and the mean silhouette of each cluster, which the recorded code did not compute; the clusters, the Unique file and the files written are unchanged.
  The silhouette is computed from the patristic distances of the tree FastTree returned, before weakly supported splits are collapsed.
  `collapse_all` adds the length of a collapsed branch to each of its children, which keeps the depth of every tip but lengthens the path between two tips below the collapsed split, so distances taken afterwards would depend on the support threshold.
  The Unique sequences are not a cluster and are left out of the score.
  `CDK.aln`, `DYRK.fasta` and `PKL.aln` log mean silhouettes of 0.398, 0.630 and 0.367; computed on the collapsed trees they would be 0.393, 0.596 and 0.355.
  The values are logged to three decimals, and none of the 33 values of the three runs lies within 3.9e-5 of a rounding boundary, so a difference in the last bits of a sum between platforms does not change the log.
- **Repeated runs of `add`** (one run, not replayed).
  `add` skips the hits whose identifier any profile of the tree already holds, keeps the longest domain of a sequence with several, skips with a warning a profile that has hits that it does not yet hold but whose HMM finds none of its own sequences at 0.001 or better, and writes its temporary sequence file outside the profile tree.
  It also accepts sequence names that contain `/`, which `test_sequence_names_may_contain_a_slash` in `tests/test_add.py` covers.
  It appends its records through `_share.write_fasta`, so each sequence is on one line where the recorded run wrapped it at 60 columns.
  It also writes a line ending before the first appended record when the profile file lacks one, where the recorded code would have joined that record to the last sequence line.
  Every recorded profile ends with a line ending, so the recorded run does not exercise this.
  On `add.json` the 30 recorded records, their headers and their sequences are otherwise unchanged; a second run on the result changes no file.
  The run needs HMMER, so `replay.py` does not replay it and `differences.json` has no entry for it; `test_add_appends_the_recorded_hits_once` compares it instead.
- **Acceptance threshold of `add`** (one run, not replayed; no replayed case moves).
  `add` accepts a hit when the bit score of its best domain reaches the profile's threshold, the higher of two best-domain bit scores, each counting only the sequences that the profile's HMM hits at 0.001 or better: the lowest among the profile's own sequences, and the highest among the sequences of the subfamily profiles in other directories, leaving out consensus records and the members of catch-alls.
  The recorded code compared the hit's best-domain e-value with the highest best-domain e-value among the profile's own sequences.
  hmmsearch scales a P-value by the number of sequences searched, so the recorded rule compared an e-value computed over the 8,320 sequences of `TgME49.pep` with e-values computed over the profile's own sequences, and for a profile of n sequences demanded a P-value 8,320/n times smaller than that of the profile's weakest sequence: 213 times for Aurora's 39.
  Bit scores do not depend on the number of sequences searched.
  The weakest sequence alone would admit kinases of families that the tree has no profile for, since one divergent member sets it: the weakest sequences of `STE11`, `MAPKAPK` and `CAMK1` score 126.1, 149.0 and 101.2 bits, and with that threshold the first run of the next entry appended two NEK kinases and a CK2-family kinase to `STE11` and CDPK9 to `MAPKAPK`.
  The members of other groups set `CAMK1`'s threshold to 242.9 bits, `STE11`'s to 225.9 and `MAPKAPK`'s to 194.6.
  Sister subfamilies in one directory are left out, since a member of a sister subfamily outscores the weakest sequence of 20 of the 32 subfamilies.
  The threshold also rejects TGME49_203010, an Aurora kinase, at 254.0 bits against `Aurora`, whose weakest sequence scores 250.6 and whose threshold, set by a member of `CAMK1`, is 259.9; the weakest-sequence threshold alone would have appended it.
  The table parser now also returns the bit score from the row that wins the assignment; the replay compares the profile and e-value that the recorded code returned and leaves the score out, so its cases do not move.
  `add.filter_profile` now yields each hit's bit score instead of its e-value; the replay looks up the e-value of each selected sequence in the table, so its cases do not move either.
  `add.filter_evalue` is replaced by `add.filter_score`, which replays its four recorded cases: within one profile a smaller e-value is never a lower bit score, so each recorded e-value threshold becomes the lowest bit score among the profile's rows that it keeps, or an infinite score when it keeps none, and each case selects the recorded hits.
  The profile's own sequences are the rows of the alignment from which `build` made its HMM that are also in its FASTA file, where the recorded code scored the whole FASTA file; a record that `add` appended since the last build therefore cannot lower the threshold.
  A profile of at most ten sequences scores each of them against an HMM built from its alignment without that row.
  In the built tree every alignment holds exactly the sequences of its FASTA file, so the first of these two changes moves nothing on `add.json`.
  The one profile of ten or fewer sequences that a hit reaches there is `CMGC-Unclass`, with nine: under the assignment of the next entry, TGME49_204280 reaches its threshold, 184.0 bits, which a member of another group sets; its weakest sequence scores 117.1 bits held out and 171.1 against its own HMM, so holding sequences out does not decide the outcome.
  The hit is borderline (see "Borderline hits of `add`") and is not written.
- **Assignment of `add` hits** (one run, not replayed; no replayed case moves).
  `add` assigns each hit only among the profiles that have a FASTA file, where the recorded code assigned it among every profile of `SomeKinases_all.hmm` and so added no hit whose best match was a group-level profile: 55 of the 118 sequences of `TgME49.pep` that pass `-E 1e-20`.
  A hit whose best profile is a catch-all, a profile named with the suffix `-Unclassified`, `-Unclass` or `-Unique`, goes instead to the subfamily profile of smallest e-value, among those in the catch-all's directory or below it, whose threshold its bit score reaches, if there is one, and back to the catch-all if that subfamily then rejects it for its length (`-L`).
  With the threshold of the previous entry and the rule of "Borderline hits of `add`", the first run on `add.json` appends 21 records to 9 files, where the recorded run appended 30 to 12; `REMOVED_RECORDS` in `tests/test_integration.py` lists the recorded records that it does not append, by file, and it appends none that the recorded run did not.
  The four hits that match `CAMK-Unclass` best and reach `CDPK`'s threshold are passed to `CDPK`, and all four are borderline; one of them, TGME49_291050, an SNF1-family kinase, is the record that the recorded run appended to `CAMK-Unclass`.
  The recorded records that the run does not append are those that the threshold rejects, TGME49_209050 and TGME49_239130, TKL kinases, and TGME49_244620, a NEK kinase, from `STE11`; TGME49_249260, a CDK, and TGME49_254630, a CMGC kinase, from `MAPK`; and TGME49_218400, a NEK kinase, and TGME49_272200, an AGC kinase, from `STE-Unclassified`.
  It also does not append TGME49_255710, which the recorded run appended to `CDPK` and which is now borderline.
  Alternating `build` and `add` on that tree adds 21 records, then 4, then none, where assigning among every profile added records in ten rounds and left TGME49_244620 in both `STE11` and `STE-Unclassified`.
  `scan` still assigns among every profile, so none of its runs moves; `parse_scanned_hits` is replaced by `parse_scanned_rows`, which reads every row of the table, and `best_hits`, which selects as `parse_scanned_hits` did and which `add` applies to the rows of its own profiles, and the replay composes the two, so its cases do not move.
- **Borderline hits of `add`** (one run, not replayed; no replayed case moves).
  `add` does not write a hit that clears its profile's threshold by less than `_share.BORDERLINE_MARGIN`, 10 bits, or that does not lead the closest other profile with a FASTA file by 10 bits, not counting the catch-all that passed the hit on; it logs the hit as a `* borderline` warning instead, with both scores, where the recorded code had no such rule.
  A hit that a catch-all passed on and receives back, because the subfamily rejected it for its length, is tested again as a hit of the catch-all; on `add.json` no such hit reaches that test, since every hit that `CAMK-Unclass` passes on is borderline for `CDPK` and is removed before the catch-alls are processed.
  On `add.json` seven hits are borderline: the four that `CAMK-Unclass` passes to `CDPK`, TGME49_237890 and TGME49_255710 for `CDPK`, and TGME49_204280 for `CMGC-Unclass`.
  The recorded run appended TGME49_255710 to `CDPK` and TGME49_291050 to `CAMK-Unclass`; `CDPK` would reject TGME49_237890 for its length in any case, since the rule is applied before `-L`.
  The functions that the replay runs (`filter_profile`, `filter_score`, `parse_scanned_rows`, `best_hits`) are unchanged.
- **`refine`** (one run, not replayed; no replayed case moves).
  `refine` holds out each sequence of a profile in turn, builds an HMM from the profile's alignment without its row, and compares the best-domain bit score of the held-out sequence with the best bit score of the unclassified set against the same HMM; the lowest-scoring sequence that scores no higher is dropped and joins the unclassified set, the rest are tested again, and without `-n` the dropped records move to the unclassified file.
  The recorded code compared an e-value computed over a one-sequence file with e-values computed over the unclassified set, took the worst unclassified e-value as the threshold, tested only the seeds left after dropping two dictionary entries, stopped at the first seed it kept, raised `KeyError` when the held-out HMM missed its seed, and changed no file without `-n` while logging that it had.
  The held-out HMM is built from the rows of the alignment of the last build with `hmmbuild` alone, where the recorded code realigned the remaining seeds with MAFFT, so `refine` refuses a profile whose alignment does not hold exactly the sequences of its FASTA file, where the recorded code required only its `.hmm` file.
  It also skips catch-all profiles, leaves the consensus records of the unclassified set, those whose identifier ends with `consensus`, out of the comparison, accepts a profile name for `-x`, and refuses a sequence that is in both a profile and the unclassified set, before it changes any file.
  `Testfam-Unclassified.fasta`, the unclassified set of `refine.json`, holds no consensus record.
  A drop whose held-out score is within `_share.BORDERLINE_MARGIN` of the best unclassified score is marked `-- borderline` in its log line and still made; the recorded run drops nothing.
  `_share.get_hmm_evalues`, which only the recorded `refine` called, is removed; `test_member_evalues_reproduce_recorded_values` replays its `hmmer.json` cases through `_share.get_hmm_hits`, and they do not move.
  On `refine.json` the run tests the same two profiles, `ULK` and `Aurora`, keeps every sequence and changes no file, as recorded; for each profile it logs the kept sequence of lowest held-out score, `13031_Tt` at 147.1 bits and `GK043_Gl` at 152.8, against the best unclassified score, 42.2 and 35.5 bits, where the recorded run logged one seed of each with e-values.
