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
fammer calls the FastTree program `fasttree`, whereas bioconda installs it as `FastTree`; the case-insensitive default macOS file system resolves one name to the other, but on Linux a symbolic link named `fasttree` is needed.

### The MAFFT wrapper

`build.align_profiles` aligns the family consensus sequences of every directory with `mafft --globalgenafpair`.
MAFFT stopped implementing that option in 2012 and stopped accepting it at the turn of 2014, as its source releases show:

| MAFFT | `--globalgenafpair` |
|---|---|
| 6.864 (November 2011) and earlier | global pairwise alignment with generalized affine gap costs |
| 6.903 (May 2012) to 7.123 (October 2013) | accepted, but the pairwise step ignores it and runs the same global algorithm as `--globalpair` |
| 7.130 (December 2013) | prints `No longer supported!` and exits without an alignment |
| 7.141 (February 2014) and later, including every packaged build | exits with status 1 and `--globalgenaf is no longer supported!` |

The unmodified code therefore cannot build any profile tree with a MAFFT that can be installed today.
The capture puts `tool-shims/` first on `PATH`; its `mafft` script replaces `--globalgenafpair` with `--globalpair` and runs the real MAFFT, and every other MAFFT call passes through unchanged.
For the MAFFT releases contemporary with fammer's development in 2012 and 2013, this replacement computes what the original option computed.
The current code still passes the option, so `fammer build` needs the wrapper on any installable MAFFT; a change that replaces the option with something other than `--globalpair` must regenerate the references with a matching wrapper.

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
- `share-scan-add.json`: `_share.parse_scanned_hits` on the example in its docstring, on a synthetic table with ties, and on `tool-outputs/hmmsearch/twenty.fasta.SomeKinases_all.tbl`; `scan.write_summary` and `scan.write_table` output; `add.filter_profile` and `add.filter_evalue`; `_share.write_fasta`.
- `update_fasta.json`: `ungap_and_unique` on synthetic record sets and on `profiles/ROPK/PKL.aln`, and `is_same_aln_fasta` against several FASTA variants of one alignment.
- `tmalign.json`: `read_tmalign_as_seqrec_pair` on each file in `tool-outputs/tmalign/`, as written and with its final CPU-time line removed; `tmscore_from_description`, `clean_and_dedupe_seqs`, `purge_seqs` and `mst_pairs`.
- `cluster.json`: `geom_mean`, `cluster.write_fasta`, and `load_tree` on input that is neither `.aln` nor `.fasta`.
- `hmmer.json`: `_share.get_hmm_evalues`, `add.get_hmm_len` and `add.walk_profiles`, which need HMMER or built profiles and so run in the end-to-end workspace after `fammer build`.

### `py27/e2e/`

Command-line runs of `fammer.py` and `tmalign.py` in a scratch directory holding a copy of the tracked files of `profiles/SomeKinases` and the sequence files of `test/`.
Each run records its arguments, exit status, standard output, standard error, and the files it created, modified or deleted.

- `build.json`: `fammer build SomeKinases --hmmer --no-pdb --tree`, with a manifest of the 213 files written (size and SHA-256; `.hmm` files hashed without their `DATE` line; HMMER's binary `.h3?` files by size only, because they embed the date).
  The `.aln`, `.cons.seq`, `.families.fa` and `.nwk` files are kept in full under `e2e/build/`.
- `build-clean.json`: `fammer build ... --clean` on the built tree, which rebuilds nothing and deletes the intermediate files listed in the task graph.
- `scan.json`: `fammer scan` on `test/twenty.fasta` and `test/TgME49.pep`, with the summary and table outputs, the `-n`, `-x`, `-o` and `-O` options, and a stricter E-value.
- `update-fasta.json`: `fammer update-fasta` on the built tree and on one alignment, each before and after the last record of an alignment (`CDK.aln`, then `DYRK.aln`) was deleted.
- `add.json`: `fammer add SomeKinases TgME49.pep` on the built tree, with the text appended to each FASTA file.
- `refine.json`: `fammer refine SomeKinases -n -u SomeKinases/Testfam-Unclassified.fasta`.
- `pdb.json`: `fammer build` with structures enabled, and `tmalign.py` on two and on three structures.
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
- **The structure-alignment path fails.**
  TM-align 20240303 ends its output with a CPU-time line, so `read_tmalign_as_seqrec_pair` fails its assertion on all three recorded pairs; with that line removed it parses them.
  `mst_pairs` raises `TypeError` because Biopython 1.67 and later cannot hash a `SeqRecord`.
  Its `identity_hash_emulation` case is not this environment's behavior: it gives `SeqRecord` the identity hash of Biopython before 1.67 and shows that the function then returns all six pairs of four structures, not a three-edge spanning tree.
- **`cluster` fails on identifiers containing `:`.**
  Newick reserves `:`, so tip names such as `TAIR|locus:2011201` are read back truncated, and the `CDK.aln` run raises `KeyError` after writing its first clusters.
- **`cluster.geom_mean` underflows.**
  The product of 249 or more branch lengths of 0.05 is zero, as is any product containing a zero-length branch.
  FastTree gives 9 of the 631 tips of `PKL.aln` a length of zero, and the product of the other 622 also underflows to zero, so the threshold is zero and 418 internal clades become candidates to cut.

## Interpreting a mismatch

- **Consensus sequences.**
  `fammer build` writes biofrills consensus sequences into `.cons.seq` files and aligns them, so the alignments and profiles above every family depend on biofrills.
  The Python 3 port of biofrills replaces its dictionary-order tie-breaking with an explicit rule, which changes some consensus sequences.
  Compare the `.cons.seq` files under `e2e/build/` first: a difference there comes from biofrills, and differences in the alignments that consume them follow from it.
- **External programs.**
  End-to-end results and `tool-outputs/` depend on the program versions above; a different MAFFT, HMMER or FastTree version is expected to change them.
- **Cluster numbering.**
  `cluster` numbers clusters of equal size in the dictionary order of `Bio.Phylo` clade objects, which Python 2 hashes by identity, so the numbering and the colors in its phyloXML output varied between runs.
  `e2e/cluster.json` therefore records the clusters as membership lists ordered by size and content, not by output file name, and does not record the phyloXML file.
- **File-system enumeration order.**
  Several code paths visit files in the order the file system lists them: `taskify_subdirs` lists PDB files with an unsorted `glob`, `build.all_hmm` concatenates profiles in the order `find` returns them, and `add`, `update-fasta` and `refine` walk the tree with `os.walk` or `glob`.
  These fixtures follow the order of the case-insensitive APFS volume that ran the capture.
  On another file system, the order of the `.pdb` dependencies in `build-*.json`, the SHA-256 of `SomeKinases_all.hmm`, the row order of the hmmsearch table and its checksum in `scan.json`, and the order of the log lines of `add`, `refine` and the structure build all differ while their content does not.
  Compare those lists as sets, the table rows and log lines as multisets, and treat the two checksums as specific to this file system.
- **Python 2 dictionary order.**
  `scan -o` and `scan -O` write sequences in the iteration order of a dictionary, which is neither the order of the hmmsearch table nor sorted; compare those files as sets of records.
  `refine` drops two entries of a dictionary before sorting the rest, so which seeds it evaluates depends on Python 2 string hashing; the capture runs with `PYTHONHASHSEED=0`.
- **Standard error.**
  Remove the `[@line]` part of each log line and reduce a traceback to its final exception line before comparing.
- **Not captured.**
  A successful structure alignment, which no current toolchain can produce from the unmodified code.
  `--mapgaps` end to end, because the MAPGAPS programs are not packaged; its task graph is in `build-*.json`, and `mapgaps/` holds three of the profile files it writes before running `run_map`.
  `scan -A` and `build --align-script`, which have no effect.

## Deliberate differences

`replay.py` runs the recorded cases again with the current sources of this repository and writes each case whose outcome differs from the recorded one to `differences.json`.
It replays every case under `library/` except those of `hmmer.json`, which need HMMER and built profiles, and the runs of `e2e/update-fasta.json`.
The cases, and the way an outcome is compared with a recorded one, are defined in `tests/reference_cases.py`, which `tests/test_reference.py` uses as well.
A library entry gives the fixture file relative to `py27/`, the section, the case name or position, and the recorded and current outcomes.
Values are compared exactly, and an exception by its type only; a recorded built-in exception type also accepts its subclasses, so the `FileNotFoundError` of Python 3 matches a recorded `OSError`.
The cases of `library/mapgaps.json` form three sections: `consensus` compares the recorded consensus sequence with the one biofrills now computes, `mg_aln2cma` runs the function, and `biocma.write` runs biocma on the alignment with the recorded consensus prepended; the last two compare the SHA-256 of the written profile with that of the file kept under `library/mapgaps/`.
An end-to-end case is replayed in-process on a copy of `profiles/SomeKinases` with the alignments under `e2e/build/`, so it needs no external program; its entry gives the fields that differ among the log lines (without line numbers) and the lists of created, modified and deleted files, and, for each file whose contents differ, the record counts and checksums of both versions with the identifiers found in only one of them.
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
