import logging
import os
import shutil
import subprocess
from os.path import exists, isfile

logger = logging.getLogger(__name__)

# A Git checkout writes files in index order and gives each the time at which
# it was written, so a tracked target can be written about a millisecond before
# the dependency it was built from (profiles/ROPK/PKL.aln before PKL.fasta).
# Within this many seconds, modification times say nothing about which file
# changed last, and a task's matches_dep decides instead. A fammer command
# takes about 0.2 seconds to start, so a dependency that update-fasta rewrites
# after an edit of the target is newer by more than this.
CHECKOUT_SKEW = 0.1


def ext(path, extension):
    """Replace the extension on a path (or iterable of paths)."""
    if isinstance(path, Task):
        return ext(path.target, extension)
    if isinstance(path, str):
        return path.rsplit(".", 1)[0] + "." + extension
    return [ext(p, extension) for p in path]


def noext(path):
    """Replace the extension on a path (or iterable of paths)."""
    if isinstance(path, Task):
        return noext(path.target)
    if isinstance(path, str):
        return path.rsplit(".", 1)[0]
    return [noext(p) for p in path]


def sh(cmd):
    """Execute a command through the shell. Print any error messages."""
    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    except OSError:
        # NB: Not sure what could trigger this; check_output eats most errors
        logger.warning("*** Failed command: %s", str(cmd))
        raise
    except subprocess.CalledProcessError as exc:
        output = exc.output.decode("utf-8", "replace")
        raise RuntimeError(
            f"Failed command (returned {exc.returncode}):\n{cmd}\n\n{output}"
        ) from exc


def which(exenames):
    """Returns the full path to an executable under any of the given names.

    The names are tried in order, and None is returned if none is found.

    Example::

        >>> which(['raxml', 'raxmlHPC', 'raxmlHPC-PTHREADS'])
        '/usr/local/bin/raxmlHPC'
    """
    for exe in exenames:
        path = shutil.which(exe)
        if path:
            return path
    return None


def is_empty(path):
    """True if a file holds no more than a single byte (has size < 2).

    A single byte is ignored because shell and I/O operations sometimes create
    a file containing only a space or an end-of-line character.
    """
    return os.stat(path).st_size < 2


def ends_with_newline(path):
    """True if the last byte of a file is a newline, or the file is empty.

    Text appended to a file whose last line has no ending would continue that
    line.
    """
    with open(path, "rb") as handle:
        if handle.seek(0, os.SEEK_END) == 0:
            return True
        handle.seek(-1, os.SEEK_END)
        return handle.read(1) == b"\n"


class Task:
    """Solve the cache-invalidation problem.

    :Parameters:
        `target`: string
            Filename or directory to build.
        `action`: callable
            Function to call to build the target.
            Implicitly uses self.depends as inputs and writes to self.target,
            though this isn't checked.
        `depends`: iterable of Task instances
            Collection of dependency Tasks or strings required to be up-to-date
            to build self.target. Converted to a list automatically, if needed.
        `cleans`: iterable of strings
            File names to be removed later. Converted to a list automatically.
        `matches_dep`: callable
            Optional function of (target, dependency) file names, true if the
            target still corresponds to the dependency. A target older than a
            dependency by less than CHECKOUT_SKEW seconds is up to date if this
            holds; without it, any older target is out of date.
    """

    @staticmethod
    def _enlist(x):
        if not x:
            return []
        if isinstance(x, Task) or isinstance(x, str):
            return [x]
        return list(x)

    def __init__(
        self,
        target,
        action=lambda x: None,
        depends=None,
        cleans=None,
        kwargs=None,
        matches_dep=None,
    ):
        # Filename or directory to build
        assert callable(action)
        self.target = target
        self.action = action
        self.depends = self._enlist(depends)
        self.kwargs = kwargs or {}
        self.cleans = self._enlist(cleans)
        self.matches_dep = matches_dep

    def build(self):
        for dep in self.depends:
            if isinstance(dep, Task) and not dep.is_up_to_date():
                dep.build()
        if exists(self.target) and all(
            self.is_current_for(dep) for dep in self.depends
        ):
            logger.info("%s: Up to date.", self.target)
        else:
            logger.info("%s: Building...", self.target)
            try:
                self.action(self, **self.kwargs)
            except:
                logger.warning("*** Failed action @ %s", self.target)
                raise

    def is_current_for(self, dep):
        """True if the target need not be rebuilt on account of one dependency.

        A target at least as new as the dependency is current; one older by
        less than CHECKOUT_SKEW seconds is current if matches_dep holds.
        """
        if not exists(self.target):
            return False
        dep_name = str(dep)
        lag = os.stat(dep_name).st_mtime - os.stat(self.target).st_mtime
        if lag <= 0:
            return True
        return (
            lag < CHECKOUT_SKEW
            and self.matches_dep is not None
            and self.matches_dep(self.target, dep_name)
        )

    def is_up_to_date(self):
        """True if the target does not need to be rebuilt."""
        return exists(self.target) and all(
            (dep.is_up_to_date() if isinstance(dep, Task) else exists(dep))
            and self.is_current_for(dep)
            for dep in self.depends
        )

    def clean(self):
        """Recursively delete auxilliary files."""
        for fname in self.cleans:
            if isfile(fname):
                logger.info("rm %s", fname)
                os.remove(fname)
        for dep in self.depends:
            if isinstance(dep, Task):
                dep.clean()

    def __str__(self):
        return self.target
