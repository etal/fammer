import logging
import os
import subprocess
from os.path import exists, isfile


def ext(path, extension):
    """Replace the extension on a path (or iterable of paths)."""
    if isinstance(path, Task):
        return ext(path.target, extension)
    if isinstance(path, basestring):
        return path.rsplit('.', 1)[0] + '.' + extension
    return [ext(p, extension) for p in path]


def noext(path):
    """Replace the extension on a path (or iterable of paths)."""
    if isinstance(path, Task):
        return noext(path.target)
    if isinstance(path, basestring):
        return path.rsplit('.', 1)[0]
    return [noext(p) for p in path]


def sh(cmd):
    """Execute a command through the shell. Print any error messages."""
    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    except OSError:
        # NB: Not sure what could trigger this; check_output eats most errors
        logging.warning("*** Failed command: %s", str(cmd))
        raise
    except subprocess.CalledProcessError, exc:
        raise RuntimeError("Failed command (returned %s):\n%s\n\n%s"
                           % (exc.returncode, str(cmd), exc.output))


def which(exenames):
    """Returns the full path to an executable under any of the given names.

    Example::

        >>> which(['raxml', 'raxmlHPC', 'raxmlHPC-PTHREADS'])
        '/usr/local/bin/raxmlHPC'
    """
    for exe in exenames:
        try:
            path = subprocess.check_output(['which', exe])
            return path.strip()
        except subprocess.CalledProcessError:
            pass


def is_empty(path):
    """True if a file contains more than a single byte (has size >1).

    We ignore a single byte because shell and I/O operations sometimes create a
    file containing a space or EOL character, and that's meaningless.
    """
    return os.stat(path).st_size < 2


class Task(object):
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
    """

    @staticmethod
    def _enlist(x):
        if not x:
            return []
        if isinstance(x, Task) or isinstance(x, basestring):
            return [x]
        return list(x)

    def __init__(self, target, action=lambda x: None, depends=None, cleans=None,
            kwargs=None):
        # Filename or directory to build
        assert callable(action)
        self.target = target
        self.action = action
        self.depends = self._enlist(depends)
        self.kwargs = kwargs or {}
        self.cleans = self._enlist(cleans)

    def build(self):
        for dep in self.depends:
            if isinstance(dep, Task) and not dep.is_up_to_date():
                dep.build()
        if exists(self.target) and all(self.is_newer_than(dep)
                for dep in self.depends):
            logging.info("%s: Up to date.", self.target)
        else:
            logging.info("%s: Building...", self.target)
            try:
                self.action(self, **self.kwargs)
            except:
                logging.warning("*** Failed action @ %s", self.target)
                raise

    def is_newer_than(self, dep):
        if not exists(self.target):
            return False
        return (os.stat(self.target).st_mtime
                >= os.stat(str(dep)).st_mtime)

    def is_up_to_date(self):
        """True if the target does not need to be rebuilt."""
        return exists(self.target) and all(
                (dep.is_up_to_date()
                    if isinstance(dep, Task)
                    else exists(dep))
                and self.is_newer_than(dep)
                for dep in self.depends)

    def clean(self):
        """Recursively delete auxilliary files."""
        for fname in self.cleans:
            if isfile(fname):
                logging.info('rm %s', fname)
                os.remove(fname)
        for dep in self.depends:
            if isinstance(dep, Task):
                dep.clean()

    def __cmp__(self, other):
        """Enable sorting alphabetically by target name."""
        return cmp(self.target, str(other))

    def __str__(self):
        return self.target

