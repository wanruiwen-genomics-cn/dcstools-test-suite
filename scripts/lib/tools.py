"""Locate and run the external tools the comparison scripts depend on.

Resolution order for each tool:

1. an explicit path passed on the command line
2. the matching environment variable (SAMTOOLS, BCFTOOLS, DCS_HOME, ...)
3. ``env/tools.conf`` at the repo root, ``KEY=/path/to/binary`` per line
4. PATH

On this cluster samtools is not on PATH and bcftools lives outside PATH, so
``env/tools.conf`` is the intended place to pin them per test machine.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TOOLS_CONF = os.path.join(REPO_ROOT, "env", "tools.conf")

_conf_cache = None


class ToolNotFound(RuntimeError):
    """Raised when a required external binary cannot be located."""


class CommandFailed(RuntimeError):
    """Raised when an external command exits non-zero."""

    def __init__(self, cmd, returncode, stderr=""):
        self.cmd = cmd
        self.returncode = returncode
        self.stderr = stderr or ""
        tail = self.stderr.strip().splitlines()[-5:]
        super().__init__(
            "command failed (exit %d): %s%s"
            % (returncode, " ".join(cmd), ("\n  " + "\n  ".join(tail)) if tail else "")
        )


def _load_conf():
    global _conf_cache
    if _conf_cache is not None:
        return _conf_cache
    conf = {}
    if os.path.isfile(TOOLS_CONF):
        with open(TOOLS_CONF) as handle:
            for raw in handle:
                line = raw.split("#", 1)[0].strip()
                if not line or "=" not in line:
                    continue
                key, value = line.split("=", 1)
                conf[key.strip().upper()] = os.path.expandvars(value.strip().strip("\"'"))
    _conf_cache = conf
    return conf


def resolve(name, env_var=None, explicit=None):
    """Return an executable path for ``name`` or raise ToolNotFound."""
    env_var = env_var or name.upper()
    candidates = [
        (explicit, "--%s-bin" % name),
        (os.environ.get(env_var), "$%s" % env_var),
        (_load_conf().get(env_var), "env/tools.conf:%s" % env_var),
    ]
    for candidate, _origin in candidates:
        if not candidate:
            continue
        expanded = os.path.expanduser(os.path.expandvars(candidate))
        if os.path.isfile(expanded) and os.access(expanded, os.X_OK):
            return expanded
        found = shutil.which(expanded)
        if found:
            return found
    found = shutil.which(name)
    if found:
        return found
    raise ToolNotFound(
        "%s not found. Set $%s, add %s=<path> to %s, or put it on PATH."
        % (name, env_var, env_var, TOOLS_CONF)
    )


def run(cmd, check=True, timeout=None, env=None):
    """Run ``cmd`` and capture its output as text."""
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        timeout=timeout,
        env=env,
    )
    if check and proc.returncode != 0:
        raise CommandFailed(cmd, proc.returncode, proc.stderr)
    return proc


def stream_lines(cmd, max_lines=None):
    """Yield stdout lines of ``cmd`` without buffering the whole output.

    stderr goes to a temp file so a chatty tool cannot deadlock on a full pipe.
    Stopping early (``max_lines`` or the caller breaking out) kills the child
    instead of waiting for it.
    """
    stderr_file = tempfile.TemporaryFile(mode="w+")
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=stderr_file,
        universal_newlines=True,
        bufsize=1 << 20,
    )
    drained = False
    try:
        for index, line in enumerate(proc.stdout):
            if max_lines is not None and index >= max_lines:
                break
            yield line
        else:
            drained = True
    finally:
        try:
            if drained:
                proc.stdout.close()
                returncode = proc.wait()
                if returncode != 0:
                    stderr_file.seek(0)
                    raise CommandFailed(cmd, returncode, stderr_file.read())
            else:
                proc.stdout.close()
                proc.kill()
                proc.wait()
        finally:
            stderr_file.close()


def samtools(explicit=None):
    return resolve("samtools", "SAMTOOLS", explicit)


def bcftools(explicit=None):
    return resolve("bcftools", "BCFTOOLS", explicit)
