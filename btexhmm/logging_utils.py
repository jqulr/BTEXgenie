import io
import shlex
import subprocess
import sys
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from pathlib import Path


VISIBLE_PREFIXES = ("[*]", "[info]", "[warn]", "[err]")


class FilteredLogStream(io.TextIOBase):
    def __init__(self, log_fh, terminal_stream):
        self.log_fh = log_fh
        self.terminal_stream = terminal_stream
        self._buffer = ""

    def write(self, s):
        if not s:
            return 0
        self._buffer += s
        while True:
            idx = self._buffer.find("\n")
            if idx == -1:
                break
            line = self._buffer[: idx + 1]
            self._buffer = self._buffer[idx + 1 :]
            self._emit(line)
        return len(s)

    def flush(self):
        if self._buffer:
            self._emit(self._buffer)
            self._buffer = ""
        self.log_fh.flush()
        self.terminal_stream.flush()

    def _emit(self, line):
        self.log_fh.write(line)
        stripped = line.lstrip()
        if stripped.startswith(VISIBLE_PREFIXES):
            self.terminal_stream.write(line)


@contextmanager
def command_logger(log_path: Path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as log_fh:
        out_stream = FilteredLogStream(log_fh, sys.__stdout__)
        err_stream = FilteredLogStream(log_fh, sys.__stderr__)
        with redirect_stdout(out_stream), redirect_stderr(err_stream):
            yield log_path
        out_stream.flush()
        err_stream.flush()


def run_logged_command(cmd):
    print(f"[cmd] {shlex.join([str(part) for part in cmd])}", flush=True)
    proc = subprocess.Popen(
        [str(part) for part in cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        sys.stdout.write(line)
    rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)
