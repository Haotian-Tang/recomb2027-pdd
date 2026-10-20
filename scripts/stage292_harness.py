#!/usr/bin/env python3
"""Stage292 harness: hashes, serialisation and isolated forked job execution.

House pattern (Stage44): one forked child per job with RLIMIT_AS and a hard
wall clock; every run is retained, including timeouts and errors.
"""
import hashlib
import json
import os
import resource
import signal
import sys
import tempfile
import time
import traceback
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with open(ROOT / path if not os.path.isabs(str(path)) else path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set, frozenset)):
        return [encode(v) for v in x]
    return x


def isolated_run(function, payload, wall_seconds, memory_bytes):
    """Run function(payload) in a forked child; always return a record."""
    with tempfile.TemporaryDirectory(prefix='stage292_') as temporary:
        path = Path(temporary) / 'result.json'
        start = time.perf_counter()
        pid = os.fork()
        if pid == 0:
            try:
                resource.setrlimit(resource.RLIMIT_AS, (memory_bytes, memory_bytes))
                answer = dict(status='ok', value=function(payload))
            except BaseException as exc:                     # retained, never dropped
                answer = dict(status='error', error_type=type(exc).__name__,
                              error=repr(exc), traceback=traceback.format_exc())
            try:
                path.write_text(json.dumps(encode(answer)))
            finally:
                os._exit(0)
        timed_out = False
        while True:
            waited, _, usage = os.wait4(pid, os.WNOHANG)
            if waited:
                break
            if time.perf_counter() - start >= wall_seconds:
                timed_out = True
                os.kill(pid, signal.SIGKILL)
                os.wait4(pid, 0)
                break
            time.sleep(0.005)
        elapsed = time.perf_counter() - start
        if timed_out:
            return dict(status='timeout', wall_seconds=elapsed,
                        wall_limit=wall_seconds)
        if not path.exists():
            return dict(status='died', wall_seconds=elapsed)
        answer = json.loads(path.read_text())
        answer['wall_seconds'] = elapsed
        return answer

