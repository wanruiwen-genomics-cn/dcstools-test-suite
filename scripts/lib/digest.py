"""Order-sensitive and order-insensitive streaming digests.

Both are computed in a single pass with constant memory, which matters for 30x
WGS BAMs.

The order-insensitive digest is the sum of the per-record MD5 values modulo
2**128. Two files agree on it exactly when they contain the same multiset of
records, regardless of the order they appear in. That distinction is the useful
one for dcstools: multi-threaded sort/markdup can legitimately reorder records
that share a coordinate, so "same records, different order" deserves a
different verdict than "different records".
"""

from __future__ import annotations

import hashlib

_MOD = 1 << 128


class StreamDigest(object):
    def __init__(self):
        self._ordered = hashlib.md5()
        self._unordered = 0
        self.count = 0

    def update(self, canonical):
        data = canonical.encode("utf-8", "surrogateescape")
        self._ordered.update(data)
        self._ordered.update(b"\n")
        self._unordered = (
            self._unordered + int.from_bytes(hashlib.md5(data).digest(), "big")
        ) % _MOD
        self.count += 1

    @property
    def ordered(self):
        return self._ordered.hexdigest()

    @property
    def unordered(self):
        return "%032x" % self._unordered

    def as_dict(self):
        return {"records": self.count, "ordered": self.ordered, "unordered": self.unordered}
