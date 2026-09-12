"""Generate a deterministic placeholder figure for celltypist.

This produces a single PNG that is then read + re-saved via Pillow
(matching worker C §0.5.5 "PNG hash 稳定" requirement: png::writePNG(readPNG(f), f)
on the R side, Pillow Image.save on the Python side).

Output:
    resources/built-in-plugins/celltypist/tools/figures/annotation_placeholder.png

Run:
    /Applications/anaconda3/envs/biof3-py-runtime/bin/python3 tools/figures/generate_placeholder_figure.py

Idempotent: re-running produces byte-identical output (modulo the read+save
hash-rewrite loop), so the file's sha256 stays stable for Result Studio
artifact-lineage.
"""

from __future__ import annotations

import io
import struct
import zlib
from pathlib import Path

OUT_PATH = Path(__file__).parent / "annotation_placeholder.png"


def _png_chunk(tag: bytes, data: bytes) -> bytes:
    length = struct.pack(">I", len(data))
    crc = struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)
    return length + tag + data + crc


def _build_minimal_png(width: int = 256, height: int = 64) -> bytes:
    """Build a minimal valid PNG: 256×64 transparent strip with text-style grey.

    Deterministic, no timestamps, no compression metadata that varies.
    """
    sig = b"\x89PNG\r\n\x1a\n"
    ihdr_data = struct.pack(">IIBBBBB", width, height, 8, 6, 0, 0, 0)  # RGBA, 8-bit, no interlace
    ihdr = _png_chunk(b"IHDR", ihdr_data)

    # Build scanlines: filter byte 0 (None) + RGBA pixels
    raw = bytearray()
    for _ in range(height):
        raw.append(0)
        for x in range(width):
            # violet gradient (celltypist brand color theme)
            intensity = (x % 64) * 3 + 128
            r = (intensity + ((x // 8) % 16)) % 256
            g = (intensity // 2 + 64) % 256
            b = (intensity + 48) % 256
            a = 255
            raw.extend((r, g, b, a))
    idat = _png_chunk(b"IDAT", zlib.compress(bytes(raw), level=9))
    iend = _png_chunk(b"IEND", b"")
    return sig + ihdr + idat + iend


def _canonicalize_png(buf: bytes) -> bytes:
    """Read + save via Pillow so the final PNG hash matches Pillow's encoder
    (parallels R png::writePNG(png::readPNG(f), f))."""
    from PIL import Image

    img = Image.open(io.BytesIO(buf))
    img.load()
    out = io.BytesIO()
    img.save(out, format="PNG", optimize=False)
    return out.getvalue()


def main() -> None:
    raw = _build_minimal_png()
    canonical = _canonicalize_png(raw)
    OUT_PATH.write_bytes(canonical)
    import hashlib

    print(
        f"[placeholder figure] wrote {OUT_PATH} ({len(canonical)} bytes, "
        f"sha256-prefix={hashlib.sha256(canonical).hexdigest()[:12]})"
    )


if __name__ == "__main__":
    main()
