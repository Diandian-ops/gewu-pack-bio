"""Generate a deterministic placeholder figure for cellchat.

This produces a single 256×64 PNG that is then read + re-saved via Pillow
(matching worker C §0.5.5 "PNG hash 稳定" requirement: png::writePNG(readPNG(f), f)
on the R side, Pillow Image.save on the Python side).

Run:
    python3 resources/built-in-plugins/cellchat/tools/figures/generate_placeholder_figure.py

Idempotent: re-running produces byte-identical output (modulo the read+save
hash-rewrite loop), so the file's sha256 stays stable for Result Studio
artifact-lineage.
"""

from __future__ import annotations

import io
import struct
import zlib
from pathlib import Path

OUT_PATH = Path(__file__).parent / "heatmap_placeholder.png"


def _png_chunk(tag: bytes, data: bytes) -> bytes:
    length = struct.pack(">I", len(data))
    crc = struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)
    return length + tag + data + crc


def _build_minimal_png(width: int = 256, height: int = 64) -> bytes:
    """Build a minimal valid PNG: 256×64 RGBA strip with deterministic gradient.

    No timestamps, no compression metadata that varies; transparent padding
    so a downstream Pillow Image.open + Image.save round-trip stays stable.
    """
    sig = b"\x89PNG\r\n\x1a\n"
    ihdr_data = struct.pack(">IIBBBBB", width, height, 8, 6, 0, 0, 0)  # RGBA, 8-bit, no interlace
    ihdr = _png_chunk(b"IHDR", ihdr_data)

    raw = bytearray()
    for _ in range(height):
        raw.append(0)
        for x in range(width):
            intensity = (x % 64) * 3 + 128
            r = g = b = (intensity + ((x // 8) % 16)) % 256
            a = 255
            raw.extend((r, g, b, a))

    idat = _png_chunk(b"IDAT", zlib.compress(bytes(raw), level=9))
    iend = _png_chunk(b"IEND", b"")
    return sig + ihdr + idat + iend


def main() -> Path:
    png_bytes = _build_minimal_png()
    OUT_PATH.write_bytes(png_bytes)
    # Pillow read+save rewrite (R png::writePNG(png::readPNG(f), f) equivalent)
    try:
        from PIL import Image

        with Image.open(OUT_PATH) as img:
            img.load()
            target = img.convert("RGBA")
            target.save(OUT_PATH, format="PNG", optimize=False)
    except ImportError:
        # Pillow not available; keep raw bytes
        pass

    size = OUT_PATH.stat().st_size
    print(f"[cellchat-placeholder] Wrote {OUT_PATH} ({size} bytes)")
    return OUT_PATH


if __name__ == "__main__":
    main()