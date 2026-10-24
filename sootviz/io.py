"""Loading soot-foil images and giving them a physical scale."""
from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Foil:
    """A soot-foil record with a physical scale.

    image : 2-D float array, square pixels, row 0 at the top
    s     : pixels per physical unit
    width, height : physical extent
    unit  : label only ("mm", "m", ...)
    """
    image: np.ndarray
    s: float
    width: float
    height: float
    unit: str = "mm"
    name: str = "foil"

    @property
    def extent(self):
        return (0.0, self.width, self.height, 0.0)

    def px(self, length: float) -> int:
        """Physical length -> pixels, at least one."""
        return max(1, int(round(length * self.s)))

    def __repr__(self):
        ny, nx = self.image.shape
        return (f"Foil({self.name!r}, {nx}x{ny} px, {self.width}x{self.height} "
                f"{self.unit}, {self.s:.1f} px/{self.unit})")


def _square(arr: np.ndarray, width: float, height: float):
    """Resample so that pixels are square; angles are then geometric."""
    from PIL import Image
    ny, nx = arr.shape
    sx, sy = nx / width, ny / height
    if abs(sx - sy) / max(sx, sy) < 0.005:
        return arr, sx
    s = max(sx, sy)
    im = Image.fromarray(arr.astype(np.float32), mode="F")
    out = im.resize((int(round(width * s)), int(round(height * s))),
                    Image.LANCZOS)
    return np.asarray(out, dtype=float), s


def load_foil(path, width: float, height: float, unit: str = "mm",
              pdf_page: int | None = None, pdf_index: int = 0,
              invert: bool = False, name: str | None = None) -> Foil:
    """Load a soot foil and attach its physical extent.

    path       image file, or a PDF from which a raster is extracted
    width,
    height     physical size of the IMAGE (not of the figure around it)
    pdf_page   1-based page to pull rasters from (requires poppler's pdfimages)
    pdf_index  which raster on that page, ordered as the page stores them
    invert     negate the image (use when the published figure is inverted)
    """
    from PIL import Image
    path = Path(path)
    if path.suffix.lower() == ".pdf":
        arr = _from_pdf(path, pdf_page, pdf_index)
    else:
        arr = np.asarray(Image.open(path).convert("L"), dtype=float)
    if invert:
        arr = arr.max() - arr
    arr, s = _square(arr, width, height)
    return Foil(arr, s, width, height, unit, name or path.stem)


def _from_pdf(path: Path, page: int | None, index: int) -> np.ndarray:
    from PIL import Image
    with tempfile.TemporaryDirectory() as td:
        cmd = ["pdfimages", "-png"]
        if page is not None:
            cmd += ["-f", str(page), "-l", str(page)]
        cmd += [str(path), f"{td}/img"]
        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except FileNotFoundError:  # pragma: no cover
            raise RuntimeError("pdfimages not found: install poppler-utils, "
                               "or export the figure as PNG first")
        files = sorted(Path(td).glob("img*.png"))
        if not files:
            raise RuntimeError(f"no raster images found in {path}"
                               + (f" page {page}" if page else ""))
        rasters = [f for f in files if min(Image.open(f).size) > 64]
        if index >= len(rasters):
            raise IndexError(f"page has {len(rasters)} rasters, asked for "
                             f"index {index}")
        return np.asarray(Image.open(rasters[index]).convert("L"), dtype=float)


def list_pdf_rasters(path, page: int | None = None):
    """Report the rasters on a PDF page, to choose `pdf_index`."""
    from PIL import Image
    with tempfile.TemporaryDirectory() as td:
        cmd = ["pdfimages", "-png"]
        if page is not None:
            cmd += ["-f", str(page), "-l", str(page)]
        cmd += [str(path), f"{td}/img"]
        subprocess.run(cmd, check=True, capture_output=True)
        out = []
        for i, f in enumerate(sorted(Path(td).glob("img*.png"))):
            im = Image.open(f)
            if min(im.size) > 64:
                out.append((len(out), im.size, im.size[0] / im.size[1]))
        return out


def synth_foil(width=30.0, height=20.0, s=25.0, lam=2.5, aspect=1.8,
               skew=0.0, jitter=0.12, noise=0.25, contrast=0.6,
               vignette=0.5, unit="mm", seed=0) -> Foil:
    """A synthetic foil with a KNOWN cell width, for validating the pipeline.

    lam      cell width (the quantity the pipeline should recover)
    aspect   cell length / cell width
    skew     degrees added to one family, so the lattice is not symmetric
    jitter   random displacement of the lattice, in units of lam
    noise    white noise amplitude relative to the track contrast
    contrast track darkness relative to the background
    vignette large-scale illumination variation (0 = none)
    """
    rng = np.random.default_rng(seed)
    nx, ny = int(width * s), int(height * s)
    X, Y = np.meshgrid(np.arange(nx) / s, np.arange(ny) / s)
    L = aspect * lam
    a1 = np.degrees(np.arctan(lam / L)) + skew
    a2 = -np.degrees(np.arctan(lam / L))
    img = np.ones((ny, nx))
    for ang in (a1, a2):
        t = np.tan(np.radians(ang))
        phase = (Y - t * X) / lam
        if jitter:
            w = rng.normal(0, jitter, size=(8, 8))
            from scipy import ndimage as ndi
            big = ndi.zoom(w, (ny / 8, nx / 8), order=3)[:ny, :nx]
            phase = phase + big
        d = np.abs(((phase + 0.5) % 1.0) - 0.5)          # distance to a track
        img -= contrast * np.exp(-(d * lam / (0.05 * lam)) ** 2 / 2)
    img -= img.min()
    img /= img.max() or 1.0
    if vignette:
        img *= 1.0 - vignette * (0.5 + 0.5 * np.cos(2 * np.pi * X / width))
    img += noise * contrast * rng.normal(size=img.shape)
    rng_ = img.max() - img.min()
    img = 40 + 180 * (img - img.min()) / (rng_ if rng_ > 0 else 1.0)
    return Foil(img, s, width, height, unit, f"synthetic(lam={lam})")
