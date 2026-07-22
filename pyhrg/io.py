"""
pyhrg.io — reading CHM rasters and exporting crowns and tree tops.

All disk access lives here; :mod:`pyhrg.hrg`, :mod:`pyhrg.chm` and
:mod:`pyhrg.treetops` work on plain arrays. rasterio and fiona are
imported lazily so that the array-based API stays usable without them.

Copyright (C) 2025 Igor Pawelec
Licence: GPLv3 — see LICENSE.
"""

import os
import numpy as np


# ── Lazy imports ──────────────────────────────────────────────────────

__all__ = [
    "read_chm", "save_crowns_raster", "save_segments", "save_tree_tops",
]


def _ensure_rasterio():
    try:
        import rasterio
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "rasterio is required for reading and writing rasters.\n"
            "  conda install -c conda-forge rasterio"
        ) from e
    return rasterio


def read_chm(path, window=None, band=1):
    """
    Read a CHM raster.

    Parameters
    ----------
    path : str or Path
        GeoTIFF (or any GDAL-readable raster).
    window : tuple, optional
        ``(col_off, row_off, width, height)`` to read a sub-region only.
        The returned transform is shifted accordingly, so exported crowns
        stay georeferenced. Use this to work tile-by-tile on rasters that
        would not fit in memory.
    band : int
        Band index, 1-based.

    Returns
    -------
    chm : ndarray, shape (rows, cols), dtype float32
    transform : affine.Affine
    crs : rasterio.crs.CRS

    Notes
    -----
    The raster nodata value is *not* converted to NaN — the growing works
    on raw values, and a nodata of -9999 would simply fall below any sane
    ``mask_thresh``. Convert explicitly if your nodata is a plausible
    height.
    """
    rasterio = _ensure_rasterio()
    with rasterio.open(path) as src:
        if window is not None:
            from rasterio.windows import Window
            if isinstance(window, (list, tuple)) and len(window) == 4:
                window = Window(*window)
            chm = src.read(band, window=window)
            transform = src.window_transform(window)
        else:
            chm = src.read(band)
            transform = src.transform
        crs = src.crs
    return chm.astype(np.float32, copy=False), transform, crs


def _ensure_fiona():
    try:
        import fiona
        return fiona
    except ImportError as e:
        raise ImportError(
            "fiona is required for vector I/O.\n"
            "Install with:  conda install -c conda-forge fiona"
        ) from e
    except OSError as e:
        raise OSError(
            "fiona found but failed to load (DLL/shared library error).\n"
            "Fix:  conda install -c conda-forge fiona --force-reinstall\n"
            f"Original error: {e}"
        ) from e


def _ensure_rasterio_features():
    try:
        from rasterio.features import shapes
        return shapes
    except ImportError as e:
        raise ImportError(
            "rasterio is required for raster vectorization.\n"
            "Install with:  conda install -c conda-forge rasterio"
        ) from e


def _ensure_morphology():
    from skimage.morphology import closing, disk
    return closing, disk


# ── Public API ────────────────────────────────────────────────────────

def save_segments(segments: np.ndarray,
                  out_path: str,
                  fname: str,
                  transform,
                  crs_wkt: str,
                  chm_array: np.ndarray,
                  closing_radius: int = 0,
                  driver: str = "ESRI Shapefile") -> None:
    """
    Save crown segments as vector file + raw raster dump.

    Uses single-pass vectorization with precomputed per-crown statistics.
    Much faster than per-segment iteration for large numbers of crowns.

    Outputs
    -------
    - RAW raster: {fname}.bin + {fname}.vrt
    - Vector file with attributes: id, max_height, area_m2, crown_diameter
    """
    fiona = _ensure_fiona()
    shapes = _ensure_rasterio_features()

    # --- 1) RAW dump (.bin) + VRT ---
    bin_path = os.path.join(out_path, f"{fname}.bin")
    vrt_path = os.path.join(out_path, f"{fname}.vrt")
    segments.astype(np.int32).tofile(bin_path)

    rows, cols = segments.shape
    pixel_area = abs(transform.a * transform.e)
    gtx = (transform.c, transform.a, 0, transform.f, 0, -transform.e)
    with open(vrt_path, "w") as f:
        f.write(f'<VRTDataset rasterXSize="{cols}" rasterYSize="{rows}">\n')
        f.write(f'  <SRS>{crs_wkt}</SRS>\n')
        f.write(f'  <GeoTransform>{",".join(map(str, gtx))}</GeoTransform>\n')
        f.write('  <VRTRasterBand dataType="Int32" band="1">\n')
        f.write(f'    <SourceFilename relativeToVRT="1">{fname}.bin</SourceFilename>\n')
        f.write('    <ImageOffset>0</ImageOffset>\n')
        f.write('    <PixelOffset>4</PixelOffset>\n')
        f.write(f'    <LineOffset>{4*cols}</LineOffset>\n')
        f.write('  </VRTRasterBand>\n')
        f.write('</VRTDataset>\n')

    # --- 2) Optional morphological closing on the FULL raster ---
    seg_data = segments.astype(np.int32)
    if closing_radius > 0:
        morph_closing, disk = _ensure_morphology()
        # Close each segment mask — but do it efficiently via label dilation
        # For small numbers of crowns, per-label closing is acceptable
        # For large numbers, we close the binary mask then re-label
        closed = np.zeros_like(seg_data)
        seg_ids = np.unique(seg_data)
        seg_ids = seg_ids[seg_ids != 0]
        selem = disk(closing_radius)
        for sid in seg_ids:
            m = morph_closing(seg_data == sid, selem)
            closed[m & (closed == 0)] = sid
        seg_data = closed

    # --- 3) Precompute per-crown stats (one pass) ---
    max_id = int(seg_data.max()) + 1 if seg_data.max() > 0 else 1
    pixel_counts = np.bincount(seg_data.ravel(), minlength=max_id)

    # Max height per crown — use np.maximum.at (in-place, no copy)
    max_heights = np.full(max_id, -np.inf, dtype=np.float32)
    flat_seg = seg_data.ravel()
    flat_chm = chm_array.ravel()
    if flat_chm.dtype != np.float32:
        flat_chm = flat_chm.astype(np.float32)
    np.maximum.at(max_heights, flat_seg, flat_chm)
    max_heights[0] = 0.0  # background

    # --- 4) Single-pass vectorization ---
    ext_map = {"ESRI Shapefile": ".shp", "GPKG": ".gpkg", "GeoJSON": ".geojson"}
    ext = ext_map.get(driver, ".shp")
    vec_path = os.path.join(out_path, f"{fname}{ext}")

    schema = {
        'geometry': 'Polygon',
        'properties': {
            'id': 'int',
            'max_height': 'float',
            'area_m2': 'float',
            'crown_diameter': 'float'
        }
    }

    seg_mask = seg_data > 0
    with fiona.open(
        vec_path,
        'w',
        driver=driver,
        crs_wkt=crs_wkt,
        schema=schema
    ) as dst:
        for geom, val in shapes(seg_data, mask=seg_mask, transform=transform):
            seg_id = int(val)
            if seg_id <= 0:
                continue

            n_px = int(pixel_counts[seg_id])
            area = n_px * pixel_area
            diam = 2.0 * np.sqrt(area / np.pi)
            max_h = float(max_heights[seg_id])

            dst.write({
                'geometry': geom,
                'properties': {
                    'id': seg_id,
                    'max_height': round(max_h, 2),
                    'area_m2':    round(area, 2),
                    'crown_diameter': round(diam, 2)
                }
            })


def save_tree_tops(corrected_tops: np.ndarray,
                   out_path: str,
                   fname: str,
                   transform,
                   crs_wkt: str,
                   chm: np.ndarray,
                   driver: str = "ESRI Shapefile") -> None:
    """
    Save corrected tree tops as point vector file.

    Parameters
    ----------
    corrected_tops : ndarray (n, 2)
        Tree top positions as (row, col) pixel coordinates.
    out_path, fname, transform, crs_wkt, chm : see save_segments.
    driver : str
        Fiona driver: "ESRI Shapefile" (default), "GPKG", "GeoJSON".
    """
    fiona = _ensure_fiona()

    ext_map = {"ESRI Shapefile": ".shp", "GPKG": ".gpkg", "GeoJSON": ".geojson"}
    ext = ext_map.get(driver, ".shp")
    vec_path = os.path.join(out_path, fname + f"_treetops{ext}")

    schema = {
        'geometry': 'Point',
        'properties': {
            'id': 'int',
            'height': 'float'
        }
    }

    coords = np.array(corrected_tops, dtype=float)
    rows = coords[:, 0].astype(int)
    cols = coords[:, 1].astype(int)
    heights = chm[rows, cols]
    heights = np.round(heights, 2)

    with fiona.open(
        vec_path,
        'w',
        driver=driver,
        crs_wkt=crs_wkt,
        schema=schema
    ) as dst:
        for idx, (r, c) in enumerate(coords):
            # +0.5 puts the point at the pixel's centre. An affine transform
            # maps grid coordinates where whole numbers land on pixel
            # *corners*, while a tree top from center_of_mass is an array
            # index, and array indices refer to pixel centres. Without the
            # shift every exported point sat half a pixel up and to the left
            # -- 0.25 m on the 0.5 m test rasters, systematic, and enough to
            # matter when comparing against field-measured stems. This is
            # what rasterio.transform.xy does.
            x, y = transform * (float(c) + 0.5, float(r) + 0.5)
            geom = {"type": "Point", "coordinates": (x, y)}
            dst.write({
                'geometry': geom,
                'properties': {
                    'id': idx,
                    'height': float(heights[idx])
                }
            })


def save_crowns_raster(crowns: np.ndarray,
                       output_path: str,
                       transform,
                       crs_wkt: str,
                       nodata: int = 0,
                       compress: str = "deflate") -> None:
    """
    Save crown label raster as GeoTIFF.

    Parameters
    ----------
    crowns : ndarray (int32)
        Crown label raster. 0 = background.
    output_path : str
        Output GeoTIFF path.
    transform : affine.Affine
        Geotransform.
    crs_wkt : str
        CRS as WKT string.
    nodata : int
        NoData value. Default 0.
    compress : str
        Compression: "deflate", "lzw", "zstd", "none". Default "deflate".
    """
    try:
        import rasterio
        from rasterio.crs import CRS
    except ImportError as e:
        raise ImportError(
            "rasterio required. Install: conda install -c conda-forge rasterio"
        ) from e

    rows, cols = crowns.shape
    profile = {
        'driver': 'GTiff',
        'dtype': 'int32',
        'width': cols,
        'height': rows,
        'count': 1,
        'crs': CRS.from_wkt(crs_wkt) if isinstance(crs_wkt, str) else crs_wkt,
        'transform': transform,
        'nodata': nodata,
    }
    if compress and compress.lower() != "none":
        profile['compress'] = compress

    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(crowns.astype(np.int32), 1)
