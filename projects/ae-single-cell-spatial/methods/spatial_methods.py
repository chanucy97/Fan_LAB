"""Reusable spatial scoring and boundary-distance methods.

Adapted from the study's spatial program-scoring/directional-extent workflow.
Gene lists, masks, region labels, calibration and all measurements are supplied
by the caller. This module contains no study observations or saved results.
"""
import numpy as np
from scipy import ndimage as ndi


def program_score(counts, total_counts, reference_mask, scale_quantile=0.99):
    """Score a caller-selected gene matrix (bins x genes), without gene labels.

    Uses log1p(CP10k), each gene's positive-reference quantile, clipping to
    [0, 1], then an unweighted mean. Returns score and detected-gene support.
    The caller must select the gene set before calling this function.
    """
    x = np.asarray(counts, dtype=float)
    totals = np.asarray(total_counts, dtype=float)
    reference = np.asarray(reference_mask, dtype=bool)
    if x.ndim != 2 or x.shape[1] == 0 or totals.shape != (len(x),) or reference.shape != (len(x),):
        raise ValueError('Expected bins-by-genes counts and aligned total/reference vectors')
    if not np.all(np.isfinite(x)) or np.any(x < 0) or not np.all(np.isfinite(totals)) or np.any(totals < 0):
        raise ValueError('Counts and totals must be finite and non-negative')
    if np.any(x.sum(axis=1) > totals + 1e-8) or not reference.any() or not 0 < scale_quantile <= 1:
        raise ValueError('Check library totals, reference mask and quantile')
    normalized = np.log1p(1e4 * x / np.maximum(totals[:, None], 1.0))
    scaled = np.zeros_like(normalized)
    for j in range(x.shape[1]):
        positive = normalized[reference & (normalized[:, j] > 0), j]
        scale = float(np.quantile(positive, scale_quantile)) if positive.size else 1.0
        if not np.isfinite(scale) or scale <= 0:
            scale = 1.0
        scaled[:, j] = np.clip(normalized[:, j] / scale, 0.0, 1.0)
    return scaled.mean(axis=1), (x > 0).sum(axis=1)


def directional_extent(score, support, valid, region_labels, bin_um,
                       score_quantile=0.80, min_detected_genes=2,
                       min_component_bins=6, seed_band_um=75.0):
    """Measure boundary-connected signal in eight centroid-referenced sectors.

    Region labels are arbitrary caller-provided positive integers; zero denotes
    exterior. Distance is nearest bin-centre distance minus half a bin. Regions
    compete for the nearest exterior territory. Each sector is analysed using
    eight-neighbour components, a minimum size and a boundary seed-band rule.
    Pad the grid with an invalid frame to represent the analysis-field edge.
    The returned list contains freshly computed values, never saved study data.
    """
    score = np.asarray(score, dtype=float)
    support = np.asarray(support, dtype=float)
    valid = np.asarray(valid, dtype=bool)
    labels = np.asarray(region_labels)
    if score.ndim != 2 or any(a.shape != score.shape for a in (support, valid, labels)):
        raise ValueError('All inputs must be aligned two-dimensional grids')
    if not np.isfinite(bin_um) or bin_um <= 0 or not 0 <= score_quantile <= 1:
        raise ValueError('Invalid spatial calibration or score quantile')
    if np.any(labels < 0) or np.any(labels != labels.astype(int)) or not np.any(labels > 0):
        raise ValueError('Region labels must be non-negative integers with a positive region')
    if not np.all(np.isfinite(support[valid])) or np.any(support[valid] < 0):
        raise ValueError('Detected-gene support must be finite and non-negative')
    if min_detected_genes < 0 or min_component_bins < 1 or seed_band_um < 0:
        raise ValueError('Invalid support, component or seed-band parameter')
    if valid[0].any() or valid[-1].any() or valid[:, 0].any() or valid[:, -1].any():
        raise ValueError('Pad the grid with an invalid frame to represent field-edge censoring')
    reference = score[valid & np.isfinite(score)]
    if reference.size == 0:
        raise ValueError('No finite scores in the analysis mask')
    threshold = float(np.quantile(reference, score_quantile))
    union = labels > 0
    _, indices = ndi.distance_transform_edt(~union, return_indices=True)
    nearest = labels[indices[0], indices[1]]
    structure = np.ones((3, 3), dtype=np.uint8)
    edge_adjacent = ndi.binary_dilation(~valid, structure=structure)
    yy, xx = np.indices(score.shape)
    directions = ('E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE')
    candidate = valid & np.isfinite(score) & (score >= threshold) & (support >= min_detected_genes) & ~union
    rows = []
    for code in np.unique(labels[union]):
        mask = labels == code
        cy, cx = ndi.center_of_mass(mask)
        angle = (np.degrees(np.arctan2(cy - yy, xx - cx)) + 360.0) % 360.0
        wedges = np.floor(((angle + 22.5) % 360.0) / 45.0).astype(int)
        distance = np.maximum(ndi.distance_transform_edt(~mask, sampling=bin_um) - bin_um / 2.0, 0.0)
        territory = valid & ~union & (nearest == code)
        for wedge, direction in enumerate(directions):
            sector = territory & (wedges == wedge)
            components, number = ndi.label(candidate & sector, structure=structure)
            keep = np.zeros_like(valid)
            for component in range(1, number + 1):
                member = components == component
                if member.sum() >= min_component_bins and np.any(member & (distance <= seed_band_um)):
                    keep |= member
            width = float(distance[keep].max()) if keep.any() else 0.0
            envelope = sector & (distance <= width) if width > 0 else np.zeros_like(valid)
            rows.append({'region_label': int(code), 'direction': direction,
                         'score_threshold': threshold, 'width_um': width,
                         'connected_area_um2': float(keep.sum() * bin_um**2),
                         'envelope_area_um2': float(envelope.sum() * bin_um**2),
                         'right_censored': bool(np.any(keep & edge_adjacent))})
    return rows
