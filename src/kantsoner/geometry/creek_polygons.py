"""Helpers that prepare FKB creek lines before they are buffered into the
water: touching segments are joined into runs, and each run is trimmed back
at its ends by `CREEK_TRIM_M`.
"""

from shapely.ops import linemerge, substring, unary_union

# Taken off each end of a creek line before it is buffered.
CREEK_TRIM_M = 5.0


def _dissolve(geoms):
    """Touching creek segments joined into continuous runs.

    FKB-Vann ships a creek as a chain of short segments -- a median of 29 m
    here -- so trimming each one back would open a gap at every internal
    joint and shatter the creek rather than shorten it. Dissolving first
    means only the true ends of a creek are ever pulled back.
    """
    merged = linemerge(unary_union([g for g in geoms if g is not None]))
    if merged.is_empty:
        return []
    if merged.geom_type == "LineString":
        return [merged]
    return [g for g in merged.geoms if g.geom_type == "LineString"]


def _trim(line, trim):
    """Shorten a line by `trim` at both ends, leaving at least 2 m."""
    if line is None or line.is_empty or line.geom_type != "LineString":
        return line
    take = min(trim, max((line.length - 2.0) / 2.0, 0.0))
    if take <= 0:
        return line
    return substring(line, take, line.length - take)
