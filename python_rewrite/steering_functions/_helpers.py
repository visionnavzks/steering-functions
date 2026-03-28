"""Shared helpers to reduce repetitive patterns in steering function modules."""

from steering_functions.configuration import Configuration
from steering_functions.utilities import global_frame_change, PI

_INF = float("inf")


# ---------------------------------------------------------------------------
# Direction-sign helpers  (replace 4-branch  if/elif/elif/else  blocks)
# ---------------------------------------------------------------------------

def _direction_sign(left, forward):
    """Return +1 when *left* == *forward*, else -1.

    Used as a multiplier for delta-y and alpha in tangent calculations.
    """
    return 1 if (left == forward) else -1


def _theta_pi(forward):
    """Return the theta offset: 0.0 when *forward*, PI when backward."""
    return 0.0 if forward else PI


def _make_tangent_point(xc, yc, theta, dx, dy, theta_offset, kappa):
    """Convenience: global_frame_change + Configuration in one call."""
    x, y = global_frame_change(xc, yc, theta, dx, dy)
    return Configuration(x, y, theta + theta_offset, kappa)


# ---------------------------------------------------------------------------
# Generic path-family evaluation  (replaces the 100+ line dispatcher blocks)
# ---------------------------------------------------------------------------

def _evaluate_families(h, c1, c2, families, arrays):
    """Evaluate path families and store results in *arrays*.

    Parameters
    ----------
    h : helper object with ``*_exists`` / ``*_path`` methods
    c1, c2 : start / end circles
    families : sequence of ``(type_id, family_name, nullable, field_names)``
        *type_id*    – integer enum value (e.g. ``tp.TT``)
        *family_name* – string like ``'TT'``, used to look up ``h.TT_exists``
                         and ``h.TT_path``
        *nullable*    – True if ``*_path`` may return ``None``
        *field_names* – tuple of array keys to unpack from the result tuple
                        (positions 1, 2, … after the leading length element)
    arrays : dict  ``{'length': [...], 'cstart': [...], ... }``
    """
    for type_id, name, nullable, fields in families:
        if not getattr(h, f'{name}_exists')(c1, c2):
            continue
        r = getattr(h, f'{name}_path')(c1, c2)
        if nullable and r is None:
            continue
        arrays['length'][type_id] = r[0]
        for i, field in enumerate(fields):
            arrays[field][type_id] = r[i + 1]


def _make_path_arrays(n=18):
    """Allocate the standard arrays used by ``_evaluate_families``."""
    return {
        'length': [_INF] * n,
        'qi1': [None] * n, 'qi2': [None] * n,
        'qi3': [None] * n, 'qi4': [None] * n,
        'cstart': [None] * n, 'cend': [None] * n,
        'ci1': [None] * n, 'ci2': [None] * n,
    }


def _select_best_path(arrays, n, c1_start, c2_start, kappa, sigma,
                       path_cls, type_cls):
    """Pick the shortest path from *arrays* and return a path object."""
    L = arrays['length']
    best = min(range(n), key=lambda i: L[i])
    return path_cls(
        c1_start, c2_start, type_cls(best),
        kappa, sigma,
        arrays['qi1'][best], arrays['qi2'][best],
        arrays['qi3'][best], arrays['qi4'][best],
        arrays['cstart'][best], arrays['cend'][best],
        arrays['ci1'][best], arrays['ci2'][best],
        L[best],
    )


# ---------------------------------------------------------------------------
# Table-driven get_controls  (replaces the 18-branch elif chains)
# ---------------------------------------------------------------------------

def _dispatch_controls(path, table, control_fns):
    """Build the control list for *path* using a lookup *table*.

    Parameters
    ----------
    path : a path object with attributes ``type``, ``cstart``, ``cend``,
           ``qi1``–``qi4``, ``ci1``–``ci2``, ``start``, ``end``
    table : dict  ``{ type_id: [ (fn_key, *args), ... ] }``
        Each entry is a list of steps.  A step is a tuple whose first element
        selects the function via *control_fns*, and remaining elements are
        attribute names on *path* (resolved with ``getattr``).
        The last element for turn functions is the boolean *order* flag
        (not an attribute name).
    control_fns : dict  ``{ fn_key: callable }``
    """
    controls = []
    steps = table.get(path.type)
    if steps is None:
        return controls
    for step in steps:
        fn_key = step[0]
        fn = control_fns[fn_key]
        if fn_key == 'empty':
            fn(controls)
        elif fn_key == 'straight':
            fn(getattr(path, step[1]), getattr(path, step[2]), controls)
        else:
            # turn function: (circle_attr, config_attr, order_bool)
            fn(getattr(path, step[1]), getattr(path, step[2]),
               step[3], controls)
    return controls
