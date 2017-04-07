# External dependencies
from __future__ import division, absolute_import, print_function
import numpy as np

# Internal dependencies
from svgpathtools import *


def directional_field(curve, tvals=3, asize=1e-2,
                      colored=False):
    """Creates/Displays a directional_field for the given curve.

    Parameters
    ----------

    curve: `Path` or any path segment object
        the curve to find the directional field for.

    tvals: iterable or `int`
        if  tvals is a `int`, then will determine the number of arrows per
        segment.  Otherwise, an arrow will be placed at `curve.point(t)` for
        each `t` in `tvals`.

    asize: float-like or int-like
        relative size of the arrowheads (relative to `curve.length()`).

    colored: `bool`
        if true (defaults `False`), will return a 3-tuple of (arrowhead paths,
        tvals, colors).

    Returns
    -------

    A path composed of arrowheads, unless `colored` is true, see above for that
    case.

    Example Usages
    -------------_

    >>> p = Path(Line(0,1), Line(1,1+1j), Line(1+1j,1j), Line(1j,1))
    >>> disvg([p, directional_field(p)])

    >>> p = Path(Line(0,1), Line(1,1+1j), Line(1+1j,1j), Line(1j,1))
    >>> arrows, tvals, colors  = directional_field(p, colored=True)
    >>> disvg([p]+arrows, ['green']+colors)

    """

    if isinstance(tvals, int):
        if isinstance(curve, Path):
            tvals = [curve.t2T(k, t) for k in range(tvals)
                                      for t in np.linspace(0, 1, tvals)]
        else:
            tvals = np.linspace(0, 1, tvals)

    size = asize * curve.length()

    arrows = []
    for t in tvals:
        pt = curve.point(t)
        ut = curve.unit_tangent(t)
        un = curve.normal(t)
        l1 = Line(pt, pt + size*(un - ut)/2).reversed()
        l2 = Line(pt, pt + size*(-un - ut)/2)
        if colored:
            arrows.append(Path(l1, l2))
        else:
            arrows += [l1, l2]

    if colored:
        colors = [(int(255*t), 0, 0) for t in tvals]
        return arrows, tvals, colors
    else:
        return Path(*arrows)



p = Path(Line(0,1), Line(1,1+1j), Line(1+1j,1j), Line(1j,0))
disvg([p, directional_field(p)])
