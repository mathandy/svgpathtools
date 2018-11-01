"""This submodule contains the class definitions of the the main five classes
svgpathtools is built around: Path, Line, QuadraticBezier, CubicBezier, and Arc."""

# External dependencies
from __future__ import division, absolute_import, print_function

# jps added tan, atan:
from math import sqrt, cos, sin, acos, degrees, radians, log, pi, tan, atan2
from cmath import exp, sqrt as csqrt, phase
from warnings import warn
from operator import itemgetter
import numpy as np
try:
    from scipy.integrate import quad
    _quad_available = True
except ImportError:   # jps added 'ImportError'
    _quad_available = False

# Internal dependencies
from .bezier    import bezier_intersections, bezier_bounding_box, split_bezier
from .bezier    import bezier_by_line_intersections, polynomial2bezier, bezier2polynomial
from .misctools import BugException
from .polytools import rational_limit, polyroots, polyroots01, imag, real


# Default Parameters  ##########################################################

# path segment .length() parameters for arc length computation
LENGTH_MIN_DEPTH  = 5
LENGTH_ERROR      = 1e-12
USE_SCIPY_QUAD    = True  # for elliptic Arc segment arc length computation

# path segment .ilength() parameters for inverse arc length computation
ILENGTH_MIN_DEPTH = 5
ILENGTH_ERROR     = 1e-12
ILENGTH_S_TOL     = 1e-12
ILENGTH_MAXITS    = 10000

# compatibility / implementation related warnings and parameters
CLOSED_WARNING_ON = True

# d-string printing defaults:
SUBPATH_TO_SUBPATH_SPACE = ' '
COMMAND_TO_NUMBER_SPACE = ''

_NotImplemented4ArcException = \
    Exception("This method has not yet been implemented for Arc objects.")
# _NotImplemented4QuadraticException = \
#     Exception("This method has not yet been implemented for QuadraticBezier "
#               "objects.")
_is_smooth_from_warning = \
    ("The name of this method is somewhat misleading (yet kept for "
     "compatibility with scripts created using svg.path 2.0).  This method "
     "is meant only for d-string creation and should NOT be used to check "
     "for kinks.  To check a segment for differentiability, use the "
     "joins_smoothly_with() method instead or the kinks() function (in "
     "smoothing.py).\nTo turn off this warning, set "
     "warning_on=False.")


# Miscellaneous  ###############################################################


def bbox2path(xmin, xmax, ymin, ymax):
    """Converts a bounding box 4-tuple to a Path object."""
    b = Line(xmin + 1j * ymin, xmax + 1j * ymin)
    t = Line(xmin + 1j * ymax, xmax + 1j * ymax)
    r = Line(xmax + 1j * ymin, xmax + 1j * ymax)
    l = Line(xmin + 1j * ymin, xmin + 1j * ymax)
    return Path(b, r, t.reversed(), l.reversed())


def is_smooth_join(coming_from, going_to, wrt_parameterization=False):
    """
    Checks if coming_from path / subpath / segment joins smoothly with
    going_to path / subpath / segment.  By default, this only checks that
    going_to starts moving (at t=0) in the same direction (and from the same
    positive) as coming_from stopped moving (at t=1).  To check if the tangent
    magnitudes also match, set wrt_parameterization=True.
    """
    if wrt_parameterization:
        return going_to.start == coming_from.end and going_to.derivative(
            0) == coming_from.derivative(1)
    else:
        return going_to.start == coming_from.end and going_to.unit_tangent(
            0) == coming_from.unit_tangent(1)


def segment_iterator_of(thing, back_to_front=False):
    if isinstance(thing, Segment):
        return [thing]

    if isinstance(thing, Subpath):
        if back_to_front:
            return reversed(thing)
        return thing

    if isinstance(thing, Path):
        return thing.segment_iterator(back_to_front=back_to_front)

    assert False


# Conversion  ###################################################################


def bezier_segment(*bpoints):
    """ (publicly-minded function ?) """
    if len(bpoints) == 2:
        return Line(*bpoints)
    elif len(bpoints) == 4:
        return CubicBezier(*bpoints)
    elif len(bpoints) == 3:
        return QuadraticBezier(*bpoints)
    else:
        assert len(bpoints) in (2, 3, 4)


def poly2bez(poly, return_bpoints=False):
    """Converts a cubic or lower order Polynomial object (or a sequence of
    coefficients) to a CubicBezier, QuadraticBezier, or Line object as
    appropriate.  If return_bpoints=True then this will instead only return
    the control points of the corresponding Bezier curve.
    Note: The inverse operation is available as a method of CubicBezier,
    QuadraticBezier and Line objects."""
    bpoints = polynomial2bezier(poly)
    if return_bpoints:
        return bpoints
    else:
        return bezier_segment(*bpoints)


def bez2poly(bez, numpy_ordering=True, return_poly1d=False):
    """Converts a Bezier object or tuple of Bezier control points to a tuple
    of coefficients of the expanded polynomial.
    return_poly1d : returns a numpy.poly1d object.  This makes computations
    of derivatives / anti - derivatives and many other operations quite quick.
    numpy_ordering : By default (to accommodate numpy) the coefficients will
    be output in reverse standard order.
    Note:  This function is redundant thanks to the .poly() method included
    with all bezier segment classes."""
    if isinstance(bez, BezierSegment):
        bez = bez.bpoints()
    return bezier2polynomial(bez,
                             numpy_ordering=numpy_ordering,
                             return_poly1d=return_poly1d)


# Geometric  ####################################################################


def rotate(curve, degs, origin=None):
    """Returns curve rotated by `degs` degrees (CCW) around the point `origin`
    (a complex number).  By default origin is either `curve.point(0.5)`, or in
    the case that curve is an Arc object, `origin` defaults to `curve.center`.
    """
    def transform(z):
        return exp(1j * radians(degs)) * (z - origin) + origin

    if origin is None:
        if isinstance(curve, Arc):
            origin = curve.center
        else:
            origin = curve.point(0.5)

    if isinstance(curve, Path):
        return Path(*[rotate(subpath, degs, origin=origin) for subpath in curve.subpath_iterator()])
    elif isinstance(curve, Subpath):
        to_return = Subpath(*[rotate(seg, degs, origin=origin) for seg in curve])
        if curve.Z():
            to_return.set_Z()
        return to_return
    elif isinstance(curve, BezierSegment):
        return bezier_segment(*[transform(bpt) for bpt in curve.bpoints()])
    elif isinstance(curve, Arc):
        new_start = transform(curve.start)
        new_end = transform(curve.end)
        new_rotation = curve.rotation + degs
        return Arc(new_start, radius=curve.radius, rotation=new_rotation,
                   large_arc=curve.large_arc, sweep=curve.sweep, end=new_end)
    else:
        raise TypeError("Input `curve` should be a Path, Line, "
                        "QuadraticBezier, CubicBezier, or Arc object.")


def translate(curve, z0):
    """Shifts the curve by the complex quantity z such that
    translate(curve, z0).point(t) = curve.point(t) + z0"""
    if isinstance(curve, Path):
        return Path(*[translate(subpath, z0) for subpath in curve.subpath_iterator()])
    elif isinstance(curve, Subpath):
        to_return = Subpath(*[translate(seg, z0) for seg in curve])
        if curve.Z():
            to_return.set_Z()
        return to_return
    elif isinstance(curve, BezierSegment):
        return bezier_segment(*[bpt + z0 for bpt in curve.bpoints()])
    elif isinstance(curve, Arc):
        new_start = curve.start + z0
        new_end = curve.end + z0
        return Arc(new_start, radius=curve.radius, rotation=curve.rotation,
                   large_arc=curve.large_arc, sweep=curve.sweep, end=new_end)
    else:
        raise TypeError("Input `curve` should be a Path, Line, "
                        "QuadraticBezier, CubicBezier, or Arc object.")


def bezier_unit_tangent(seg, t):
    """Returns the unit tangent of the segment at t.

    Notes
    - - - - -
    If you receive a RuntimeWarning, try the following:
    >>> import numpy
    >>> old_numpy_error_settings = numpy.seterr(invalid='raise')
    This can be undone with:
    >>> numpy.seterr(**old_numpy_error_settings)
    """
    def compute_error_message(place):
        bef = seg.poly().deriv()(t - 1e-4)
        aft = seg.poly().deriv()(t + 1e-4)
        mes = ("thrown at %s in bezier_unit_tangent:" % place +
               "unit tangent appears to not be well-defined at " +
               "t = {}, \n".format(t) +
               "seg.poly().deriv()(t - 1e-4) = {}\n".format(bef) +
               "seg.poly().deriv()(t + 1e-4) = {}".format(aft))
        return mes

    assert 0 <= t <= 1
    dseg = seg.derivative(t)

    # Note: dseg might be numpy value, use np.seterr(invalid='raise')
    try:
        unit_tangent = dseg / abs(dseg)
    except (ZeroDivisionError, FloatingPointError):
        # [Author 1] a previous author suggested this
        # but---alas!---csqrt(x^2) is not equal to x for all x in C
        # (i.e., the square root might be incorrect):

        if False:  #
            # [Prev Author] This may be a removable singularity, if so we just need
            # to compute the limit.
            # Note: limit{{dseg / abs(dseg)} = sqrt(limit{dseg**2 / abs(dseg)**2})
            dseg_poly = seg.poly().deriv()
            dseg_abs_squared_poly = (real(dseg_poly)**2 +
                                     imag(dseg_poly)**2)
            try:
                return csqrt(rational_limit(dseg_poly**2, dseg_abs_squared_poly, t))
            except ValueError:
                raise ValueError(compute_error_message("pt A"))

        # [Author 1]...so here is another hopeful (partial) fix:
        if True:
            bpoints = list(seg.bpoints)
            if len(bpoints) > 2 and t == 0 and np.isclose(bpoints[0], bpoints[1]):
                try:
                    dif = bpoints[2] - bpoints[0]
                    unit_tangent = dif / abs(dif)
                except (ZeroDivisionError, FloatingPointError):
                    raise ValueError(compute_error_message("pt B"))

            elif len(bpoints) > 2 and t == 0 and np.isclose(bpoints[-2], bpoints[-1]):
                try:
                    dif = bpoints[-1] - bpoints[-3]
                    unit_tangent = dif / abs(dif)
                except (ZeroDivisionError, FloatingPointError):
                    raise ValueError(compute_error_message("pt C"))

            else:
                raise ValueError(compute_error_message("pt D"))

    return unit_tangent


def segment_curvature(self, t, use_inf=False):
    """returns the curvature of the segment at t.

    Notes
    - - - - -
    If you receive a RuntimeWarning, run command
    >>> old = np.seterr(invalid='raise')
    This can be undone with
    >>> np.seterr(**old)
    """
    dz = self.derivative(t)
    ddz = self.derivative(t, n=2)
    dx, dy = dz.real, dz.imag
    ddx, ddy = ddz.real, ddz.imag
    old_np_seterr = np.seterr(invalid='raise')
    try:
        kappa = abs(dx * ddy - dy * ddx) / sqrt(dx * dx + dy * dy)**3
    except (ZeroDivisionError, FloatingPointError):
        # tangent vector is zero at t, use polytools to find limit
        p = self.poly()
        dp = p.deriv()
        ddp = dp.deriv()
        dx, dy = real(dp), imag(dp)
        ddx, ddy = real(ddp), imag(ddp)
        f2 = (dx * ddy - dy * ddx)**2
        g2 = (dx * dx + dy * dy)**3
        lim2 = rational_limit(f2, g2, t)
        if lim2 < 0:  # impossible, must be numerical error
            return 0
        kappa = sqrt(lim2)
    finally:
        np.seterr(**old_np_seterr)
    return kappa


def bezier_radialrange(seg, origin, return_all_global_extrema=False):
    """returns the tuples (d_min, t_min) and (d_max, t_max) which minimize and
    maximize, respectively, the distance d = |self.point(t) - origin|.
    return_all_global_extrema:  Multiple such t_min or t_max values can exist.
    By default, this will only return one. Set return_all_global_extrema=True
    to return all such global extrema."""

    def _radius(tau):
        return abs(seg.point(tau) - origin)

    shifted_seg_poly = seg.poly() - origin
    r_squared = real(shifted_seg_poly)**2 + imag(shifted_seg_poly)**2
    extremizers = [0, 1] + polyroots01(r_squared.deriv())
    extrema = [(_radius(t), t) for t in extremizers]

    if return_all_global_extrema:
        raise NotImplementedError
    else:
        seg_global_min = min(extrema, key=itemgetter(0))
        seg_global_max = max(extrema, key=itemgetter(0))
        return seg_global_min, seg_global_max


def closest_point_in_path(pt, path):
    """returns (|path.seg.point(t) - pt|, t, seg_idx) where t and seg_idx
    minimize the distance between pt and curve path[idx].point(t) for 0<=t<=1
    and any seg_idx.
    Warning:  Multiple such global minima can exist.  This will only return
    one."""
    return path.radialrange(pt)[0]


def farthest_point_in_path(pt, path):
    """returns (|path.seg.point(t)-pt|, t, seg_idx) where t and seg_idx
    maximize the distance between pt and curve path[idx].point(t) for 0<=t<=1
    and any seg_idx.
    : rtype : object
    : param pt:
    : param path:
    Warning:  Multiple such global maxima can exist.  This will only return
    one."""
    return path.radialrange(pt)[1]


def segment_length(curve, start, end, start_point, end_point,
                   error=LENGTH_ERROR, min_depth=LENGTH_MIN_DEPTH, depth=0):
    """Recursively approximates the length by straight lines"""
    mid = (start + end) / 2
    mid_point = curve.point(mid)
    length = abs(end_point - start_point)
    first_half = abs(mid_point - start_point)
    second_half = abs(end_point - mid_point)

    length2 = first_half + second_half
    if (length2 - length > error) or (depth < min_depth):
        # Calculate the length of each segment:
        depth += 1
        return (segment_length(curve, start, mid, start_point, mid_point,
                               error, min_depth, depth) +
                segment_length(curve, mid, end, mid_point, end_point,
                               error, min_depth, depth))
    # This is accurate enough.
    return length2


def inv_arclength(curve, s, s_tol=ILENGTH_S_TOL, maxits=ILENGTH_MAXITS,
                  error=ILENGTH_ERROR, min_depth=ILENGTH_MIN_DEPTH):
    """INPUT: curve should be a CubicBezier, Line, of Path of CubicBezier
    and / or Line objects.
    OUTPUT: Returns a float, t, such that the arc length of curve from 0 to
    t is approximately s.
    s_tol - exit when |s(t) - s| < s_tol where
        s(t) = seg.length(0, t, error, min_depth) and seg is either curve or,
        if curve is a Path object, then seg is a segment in curve.
    error - used to compute lengths of cubics and arcs
    min_depth - used to compute lengths of cubics and arcs
    Note:  This function is not designed to be efficient, but if it's slower
    than you need, make sure you have scipy installed."""

    curve_length = curve.length(error=error, min_depth=min_depth)
    assert curve_length > 0
    if not 0 <= s <= curve_length:
        raise ValueError("s is not in interval [0, curve.length()].")

    if s == 0:
        return 0
    if s == curve_length:
        return 1

    if isinstance(curve, Path):
        seg_lengths = [seg.length(error=error, min_depth=min_depth) for seg in curve]
        lsum = 0
        # Find which segment the point we search for is located on
        for k, len_k in enumerate(seg_lengths):
            if lsum <= s <= lsum + len_k:
                t = inv_arclength(curve[k], s - lsum, s_tol=s_tol, maxits=maxits, error=error, min_depth=min_depth)
                return curve.t2T(k, t)
            lsum += len_k
        return 1

    elif isinstance(curve, Line):
        return s / curve.length(error=error, min_depth=min_depth)

    elif (isinstance(curve, QuadraticBezier) or
          isinstance(curve, CubicBezier) or
          isinstance(curve, Arc)):
        t_upper = 1
        t_lower = 0
        iteration = 0
        while iteration < maxits:
            iteration += 1
            t = (t_lower + t_upper) / 2
            s_t = curve.length(t1=t, error=error, min_depth=min_depth)
            if abs(s_t - s) < s_tol:
                return t
            elif s_t < s:  # t too small
                t_lower = t
            else:  # s < s_t, t too big
                t_upper = t
            if t_upper == t_lower:
                warn("t is as close as a float can be to the correct value, "
                     "but |s(t) - s| = {} > s_tol".format(abs(s_t - s)))
                return t
        raise Exception("Maximum iterations reached with s(t) - s = {}."
                        "".format(s_t - s))
    else:
        raise TypeError("First argument must be a Line, QuadraticBezier, "
                        "CubicBezier, Arc, or Path object.")


def multisplit_of(croppable, ts):
    """takes a possibly empty list of t-values t1, ..., tn such that
    0 < t1 < ... < tn < 1 and returns a list of the same type as 'croppable'
    whose union is 'croppable' and whose endpoints are self.point(0),
    self.point(t1), ..., self.point(tn); 'croppable' can be a segment, a
    subpath or a path!"""
    assert all(x < y for x, y in zip(ts, ts[1:]))
    assert len(ts) == 0 or (ts[0] > 0 and ts[-1] < 1)
    ts.insert(0, 0)
    ts.append(1)
    pieces = []
    for t0, t1 in zip(ts, ts[1:]):
        pieces.append(croppable.cropped(t0, t1))
    return pieces


def solve_complex_linear_congruence(a, z, b, w):
    # solves a + l * z = b + m * w for real l, m
    M = np.matrix(
        [[real(z), -real(w)],
         [imag(z), -imag(w)]])
    try:
        lm = M.I * np.matrix([[real(b - a)], [imag(b - a)]])
    except np.linalg.linalg.LinAlgError:
        raise ValueError
    return lm.item(0, 0), lm.item(1, 0)


def offset_intersection(a, b, c, offset_amount):
    normal_a_b = (b - a) * (-1j) / abs(b - a)
    normab_b_c = (c - b) * (-1j) / abs(c - b)
    A = a + normal_a_b * offset_amount
    C = c + normab_b_c * offset_amount

    # first possible system:
    # A + l * (b - a) = C + m * (b - c)

    # second possible system:
    # A + l * (b - a) = b + m * (sum_of_normals)

    sum_of_normals = normal_a_b + normab_b_c

    first_system_determinant = imag((b - a) * (b - c).conjugate())
    secnd_system_determinant = imag((b - a) * sum_of_normals.conjugate())

    try:
        if (abs(first_system_determinant) > abs(secnd_system_determinant)):
            l, m = solve_complex_linear_congruence(A, b - a, C, b - c)
            return A + l * (b - a)
        else:
            l, m = solve_complex_linear_congruence(A, b - a, b, sum_of_normals)
            return A + l * (b - a)
    except ValueError:
        raise


def divergence_of_offset(p1, p2, putatitve_offset, safety=10, early_return_threshold=None):
    safety = int(max(1, safety))
    t_min  = 0
    max_distance = 0
    # t = 0, 1 are not tested (not useful):
    for t in [(x / (safety + 1)) for x in range(1, safety + 1)]:
        extruded = p1.point(t) + p1.normal(t) * putatitve_offset
        if t_min > 0:
            __, p2 = p2.split(t_min)  # keep the sequence of closest points in monotone order
        (d_min, t_min) = closest_point_in_path(extruded, p2)
        assert d_min == abs(p2.point(t_min) - extruded)
        max_distance = max(max_distance, d_min)
        if early_return_threshold is not None and max_distance > early_return_threshold:
            return max_distance
    return max_distance


def segment_offset_as_subpath(segment, amount, two_sided=False, quality=0.01, safety=10):
    segs = [segment, segment.reversed()] if two_sided else [segment]
    naive = []
    max_divergence = 0
    for s in segs:
        naive.append(s.naive_offset(amount))
        divergence = divergence_of_offset(s, naive[-1], amount, early_return_threshold=quality, safety=safety)
        max_divergence = max(divergence, max_divergence)
        if max_divergence > quality:
            break
    if max_divergence <= quality:
        way_out  = Subpath(naive[0])
        way_in   = Subpath(naive[1]) if two_sided else Subpath()
        skeleton = Subpath(segment)
    else:
        s1, s2 = segment.split(0.5)
        wo1, sk1, wi1 = segment_offset_as_subpath(s1, amount, two_sided, quality=quality, safety=safety)
        wo2, sk2, wi2 = segment_offset_as_subpath(s2, amount, two_sided, quality=quality, safety=safety)
        way_out  = wo1.extend(wo2)
        way_in   = wi2.extend(wi1)
        skeleton = sk1.extend(sk2)
    return way_out, skeleton, way_in


def compute_offset_joining_subpath(seg1, off1, seg2, off2, offset_amount, join='miter', miter_limit=4):
    """returns a triple of the form 'subpath, t1, t2' where 'subpath' is the connecting
    subpath, t1 and t2 are the new end and start points of off1, off2, respectively"""
    assert all(isinstance(x, Segment) for x in [seg1, off1, seg2, off2])
    assert seg1.end == seg2.start
    assert join in ['miter', 'round', 'bevel']

    if join == 'bevel':
        join = 'miter'
        miter_limit = 1

    assert join in ['miter', 'round']

    base_corner = seg1.end

    n1 = seg1.normal(1)
    n2 = seg2.normal(0)

    theoretical_start = seg1.end   + n1 * offset_amount
    theoretical_end   = seg2.start + n2 * offset_amount

    if not np.isclose(theoretical_start, off1.end):
        print("seg1:", seg1)
        print("off1:", off1)

    assert np.isclose(theoretical_start, off1.end)  # also checks offset_amount!
    assert np.isclose(theoretical_end, off2.start)

    tangent1_base = seg1.unit_tangent(1)
    tangent2_base = seg2.unit_tangent(0)

    tangent1_offset = off1.unit_tangent(1)
    tangent2_offset = off2.unit_tangent(0)

    assert np.isclose(tangent1_base, tangent1_offset)
    assert np.isclose(tangent2_base, tangent2_offset)

    # note: real(w * z.conjugate()) is the dot product of w, z as vectors

    if offset_amount * real(n1 * tangent2_base.conjugate()) > 0:
        # acute case
        intersections = off1.intersect(off2)
        if len(intersections) > 0:
            a1, a2 = intersections[-1]
            t1, t2 = a1.t, a2.t
            assert 0 <= t1 <= 1 and 0 <= t2 <= 1
            assert np.isclose(off1.point(t1), off2.point(t2))
            if t2 > 0 and t1 < 1:
                return None, t1, t2
        return Subpath(Line(off1.end, off2.start)), 1, 0

    if offset_amount * real(n1 * tangent2_base.conjugate()) < 0:
        # obtuse case
        if join == 'miter':
            a = base_corner - tangent1_base
            b = base_corner
            c = base_corner + tangent2_base

            apex = offset_intersection(a, b, c, offset_amount)
            z1 = apex - off1.end
            z2 = apex - off2.start
            assert real(z1 * tangent1_base.conjugate()) > 0
            assert real(z2 * tangent2_base.conjugate()) < 0
            assert np.isclose(imag(z1 * tangent1_base.conjugate()), 0)
            assert np.isclose(imag(z2 * tangent2_base.conjugate()), 0)

            miter = abs(apex - base_corner) / abs(offset_amount)
            if miter > miter_limit:
                return Subpath(Line(off1.end, off2.start)), 1, 0
            else:
                # the only case we actually need a path instead of a segment:
                return Subpath(Line(off1.end, apex),
                               Line(apex, off2.start)), 1, 0

        assert join == 'round'

        r1 = abs(off1.end - base_corner)
        r2 = abs(off2.start - base_corner)

        assert np.isclose(r1, r2)

        sweep = 1 if imag(tangent1_base * tangent2_base.conjugate()) < 0 else 0

        return Subpath(Arc(off1.end, r1 + 1j * r1, 0, 0, sweep, off2.start)), 1, 0

    assert real(n1 * tangent2_base.conjugate()) == 0
    assert np.isclose(off1.end, off2.start)
    assert False


def join_offset_segments_into_subpath(skeleton, offsets, putative_amount, join, miter_limit):
    """internal function, assumes all the following:"""
    assert isinstance(skeleton, Subpath) and isinstance(offsets, Path)
    assert len(skeleton) == offsets.num_segments()
    assert skeleton.is_bezier()

    if len(skeleton) == 0:
        return Subpath()

    assert all(isinstance(thing, Segment) for thing in skeleton)
    assert all(isinstance(thing, Subpath) for thing in offsets)

    segment_pairs = [(u, v) for u, v in zip(skeleton,
                                            offsets.segment_iterator())]

    if skeleton.Z():
        segment_pairs.append((skeleton[0], offsets[0][0]))

    successive_pairs = zip(segment_pairs, segment_pairs[1:])

    to_return = Subpath(offsets[0][0])

    for index, (u, v) in enumerate(successive_pairs):
        is_loop_around_iteration = skeleton.Z() and index == len(segment_pairs) - 2
        assert isinstance(u[0], Segment)
        assert isinstance(v[0], Segment)
        assert isinstance(v[1], Segment)
        assert isinstance(u[1], Segment)

        seg1 = u[0]
        seg2 = v[0]
        off1 = u[1] if index == 0 else to_return[-1]
        off2 = v[1] if not is_loop_around_iteration else to_return[0]

        assert all(isinstance(x, Segment) for x in [seg1, seg2, off1, off2])

        if np.isclose(off1.end, off2.start):
            if off1.end != off2.start:
                o2 = off2.tweaked(start=off1.end)
            else:
                o2 = off2

        else:
            p, t1, t2 = compute_offset_joining_subpath(seg1, off1, seg2, off2, putative_amount, join=join, miter_limit=miter_limit)

            assert t1 > 0 and t2 < 1
            o1 = off1 if t1 == 1 else off1.cropped(0, t1)
            o2 = off2 if t2 == 0 else off2.cropped(t2, 1)

            if p is None:
                if t1 == 1 and t2 == 0:
                    p = Subpath(Line(off1.end, off2.start))
                else:
                    assert np.isclose(o1.end, o2.start)
                    o2 = o2.tweaked(start=o1.end)

            to_return[-1] = o1  # (overwrite previous )

            if p is not None:
                assert p.start == o1.end
                assert p.end == o2.start
                to_return.extend(p)

        if not is_loop_around_iteration:
            to_return.append(o2)
        else:
            to_return[0] = o2
            to_return.set_Z(forceful=False)

    assert to_return.Z() == skeleton.Z()

    return to_return


def endcap_for_curve(p, offset_amount, cap_style):
    n     = p.normal(1)
    start = p.end + n * offset_amount
    end   = p.end - n * offset_amount
    t     = p.unit_tangent(1) * offset_amount

    if cap_style == 'square':
        l1 = Line(start, start + t)
        l2 = Line(start + t, end + t)
        l3 = Line(end + t, end)
        return Subpath(l1, l2, l3)

    if cap_style == 'round':
        mid = p.end + t
        a1 = Arc(start, offset_amount + 1j * offset_amount, 0, 0, 1, mid)
        a2 = Arc(mid,   offset_amount + 1j * offset_amount, 0, 0, 1, end)
        return Subpath(a1, a2)

    if cap_style == 'butt':
        return Subpath(Line(start, end))

    assert False


def crop_bezier(seg, t0, t1):
    """returns a cropped copy of this segment which starts at self.point(t0)
    and ends at self.point(t1)."""
    if t0 < t1:
        swap = False
    else:
        t1, t0 = t0, t1
        swap = True

    if t0 == 0:
        cropped_seg = seg.split(t1)[0]

    elif t1 == 1:
        cropped_seg = seg.split(t0)[1]

    else:
        pt1 = seg.point(t1)

        # trim off the 0 <= t < t0 part
        trimmed_seg = seg.split(t0)[1]

        # find the adjusted t1 (i.e. the t1 such that
        # trimmed_seg.point(t1) ~= pt))and trim off the t1 < t <= 1 part

        # recent note: I suspect there's a more elegant way to do this
        # than to call radialrange...
        t1_adj = bezier_radialrange(trimmed_seg, pt1)[0][1]

        # let's try this:
        assert np.isclose(t1_adj, (t1 - t0) / (1 - t0))

        print("hey don't forget to remove bezier_radialrange thing (?)")

        t1_adj = 1 if t0 == 1 else (t1 - t0) / (1 - t0)
        cropped_seg = trimmed_seg.split(t1_adj)[0]

    assert isinstance(cropped_seg, BezierSegment)
    assert type(cropped_seg) == type(seg)

    return cropped_seg if not swap else cropped_seg.reversed()


# Main Classes #################################################################

# README

# The main user-manipulated classes are:

# Line
# QuadraticBezier
# CubicBezier
# Arc
# Subpath
# Path

# The four first four types of objects ('Line', 'QuadraticBezier', 'CubicBezier'
# and 'Arc') are commonly known as "segments".

# A subpath is a list of end-to-end contiguous segments (possibly closed).

# Finally a path is an ordered list of subpaths.

# There are also superclasses to help group and categorize the various types of
# segments. These are 'Segment', from which all segments inherit, and
# 'BezierSegment', from which only 'Line', 'QuadraticBezier' and 'CubicBezier'
# inherit. The class inheritance diagram for all segments is as follows:

# - Segment:
#   - BezierSegment:
#     - Line
#     - QuadraticBezier
#     - CubicBezier
#   - Arc

# Note that only 'Arc', 'Line', 'QuadraticBezier' and 'CubicBezier' are meant to
# be instantiated.

# The superclasses 'Segment' and 'BezierSegment' are useful for identifying types
# of segments via "isinstance(...)" as well as to reduce code size (as always with
# inheritance).


class Address(object):
    def __init__(self, **kw):
        assert all(key in ['t', 'T', 'W', 'segment_index', 'subpath_index'] for key in kw)
        self._t = kw.get('t')
        self._T = kw.get('T', None)
        self._W = kw.get('W', None)
        self._segment_index = kw.get('segment_index', None)
        self._subpath_index = kw.get('subpath_index', None)

    def __repr__(self):
        return 'Address t={} T={} W={} segment_index={} subpath_index={}'.format(self.t, self.T, self.W, self.segment_index, self.subpath_index)

    # this is somehow useful because otherwise float '1.0' vs int '1' causes two
    # same addresses to be considered unequal:
    def __eq__(self, other):
        if not isinstance(other, Address):
            return NotImplemented
        return \
            self._t == other._t and \
            self._T == other._T and \
            self._W == other._W and \
            self._segment_index == other._segment_index and \
            self._subpath_index == other._subpath_index

    def complete(self):
        return self._t is not None and \
            self._T is not None and \
            self._W is not None and \
            self._segment_index is not None and \
            self._subpath_index is not None

    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, val):
        if val is not None:
            if self._t is None:
                self._t = val
            else:
                if self.t != val:
                    raise ValueError("Attempt to overwrite Address.t")

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, val):
        if val is not None:
            if self._T is None:
                self._T = val
            else:
                if self._T != val:
                    raise ValueError("Attempt to overwrite Address.T")

    @property
    def W(self):
        return self._W

    @W.setter
    def W(self, val):
        if val is not None:
            if self._W is None:
                self._W = val
            else:
                if self._W != val:
                    raise ValueError("Attempt to overwrite Address.W")

    @property
    def segment_index(self):
        return self._segment_index

    @segment_index.setter
    def segment_index(self, val):
        if val is not None:
            if self._segment_index is None:
                self._segment_index = val
            else:
                if self._segment_index != val:
                    raise ValueError("Attempt to overwrite Address.segment_index")

    @property
    def subpath_index(self):
        return self._subpath_index

    @subpath_index.setter
    def subpath_index(self, val):
        if val is not None:
            if self._subpath_index is None:
                self._subpath_index = val
            else:
                if self._subpath_index != val:
                    raise ValueError("Attempt to overwrite Address.subpath_index")


def address_pair_from_t1t2(t1, t2):
    return (Address(t=t1), Address(t=t2))


def address_pair_from_t1t2_tuple(p):
    return (Address(t=p[0]), Address(t=p[1]))


class Segment(object):
    # shared implementation functions:

    def rotated(self, degs, origin=None):
        """Returns a copy of self rotated by `degs` degrees (CCW) around the
        point `origin` (a complex number).  By default `origin` is either
        `self.point(0.5)`, or in the case that self is an Arc object,
        `origin` defaults to `self.center`."""
        return rotate(self, degs, origin=origin)

    def translated(self, z0):
        """Returns a copy of self shifted by the complex quantity `z0` such
        that self.translated(z0).point(t) = self.point(t) + z0 for any t."""
        return translate(self, z0)

    def joins_smoothly_with(self, previous, wrt_parameterization=False):
        """ see docstring for is_smooth_join """
        return is_smooth_join(previous, self, wrt_parameterization)

    def multisplit(self, ts):
        """takes a possibly empty list of t-values t1, ..., tn such that
        0 < t1 < ... < tn < 1 and returns a list of Lines whose union is this
        segment and whose endpoints are self.point(0), self.point(t1), ...,
        self.point(tn)"""
        return multisplit_of(self, ts)

    def normal(self, t=None):
        """returns the (right hand rule) unit normal vector to self at t."""
        return - 1j * self.unit_tangent(t)

    def ilength(self, s, s_tol=ILENGTH_S_TOL, maxits=ILENGTH_MAXITS,
                error=ILENGTH_ERROR, min_depth=ILENGTH_MIN_DEPTH):
        """Returns a float, t, such that self.length(0, t) is approximately s.
        See the inv_arclength() docstring for more details."""
        return inv_arclength(self, s, s_tol=s_tol, maxits=maxits, error=error,
                             min_depth=min_depth)

    def curvature(self, t):
        """returns the curvature of the segment at t."""
        """(overwritten by Line, for which the value is 0)"""
        return segment_curvature(self, t)


class BezierSegment(Segment):
    # shared implementation functions:

    def offset(self, amount, two_sided=False, quality=None, safety=None):
        """Does not provide linecaps etc; use Subpath.offset for that"""
        return segment_offset_as_subpath(self, amount, two_sided, quality, safety)

    def radialrange(self, origin, return_all_global_extrema=False):
        """returns the tuples (d_min, t_min) and (d_max, t_max) which minimize
        and maximize, respectively, the distance d = |self.point(t) - origin|."""
        return bezier_radialrange(self, origin, return_all_global_extrema=return_all_global_extrema)

    def bbox(self):
        """(overwritten by Line)"""
        """returns the bounding box for the segment in the form
        (xmin, xmax, ymin, ymax)."""
        return bezier_bounding_box(self)

    def unit_tangent(self, t):
        """(overwritten by Line)"""
        """returns the unit tangent vector of the segment at t (centered at
        the origin and expressed as a complex number).  If the tangent
        vector's magnitude is zero, this method will find the limit of
        self.derivative(tau) / abs(self.derivative(tau)) as tau approaches t."""
        return bezier_unit_tangent(self, t)

    def __eq__(self, other):
        if type(self) != type(other):
            return NotImplemented
        return all(x == y for x, y in zip(self.bpoints(), other.bpoints()))

    def __ne__(self, other):
        if type(self) != type(other):
            return NotImplemented
        return not self == other

    def __getitem__(self, item):
        return self.bpoints()[item]

    def __len__(self):
        return len(self.bpoints())

    def intersect(self, other_seg, tol=1e-12):
        """Finds the intersections of two segments. Returns a list of tuples
        (I1, I2) of Intersection objects, such that self.point(I1.t) ==
        other_seg.point(I2.t).

        Note: This will fail if the two segments coincide for more than a
        finite collection of points."""
        assert self != other_seg

        if isinstance(other_seg, Line):
            t1t2s = bezier_by_line_intersections(self, other_seg)
            return [address_pair_from_t1t2(t1, t2) for t1, t2 in t1t2s]

        elif isinstance(other_seg, BezierSegment):
            t1t2s = bezier_intersections(self, other_seg,
                                         longer_length=max(self.length(), other_seg.length()),
                                         tol=tol, tol_deC=tol)
            return [address_pair_from_t1t2(t1, t2) for t1, t2 in t1t2s]

        elif isinstance(other_seg, Arc):
            return [(I1, I2) for I2, I1 in other_seg.intersect(self)]

        elif isinstance(other_seg, Path):
            raise TypeError(
                "other_seg must be a path segment, not a Path object, use "
                "Path.intersect() or Subpath.intersect().")

        else:
            raise TypeError("other_seg must be a path segment.")

    def split(self, t):
        """returns two segments of same type whose union is this segment and
        which join at self.point(t). (Overwritten by Line.)"""

        bpoints1, bpoints2 = split_bezier(self.bpoints(), t)

        if isinstance(self, QuadraticBezier):
            return QuadraticBezier(*bpoints1), QuadraticBezier(*bpoints2)

        elif isinstance(self, CubicBezier):
            return CubicBezier(*bpoints1), CubicBezier(*bpoints2)

        raise BugException

    def cropped(self, t0, t1):
        """returns a cropped copy of the segment which starts at self.point(t0)
        and ends at self.point(t1). Allows t1 >= t0. (Overwritten by Line.)"""
        return crop_bezier(self, t0, t1)


class Line(BezierSegment):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __repr__(self):
        return 'Line(start=%s, end=%s)' % (self.start, self.end)

    def tweaked(self, start=None, end=None):
        start = start if start is not None else self.start
        end = end if end is not None else self.end
        return Line(start, end)

    def point(self, t):
        """returns the coordinates of the Bezier curve evaluated at t."""
        if t == 1:
            # (this special case mainly here because otherwise we would
            # not have self.end == self.point(1), due to rounding issues)
            return self.end
        return self.start + t * (self.end - self.start)

    def length(self, t0=0, t1=1, error=None, min_depth=None):
        """returns the length of the line segment between t0 and t1."""
        return abs(self.end - self.start) * (t1 - t0)

    def bpoints(self):
        """returns the Bezier control points of the segment."""
        return self.start, self.end

    def poly(self, return_coeffs=False):
        """returns the line as a Polynomial object."""
        p = self.bpoints()
        coeffs = ([p[1] - p[0], p[0]])
        if return_coeffs:
            return coeffs
        else:
            return np.poly1d(coeffs)

    def derivative(self, t=None, n=1):
        """returns the nth derivative of the segment at t."""
        assert self.end != self.start
        if n == 1:
            return self.end - self.start
        elif n > 1:
            return 0
        else:
            raise ValueError("n should be a positive integer.")

    def unit_tangent(self, t=None):
        """returns the unit tangent of the segment at t."""
        assert self.end != self.start
        dseg = self.end - self.start
        return dseg / abs(dseg)

    def curvature(self, t):
        """returns the curvature of the line, which is always zero."""
        return 0

    def reversed(self):
        """returns a copy of the Line object with its orientation reversed"""
        return Line(self.end, self.start)

    def bbox(self):
        """returns the bounding box for the segment in the form
        (xmin, xmax, ymin, ymax)."""
        xmin = min(self.start.real, self.end.real)
        xmax = max(self.start.real, self.end.real)
        ymin = min(self.start.imag, self.end.imag)
        ymax = max(self.start.imag, self.end.imag)
        return xmin, xmax, ymin, ymax

    def cropped(self, t0, t1):
        """returns a cropped copy of this segment which starts at self.point(t0)
        and ends at self.point(t1). Allows t1 >= t0."""
        return Line(self.point(t0), self.point(t1))

    def split(self, t):
        """returns two Lines, whose union is this segment and which join at
        self.point(t)."""
        pt = self.point(t)
        return Line(self.start, pt), Line(pt, self.end)

    def naive_offset(self, amount):
        """performs a one-sided offset of the line by amount"""
        n = self.normal()
        return Line(self.start + amount * n, self.end + amount * n)


class QuadraticBezier(BezierSegment):
    # (jpsteinb: I don't understand this pickle stuff; in meanwhile, I leave
    # this here (and ditto for CubicBezier):)
    # For compatibility with old pickle files.
    _length_info = {'length': None, 'bpoints': None}

    def __init__(self, start, control, end):
        self.start = start
        self.end = end
        self.control = control
        # used to know if self._length needs to be updated:
        self._length_info = {'length': None, 'bpoints': None}

    def __repr__(self):
        return 'QuadraticBezier(start=%s, control=%s, end=%s)' % \
               (self.start, self.control, self.end)

    def tweaked(self, start=None, control=None, end=None):
        start = start if start is not None else self.start
        control = control if control is not None else self.control
        end = end if end is not None else self.end
        return QuadraticBezier(start, control, end)

    def is_smooth_from(self, previous, warning_on=True):
        """[Warning: The name of this method is somewhat misleading (yet kept
        for compatibility with scripts created using svg.path 2.0).  This
        method is meant only for d string creation and should not be used to
        check for kinks.  To check a segment for differentiability, use the
        joins_smoothly_with() method instead.]"""
        if warning_on:
            warn(_is_smooth_from_warning)
        if isinstance(previous, QuadraticBezier):
            return (self.start == previous.end and
                    (self.control - self.start) == (
                        previous.end - previous.control))
        else:
            return self.control == self.start

    def point(self, t):
        """returns the coordinates of the Bezier curve evaluated at t."""
        return (1 - t)**2 * self.start + 2 * (1 - t) * t * self.control + t**2 * self.end

    def length(self, t0=0, t1=1, error=None, min_depth=None):
        if t0 == 1 and t1 == 0:
            if self._length_info['bpoints'] == self.bpoints():
                return self._length_info['length']
        a = self.start - 2 * self.control + self.end
        b = 2 * (self.control - self.start)
        a_dot_b = a.real * b.real + a.imag * b.imag

        if abs(a) < 1e-12:
            s = abs(b) * (t1 - t0)
        elif abs(a_dot_b + abs(a) * abs(b)) < 1e-12:
            tstar = abs(b) / (2 * abs(a))
            if t1 < tstar:
                return abs(a) * (t0**2 - t1**2) - abs(b) * (t0 - t1)
            elif tstar < t0:
                return abs(a) * (t1**2 - t0**2) - abs(b) * (t1 - t0)
            else:
                return abs(a) * (t1**2 + t0**2) - abs(b) * (t1 + t0) + \
                    abs(b)**2 / (2 * abs(a))
        else:
            c2 = 4 * (a.real**2 + a.imag**2)
            c1 = 4 * a_dot_b
            c0 = b.real**2 + b.imag**2

            beta = c1 / (2 * c2)
            gamma = c0 / c2 - beta**2

            dq1_mag = sqrt(c2 * t1**2 + c1 * t1 + c0)
            dq0_mag = sqrt(c2 * t0**2 + c1 * t0 + c0)
            logarand = (sqrt(c2) * (t1 + beta) + dq1_mag) / \
                       (sqrt(c2) * (t0 + beta) + dq0_mag)

            s = (t1 + beta) * dq1_mag - (t0 + beta) * dq0_mag + \
                gamma * sqrt(c2) * log(logarand)
            s /= 2

        if t0 == 1 and t1 == 0:
            self._length_info['length'] = s
            self._length_info['bpoints'] = self.bpoints()
            return self._length_info['length']
        else:
            return s

    def bpoints(self):
        """returns the Bezier control points of the segment."""
        return self.start, self.control, self.end

    def poly(self, return_coeffs=False):
        """returns the quadratic as a Polynomial object."""
        p = self.bpoints()
        coeffs = (p[0] - 2 * p[1] + p[2], 2 * (p[1] - p[0]), p[0])
        if return_coeffs:
            return coeffs
        else:
            return np.poly1d(coeffs)

    def derivative(self, t, n=1):
        """returns the nth derivative of the segment at t.
        Note: Bezier curves can have points where their derivative vanishes.
        If you are interested in the tangent direction, use the unit_tangent()
        method instead."""
        p = self.bpoints()
        if n == 1:
            return 2 * ((p[1] - p[0]) * (1 - t) + (p[2] - p[1]) * t)
        elif n == 2:
            return 2 * (p[2] - 2 * p[1] + p[0])
        elif n > 2:
            return 0
        else:
            raise ValueError("n should be a positive integer.")

    def reversed(self):
        """returns a copy of the QuadraticBezier object with its orientation
        reversed."""
        new_quad = QuadraticBezier(self.end, self.control, self.start)
        if self._length_info['length']:
            new_quad._length_info = self._length_info
            new_quad._length_info['bpoints'] = (
                self.end, self.control, self.start)
        return new_quad

    def naive_offset(self, amount):
        start   = self.start + self.normal(0) * amount
        end     = self.end   + self.normal(1) * amount
        control = offset_intersection(self.start, self.control, self.end, amount)
        return QuadraticBezier(start, control, end)


class CubicBezier(BezierSegment):
    # For compatibility with old pickle files.
    _length_info = {'length': None, 'bpoints': None, 'error': None,
                    'min_depth': None}

    def __init__(self, start, control1, control2, end):
        self.start = start
        self.control1 = control1
        self.control2 = control2
        self.end = end

        # used to know if self._length needs to be updated
        self._length_info = {'length': None, 'bpoints': None, 'error': None,
                             'min_depth': None}

    def __repr__(self):
        return 'CubicBezier(start=%s, control1=%s, control2=%s, end=%s)' % \
               (self.start, self.control1, self.control2, self.end)

    def tweaked(self, start=None, control1=None, control2=None, end=None):
        start = start if start is not None else self.start
        control1 = control1 if control1 is not None else self.control1
        control2 = control2 if control2 is not None else self.control2
        end = end if end is not None else self.end
        return CubicBezier(start, control1, control2, end)

    def is_smooth_from(self, previous, warning_on=True):
        """[Warning: The name of this method is somewhat misleading (yet kept
        for compatibility with scripts created using svg.path 2.0).  This
        method is meant only for d string creation and should not be used to
        check for kinks.  To check a segment for differentiability, use the
        joins_smoothly_with() method instead.]"""
        if warning_on:
            warn(_is_smooth_from_warning)
        if isinstance(previous, CubicBezier):
            return (self.start == previous.end and
                    (self.control1 - self.start) == (
                        previous.end - previous.control2))
        else:
            return self.control1 == self.start

    def point(self, t):
        """Evaluate the cubic Bezier curve at t using Horner's rule."""
        # algebraically equivalent to
        # P0 * (1 - t)**3 + 3 * P1 * t * (1 - t)**2 + 3 * P2 * (1 - t) * t**2 + P3 * t**3
        # for (P0, P1, P2, P3) = self.bpoints()
        return self.start + t * (
            3 * (self.control1 - self.start) + t * (
                3 * (self.start + self.control2) - 6 * self.control1 + t * (
                    - self.start + 3 * (self.control1 - self.control2) + self.end)))

    def length(self, t0=0, t1=1, error=LENGTH_ERROR, min_depth=LENGTH_MIN_DEPTH):
        """Calculate the length of the path up to a certain position"""
        if t0 == 0 and t1 == 1:
            if self._length_info['bpoints'] == self.bpoints() \
                    and self._length_info['error'] <= error \
                    and self._length_info['min_depth'] >= min_depth:
                return self._length_info['length']

        # using scipy.integrate.quad is quick
        if _quad_available:
            s = quad(lambda tau: abs(self.derivative(tau)), t0, t1, epsabs=error, limit=1000)[0]
        else:
            s = segment_length(self, t0, t1, self.point(t0), self.point(t1),
                               error, min_depth, 0)

        if t0 == 0 and t1 == 1:
            self._length_info['length'] = s
            self._length_info['bpoints'] = self.bpoints()
            self._length_info['error'] = error
            self._length_info['min_depth'] = min_depth
            return self._length_info['length']
        else:
            return s

    def bpoints(self):
        """returns the Bezier control points of the segment."""
        return self.start, self.control1, self.control2, self.end

    def poly(self, return_coeffs=False):
        """Returns a the cubic as a Polynomial object."""
        p = self.bpoints()
        coeffs = (- p[0] + 3 * (p[1] - p[2]) + p[3],
                  3 * (p[0] - 2 * p[1] + p[2]),
                  3 * (- p[0] + p[1]),
                  p[0])
        if return_coeffs:
            return coeffs
        else:
            return np.poly1d(coeffs)

    def derivative(self, t, n=1):
        """returns the nth derivative of the segment at t.
        Note: Bezier curves can have points where their derivative vanishes.
        If you are interested in the tangent direction, use the unit_tangent()
        method instead."""
        p = self.bpoints()

        if n == 1:
            return \
                3 * (p[1] - p[0]) * (1 - t)**2 + \
                6 * (p[2] - p[1]) * (1 - t) * t + \
                3 * (p[3] - p[2]) * t**2

        elif n == 2:
            return \
                6 * (1 - t) * (p[2] - 2 * p[1] + p[0]) + \
                6 * t * (p[3] - 2 * p[2] + p[1])

        elif n == 3:
            return 6 * (p[3] - 3 * (p[2] - p[1]) - p[0])

        elif n > 3:
            return 0

        else:
            raise ValueError("n should be a positive integer.")

    def reversed(self):
        """returns a copy of the CubicBezier object with its orientation
        reversed."""
        new_cub = CubicBezier(self.end, self.control2, self.control1,
                              self.start)
        if self._length_info['length']:
            new_cub._length_info = self._length_info
            new_cub._length_info['bpoints'] = (
                self.end, self.control2, self.control1, self.start)
        return new_cub

    def naive_offset(self, amount):
        start    = self.start + self.normal(0) * amount
        end      = self.end   + self.normal(1) * amount
        control1 = offset_intersection(self.start, self.control1, self.control2, amount)
        control2 = offset_intersection(self.control1, self.control2, self.end, amount)
        return CubicBezier(start, control1, control2, end)


class Arc(BezierSegment):
    def __init__(self, start, radius, rotation, large_arc, sweep, end,
                 autoscale_radius=True):
        """
        This should be thought of as a part of an ellipse connecting two
        points on that ellipse, start and end.
        Parameters
        - - - - - - - - - -
        start : complex
            The start point of the curve. Note: `start` and `end` cannot be the
            same.  To make a full ellipse or circle, use two `Arc` objects.
        radius : complex
            rx + 1j * ry, where rx and ry are the radii of the ellipse (also
            known as its semi - major and semi - minor axes, or vice - versa or if
            rx < ry).
            Note: If rx = 0 or ry = 0 then this arc is treated as a
            straight line segment joining the endpoints.
            Note: If rx or ry has a negative sign, the sign is dropped; the
            absolute value is used instead.
            Note:  If no such ellipse exists, the radius will be scaled so
            that one does (unless autoscale_radius is set to False).
        rotation : float
            This is the CCW angle (in degrees) from the positive x - axis of the
            current coordinate system to the x - axis of the ellipse.
        large_arc : bool
            Given two points on an ellipse, there are two elliptical arcs
            connecting those points, the first going the short way around the
            ellipse, and the second going the long way around the ellipse.  If
            `large_arc == False`, the shorter elliptical arc will be used.  If
            `large_arc == True`, then longer elliptical will be used.
            In other words, `large_arc` should be 0 for arcs spanning less than
            or equal to 180 degrees and 1 for arcs spanning greater than 180
            degrees.
        sweep : bool
            For any acceptable parameters `start`, `end`, `rotation`, and
            `radius`, there are two ellipses with the given major and minor
            axes (radii) which connect `start` and `end`.  One which connects
            them in a CCW fashion and one which connected them in a CW
            fashion.  If `sweep == True`, the CCW ellipse will be used.  If
            `sweep == False`, the CW ellipse will be used.  See note on curve
            orientation below.
        end : complex
            The end point of the curve. Note: `start` and `end` cannot be the
            same.  To make a full ellipse or circle, use two `Arc` objects.
        autoscale_radius : bool
            If `autoscale_radius == True`, then will also scale `self.radius`
            in the case that no ellipse exists with the input parameters
            (see inline comments for further explanation).

        Derived Parameters / Attributes
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.theta : float
            This is the phase (in degrees) of self.u1transform(self.start).
            It is $\theta_1$ in the official documentation and ranges from
            - 180 to 180.
        self.delta : float
            This is the angular distance (in degrees) between the start and
            end of the arc after the arc has been sent to the unit circle
            through self.u1transform().
            It is $\Delta \theta$ in the official documentation and ranges from
            - 360 to 360; being positive when the arc travels CCW and negative
            otherwise (i.e. is positive / negative when sweep == True / False).
        self.center : complex
            This is the center of the arc's ellipse.
        self.phi : float
            The arc's rotation in radians, i.e. `radians(self.rotation)`.
        self.rot_matrix : complex
            Equal to `exp(1j * self.phi)` which is also equal to
            `cos(self.phi) + 1j * sin(self.phi)`.


        Note on curve orientation (CW vs CCW)
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        The notions of clockwise (CW) and counter - clockwise (CCW) are reversed
        in some sense when viewing SVGs (as the y coordinate starts at the top
        of the image and increases towards the bottom).
        """
        assert start != end
        assert radius.real != 0 and radius.imag != 0

        self.start = start
        self.radius = abs(radius.real) + 1j * abs(radius.imag)
        self.rotation = rotation
        self.large_arc = bool(large_arc)
        self.sweep = bool(sweep)
        self.end = end
        self.autoscale_radius = autoscale_radius

        # Convenience parameters
        self.phi = radians(self.rotation)
        self.rot_matrix = exp(1j * self.phi)

        # Derive derived parameters
        self._parameterize()

    def __repr__(self):
        params = (self.start, self.radius, self.rotation,
                  self.large_arc, self.sweep, self.end)
        return ("Arc({}start={}, radius={}, rotation={}, "
                "large_arc={}, sweep={}, end={}{})".format(*params))

    def tweaked(self, start=None, radius=None, rotation=None, large_arc=None, sweep=None, end=None, autoscale_radius=None):
        start = start if start is not None else self.start
        radius = radius if radius is not None else self.radius
        rotation = rotation if rotation is not None else self.rotation
        large_arc = large_arc if large_arc is not None else self.large_arc
        sweep = sweep if sweep is not None else self.sweep
        end = end if end is not None else self.end
        return Arc(start, radius, rotation, large_arc, sweep, end, autoscale_radius)

    def __eq__(self, other):
        if not isinstance(other, Arc):
            return NotImplemented
        return self.start == other.start and self.end == other.end \
            and self.radius == other.radius \
            and self.rotation == other.rotation \
            and self.large_arc == other.large_arc and self.sweep == other.sweep

    def __ne__(self, other):
        if not isinstance(other, Arc):
            return NotImplemented
        return not self == other

    def _parameterize(self):
        # See http://www.w3.org/TR/SVG/implnote.html#ArcImplementationNotes
        # my notation roughly follows theirs
        rx = self.radius.real
        ry = self.radius.imag
        rx_sqd = rx * rx
        ry_sqd = ry * ry

        # Transform z - > z' = x' + 1j * y'
        # = self.rot_matrix**(-1) * (z - (end + start) / 2)
        # coordinates.  This translates the ellipse so that the midpoint
        # between self.end and self.start lies on the origin and rotates
        # the ellipse so that the its axes align with the xy - coordinate axes.
        # Note:  This sends self.end to - self.start
        zp1 = (1 / self.rot_matrix) * (self.start - self.end) / 2
        x1p, y1p = zp1.real, zp1.imag
        x1p_sqd = x1p * x1p
        y1p_sqd = y1p * y1p

        # Correct out of range radii
        # Note: an ellipse going through start and end with radius and phi
        # exists if and only if radius_check is true
        radius_check = (x1p_sqd / rx_sqd) + (y1p_sqd / ry_sqd)
        if radius_check > 1:
            if self.autoscale_radius:
                rx *= sqrt(radius_check)
                ry *= sqrt(radius_check)
                self.radius = rx + 1j * ry
                rx_sqd = rx * rx
                ry_sqd = ry * ry
            else:
                raise ValueError("No such elliptic arc exists.")

        # Compute c'=(c_x', c_y'), the center of the ellipse in (x', y') coords
        # Noting that, in our new coord system, (x_2', y_2') = (- x_1', -x_2')
        # and our ellipse is cut out by of the plane by the algebraic equation
        # (x'-c_x')**2 / r_x**2 + (y'-c_y')**2 / r_y**2 = 1,
        # we can find c' by solving the system of two quadratics given by
        # plugging our transformed endpoints (x_1', y_1') and (x_2', y_2')
        tmp = rx_sqd * y1p_sqd + ry_sqd * x1p_sqd
        radicand = (rx_sqd * ry_sqd - tmp) / tmp
        try:
            radical = sqrt(radicand)
        except ValueError:
            radical = 0
        if self.large_arc == self.sweep:
            cp = - radical * (rx * y1p / ry - 1j * ry * x1p / rx)
        else:
            cp = radical * (rx * y1p / ry - 1j * ry * x1p / rx)

        # The center in (x, y) coordinates is easy to find knowing c'
        self.center = exp(1j * self.phi) * cp + (self.start + self.end) / 2

        # Now we do a second transformation, from (x', y') to (u_x, u_y)
        # coordinates, which is a translation moving the center of the
        # ellipse to the origin and a dilation stretching the ellipse to be
        # the unit circle
        u1 = (x1p - cp.real) / rx + 1j * (y1p - cp.imag) / ry  # transformed start
        u2 = (- x1p - cp.real) / rx + 1j * (- y1p - cp.imag) / ry  # transformed end

        # Now compute theta and delta (we'll define them as we go)
        # delta is the angular distance of the arc (w.r.t the circle)
        # theta is the angle between the positive x'-axis and the start point
        # on the circle
        if u1.imag > 0:
            self.theta = degrees(acos(u1.real))
        elif u1.imag < 0:
            self.theta = - degrees(acos(u1.real))
        else:
            if u1.real > 0:  # start is on pos u_x axis
                self.theta = 0
            else:  # start is on neg u_x axis
                # Note: This behavior disagrees with behavior documented in
                # http://www.w3.org / TR / SVG / implnote.html  #ArcImplementationNotes
                # where theta is set to 0 in this case.
                self.theta = 180

        det_uv = u1.real * u2.imag - u1.imag * u2.real

        acosand = u1.real * u2.real + u1.imag * u2.imag
        if acosand > 1 or acosand < -1:
            acosand = round(acosand)
        if det_uv > 0:
            self.delta = degrees(acos(acosand))
        elif det_uv < 0:
            self.delta = - degrees(acos(acosand))
        else:
            if u1.real * u2.real + u1.imag * u2.imag > 0:
                # u1 == u2
                self.delta = 0
            else:
                # u1 == - u2
                # Note: This behavior disagrees with behavior documented in
                # http://www.w3.org/TR/SVG/implnote.html  #ArcImplementationNotes
                # where delta is set to 0 in this case.
                self.delta = 180

        if not self.sweep and self.delta >= 0:
            self.delta -= 360
        elif self.large_arc and self.delta <= 0:
            self.delta += 360

    def point(self, t):
        if t == 0:
            return self.start
        if t == 1:
            return self.end
        angle = radians(self.theta + t * self.delta)
        cosphi = self.rot_matrix.real
        sinphi = self.rot_matrix.imag
        rx = self.radius.real
        ry = self.radius.imag

        # z = self.rot_matrix * (rx * cos(angle) + 1j * ry * sin(angle)) + self.center
        x = rx * cosphi * cos(angle) - ry * sinphi * sin(angle) + self.center.real
        y = rx * sinphi * cos(angle) + ry * cosphi * sin(angle) + self.center.imag
        return complex(x, y)

    def centeriso(self, z):
        """This is an isometry that translates and rotates self so that it
        is centered on the origin and has its axes aligned with the xy axes."""
        return (1 / self.rot_matrix) * (z - self.center)

    def icenteriso(self, zeta):
        """This is an isometry, the inverse of standardiso()."""
        return self.rot_matrix * zeta + self.center

    def u1transform(self, z):
        """This is an affine transformation (same as used in
        self._parameterize()) that sends self to the unit circle."""
        zeta = (1 / self.rot_matrix) * (z - self.center)  # same as centeriso(z)
        x, y = real(zeta), imag(zeta)
        return x / self.radius.real + 1j * y / self.radius.imag

    def iu1transform(self, zeta):
        """This is an affine transformation, the inverse of
        self.u1transform()."""
        x = real(zeta)
        y = imag(zeta)
        z = x * self.radius.real + y * self.radius.imag
        return self.rot_matrix * z + self.center

    def length(self, t0=0, t1=1, error=LENGTH_ERROR, min_depth=LENGTH_MIN_DEPTH):
        """The length of an elliptical large_arc segment requires numerical
        integration, and in that case it's simpler to just do a geometric
        approximation, as for cubic bezier curves."""
        assert 0 <= t0 <= 1 and 0 <= t1 <= 1
        if _quad_available:
            return quad(lambda tau: abs(self.derivative(tau)), t0, t1,
                        epsabs=error, limit=1000)[0]
        else:
            return segment_length(self, t0, t1, self.point(t0), self.point(t1),
                                  error, min_depth, 0)

    def derivative(self, t, n=1):
        """returns the nth derivative of the segment at t."""
        angle = radians(self.theta + t * self.delta)
        phi = radians(self.rotation)
        rx = self.radius.real
        ry = self.radius.imag
        k = (self.delta * 2 * pi / 360)**n  # ((d / dt)angle)**n

        if n % 4 == 0 and n > 0:
            return rx * cos(phi) * cos(angle) - ry * sin(phi) * sin(angle) + 1j * (
                rx * sin(phi) * cos(angle) + ry * cos(phi) * sin(angle))
        elif n % 4 == 1:
            return k * (- rx * cos(phi) * sin(angle) - ry * sin(phi) * cos(angle) + 1j * (
                - rx * sin(phi) * sin(angle) + ry * cos(phi) * cos(angle)))
        elif n % 4 == 2:
            return k * (- rx * cos(phi) * cos(angle) + ry * sin(phi) * sin(angle) + 1j * (
                - rx * sin(phi) * cos(angle) - ry * cos(phi) * sin(angle)))
        elif n % 4 == 3:
            return k * (rx * cos(phi) * sin(angle) + ry * sin(phi) * cos(angle) + 1j * (
                rx * sin(phi) * sin(angle) - ry * cos(phi) * cos(angle)))
        else:
            raise ValueError("n should be a positive integer.")

    def unit_tangent(self, t):
        """returns the unit tangent vector of the segment at t (centered at
        the origin and expressed as a complex number)."""
        dseg = self.derivative(t)
        return dseg / abs(dseg)

    def reversed(self):
        """returns a copy of the Arc object with its orientation reversed."""
        return Arc(self.end, self.radius, self.rotation, self.large_arc,
                   not self.sweep, self.start)

    def phase2t(self, psi):
        """Given phase -pi < psi <= pi,
        returns the t value such that
        exp(1j * psi) = self.u1transform(self.point(t)).
        """
        def _deg(rads, domain_lower_limit):
            # Convert rads to degrees in [0, 360) domain
            degs = degrees(rads % (2 * pi))

            # Convert to [domain_lower_limit, domain_lower_limit + 360) domain
            k = domain_lower_limit // 360
            degs += k * 360
            if degs < domain_lower_limit:
                degs += 360
            return degs

        if self.delta > 0:
            degs = _deg(psi, domain_lower_limit=self.theta)
        else:
            degs = _deg(psi, domain_lower_limit=self.theta)
        return (degs - self.theta) / self.delta

    def t2lambda(self, t):
        p = self.point(t) - self.center
        q = (1 / self.rot_matrix) * p
        return atan2(imag(q), real(q))

    def lambda2eta(self, lamda):
        return atan2(sin(lamda) / self.radius.imag, cos(lamda) / self.radius.real)

    def Maisonobe_E(self, eta):
        """(this function is unused; see next function)"""
        self.center + self.rot_matrix * (self.rx * cos(eta) + 1j * self.ry * sin(eta))

    def Maisonobe_E_prime(self, eta):
        """see paper 'Drawing an elliptical arc using polylines, quadratic
        or cubic Bezier curves' by L. Maisonobe, 2003, sections 2.2.1 and 3.4.1"""
        return self.rot_matrix * (-self.radius.real * sin(eta) + 1j * self.radius.imag * cos(eta))

    def Maisonobe_cubic_interpolation(self, t1, t2):
        """see paper 'Drawing an elliptical arc using polylines, quadratic
        or cubic Bezier curves' by L. Maisonobe, 2003, sections 2.2.1 and 3.4.1"""
        assert 0 <= t1 < t2 <= 1
        start = self.point(t1)
        end   = self.point(t2)
        eta1  = self.lambda2eta(self.t2lambda(t1))
        eta2  = self.lambda2eta(self.t2lambda(t2))
        discr = 4 + 3 * tan((eta2 - eta1) * 0.5)**2
        alpha = sin(eta2 - eta1) * (sqrt(discr) - 1) / 3
        control1 = start + alpha * self.Maisonobe_E_prime(eta1)
        control2 = end   - alpha * self.Maisonobe_E_prime(eta2)
        return CubicBezier(start, control1, control2, end)

    def converted_to_bezier_subpath(self, quality=0.01, safety=5):
        assert quality > 0
        safety = int(min(4, safety))
        other = self.Maisonobe_cubic_interpolation(0, 1)
        assert other.start == self.start
        assert other.end == self.end
        divergence = divergence_of_offset(self, other, 0, safety=safety, early_return_threshold=quality)
        if divergence <= quality:
            return Subpath(other), [0, 1]
        first_half, secnd_half = self.split(0.5)
        assert first_half.end == secnd_half.start
        P1, ts1 = first_half.converted_to_bezier_subpath(quality, safety)
        P2, ts2 = secnd_half.converted_to_bezier_subpath(quality, safety)
        assert P1.end == P2.start
        allts = [0.5 * t for t in ts1] + [0.5 + 0.5 * t for t in ts2[1:]]
        return P1.extend(P2), allts

    def offset(self, amount, two_sided=False, quality=0.01, safety=5):
        cpath, _ = self.converted_to_bezier_subpath(quality * 0.5, safety)
        assert isinstance(cpath, Subpath)
        return cpath.offset(amount, two_sided, quality * 0.5, safety)

    def intersect(self, other_seg, tol=1e-12):
        """
        NOT FULLY IMPLEMENTED. Supposed to return a list of tuples (t1, t2) such
        that self.point(t1) == other_seg.point(t2).

        Note: This will fail if the two segments coincide for more than a
        finite collection of points.

        Note: Arc related intersections are only partially supported, i.e. are
        only half-heartedly implemented and not well tested. Please feel free
        to let me know if you're interested in such a feature -- or even better
        please submit an implementation if you want to code one."""
        if isinstance(other_seg, BezierSegment):
            u1poly = self.u1transform(other_seg.poly())
            u1poly_mag2 = real(u1poly)**2 + imag(u1poly)**2
            t2s = polyroots01(u1poly_mag2 - 1)
            t1s = [self.phase2t(phase(u1poly(t2))) for t2 in t2s]
            return [address_pair_from_t1t2(t1, t2) for t1, t2 in zip(t1s, t2s)]

        elif isinstance(other_seg, Arc):
            assert other_seg != self
            # This could be made explicit to increase efficiency
            longer_length = max(self.length(), other_seg.length())
            inters = bezier_intersections(self, other_seg,
                                          longer_length=longer_length,
                                          tol=tol, tol_deC=tol)

            # ad hoc fix for redundant solutions
            if len(inters) > 2:
                def keyfcn(tpair):
                    t1, t2 = tpair
                    return abs(self.point(t1) - other_seg.point(t2))
                inters.sort(key=keyfcn)
                for idx in range(1, len(inters) - 1):
                    if (abs(inters[idx][0] - inters[idx + 1][0]) < abs(inters[idx][0] - inters[0][0])):
                        return [address_pair_from_t1t2_tuple(inters[0]),
                                address_pair_from_t1t2_tuple(inters[idx])]
                else:
                    return [address_pair_from_t1t2_tuple(inters[0]),
                            address_pair_from_t1t2_tuple(inters[-1])]
            return [address_pair_from_t1t2_tuple(inters[0]),
                    address_pair_from_t1t2_tuple(inters[1])]

        else:
            raise TypeError("other_seg should be a Arc, Line, "
                            "QuadraticBezier, or CubicBezier object.")

    def bbox(self):
        """returns a bounding box for the segment in the form
        (xmin, xmax, ymin, ymax)."""
        # a(t) = radians(self.theta + self.delta * t)
        #      = (2 * pi / 360) * (self.theta + self.delta * t)
        # x'=0: ~~~~~~~~~
        # - rx * cos(phi) * sin(a(t)) = ry * sin(phi) * cos(a(t))
        # - (rx / ry) * cot(phi) * tan(a(t)) = 1
        # a(t) = arctan(- (ry / rx)tan(phi)) + pi * k === atan_x
        # y'=0: ~~~~~~~~~~
        # rx * sin(phi) * sin(a(t)) = ry * cos(phi) * cos(a(t))
        # (rx / ry) * tan(phi) * tan(a(t)) = 1
        # a(t) = arctan((ry / rx) * cot(phi))
        # atanres = arctan((ry / rx) * cot(phi)) === atan_y
        # ~~~~~~~~
        # (2 * pi / 360) * (self.theta + self.delta * t) = atanres + pi * k
        # Therfore, for both x' and y', we have...
        # t = ((atan_{x / y} + pi * k) * (360 / (2 * pi)) - self.theta) / self.delta
        # for all k s.t. 0 < t < 1
        from math import atan, tan

        if cos(self.phi) == 0:
            atan_x = pi / 2
            atan_y = 0
        elif sin(self.phi) == 0:
            atan_x = 0
            atan_y = pi / 2
        else:
            rx, ry = self.radius.real, self.radius.imag
            atan_x = atan(- (ry / rx) * tan(self.phi))
            atan_y = atan((ry / rx) / tan(self.phi))

        def angle_inv(ang, k):  # inverse of angle from Arc.derivative()
            return ((ang + pi * k) * (360 / (2 * pi)) - self.theta) / self.delta

        xtrema = [self.start.real, self.end.real]
        ytrema = [self.start.imag, self.end.imag]

        for k in range(- 4, 5):
            tx = angle_inv(atan_x, k)
            ty = angle_inv(atan_y, k)
            if 0 <= tx <= 1:
                xtrema.append(self.point(tx).real)
            if 0 <= ty <= 1:
                ytrema.append(self.point(ty).imag)
        return min(xtrema), max(xtrema), min(ytrema), max(ytrema)

    def split(self, t):
        """returns two segments, whose union is this segment and which join
        at self.point(t)."""
        return self.cropped(0, t), self.cropped(t, 1)

    def cropped(self, t0, t1):
        """returns a cropped copy of this segment which starts at self.point(t0)
        and ends at self.point(t1). Allows t1 > t0 but not t0 == t1."""
        if abs(self.delta * (t1 - t0)) <= 180:
            new_large_arc = 0
        else:
            new_large_arc = 1
        return Arc(self.point(t0), radius=self.radius, rotation=self.rotation,
                   large_arc=new_large_arc, sweep=self.sweep,
                   end=self.point(t1), autoscale_radius=self.autoscale_radius)

    def radialrange(self, origin, return_all_global_extrema=False):
        """returns the tuples (d_min, t_min) and (d_max, t_max) which minimize
        and maximize, respectively, the distance,
        d = |self.point(t) - origin|."""

        u1orig = self.u1transform(origin)

        # Transform to a coordinate system where the ellipse is centered
        # at the origin and its axes are horizontal / vertical
        zeta0 = self.centeriso(origin)
        a, b = self.radius.real, self.radius.imag
        x0, y0 = zeta0.real, zeta0.imag

        # Find t s.t. z'(t)
        a2mb2 = (a**2 - b**2)
        if u1orig.imag:  # x != x0

            coeffs = [a2mb2**2,
                      2 * a2mb2 * b**2 * y0,
                      (- a**4 + (2 * a**2 - b**2 + y0**2) * b**2 + x0**2) * b**2,
                      - 2 * a2mb2 * b**4 * y0,
                      - b**6 * y0**2]
            ys = polyroots(coeffs, realroots=True,
                           condition=lambda r: - b <= r <= b)
            xs = (a * sqrt(1 - y**2 / b**2) for y in ys)

        else:  # This case is very similar, see notes and assume instead y0!=y
            b2ma2 = (b**2 - a**2)
            coeffs = [b2ma2**2,
                      2 * b2ma2 * a**2 * x0,
                      (- b**4 + (2 * b**2 - a**2 + x0**2) * a**2 + y0**2) * a**2,
                      - 2 * b2ma2 * a**4 * x0,
                      - a**6 * x0**2]
            xs = polyroots(coeffs, realroots=True,
                           condition=lambda r: - a <= r <= a)
            ys = (b * sqrt(1 - x**2 / a**2) for x in xs)

        raise _NotImplemented4ArcException


class Subpath():
    """
    A subpath is a sequence of contiguous path segments; may or may not be
    closed via a Z value; see functions hard_stitch, soft_stitch, soft_unstitch.
    Can mostly be used autonomously from Path (e.g., has its own d() function)."""

    def __init__(self, *segments):
        self._start = None
        self._end = None
        self._z = False
        self._segments = []
        for s in segments:
            self.append(s)
        self.basic_reset()

    def check_health(self):
        assert all(isinstance(thing, Segment) for thing in self)
        for s, t in zip(self, self[1:]):
            assert s.end == t.start
        if self._z:
            assert len(self) > 0
            assert self[0].start == self[-1].end
        if len(self) > 0:
            assert self[0].start == self._start
            assert self[-1].end == self._end
        else:
            assert self._z is False
            assert self._start is None
            assert self._end is None

    def basic_reset(self):
        self._length = None
        self._lengths = None
        if len(self) > 0:
            self._start = self._segments[0].start
            self._end = self._segments[-1].end
            if self._z and self._start != self._end:
                self._z = False
        else:
            self._start = None
            self._end = None
            self._z = False  # debatable (but clean)

    def __len__(self):
        return len(self._segments)

    def __getitem__(self, index):
        return self._segments[index]

    def mod_index(self, index, use_Z):
        return index if not use_Z or not self._z else index % len(self)

    def prev_segment(self, index, use_Z=False):
        assert 0 <= index <= len(self) - 1
        prev_index = self.mod_index(index - 1, use_Z)
        return self[prev_index] if prev_index >= 0 else None

    def next_segment(self, index, use_Z=False):
        assert 0 <= index <= len(self) - 1
        next_index = self.mod_index(index + 1, use_Z)
        return self[next_index] if next_index <= len(self) - 1 else None

    def __setitem__(self, index, value):
        assert isinstance(value, Segment)
        assert -len(self) <= index <= len(self) - 1
        if len(self) > 0:
            index = index % len(self)
            prev = self.prev_segment(index, use_Z=True)
            next = self.next_segment(index, use_Z=True)
        else:
            prev = next = None
        if prev and prev.end != value.start:
            raise ValueError("subpath segment insertion results in .start discontinuity")
        if next and next.start != value.end:
            raise ValueError("subpath segment insertion results in .end discontinuity")
        self._segments[index] = value
        self.basic_reset()
        return value

    def __delitem__(self, index):
        assert -len(self) <= index <= len(self) - 1
        if len(self) > 0:
            index = index % len(self)
            prev = self.prev_segment(index, use_Z=True)
            next = self.next_segment(index, use_Z=True)
        else:
            prev = next = None
        if prev and next and prev.end != next.start:
            raise ValueError("subpath segment deletion results in discontinuity")
        to_return = self._segments[index]
        del self._segments[index]
        self.basic_reset()
        return to_return

    def __iter__(self):
        return self._segments.__iter__()

    def append(self, value):
        self.splice(len(self), len(self), value)
        return value  # [sic]

    def extend(self, subpath):
        self.splice(len(self), len(self), subpath)
        return self  # [sic]

    def prepend(self, value):
        self.insert(0, value)

    def insert(self, index, value):
        self.splice(index, index, value)
        return value  # [sic]

    def splice(self, start_index, end_index, value):
        """replaces segments of indices start_index, ..., end_index - 1 with the
        (segments in) value, which may be a segment, a subpath, or a naively
        continuous path, as long as no discontinuity is induced; if end_index ==
        start_index, no segments are replaced, and insertion occurs right before
        the segment at start_index == end_index, so that the first segment of new
        inserted material has index start_index. Can be used with an empty value
        (of type Path or Subpath), in which case the effect is only to delete the
        segments start_index, ..., end_index-1."""

        # assertions / checking
        assert isinstance(value, Segment) or isinstance(value, Subpath)
        assert not isinstance(value, Path) or value.is_naively_continuous()
        assert 0 <= start_index <= end_index <= len(self)

        if len(self) > 0:
            prev = self.prev_segment(start_index % len(self), use_Z=True)
            next = self.next_segment((end_index - 1) % len(self), use_Z=True)
        else:
            prev = next = None

        is_empty_insert = not isinstance(value, Segment) and len(value) == 0

        if is_empty_insert:
            if prev is not None and next is not None and prev.start != next.end:
                raise ValueError("subpath splice results jumpcut discontinuity")
        else:
            if prev and prev.end != value.start:
                print("prev:", prev)
                print("value:", value)
                raise ValueError("subpath splice results in .start discontinuity")
            if next and next.start != value.end:
                raise ValueError("subpath splice results in .end discontinuity")

        # delete
        for i in range(start_index, end_index):
            del self._segments[i]

        # insert
        for seg in segment_iterator_of(value, back_to_front=True):
            self._segments.insert(start_index, seg)

        self.basic_reset()
        self.check_health()

        return self

    def reversed(self):
        to_return = Subpath(*[seg.reversed() for seg in reversed(self)])
        if self._z:
            to_return.set_Z()
        return to_return

    def __repr__(self):
        return "Subpath({})".format(
            ",\n     ".join(repr(x) for x in self._segments))

    def __eq__(self, other):
        if not isinstance(other, Subpath):
            return NotImplemented
        if len(self) != len(other):
            return False
        if any(s != o for s, o in zip(self._segments, other._segments)):
            return False
        return True

    def __ne__(self, other):
        if not isinstance(other, Subpath):
            return NotImplemented
        return not self == other

    def _calc_lengths(self, error=LENGTH_ERROR, min_depth=LENGTH_MIN_DEPTH):
        if self._length is not None:
            return

        lengths = [each.length(error=error, min_depth=min_depth) for each in
                   self._segments]

        self._length = sum(lengths)
        self._lengths = [each / self._length for each in lengths]

    def point(self, T):
        index, t = self.T2t(T)
        return self._segments[index].point(t)

    def length(self, T0=0, T1=1, error=LENGTH_ERROR, min_depth=LENGTH_MIN_DEPTH):
        self._calc_lengths(error=error, min_depth=min_depth)
        if T0 == 0 and T1 == 1:
            return self._length
        else:
            if len(self) == 1:
                return self[0].length(t0=T0, t1=T1)
            idx0, t0 = self.T2t(T0)
            idx1, t1 = self.T2t(T1)
            if idx0 == idx1:
                return self[idx0].length(t0=t0, t1=t1)
            return (self[idx0].length(t0=t0) +
                    sum(self[idx].length() for idx in range(idx0 + 1, idx1)) +
                    self[idx1].length(t1=t1))

    def ilength(self, s, s_tol=ILENGTH_S_TOL, maxits=ILENGTH_MAXITS,
                error=ILENGTH_ERROR, min_depth=ILENGTH_MIN_DEPTH):
        """Returns a float, t, such that self.length(0, t) is approximately s.
        See the inv_arclength() docstring for more details."""
        return inv_arclength(self, s, s_tol=s_tol, maxits=maxits, error=error,
                             min_depth=min_depth)

    @property
    def start(self):
        if not self._start:
            self._start = self._segments[0].start
        return self._start

    @start.setter
    def start(self, pt):
        self._start = pt
        self._segments[0].start = pt
        assert self._end is not None
        if self._z and self._end != self._start:
            self._z = False

    @property
    def end(self):
        if not self._end:
            self._end = self._segments[-1].end
        return self._end

    @end.setter
    def end(self, pt):
        self._end = pt
        self._segments[-1].end = pt
        assert self._start is not None
        if self._z and self._end != self._start:
            self._z = False

    def Z(self):
        return self._z

    def is_bezier(self):
        return all(isinstance(seg, BezierSegment) for seg in self)

    def d(self, command_number_spacing=COMMAND_TO_NUMBER_SPACE, useSandT=False):
        """
        Returns a path d-string for the path object. For an
        explanation of useSandT see the notes in the README.
        """
        if len(self) == 0:
            assert self._start is None and self._end is None and self._z is False
            return ''

        assert self._start is not None and self._end is not None
        assert self._start == self._segments[0].start
        assert self._end == self._segments[-1].end

        parts = []
        parts.append('M')
        parts.append('{},{}'.format(self._start.real, self._start.imag))
        previous_segment = None

        for index, segment in enumerate(self):
            seg_start = segment.start

            assert previous_segment is None or previous_segment.end == seg_start

            if isinstance(segment, Line):
                args = segment.end.real, segment.end.imag
                if index < len(self) - 1 or not self._z:
                    parts.append('L')
                    parts.append('{},{}'.format(*args))

            elif isinstance(segment, CubicBezier):
                if useSandT and segment.is_smooth_from(previous_segment,
                                                       warning_on=False):
                    args = (segment.control2.real, segment.control2.imag,
                            segment.end.real, segment.end.imag)
                    parts.append('S')
                    parts.append('{},{} {},{}'.format(*args))
                else:
                    args = (segment.control1.real, segment.control1.imag,
                            segment.control2.real, segment.control2.imag,
                            segment.end.real, segment.end.imag)
                    parts.append('C')
                    parts.append('{},{} {},{} {},{}'.format(*args))

            elif isinstance(segment, QuadraticBezier):
                if useSandT and segment.is_smooth_from(previous_segment,
                                                       warning_on=False):
                    args = segment.end.real, segment.end.imag
                    parts.append('T')
                    parts.append('{},{}'.format(*args))
                else:
                    args = (segment.control.real, segment.control.imag,
                            segment.end.real, segment.end.imag)
                    parts.append('Q')
                    parts.append('{},{} {},{}'.format(*args))

            elif isinstance(segment, Arc):
                args = (segment.radius.real, segment.radius.imag,
                        segment.rotation, int(segment.large_arc),
                        int(segment.sweep), segment.end.real, segment.end.imag)
                parts.append('A')
                parts.append('{},{} {} {:d},{:d} {},{}'.format(* args))

            previous_segment = segment

        if self._z:
            assert previous_segment.end == self._start
            parts.append('Z')

        return command_number_spacing.join(parts)

    def joins_smoothly_with(self, previous, wrt_parameterization=False):
        """ (see docstring for is_smooth_join) """
        return is_smooth_join(previous, self)

    def T2address(self, T=None, address=None):
        assert len(self) > 0
        a = address if address else Address()
        a.T = T  # please note this calls a setter in Address---does not overwrite
        if a.T is None:
            raise ValueError("no source of value for T in Subpath.T2address")
        a.segment_index, a.t = self.T2t(a.T)
        return a

    def T2t(self, T):
        """returns the segment index and segment parameter, t, corresponding to
        the path parameter T. This is the inverse of the Subpath.t2T() method."""
        if len(self) == 0:
            raise ValueError(".T2t() called on empty subpath")

        if T == 1:
            return len(self) - 1, 1

        if T == 0:
            return 0, 0

        self._calc_lengths()

        T0 = 0
        for seg_idx, seg_length in enumerate(self._lengths):
            T1 = T0 + seg_length  # the T-value the current segment ends on
            if T1 >= T:
                return seg_idx, (T - T0) / seg_length
            T0 = T1

        assert 0 <= T <= 1
        raise BugException

    def index_of_segment(self, seg):
        indices = [index for index, s in enumerate(self) if s is seg]
        if len(indices) == 1:
            return indices[0]
        elif len(indices) == 0:
            return -1
        else:
            raise ValueError("Subpath.index_of_segment called on multiply- \
                              appearing segment")

    def segment_index_factory(self, thing):
        if isinstance(thing, int):
            segment_index = thing
            if not 0 <= segment_index <= len(self) - 1:
                raise ValueError("segment_index out of range")
        elif isinstance(thing, Segment):
            segment_index = self.index_of_segment(thing)
            if segment_index == -1:
                raise ValueError("segment does not exist in Subpath")
        elif thing is None:
            if len(self) == 1:
                segment_index = 0
            else:
                raise ValueError("no default segment in Subpath of length != 1")
        return segment_index

    def t2T(self, t, segment_index_or_segment):
        """
        Returns a value T such that self.point(T) = self[segment_index](t),
        where 'segment_index' is an int such that:

        segment_index == segment_index_or_segment if the latter is an int in
        the appropriate range, i.e., between 0 and len(self) - 1

        self[segment_index] == segment_index_or_segment if the latter is a
        segment that appears exactly once in the path

        segment_index == 0 if segment_index_or_segment is None and len(self)
        == 1

        Throws a ValueError in all other cases, i.e., if unable to determine
        the correct segment_index from segment_index_or_segment.
        """
        segment_index = self.segment_index_factory(segment_index_or_segment)

        self._calc_lengths()

        T_seg_start = sum(self._lengths[:segment_index])
        T_seg_end = T_seg_start + self._lengths[segment_index]
        assert T_seg_end == sum(self._lengths[: segment_index + 1])
        if t == 1:  # (do not remove! numerically better)
            return T_seg_end
        return (T_seg_end - T_seg_start) * t + T_seg_start

    def t2address(self, t=None, address=None, segment_index=None, segment=None):
        """
        Takes a mandatory value t, 0 <= t <= 1, and an integer segment_index
        (mandatory unless len(self) == 1, and possibly supplied indirectly via
        the 'segment' keyword argument---see previous function---which is
        otherwise not needed), and returns an address a such that

           -- a.t == t
           -- a.segment_index == segment_index
           -- self.point(a.T) == self[segment_index].point(t)
           -- a.subpath_index == None
           -- a.W == None

        Note: t, segment_index are kept as keywords parameters mainly for
        syntactic consistency with Path.T2address.
        """
        a = address if address is not None else Address()

        a.t = t  # (uses a setter, does not overwrite non-None a.t)
        if a.t is None:
            raise ValueError("t value required in Subpath.t2address")

        a.segment_index = segment_index  # (does not overwrite non-None value)
        a.segment_index = self.segment_index_factory(segment or a.segment_index)
        if a.segment_index is None:
            raise ValueError("segment_index required in Subpath.t2address")

        a.T = self.t2T(a.t, a.segment_index)  # (ditto)

        return a

    def derivative(self, T, n=1):
        """returns the tangent vector of the Path at T (centered at the origin
        and expressed as a complex number).
        Note: Bezier curves can have points where their derivative vanishes.
        The unit_tangent() method, by contrast, attempts to compute the direction
        of the derivative at those points as well."""
        seg_idx, t = self.T2t(T)
        seg = self._segments[seg_idx]
        if self._length:
            seg_length = self._lengths[seg_idx] * self._length
        else:
            seg_length = seg.length()
        return seg.derivative(t, n=n) / seg_length**n

    def unit_tangent(self, T):
        """returns the unit tangent vector of the Path at T (centered at the
        origin and expressed as a complex number).  If the tangent vector's
        magnitude is zero, this method will attempt to find the limit of
        self.derivative(tau) / abs(self.derivative(tau)) as tau approaches T.
        See the function bezier_unit_tangent for more details."""
        seg_idx, t = self.T2t(T)
        return self._segments[seg_idx].unit_tangent(t)

    def normal(self, T):
        """returns the (right hand rule) unit normal vector to self at t."""
        return - 1j * self.unit_tangent(T)

    def multisplit(self, Ts):
        """takes a possibly empty list of t-values t1, ..., tn such that
        0 < t1 < ... < tn < 1 and returns a list of Subpaths whose union
        is self and whose endpoints are self.point(0), self.point(t1), ...,
        self.point(tn)"""
        return multisplit_of(self)

    def rotated(self, degs, origin=None):
        """
        Returns a copy of self rotated by `degs` degrees (CCW) around the
        point `origin` (a complex number).  By default `origin` is either
        `self.point(0.5)`, or in the case that self is an Arc object,
        `origin` defaults to `self.center`.
        """
        return rotate(self, degs, origin=origin)

    def translated(self, z0):
        """Returns a copy of self shifted by the complex quantity `z0` such
        that self.translated(z0).point(t) = self.point(t) + z0 for any t."""
        return translate(self, z0)

    def curvature(self, T):
        """returns the curvature of the subpath at T while checking for
        possible non-differentiability at T. Outputs float('inf') if not
        differentiable at T."""
        seg_idx, t = self.T2t(T)

        seg = self[seg_idx]

        if np.isclose(t, 0) and (seg_idx != 0 or self._z):
            previous_seg_in_path = self._segments[
                (seg_idx - 1) % len(self._segments)]
            if not seg.joins_smoothly_with(previous_seg_in_path):
                return float('inf')

        elif np.isclose(t, 1) and (seg_idx != len(self) - 1 or self._z):
            next_seg_in_path = self._segments[
                (seg_idx + 1) % len(self._segments)]
            if not next_seg_in_path.joins_smoothly_with(seg):
                return float('inf')

        # curvature is invariant under reparameterization, so we can
        # use the segment's own parameterization (?):
        return seg.curvature(t)

    def area(self, quality=0.01, safety=3):
        """
        Returns the area enclosed by the subpath, assuming self._end == self._start,
        after converting arcs to cubic bezier segments.

        The 'quality' and 'safety' parameters control the latter conversion. See
        the docstring for Subpath.converted_bezier_path for more details.

        Note: negative area results from CW (as opposed to CCW)
        parameterization of the Path object.
        """
        assert self._end == self._start
        cubicized = self.converted_to_bezier(quality=quality, safety=safety, reuse_segments=True)
        area_enclosed = 0
        for seg in cubicized:
            x         = real(seg.poly())
            dy        = imag(seg.poly()).deriv()
            integrand = x * dy
            integral  = integrand.integ()
            area_enclosed += integral(1) - integral(0)
        return area_enclosed

    def normalize_address(self, a, use_Z=True):
        """ expects a fully-formed address a; tweaks a like so: if a.t == 0
        and the preceding segment's t == 1 address points to the same point on the
        subpath (i.e., the previous segment exists!), change a to use t == 1
        and the segment_index of the previous segment.

        The 'use_Z' option determines whether self[len(self) - 1] counts as a
        "previous segment" of self[0], assuming self.Z() == True. """
        assert a.T is not None
        assert a.segment_index is not None
        assert a.t is not None
        if a.t == 0:
            prev_index = self.mod_index(a.segment_index - 1, use_Z=use_Z)
            prev = self[prev_index] if prev_index >= 0 else None
            if prev is not None:
                assert prev.end == self[a.segment_index].start
                b = self.t2address(t=1, segment_index=prev_index)
                # if b.T != a.T:
                #     print("a:", a, "b:", b)
                assert b.T == a.T or (b.T == 1 and a.T == 0 and self._z)
                a._segment_index = prev_index
                a._t = 1
                a._T = b.T

    def intersect(self, other_curve, justonemode=False, tol=1e-12, normalize=True):
        """
        Returns list of pairs of named Intersection objects by taking all the pairwise
        intersections of segments in self and segments in other_curve. Here
        other_curve can be either a segment, a subpath, or a path. In the latter case,
        the latter case, other_curve.intersect(self) is called, and intersections are
        swapped.

        If the two curves coincide for more than a finite set of
        points, this code will (should) fail.

        If justonemode==True, then returns the first intersection found.

        tol is used as a parameter passed to segment.intersect(segment); see
        implementations of segment.intersect

        If normalize==True, remove all 't==0' intersections (beginning-of-segment
        intersections) that can be written as 't==1' intersections (end-of-segment
        intersections). This option can be useful to avoid duplicates.

        Note: If the respective subpath is a geometric loop but not a topological
        loop (i.e., _z has not been set for the subpath), the adjacency between
        the first and last segments of the subpath is ignored by the 'normalize'
        option. Similarly, if other_curve is made up of several subpaths that are
        adjacent at their endpoints, these adjacencies are ignored by 'normalize'.)
        """
        subpath1 = self

        if isinstance(other_curve, Subpath):
            subpath2 = other_curve
        elif isinstance(other_curve, Segment):
            subpath2 = Subpath(other_curve)
        elif isinstance(other_curve, Path):
            reversed_intersections = other_curve.intersect(self, normalize=normalize)
            return [(a2, a1) for (a1, a2) in reversed_intersections]
        else:
            assert False

        # let...

        def append_new_intersection(a1, a2):
            if normalize:
                subpath1.normalize_address(a1)
                subpath2.normalize_address(a2)
                if (a1, a2) not in intersection_list:
                    intersection_list.append((a1, a2))
            else:
                intersection_list.append((a1, a2))

        # in...

        intersection_list = []
        for si1, seg1 in enumerate(subpath1):
            for si2, seg2 in enumerate(subpath2):
                for a1, a2 in seg1.intersect(seg2, tol=tol):
                    a1 = subpath1.t2address(segment_index=si1, address=a1)
                    a2 = subpath2.t2address(segment_index=si2, address=a2)
                    append_new_intersection(a1, a2)
                    if justonemode:
                        return intersection_list[0]

        return intersection_list

    def bbox(self):
        """returns a bounding box for the input Path object in the form
        (xmin, xmax, ymin, ymax)."""
        bbs = [seg.bbox() for seg in self]
        xmins, xmaxs, ymins, ymaxs = list(zip(*bbs))
        xmin = min(xmins)
        xmax = max(xmaxs)
        ymin = min(ymins)
        ymax = max(ymaxs)
        return xmin, xmax, ymin, ymax

    def point_outside(self):
        """ returns an arbitrary point outside the path's convex hull (barring
        some unforeseen integer-overflow type catastrophy)"""
        xmin, xmax, ymin, ymax = self.bbox()
        return xmin - 100 + (ymin - 100) * 1j

    def encloses(self, z):
        """returns true if and only if """
        if not self.isclosed():
            raise ValueError("Subpath.encloses called on non-closed path")
        intersections = self.intersect(Line(z, self.point_outside()))
        distinct_line_ts = {i[0].t for i in intersections}
        return bool(len(distinct_line_ts) % 2)

    def cropped(self, T0, T1, drop_small=False):
        """
        Returns a cropped copy of the subpath starting at self.point(T0) and
        ending at self.point(T1); T0 and T1 must be distinct, 0 < T0, T1 < 1.

        If T1 < T0 the crop interpreted as a wraparound crop. In that case the
        path's ._z value must be set.

        If drop_small==True, initial and final subpath segments seg such that
        np.isclose(seg.start, seg.eng) are dropped from the final answer.
        This can (theoretically) result in an empty subpath being returned."""
        assert T0 != T1 and 0 <= T0 <= 1 and 0 <= T1 <= 1
        assert len(self) > 0

        seg1_index, t_seg1 = self.T2t(T1)
        seg0_index, t_seg0 = self.T2t(T0)

        seg0 = self[seg0_index]
        seg1 = self[seg1_index]

        if T0 < T1:
            if seg1_index == seg0_index:
                new_path = Subpath(seg0.cropped(t_seg0, t_seg1))
            else:
                new_path = Subpath(seg0.cropped(t_seg0, 1))
                for index in range(seg0_index + 1, seg1_index):
                    new_path.append(self[index])
                new_path.append(seg1.cropped(0, t_seg1))

        else:
            if not self.Z():
                raise ValueError("Subpath.cropped does not support wraparound \
                                  cropping for non-topoligcally closed subpaths")
            # wraparound case necessarily involves more than one segment
            assert seg1_index <= seg0_index
            new_path = Subpath(seg0.cropped(t_seg0, 1))
            index = (seg0_index + 1) % len(self)
            while index != seg0_index:
                new_path.append(self[index])
                index = (index + 1) % len(self)
            new_path.append(seg1.cropped(0, t_seg1))
            assert len(new_path) >= 2

        if drop_small:
            if np.isclose(new_path[-1].start, new_path[-1].end):
                del new_path[-1]
            if len(new_path) > 0 and np.isclose(new_path[0].start, new_path[0].end):
                del new_path[0]

        return new_path

    def radialrange(self, origin, return_all_global_extrema=False):
        """returns the tuples (d_min, t_min, idx_min), (d_max, t_max, idx_max)
        which minimize and maximize, respectively, the distance
        d = |self[idx].point(t) - origin|."""
        if return_all_global_extrema:
            raise NotImplementedError
        else:
            global_min = (np.inf, None, None)
            global_max = (0, None, None)
            for seg_idx, seg in enumerate(self):
                seg_global_min, seg_global_max = seg.radialrange(origin)
                if seg_global_min[0] < global_min[0]:
                    global_min = seg_global_min + (seg_idx, )
                if seg_global_max[0] > global_max[0]:
                    global_max = seg_global_max + (seg_idx, )
            return global_min, global_max

    def converted_to_bezier(self, quality=0.01, safety=5, reuse_segments=True):
        # warning: reuses same segments when available, unless 'reuse_segments'
        #          is set to False
        new_subpath = Subpath()
        for s in self:
            if isinstance(s, Arc):
                cpath, _ = s.converted_to_bezier_subpath(quality, safety)
                new_subpath.extend(cpath)
            else:
                possibly_new_segment = s if reuse_segments else translate(s, 0)
                new_subpath.append(possibly_new_segment)
        if self._z:
            new_subpath.set_Z()
            print("Yes you are going through here; now new_subpath.Z():", new_subpath.Z())
        return new_subpath

    def set_Z(self, forceful=False):
        assert self._start == self[0].start
        assert self._end == self[-1].end
        if self._z:
            assert self._start == self._end
            print("warning: Z already set is Subpath.set_Z; ignoring")
        else:
            if self._start != self._end:
                if forceful:
                    self.append(Line(self._end, self._start))
                    self.basic_reset()
                else:
                    raise ValueError("self._end != self._start in Subpath.set_Z")
            self._z = True

    def offset(self, amount, two_sided=False, quality=0.01, safety=10, join='miter', miter_limit=4, cap='butt'):
        converted = self.converted_to_bezier(quality, safety, True)
        way_out   = Path()
        way_in    = Path()
        skeleton  = Subpath()

        for seg in converted:
            wo, sk, wi = seg.offset(amount, two_sided=two_sided, quality=quality, safety=safety)
            assert isinstance(wo, Subpath)
            assert isinstance(wi, Subpath)
            assert isinstance(sk, Subpath)
            way_out.append(wo)
            skeleton.extend(sk)
            way_in.prepend(wi)

        if self._z:
            skeleton.set_Z()

        # joins
        reversed_skeleton = skeleton.reversed() if two_sided else Subpath()
        assert len(reversed_skeleton) == way_in.num_segments()
        both_skeleton_offset_pairs = [(skeleton, way_out), (reversed_skeleton, way_in)]
        both_results = []

        for ske_off_pair in both_skeleton_offset_pairs:
            joined = join_offset_segments_into_subpath(ske_off_pair[0], ske_off_pair[1], amount, join, miter_limit)
            both_results.append(joined)
            assert joined.Z() == self.Z()

        way_out = both_results[0]  # Subpath
        way_in  = both_results[1]  # Subpath

        assert isinstance(way_out, Subpath)
        assert isinstance(way_in, Subpath)

        # line caps
        if two_sided and cap != 'none':
            if skeleton.Z():
                pass
            else:
                c_end = endcap_for_curve(skeleton, amount, cap)
                c_start = endcap_for_curve(skeleton.reversed(), amount, cap)
                assert c_start.start == way_in.end
                assert c_start.end   == way_out.start
                assert c_end.start == way_out.end
                assert c_end.end   == way_in.start
                way_out.extend(c_end)
                way_in.extend(c_start)
                way_out.extend(way_in)
                way_out.set_Z(forceful=False)
                way_in = Subpath()

        return way_out, skeleton, way_in

    def stroke(self, width, quality=0.01, safety=5, join='miter', miter_limit=4, cap='butt'):
        return self.offset(width / 2, two_sided=True, quality=quality, safety=safety, join=join, miter_limit=miter_limit, cap=cap)[0]


class Path():
    def __init__(self, *things):
        assert all(isinstance(s, Subpath) for s in things)
        self._subpaths = list(things)
        self.basic_reset()

    def basic_reset(self):
        self._length = None
        self._num_segments = sum([len(x) for x in self])
        self._start = self.recomputed_start()
        self._end = self.recomputed_end()
        assert (self._start is None) == (self._end is None) == (self._num_segments == 0)

    def num_segments(self):
        assert self._num_segments == sum([len(x) for x in self])
        return self._num_segments

    def first_nonempty(self):
        for index, x in enumerate(self):
            if len(x) > 0:
                return index, x
        return None

    def last_nonempty(self):
        for index, x in enumerate(reversed(self)):
            if len(x) > 0:
                return index, x
        return None

    def recomputed_start(self):
        for x in self:
            if len(x) > 0:
                return x.start
        return None

    def recomputed_end(self):
        for x in reversed(self):
            if len(x) > 0:
                return x.end
        return None

    def __getitem__(self, index):
        return self._subpaths[index]

    def __setitem__(self, index, value):
        assert isinstance(value, Subpath)
        assert value not in self
        self._subpaths[index] = value
        self.basic_reset()
        return value

    def append(self, value, even_if_empty=False):
        assert isinstance(value, Subpath)
        assert value not in self
        if len(value) > 0 or even_if_empty:
            self._subpaths.append(value)
            self.basic_reset()

    def prepend(self, value):
        self.insert(0, value)

    def insert(self, index, value):
        assert isinstance(value, Subpath)
        assert value not in self
        self._subpaths.insert(index, value)
        self.basic_reset()

    def extend(self, values):
        for v in values:
            assert isinstance(v, Subpath)
            assert v not in self
            self._subpaths.append(v)
        self.basic_reset()

    def __delitem__(self, index):
        to_return = self._subpaths[index]
        self.basic_reset()
        return to_return

    def __iter__(self):
        return self._subpaths.__iter__()

    def __contains__(self, x):
        return self._subpaths.__contains__(x)

    def reversed(self):
        """returns a copy of the Path object with its orientation reversed."""
        newsubpaths = [subpath.reversed() for subpath in self]
        return Path(*newsubpaths.reverse())

    def __len__(self):
        return len(self._subpaths)

    def __repr__(self):
        return "Path({})".format(
               ",\n     ".join(repr(x) for x in self))

    def __eq__(self, other):
        if not isinstance(other, Path):
            return NotImplemented
        if len(self) != len(other):
            return False
        if any(s != o for s, o in zip(self._subpaths, other._subpaths)):
            return False
        return True

    def __ne__(self, other):
        if not isinstance(other, Path):
            return NotImplemented
        return not self == other

    def segment_iterator(self, back_to_front=False):
        if back_to_front:
            for t in reversed(self):
                for seg in reversed(t):
                    yield seg
        else:
            for t in self:
                for seg in t:
                    yield seg

    def is_bezier(self):
        return all(s.is_bezier() for s in self)

    def is_naively_continuous(self):
        prev = None
        for s in self.segment_iterator():
            if prev and prev.end != s.start:
                return False
            prev = s
        return True

    def _calc_lengths(self, error=LENGTH_ERROR, min_depth=LENGTH_MIN_DEPTH):
        if self._length is not None:
            return

        lengths = [each.length(error=error, min_depth=min_depth) for each in
                   self]
        self._length = sum(lengths)
        self._lengths = [each / self._length for each in lengths]

    def point(self, W_or_address):
        # W: name of time parameter for Path
        # T: .......................for Subpath
        # t: .......................for Segment
        a = self.resolve_W_or_address(W_or_address)
        self.W2t(a)
        return self[a.subpath_index][a.segment_index].point(a.t)

    def length(self, W0=0, W1=1, error=LENGTH_ERROR, min_depth=LENGTH_MIN_DEPTH):
        self._calc_lengths(error=error, min_depth=min_depth)
        if W0 == 0 and W1 == 1:
            return self._length
        else:
            if len(self) == 1:
                return self[0].length(T0=W0, T1=W1)
            index_0, T0 = self.W2T(W0)
            index_1, T1 = self.W2T(W1)
            if index_0 == index_1:
                return self[index_0].length(T0=T0, T1=T1)
            return (self[index_0].length(T0=T0) +
                    sum(self[i].length() for i in range(index_0 + 1, index_1)) +
                    self[index_1].length(T1=T1))

    def ilength(self, s, s_tol=ILENGTH_S_TOL, maxits=ILENGTH_MAXITS,
                error=ILENGTH_ERROR, min_depth=ILENGTH_MIN_DEPTH):
        """Returns a float, t, such that self.length(0, t) is approximately s.
        See the inv_arclength() docstring for more details."""
        return inv_arclength(self, s, s_tol=s_tol, maxits=maxits, error=error,
                             min_depth=min_depth)

    def iscontinuous(self):
        num_nonempty = sum([1 for x in self if len(x) > 0])
        return num_nonempty == 1

    def isloop(self):
        if self.iscontinuous() and self.first_nonempty()[1].Z():
            return True
        return False

    def issetofloops(self):
        return all(s.Z() for s in [s for s in self if len(s) > 0])

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, pt):
        self.first_nonempty()[1].start = pt
        self.basic_reset()

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, pt):
        self.last_nonempty()[1].end = pt
        self.basic_reset()

    def d(self, subpath_spacing=SUBPATH_TO_SUBPATH_SPACE, command_number_spacing=COMMAND_TO_NUMBER_SPACE, useSandT=False):
        """
        Returns a path d-string for the path object. For an
        explanation of useSandT see the notes in the README.
        """
        return subpath_spacing.join(s.d(command_number_spacing=command_number_spacing, useSandT=useSandT) for s in [s for s in self if len(s) > 0])

    def joins_smoothly_with(self, previous, wrt_parameterization=False):
        """Checks if this Path object joins smoothly with previous
        path / segment.  By default, this only checks that this Path starts
        moving (at t=0) in the same direction (and from the same positive) as
        previous stopped moving (at t=1).  To check if the tangent magnitudes
        also match, set wrt_parameterization=True."""
        return is_smooth_join(previous, self)

    def subpath_index_factory(self, thing, should_exist=True):
        if isinstance(thing, int):
            if not 0 <= thing <= len(self) - 1:
                raise ValueError("bad index in Path.subpath_index_factory")
            return thing
        elif isinstance(thing, Subpath):
            for index, p in enumerate(self):
                if p is thing:
                    return index
        elif thing is None:
            if len(self) == 1:
                return 0
            else:
                raise ValueError("more than one subpath to choose from for \
                                  None option Path.subpath_index_factory()")
        if should_exist:
            raise ValueError("subpath not found in subpath_index_factory")
        return -1

    def W2T(self, W):
        if self._num_segments == 0:
            raise ValueError("can't compute a position on empty path")

        self._calc_lengths()

        W0 = 0
        for index, subpath in enumerate(self._subpaths):
            W1 = W0 + self._lengths[index]
            if W <= W1:
                return index, (W - W0) / (W1 - W0)
            W0 = W1

        assert 0 <= W <= 1
        raise BugException

    def W2address(self, W_or_address):
        """ constructs and returns a full address object from a W value, or fills in
        remaining fields of an unfinished address with a W value """

        if isinstance(W_or_address, Address):
            a = W_or_address
            W = a.W
        else:
            W = W_or_address
            a = Address(W=W)

        if W is None:
            raise ValueError("W is None in Path.W2address")

        # shortcuts:
        if W == 1:
            a.subpath_index = self.last_nonempty()[0]
            a.T = 1
            a.segment_index = len(self[a.subpath_index]) - 1
            a.t = 1
            return a
        if W == 0:
            a.subpath_index = self.first_nonempty()[0]
            a.T = 0
            a.segment_index = 0
            a.t = 0
            return a

        a.subpath_index, a.T = self.W2T(W)
        return self[a.subpath_index].T2address(a)

    def T2W(self, T, subpath_index_or_subpath):
        subpath_index = self.subpath_index_factory(subpath_index_or_subpath, should_exist=True)

        self._calc_lengths()

        W_seg_start = sum(self._lengths[:subpath_index])
        W_seg_end = W_seg_start + self._lengths[subpath_index]
        assert W_seg_end == sum(self._lengths[: subpath_index + 1])
        if T == 1:  # (do not remove! numerically important)
            return W_seg_end
        return (W_seg_end - W_seg_start) * T + W_seg_start

    def T2address(self, T=None, subpath_index=None, subpath=None, address=None):
        """
        Takes one of:

           -- a T-value and either 'subpath_index' or 'subpath'
           -- an incomplete Address address with at least its T- and subpath_index set

        (OK: actually, subpath_index/subpath are not even necessary if self has only 1
        nonempty subpath) and returns a full address matching the above fields, such that

           self.point(a.W) == self[subpath_index].point(T)

        where subpath_index is a.subpath_index or the input parameter subpath_index,
        depending on which is defined.

        Note: The address returned is the same object as the provided address, if
        provided; hence this function has (intentional) side effects!
        """
        a = address if address else Address()
        a.T = T
        if a.T is None:
            raise ValueError("T value missing in Path.T2address")

        # subpath_index extraction---nightmare of three-way argument supply:
        a.subpath_index = subpath_index
        a.subpath_index = self.subpath_index_factory(subpath)
        if a.subpath_index is None:
            raise ValueError("subpath_index in Path.T2address")

        if not 0 <= a.subpath_index <= len(self) - 1:
            raise ValueError("bad subpath_index in Path.T2address")

        # W
        a.W = self.T2W(a.T, a.subpath_index)

        # automatically completes / checks missing fields of a, if any:
        return self[a.subpath_index].T2address(a)

    def W2t(self, W):
        """ returns a subpath_index, segment_index, t triple """
        subpath_index, T = self.W2T(W)
        return subpath_index, self[subpath_index].T2t(T)

    def t2W(self, t, segment_index, subpath_index):
        return self.T2W(self[subpath_index].t2T(t, segment_index), subpath_index)

    def resolve_W_or_address(self, thing):
        if isinstance(thing, Address):
            return thing
        return self.W2address(thing)

    def derivative(self, W_or_address, n=1):
        """
        Given an address a or a value W, which is resolved to the default address a,
        returns the n-th derivative of the path's W-parameterization at a.
        """
        a = self.resolve_W_or_address(W_or_address)
        self._calc_lengths()
        return self[a.subpath_index].derivative(a, n=n) / self._lengths[a.subpath_index]**n

    def unit_tangent(self, W_or_address):
        """
        Given an address a or a value W, which is resolved to the default address a,
        returns self[a.subpath_index][a.segment_index].unit_tangent(a.t).
        """
        a = self.resolve_W_or_address(W_or_address)
        return self[a.subpath_index][a.segment_index].unit_tangent(a.t)

    def normal(self, W_or_address):
        """
        Given an address a or a value W, which is resolved to the default address a,
        or given an address a, returns the right-hand rule unit normal vector to self
        at a.
        """
        return - 1j * self.unit_tangent(W_or_address)

    def curvature(self, W_or_address):
        """
        Given an address a or a value W, which is resolved to the default address a,
        the curvature of the path at a, outputting float('inf') if the path is not
        differentiable at a.
        """
        a = self.resolve_W_or_address(W_or_address)
        return self[a.subpath_index][a.segment_index].curvature(a.t)

    def area(self, quality=0.01, safety=3):
        """
        Returns the directed area enclosed by the path; requires each subpath to
        be geometrically closed. (But with or without having its Z property set.)

        Negative area results from CW (as opposed to CCW) parameterization of a
        Subpath.

        Elliptical arc segments are converted to bezier paths [nb: the original
        segments remain in place]; the 'quality' and 'safety' parameters control
        this conversion, as described in Subpath.converted_to_bezier().
        """
        return sum([s.area(quality, safety) for s in self])

    def intersect(self, other_curve, normalize=True, justonemode=False, tol=1e-12):
        """
        Returns a list of pairs of addresses (a1, a2): for each intersection found,
        a1 is the intersection's adddress with respect to this path, a2 with respect
        to the other path.

        If normalize==False, makes a list of intersections between all possible pairs
        of segments in the two paths.

        If normalize==True, normalizes addresses that fall at the beginning
        of a segment (i.e., addresses with t == 0) to be included in the subpath's
        previous segment, if present (so addresses with t == 0 switch to addresses
        with t == 1, when possible), and removes duplicates. This option be helpful
        to remove duplicate intersections.

        'other_curve' can be either a path, a subpath or a segment; in the latter
        two cases, the subpath or segment is wrapper in a path resp. / subpath and
        path before computing the intersections.

        If justonemode==True, returns just the first intersection found.

        tol is used to check for redundant intersections (see comment above
        the code block where tol is used).

        Fails if the two path objects coincide for more than a finite set of points.
        """
        path1 = self
        if isinstance(other_curve, Path):
            path2 = other_curve
        elif isinstance(other_curve, Subpath):
            path2 = Path(other_curve)
        elif isinstance(other_curve, Segment):
            path2 = Path(Subpath(other_curve))
        else:
            raise ValueError("bad type for other_curve in Path.intersect()")

        assert path1 != path2

        intersection_list = []

        for ind1, sub1 in enumerate(path1):
            for ind2, sub2 in enumerate(path2):
                tweenies = sub1.intersect(sub2,
                                          normalize=normalize,
                                          justonemode=justonemode,
                                          tol=tol)
                for a1, a2 in tweenies:
                    a1.subpath_index = ind1
                    a2.subpath_index = ind2
                    a1.W = path1.T2W(a1.T, ind1)
                    a2.W = path2.T2W(a2.T, ind2)

                if justonemode and len(tweenies) > 0:
                    return tweenies[0]

                intersection_list.extend(tweenies)

        return intersection_list

    def bbox(self):
        """returns a bounding box for the input Path object in the form
        (xmin, xmax, ymin, ymax)."""
        bbs = [subpath.bbox() for subpath in self]
        xmins, xmaxs, ymins, ymaxs = list(zip(*bbs))
        xmin = min(xmins)
        xmax = max(xmaxs)
        ymin = min(ymins)
        ymax = max(ymaxs)
        return xmin, xmax, ymin, ymax

    def point_outside(self):
        """ (there are surely some integer-overflow corner cases, etc, but
        seems pretty reasonable for typical values in an svg document?) """
        xmin, xmax, ymin, ymax = self.bbox()
        return xmin - 42 + (ymin - 43) * 1j

    def encloses(self, z):
        assert self.isclosed()
        intersections = Path(Line(z, self.point_outside())).intersect(self)
        distinct = {i[0][2] for i in intersections}
        return bool(len(distinct) % 2)

    def multisplit(self, Ts):
        return multisplit_of(self, Ts)

    def cropped(self, T0, T1):
        """returns a cropped copy of the path."""
        assert T0 != T1
        if T1 == 1:
            seg1 = self[- 1]
            t_seg1 = 1
            i1 = len(self) - 1
        else:
            seg1_idx, t_seg1 = self.T2t(T1)
            seg1 = self[seg1_idx]
            if np.isclose(t_seg1, 0):
                i1 = (self.index(seg1) - 1) % len(self)
                seg1 = self[i1]
                t_seg1 = 1
            else:
                i1 = self.index(seg1)
        if T0 == 0:
            seg0 = self[0]
            t_seg0 = 0
            i0 = 0
        else:
            seg0_idx, t_seg0 = self.T2t(T0)
            seg0 = self[seg0_idx]
            if np.isclose(t_seg0, 1):
                i0 = (self.index(seg0) + 1) % len(self)
                seg0 = self[i0]
                t_seg0 = 0
            else:
                i0 = self.index(seg0)

        if T0 < T1 and i0 == i1:
            new_path = Path(seg0.cropped(t_seg0, t_seg1))
        else:
            new_path = Path(seg0.cropped(t_seg0, 1))

            # T1 < T0 must cross discontinuity case
            if T1 < T0:
                if self.isclosed():
                    raise ValueError("This path is not closed, thus T0 must "
                                     "be less than T1.")
                else:
                    for i in range(i0 + 1, len(self)):
                        new_path.append(self[i])
                    for i in range(0, i1):
                        new_path.append(self[i])

            # T0 < T1 straight - forward case
            else:
                for i in range(i0 + 1, i1):
                    new_path.append(self[i])

            if t_seg1 != 0:
                new_path.append(seg1.cropped(0, t_seg1))
        return new_path

    def radialrange(self, origin, return_all_global_extrema=False):
        """returns the tuples (d_min, t_min, idx_min), (d_max, t_max, idx_max)
        which minimize and maximize, respectively, the distance
        d = |self[idx].point(t) - origin|."""
        if return_all_global_extrema:
            raise NotImplementedError
        else:
            global_min = (np.inf, None, None)
            global_max = (0, None, None)
            for seg_idx, seg in enumerate(self):
                seg_global_min, seg_global_max = seg.radialrange(origin)
                if seg_global_min[0] < global_min[0]:
                    global_min = seg_global_min + (seg_idx, )
                if seg_global_max[0] > global_max[0]:
                    global_max = seg_global_max + (seg_idx, )
            return global_min, global_max

    def rotated(self, degs, origin=None):
        """
        Returns a copy of self rotated by `degs` degrees (CCW) around the
        point `origin` (a complex number).  By default `origin` is either
        `self.point(0.5)`, or in the case that self is an Arc object,
        `origin` defaults to `self.center`.
        """
        return rotate(self, degs, origin=origin)

    def translated(self, z0):
        """Returns a copy of self shifted by the complex quantity `z0` such
        that self.translated(z0).point(t) = self.point(t) + z0 for any t."""
        return translate(self, z0)

    def converted_to_bezier(self, quality=0.01, safety=5, reuse_segments=True):
        # warning: reuses same segments when available, unless 'reuse_segments'
        #          is set to False
        newsubpaths = []
        for s in self:
            newsubpaths.append(s.converted_to_bezier(quality=quality, safety=safety, reuse_segments=reuse_segments))
        return Path(*[newsubpaths])

    def concatenate_with(self, other_path):
        for s in other_path:
            self.append(s)
        return self

    def offset(self, amount, two_sided=False, quality=0.01, safety=5, join='miter', miter_limit=4, cap='butt'):
        skeletons = Path()
        offsets   = Path()
        for s in self:
            wo, sk, wi = s.offset(amount, two_sided, quality, safety, join, miter_limit, cap)
            assert isinstance(wo, Subpath)
            assert isinstance(wi, Subpath)
            assert isinstance(sk, Subpath)
            skeletons.append(sk)
            offsets.append(wo)
            offsets.append(wi)
        return offsets, skeletons

    def stroke(self, width, quality=0.01, safety=5, join='miter', miter_limit=4, cap='butt'):
        return self.offset(width / 2, two_sided=True, quality=quality, safety=safety, join=join, miter_limit=miter_limit, cap=cap)[0]
