import math

TAU = math.pi * 2

def mapToEllipse(P, rx, ry, cosphi, sinphi, centerx, centery):
    x, y = P
    x *= rx
    y *= ry

    xp = cosphi * x - sinphi * y
    yp = sinphi * x + cosphi * y

    return (xp + centerx, yp + centery)

def approxUnitArc(ang1, ang2):
    # If 90 degree circular arc, use a ant
    # as derived from http:#spencermortensen.com/articles/bezier-circle
    a = (0.551915024494 if (ang2 == 1.5707963267948966)
        else (-0.551915024494 if (ang2 == -1.5707963267948966)
            else 4 / 3 * math.tan(ang2 / 4)))

    x1 = math.cos(ang1)
    y1 = math.sin(ang1)
    x2 = math.cos(ang1 + ang2)
    y2 = math.sin(ang1 + ang2)

    return [
        (
            x1 - y1 * a,
            y1 + x1 * a
        ),
        (
            x2 + y2 * a,
            y2 - x2 * a
        ),
        (
            x2,
            y2
        )
    ]

def vectorAngle(ux, uy, vx, vy):
    sign = -1 if (ux * vy - uy * vx < 0) else 1

    dot = min(1, max(-1, ux * vx + uy * vy))

    return sign * math.acos(dot)

def getArcCenter(
    px,
    py,
    cx,
    cy,
    rx,
    ry,
    largeArcFlag,
    sweepFlag,
    sinphi,
    cosphi,
    pxp,
    pyp
):
    rxsq = math.pow(rx, 2)
    rysq = math.pow(ry, 2)
    pxpsq = math.pow(pxp, 2)
    pypsq = math.pow(pyp, 2)

    radicant = max(0, (rxsq * rysq) - (rxsq * pypsq) - (rysq * pxpsq))

    radicant /= (rxsq * pypsq) + (rysq * pxpsq)
    radicant = math.sqrt(radicant) * (-1 if largeArcFlag == sweepFlag else 1)

    centerxp = radicant * rx / ry * pyp
    centeryp = radicant * -ry / rx * pxp

    centerx = cosphi * centerxp - sinphi * centeryp + (px + cx) / 2
    centery = sinphi * centerxp + cosphi * centeryp + (py + cy) / 2

    vx1 = (pxp - centerxp) / rx
    vy1 = (pyp - centeryp) / ry
    vx2 = (-pxp - centerxp) / rx
    vy2 = (-pyp - centeryp) / ry

    ang1 = vectorAngle(1, 0, vx1, vy1)
    ang2 = vectorAngle(vx1, vy1, vx2, vy2)

    if (sweepFlag == 0 and ang2 > 0):
        ang2 -= TAU

    if (sweepFlag == 1 and ang2 < 0):
        ang2 += TAU

    return [ centerx, centery, ang1, ang2 ]


def arc_to_bezier(P, C, R, xAxisRotation = 0, largeArcFlag = 0, sweepFlag = 0, subdivisionRate=1):
    curves = []
    px,py = P.real,P.imag
    cx,cy = C.real,C.imag
    rx,ry = R.real,R.imag

    if (rx == 0 or ry == 0):
        return []

    sinphi = math.sin(xAxisRotation * TAU / 360)
    cosphi = math.cos(xAxisRotation * TAU / 360)

    pxp = cosphi * (px - cx) / 2 + sinphi * (py - cy) / 2
    pyp = -sinphi * (px - cx) / 2 + cosphi * (py - cy) / 2

    if (pxp == 0 and pyp == 0):
        return []

    rx = abs(rx)
    ry = abs(ry)

    lamb = pxp**2 / rx**2 + pyp**2 / ry**2

    if (lamb > 1):
        rx *= math.sqrt(lamb)
        ry *= math.sqrt(lamb)

    centerx, centery, ang1, ang2 = getArcCenter(
        px,
        py,
        cx,
        cy,
        rx,
        ry,
        largeArcFlag,
        sweepFlag,
        sinphi,
        cosphi,
        pxp,
        pyp
    )

    # If 'ang2' == 90.0000000001, then `ratio` will evaluate to
    # 1.0000000001. This causes `segments` to be greater than one, which is an
    # unecessary split, and adds extra points to the bezier curve. To alleviate
    # this issue, we round to 1.0 when the ratio is close to 1.0.
    ratio = abs(ang2) / (TAU / 4)
    if (abs(1.0 - ratio) < 0.0000001):
        ratio = 1.0

    segments = subdivisionRate*max(math.ceil(ratio), 1)

    ang2 /= segments

    for _ in range(segments):
        curves.append(approxUnitArc(ang1, ang2))
        ang1 += ang2

    newcurves = []
    for curve in curves:
        x1, y1 = mapToEllipse(curve[ 0 ], rx, ry, cosphi, sinphi, centerx, centery)
        x2, y2 = mapToEllipse(curve[ 1 ], rx, ry, cosphi, sinphi, centerx, centery)
        x, y = mapToEllipse(curve[ 2 ], rx, ry, cosphi, sinphi, centerx, centery)

        newcurves.append((x1+1j*y1, x2+1j*y2, x+1j*y))
    return newcurves 

