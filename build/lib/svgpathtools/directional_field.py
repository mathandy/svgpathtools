def directional_field(curve, tvals=np.linspace(0, 1, N), asize=1e-2, 
                      colored=False):

    size = asize * curve.length()
    arrows = []
    tvals = np.linspace(0, 1, N)
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
        return Path(arrows)