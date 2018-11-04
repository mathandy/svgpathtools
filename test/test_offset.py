import svgpathtools
import math


# fill out this stub to print/render svg path with preferred method:
def output_stub(dstring, attr):
    pass


# TEST 1 intersection accuracy

if True:
    for norm in {True, False}:
        for parser1 in {svgpathtools.parse_path, svgpathtools.parse_subpath}:
            p1 = parser1('M 0,0 0,1 1,1 1,0 Z M 2,0 2,1 3,1, 3,0 Z')
            for parser2 in {svgpathtools.parse_path, svgpathtools.parse_subpath}:
                p2 = parser2('M -1,-1 2,2 3,1')

                intersections12 = p1.intersect(p2, normalize=norm)
                intersections21 = p2.intersect(p1, normalize=norm)

                assert len(intersections12) == len(intersections21) == 3 * (1 + bool(not norm))

                for a1, a2 in intersections12:
                    assert (a2, a1) in intersections21

                for a2, a1 in intersections21:
                    assert (a1, a2) in intersections12

                distinct_a2_ts = {a2.t for a2, a1 in intersections21}

                assert distinct_a2_ts == {1 / 3, 2 / 3, 1}


# TEST 2 area computations (with double square and unit circle)

if True:
    for parser1 in {svgpathtools.parse_path, svgpathtools.parse_subpath}:
        p1 = parser1('M 0,0 0,1 1,1 1,0 Z M 2,0 2,1 3,1, 3,0 Z')
        assert p1.area() in {2, -2}
        assert p1.reversed().area() == -p1.area()

        p3 = parser1('M 6,6 a 1,1 0 0 0 -1,1 a 1,1 0 0 0 1,1 a 1,1 0 0 0 1,-1 a 1,1 0 0 0 -1,-1 Z')
        print("p3.area():", p3.area())
        print("p3.reversed().area():", p3.reversed().area())
        print("p3.reversed().area(quality=0.01):", p3.reversed().area(quality=0.01))
        print("p3.reversed().area(quality=0.001):", p3.reversed().area(quality=0.001))
        print("p3.reversed().area(quality=0.0001):", p3.reversed().area(quality=0.0001))
        print("p3.reversed().area(quality=0.00001):", p3.reversed().area(quality=0.00001))


# TEST 3 length computations

if True:
    for parser1 in {svgpathtools.parse_path, svgpathtools.parse_subpath}:
        p1 = parser1('M 0,0 0,1 1,1 1,0 Z M 2,0 2,1 3,1, 3,0 Z')
        p3 = parser1('M 6,6 a 1,1 0 0 0 -1,1 a 1,1 0 0 0 1,1 a 1,1 0 0 0 1,-1 a 1,1 0 0 0 -1,-1 Z')

        assert p1.length() == 8
        assert p1.reversed().length() == 8
        assert p1.length(0.25, 0.5) == 2
        assert p1.length(0.25, 0.75) == 4

        print("p3.length:", p3.length())
        print("tau      :", 2 * math.pi)


# TEST 4 arc -> cubic conversion

if True:
    p5 = svgpathtools.parse_subpath('M 2,1 A 1,2 45 1 0 2,0')
    arc = p5[0]
    assert isinstance(arc, svgpathtools.Arc)
    p6, _ = arc.converted_to_bezier_subpath()
    output_stub(p5.d(), {'stroke': 'black', 'stroke-width': '0.15', 'fill': 'none'})
    output_stub(p6.d(), {'stroke': 'red', 'stroke-width': '0.1', 'fill': 'none'})


# TEST 5 offsets, strokes

if True:
    p2 = svgpathtools.parse_path('M 0,0 0,2 2,2 2,1 A 1,2 45 1 0 2,0 M 5,0 l 0,2 2,0 0,-2 Z')

    p2wo, p2sk, p2wi = p2.offset(0.2, two_sided=True)

    stroke1 = p2.stroke(0.25, cap='round')
    stroke2 = p2.stroke(0.25, join='miter')
    stroke3 = p2.stroke(0.25, join='round', cap='square')

    origin = 14 - 0j
    scale = 25

    output_stub(stroke1.translated(origin).scaled(scale).d(), {'stroke': 'gray', 'stroke-width': '1.5', 'fill': 'none'})
    output_stub(p2.translated(origin).scaled(scale).d(), {'stroke': 'black', 'stroke-width': '1.5', 'fill': 'none'})
    output_stub(stroke2.translated(origin + 6j).scaled(scale).d(), {'stroke': 'pink', 'stroke-width': '1.5', 'fill': 'none'})
    output_stub(p2.translated(origin + 6j).scaled(scale).d(), {'stroke': 'black', 'stroke-width': '1.5', 'fill': 'none'})
    output_stub(stroke3.translated(origin + 12j).scaled(scale).d(), {'stroke': 'teal', 'stroke-width': '1.5', 'fill': 'none'})
    output_stub(p2.translated(origin + 12j).scaled(scale).d(), {'stroke': 'black', 'stroke-width': '1.5', 'fill': 'none'})


# TEST 6 offsets from the README

if False:
    path1 = svgpathtools.parse_path("m 288,600 c -52,-28 -42,-61 0,-97 ")
    path2 = svgpathtools.parse_path("M 151,395 C 407,485 726.17662,160 634,339").translated(300)
    path3 = svgpathtools.parse_path("m 117,695 c 237,-7 -103,-146 457,0").translated(500 + 400j)
    paths = [path1, path2, path3]

    offset_distances = [10 * k for k in range(1, 51)]
    for path in paths:
        for distance in offset_distances:
            output_stub(path.offset(distance)[0].d(), {'stroke': 'red', 'stroke-width': '2', 'fill': 'none'})
