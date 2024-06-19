import argparse

from cubicpts import tr_quadratic as quad

parser = argparse.ArgumentParser(
    description="Determine solutions of q(x,y)=xyz of height <= B, \
    where q(x,y) = (ax+1)(bx+1) + (cy+1)(dy+1)-1")
parser.add_argument("B", type=int, help="maximal height B")
parser.add_argument("a", type=int, help="parameter a", nargs='?', default=-2)
parser.add_argument("b", type=int, help="parameter b", nargs='?', default=3)
parser.add_argument("c", type=int, help="parameter c", nargs='?', default=-3)
parser.add_argument("d", type=int, help="parameter d", nargs='?', default=5)
args = parser.parse_args()

B = args.B
a = args.a
b = args.b
c = args.c
d = args.d

for (x, y, z) in quad.points(B, a, b, c, d):
    print(x, y, z)
