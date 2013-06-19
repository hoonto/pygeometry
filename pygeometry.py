# Profiling:
import cProfile

from ctypes import *
lib = cdll.LoadLibrary('./libgeometry.so')

lib.angle_box_2d.restype = None
lib.angle_contains_ray_2d.restype = c_int
lib.sphere_distance1.restype = c_double
lib.sphere_distance2.restype = c_double
lib.sphere_distance3.restype = c_double
lib.polygon_contains_point_2d.restype = c_bool
lib.polygon_contains_point_2d_2.restype = c_bool

class Geometry():

    def angle_box_2d ( dist, p1, p2, p3, p4, p5 ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        cp4 = (c_double * 2)(*p4)
        cp5 = (c_double * 2)(*p5)
        lib.angle_box_2d(c_double(dist), cp1, cp2, cp3, cp4, cp5);

    def angle_contains_ray_2d ( p1, p2, p3, p):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        cp = (c_double * 2)(*p)
        return lib.angle_contains_ray_2d(cp1, cp2, cp3, cp);

    def sphere_distance1(self, lat1, lon1, lat2, lon2, r):
        return lib.sphere_distance1(c_double(lat1), c_double(lon1), c_double(lat2), c_double(lon2), c_double(r)) 

    def sphere_distance2(self, lat1, lon1, lat2, lon2, r):
        return lib.sphere_distance2(c_double(lat1), c_double(lon1), c_double(lat2), c_double(lon2), c_double(r))

    def sphere_distance3(self, lat1, lon1, lat2, lon2, r):
        return lib.sphere_distance3(c_double(lat1), c_double(lon1), c_double(lat2), c_double(lon2), c_double(r))

    def polygon_contains_point_2d(self, v, p):
        n = len(v)
        p1 = (c_double * 2)(*p)
        pn = (c_double * len(v))(*v)
        return lib.polygon_contains_point_2d (c_int(n), pn, p1)

    def polygon_contains_point_2d_2(self, v, p):
        n = len(v)
        p1 = (c_double * 2)(*p)
        pn = (c_double * len(v))(*v)
        return lib.polygon_contains_point_2d_2 (c_int(n), pn, p1)


def test(numtests):
    g = Geometry()

    toradians = 0.0174532925;

    trange = range(0,numtests)

    for t in trange:
        g.sphere_distance1(39.7391500 * toradians, -104.9847000 * toradians, 34.0522300 * toradians, -118.2436800 * toradians, 6371)
        g.polygon_contains_point_2d([0,0,0,4,4,4,4,0], [1,1])


#test(1);

# 2000 calls to a sphere_distance1
# and 2000 calls to polygon_contains_point_2d
# = 4000 calls total
cProfile.run('test(2000)')      # ~0.062 seconds on 64bit fedora

# 20000 calls to a sphere_distance1
# and 20000 calls to polygon_contains_point_2d
# = 40000 calls total
cProfile.run('test(20000)')     # ~0.603 seconds on 64bit fedora

# 30000 calls to a sphere_distance1
# and 30000 calls to polygon_contains_point_2d
# = 60000 calls total
cProfile.run('test(30000)')     # ~0.905 seconds on 64bit fedora

# 40000 calls to a sphere_distance1
# and 40000 calls to polygon_contains_point_2d
# = 80000 calls total
cProfile.run('test(40000)')     # ~1.198 seconds on 64bit fedora

