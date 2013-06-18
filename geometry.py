from ctypes import *
lib = cdll.LoadLibrary('./libgeometry.so')

lib.py_sphere_distance1.restype = c_double
lib.py_sphere_distance2.restype = c_double
lib.py_sphere_distance3.restype = c_double
lib.py_polygon_contains_point_2d.restype = c_bool
lib.py_polygon_contains_point_2d_2.restype = c_bool

class Geometry():

    def sphere_distance1(self, lat1, lon1, lat2, lon2, r):
        return lib.py_sphere_distance1(c_double(lat1), c_double(lon1), c_double(lat2), c_double(lon2), c_double(r)) 

    def sphere_distance2(self, lat1, lon1, lat2, lon2, r):
        return lib.py_sphere_distance2(c_double(lat1), c_double(lon1), c_double(lat2), c_double(lon2), c_double(r))

    def sphere_distance3(self, lat1, lon1, lat2, lon2, r):
        return lib.py_sphere_distance3(c_double(lat1), c_double(lon1), c_double(lat2), c_double(lon2), c_double(r))

    def polygon_contains_point_2d(self, v, p):
        n = len(v)
        p1 = (c_double * 2)(*p)
        pn = (c_double * len(v))(*v)
        return lib.py_polygon_contains_point_2d (c_int(n), pn, p1)

    def polygon_contains_point_2d_2(self, v, p):
        n = len(v)
        p1 = (c_double * 2)(*p)
        pn = (c_double * len(v))(*v)
        return lib.py_polygon_contains_point_2d_2 (c_int(n), pn, p1)


    def test(self):

        toradians = 0.0174532925;

        print "sphere_distance1:", self.sphere_distance1(39.7391500 * toradians, -104.9847000 * toradians, 34.0522300 * toradians, -118.2436800 * toradians, 6371)
        print "sphere_distance2:", self.sphere_distance2(39.7391500 * toradians, -104.9847000 * toradians, 34.0522300 * toradians, -118.2436800 * toradians, 6371)
        print "sphere_distance3:", self.sphere_distance3(39.7391500 * toradians, -104.9847000 * toradians, 34.0522300 * toradians, -118.2436800 * toradians, 6371)


        print "polygon_contains_point_2d", self.polygon_contains_point_2d([0,0,0,4,4,4,4,0], [-1,-1]);
        print "polygon_contains_point_2d:", self.polygon_contains_point_2d([0,0,0,4,4,4,4,0], [1,1]);
        print "polygon_contains_point_2d:", self.polygon_contains_point_2d([0,0,0,4,4,4,4,0], [3,3]);
        print "polygon_contains_point_2d:", self.polygon_contains_point_2d([0,0,0,4,4,4,4,0], [5,5]);

g = Geometry()
g.test()
