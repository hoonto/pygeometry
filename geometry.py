from ctypes import *
lib = cdll.LoadLibrary('./libgeometry.so')

lib.py_sphere_distance1.restype = c_double
lib.py_sphere_distance2.restype = c_double
lib.py_sphere_distance3.restype = c_double
lib.py_polygon_contains_point_2d.restype = c_bool
lib.py_polygon_contains_point_2d_2.restype = c_bool

class Geometry(object):

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


g = Geometry()
print g.sphere_distance1(40.643778999999990000, -73.781993000000000000, 33.945293800000000000, -118.385562699999980000, 3963.1676)
print g.sphere_distance2(40.643778999999990000, -73.781993000000000000, 33.945293800000000000, -118.385562699999980000, 3963.1676)
print g.sphere_distance3(40.643778999999990000, -73.781993000000000000, 33.945293800000000000, -118.385562699999980000, 3963.1676)
print g.polygon_contains_point_2d([0,0,0,4,4,4,4,0], [-1,-1]);
print g.polygon_contains_point_2d([0,0,0,4,4,4,4,0], [1,1]);
print g.polygon_contains_point_2d([0,0,0,4,4,4,4,0], [3,3]);
print g.polygon_contains_point_2d([0,0,0,4,4,4,4,0], [5,5]);
print g.polygon_contains_point_2d_2([0,0,0,4,4,4,4,0], [-1,-1]);
print g.polygon_contains_point_2d_2([0,0,0,4,4,4,4,0], [3,3]);
print g.polygon_contains_point_2d_2([0,0,0,4,4,4,4,0], [5,5]);
