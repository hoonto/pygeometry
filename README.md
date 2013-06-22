PyGEOMETRY
==========

Python Geometry Library

Currently contains the following functions:

* angle_box_2d
* angle_contains_ray_2d
* angle_deg_2d
* angle_half_2d
* angle_rad_2d
* angle_rad_3d
* angle_rad_nd
* angle_turn_2d
* anglei_deg_2d
* anglei_rad_2d
* annulus_area_2d
* annulus_sector_area_2d
* annulus_sector_centroid_2d
* ball_unit_sample_2d
* ball_unit_sample_3d
* ball_unit_sample_nd
* basis_map_3d
* box_01_contains_point_2d
* box_01_contains_point_nd
* box_contains_point_2d
* box_contains_point_nd
* polygon_contains_point_2d
* polygon_contains_point_2d_2
* sphere_distance1
* sphere_distance2
* sphere_distance3

Python bridge to this excellent library: 

http://people.sc.fsu.edu/~jburkardt/c_src/geometry/geometry.html

only provides the functions described in the geometryexterns.hpp file so far.

In answer to this stackoverflow question, and for fun :)

http://stackoverflow.com/questions/17022006/fast-python-gis-library-that-supports-great-cycle-distance-and-polygon

Build Instructions
==================

clone me, then:

``` sh
gcc -c -fPIC geometry.cpp -o geometry.o
gcc -shared -Wl,-soname,libgeometry.so -o libgeometry.so  geometry.o
```

Then you can do: 

``` sh
python geometry.py
```


Notes
=====

Switched to using the C version of GEOMETRY so that we don't have to write externs



