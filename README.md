PyGEOMETRY
==========

Python Geometry Library

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



