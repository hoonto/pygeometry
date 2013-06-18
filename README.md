pygeometry
==========

Python Geometry Library

Python bridge to this excellent library: http://people.sc.fsu.edu/~jburkardt/cpp_src/geometry/geometry.html

only provides the functions described in the geometryexterns.hpp file so far.

In answer to this stackoverflow question, and for fun :)

http://stackoverflow.com/questions/17022006/fast-python-gis-library-that-supports-great-cycle-distance-and-polygon

Build Instructions
==================

clone me, then:

``` sh
g++ -c -fPIC geometry.cpp -o geometry.o
g++ -shared -Wl,-soname,libgeometry.so -o libgeometry.so  geometry.o
```

Then you can do: 

``` sh
python geometry.py
```



