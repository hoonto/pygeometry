
extern "C" {
    double py_sphere_distance1 ( double lat1, double lon1, double lat2, double lon2, double r ){
        return sphere_distance1( lat1, lon1, lat2, lon2, r );
    };

    double py_sphere_distance2 ( double lat1, double lon1, double lat2, double lon2, double r ){
        return sphere_distance2( lat1, lon1, lat2, lon2, r );
    };

    double py_sphere_distance3 ( double lat1, double lon1, double lat2, double lon2, double r ){
        return sphere_distance3( lat1, lon1, lat2, lon2, r );
    };

    bool py_polygon_contains_point_2d ( int n, double v[], double p[2] ){
        return polygon_contains_point_2d( n, v, p );
    };

    bool py_polygon_contains_point_2d_2 ( int n, double v[], double p[2] ){
        return polygon_contains_point_2d_2( n, v, p );
    };

}

