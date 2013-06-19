
extern "C" {
    double py_sphere_distance1 ( double lat1, double lon1, double lat2, double lon2, double r ){
        double result = sphere_distance1( lat1, lon1, lat2, lon2, r );
        return result;
    };
    double py_sphere_distance2 ( double lat1, double lon1, double lat2, double lon2, double r ){
        double result = sphere_distance2( lat1, lon1, lat2, lon2, r );
        return result;
    };
    double py_sphere_distance3 ( double lat1, double lon1, double lat2, double lon2, double r ){
        double result = sphere_distance3( lat1, lon1, lat2, lon2, r );
        return result;
    };
    bool py_polygon_contains_point_2d ( int n, double v[], double p[2] ){
        bool result = polygon_contains_point_2d( n, v, p );
        return result;
    };
    bool py_polygon_contains_point_2d_2 ( int n, double v[], double p[2] ){
        bool result = polygon_contains_point_2d_2( n, v, p );
        return result;
    };
}

