
extern "C" {
    double py_sphere_distance1 ( double lat1, double lon1, double lat2, double lon2, double r ){
        std::cout << "lat1: " << lat1 << std::endl;
        std::cout << "lon1: " << lon1 << std::endl;
        std::cout << "lat2: " << lat2 << std::endl;
        std::cout << "lon2: " << lon2 << std::endl;
        std::cout << "r: " << r << std::endl;
        
        double result = sphere_distance1( lat1, lon1, lat2, lon2, r );
        std::cout << result << std::endl;
        return result;
    };
}
extern "C" {
    double py_sphere_distance2 ( double lat1, double lon1, double lat2, double lon2, double r ){
        double result = sphere_distance2( lat1, lon1, lat2, lon2, r );
        std::cout << result << std::endl;
        return result;
    };
}

extern "C" {
    double py_sphere_distance3 ( double lat1, double lon1, double lat2, double lon2, double r ){
        double result = sphere_distance3( lat1, lon1, lat2, lon2, r );
        std::cout << result << std::endl;
        return result;
    };
}

extern "C" {
    bool py_polygon_contains_point_2d ( int n, double v[], double p[2] ){
        return polygon_contains_point_2d( n, v, p );
    };
}

extern "C" {
    bool py_polygon_contains_point_2d_2 ( int n, double v[], double p[2] ){
        return polygon_contains_point_2d_2( n, v, p );
    };
}

