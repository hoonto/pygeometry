# Profiling:
import cProfile

from ctypes import *
lib = cdll.LoadLibrary('./libgeometry.so')


# Return types

lib.angle_box_2d.restype = None
lib.angle_contains_ray_2d.restype = c_bool
lib.angle_deg_2d.restype = c_double
lib.angle_half_2d.restype = c_void_p 
lib.angle_rad_2d.restype = c_double
lib.angle_rad_3d.restype = c_double
lib.angle_rad_nd.restype = c_double
lib.angle_turn_2d.restype =c_double
lib.anglei_deg_2d.restype = c_double
lib.anglei_rad_2d.restype = c_double
lib.annulus_area_2d.restype = c_double
lib.annulus_sector_area_2d.restype = c_double
lib.annulus_sector_centroid_2d.restype = c_void_p # MLM: double *
lib.ball_unit_sample_2d.restype = c_void_p # MLM: double *
lib.ball_unit_sample_3d.restype = c_void_p # MLM: double *
lib.ball_unit_sample_nd.restype = c_void_p # MLM; double *
lib.basis_map_3d.restype = c_void_p #double *
lib.box_01_contains_point_2d.restype = c_bool 
lib.box_01_contains_point_nd.restype = c_bool
lib.box_contains_point_2d.restype = c_int 
lib.box_contains_point_nd.restype = c_int 
lib.box_ray_int_2d.restype = None  
lib.box_segment_clip_2d.restype = c_int 
lib.circle_arc_point_near_2d.restype =  None
lib.circle_area_2d.restype = c_double
lib.circle_dia2imp_2d.restype = None
lib.circle_exp_contains_point_2d.restype = c_int
lib.circle_exp2imp_2d.restype = None
lib.circle_imp_contains_point_2d.restype = c_int 
lib.circle_imp_line_par_int_2d.restype = None 
lib.circle_imp_point_dist_2d.restype = c_double 
lib.circle_imp_point_dist_signed_2d.restype = c_double
lib.circle_imp_point_near_2d.restype = c_double
lib.circle_imp_points_2d.restype = c_void_p # MlM: double *
lib.circle_imp_points_3d.restype = c_void_p # MLM: double *
lib.circle_imp_points_arc_2d.restype = None
lib.circle_imp_print_2d.restype = None
lib.circle_imp_print_3d.restype = None
lib.circle_imp2exp_2d.restype = None  
lib.circle_llr2imp_2d.restype = c_void_p # MLM: double *
lib.circle_lune_area_2d.restype = c_double 
lib.circle_lune_centroid_2d.restype = c_void_p # MLM; double *
lib.circle_pppr2imp_3d.restype = None  

lib.circle_ppr2imp_2d.restype = c_void_p # MLM; double *
#lib.circle_sector_area_2d.restype = c_double 
#lib.circle_sector_centroid_2d.restype = c_void_p #double *
#lib.circle_sector_contains_point_2d.restype = c_int 
#lib.circle_sector_print_2d.restype = None  
#lib.circle_triangle_area_2d.restype = c_double
#lib.circle_triple_angles_2d.restype = None 
#lib.circles_imp_int_2d.restype = None 
#lib.cone_area_3d.restype = c_double 
#lib.cone_centroid_3d.restype = c_void_p # MLM: double *
#lib.cone_volume_3d.restype = c_double
#lib.conv3d.restype = None
#lib.cos_deg.restype = c_double
#lib.cot_deg.restype = c_double
#lib.cot_rad.restype = c_double
#lib.csc_deg.restype = c_double
#lib.cube_shape_3d.restype = None
#lib.cube_size_3d.restype = None 
#lib.cylinder_point_dist_3d.restype = c_double
#lib.cylinder_point_dist_signed_3d.resytpe = c_double
#lib.cylinder_point_inside_3d.restype = c_int
#lib.cylinder_point_near_3d.restype = c_void_p # MLM: double *
#lib.cylinder_sample_3d.restype = c_void_p # MLM: double *
#lib.cylinder_volume_3d.restype = c_double
#lib.degrees_to_radians.restype = c_double
#lib.dge_det.restype = c_double
#lib.dge_fa.resytpe = c_int
#lib.dge_sl.restype = None
#lib.direction_pert_3d.restype = c_void_p # MLM: double *
#lib.direction_uniform_2d.restype = c_void_p # MLM: double *
#lib.direction_uniform_3d.restype = c_void_p # MLM: double *
#lib.direction_uniform_nd.restype = c_void_p # MLM: double *
#lib.disk_point_dist_3d.restype = c_double 
#lib.dms_to_radians.restype = c_double
#lib.dodec_shape_3d.restype = None 
#lib.dodec_size_3d.restype = None 
#lib.dual_shape_3d.restype = None 
#lib.dual_size_3d.restype = None 
#lib.ellipse_area_2d.restype = c_double
#lib.ellipse_point_dist_2d.restype = c_double
#lib.ellipse_point_near_2d.restype = c_void_p # MLM: double *
#lib.ellipse_points_2d.restype = None 
#lib.ellipse_points_arc_2d.restype = None 
#lib.enorm0_nd.restype = c_double
#lib.get_seed.restype = c_int
#lib.glob2loc_3d.restype = None 
#lib.halfplane_contains_point_2d.restype = c_int
#lib.halfspace_imp_triangle_int_3d.restype = c_int
#lib.halfspace_norm_triangle_int_3d.restype = c_int
#lib.halfspace_triangle_int_3d.restype = c_int
#lib.haversine.restype = c_double
#lib.helix_shape_3d.restype = None
#lib.hexagon_area_2d.restype = c_double
#lib.hexagon_contains_point_2d.restype = c_int
#lib.hexagon_shape_2d.restype = None
#lib.hexagon_unit_area_2d.restype = c_double
#lib.hexagon_vertices_2d.restype = None
#lib.i4_dedekind_factor.restype = c_double
#lib.i4_dedekind_sum.restype = c_double
#lib.i4_factorial2.restype = c_int
#lib.i4_gcd.restype = c_int
#lib.i4_lcm.restype = c_int
#lib.i4_max.restype = c_int
#lib.i4_min.restype = c_int
#lib.i4_modp.restype = c_int
#lib.i4_sign.restype = c_int
#lib.i4_swap.restype = None
#lib.i4_uniform.restype = c_int
#lib.i4_wrap.restype = c_int
#lib.i4col_compare.restype = c_int
#lib.i4col_find_item.restype = None
#lib.i4col_find_pair_wrap.restype = None 
#lib.i4col_sort_a.restype = None
#lib.i4col_sorted_unique_count.restype = c_int
#lib.i4col_swap.restype = None
#lib.i4mat_print.restype = None
#lib.i4mat_transpose_print.restype = None 
#lib.i4mat_transpose_print_some.restype = None
#lib.i4row_compare.restype = c_int
#lib.i4row_sort_a.restype = None
#lib.i4row_swap.restype = None
#lib.i4vec_copy.restype = None
#lib.i4vec_heap_d.restype = None
#lib.i4vec_indicator_new.restype = c_int_p
#lib.i4vec_lcm.restype = c_int
#lib.i4vec_print.restype = None
#lib.i4vec_product.restype = c_int
#lib.i4vec_reverse.restype = None
#lib.i4vec_sort_heap_a.restype = None
#lib.i4vec_sorted_unique.restype = c_int
#lib.i4vec_uniform_new.restype = c_int_p
#lib.i4vec_zero.restype = None
#lib.i4vec2_compare.restype = c_int
#lib.i4vec2_sort_a.restype = None
#lib.i4vec2_sorted_unique.restype = None
#lib.icos_shape.restype = None
#lib.icos_size.restype = None
#lib.line_exp_is_degenerate_nd.restype = c_int
#lib.line_exp_normal_2d.restype = c_void_p # MLM: double *;
#lib.line_exp_perp_2d.restype = c_void_p # double *
#lib.line_exp_point_dist_2d.restype = c_double
#lib.line_exp_point_dist_3d.restype = c_double
#lib.line_exp_point_dist_signed_2d.restype = c_double
#lib.line_exp_point_near_2d.restype = None
#lib.line_exp_point_near_3d.restype = None
#lib.line_exp2imp_2d.restype = None
#lib.line_exp2par_2d.restype = None
#lib.line_exp2par_3d.restype = None
#lib.line_imp_is_degenerate_2d.restype = c_int
#lib.line_imp_point_dist_2d.restype = c_double
#lib.line_imp_point_dist_signed_2d.restype = c_double
#lib.line_imp2exp_2d.restype = None
#lib.line_imp2par_2d.restype = None
#lib.line_par_point_dist_2d.restype = c_double
#lib.line_par_point_dist_3d.restype = c_double
#lib.line_par_point_near_2d.restype = c_void_p # MLM: double *
#lib.line_par_point_near_3d.restype = c_void_p # MLM: double *
#lib.line_par2exp_2d.restype = None
#lib.line_par2exp_3d.restype = None
#lib.line_par2imp_2d.restype = None
#lib.lines_exp_angle_3d.restype = c_double
#lib.lines_exp_angle_nd.restype = c_double
#lib.lines_exp_dist_3d.restype = c_double
#lib.lines_exp_dist_3d_2.restype = c_double
#lib.lines_exp_equal_2d.restype = c_int
#lib.lines_exp_int_2d.restype = None
#lib.lines_exp_near_3d.restype = None
#lib.lines_exp_parallel_2d.restype = c_int
#lib.lines_exp_parallel_3d.restype = c_int
#lib.lines_imp_angle_2d.restype = c_double
#lib.lines_imp_dist_2d.restype = c_double
#lib.lines_imp_int_2d.restype = None
#lib.lines_par_angle_2d.restype = c_double
#lib.lines_par_angle_3d.restype = c_double
#lib.lines_par_dist_3d.restype = c_double
#lib.lines_par_int_2d.restype = None
#lib.oc2glob_3d.restype = None 
#lib.lrline.restype = c_int
#lib.lvec_print.restype = None
#lib.minabs.restype = None
#lib.minquad.restype = c_int
#lib.octahedron_shape_3d.restype = None 
#lib.octahedron_size_3d.restype = None
#lib.parabola_ex.restype = c_int
#lib.parabola_ex2.restype = c_int
#lib.parallelogram_area_2d.restype = c_double
#lib.parallelogram_area_3d.restype = c_double
#lib.parallelogram_contains_point_2d.restype = c_int
#lib.parallelogram_contains_point_3d.restype = c_int
#lib.parallelogram_point_dist_3d.restype = c_double
#lib.parallelepiped_contains_point_3d.restype = c_int
#lib.parallelepiped_point_dist_3d.restype = c_double
#lib.perm_check.restype = c_int
#lib.perm_inv.restype = None
#lib.plane_exp_grid_3d.restype = None
#lib.plane_exp_point_dist_3d.restype = c_double
#lib.plane_exp_normal_3d.restype = None
#lib.plane_exp_pro2.restype = None
#lib.plane_exp_pro3.restype = None
#lib.plane_exp_project_3d.restype = None
#lib.plane_exp2imp_3d.restype = None
#lib.plane_exp2normal_3d.restype = None
#lib.plane_imp_is_degenerate_3d.restype = c_int
#lib.plane_imp_line_par_int_3d.restype = c_int
#lib.plane_imp_point_dist_3d.restype = c_double
#lib.plane_imp_point_dist_signed_3d.restype = c_double
#lib.plane_imp_point_near_3d.restype = None
#lib.plane_imp_segment_near_3d.restype = None
#lib.plane_imp_triangle_int_3d.restype = None
#lib.plane_imp_triangle_int_add_3d.restype = None
#lib.plane_imp_triangle_near_3d.restype = c_int
#lib.plane_imp2exp_3d.restype = None
#lib.plane_imp2normal_3d.restype = None 
#lib.plane_normal_basis_3d.restype = None
#lib.plane_normal_line_exp_int_3d.restype = c_int
#lib.plane_normal_qr_to_xyz.restype = c_void_p # MLM: double *
#lib.plane_normal_tetrahedron_intersect.restype = None
#lib.plane_normal_triangle_int_3d.restype = c_int
#lib.plane_normal_uniform_3d.restype = None
#lib.plane_normal_uniform_nd.restype = None
#lib.plane_normal_xyz_to_qr.restype = c_void_p # MLM: double *
#lib.plane_normal2exp_3d.restype = None
#lib.plane_normal2imp_3d.restype = None
#lib.planes_imp_angle_3d.restype = c_double
#lib.points_avoid_point_naive_2d.restype = c_int
#lib.points_bisect_line_imp_2d.restype = None
#lib.points_bisect_line_par_2d.restype = None
#lib.points_centroid_2d.restype = c_int
#lib.points_colin_2d.restype = c_double
#lib.points_colin_3d.restype = c_double
###lib.points_dist_2d.restype = c_double
#lib.points_dist_3d.restype = c_double
#lib.points_dist_nd.restype = c_double
#lib.points_hull_2d.restype = None 
#lib.points_plot.restype = None
#lib.points_point_near_naive_2d.restype = c_int
#lib.points_point_near_naive_3d.restype = c_int
#lib.points_point_near_naive_nd.restype = c_int
#lib.points_points_near_naive_2d.restype = c_int_p
#lib.points_points_near_naive_3d.restype = c_int_p
#lib.polar_to_xy.restype = None
#lib.polygon_1_2d.restype = c_double
#lib.polygon_angles_2d.restype = c_void_p # MLM: double *
#lib.polygon_area_2d.restype = c_double
#lib.polygon_area_2d_2.restype = c_double
#lib.polygon_area_3d.restype = c_double
#lib.polygon_area_3d_2.restype = c_double
#lib.polygon_centroid_2d.restype = c_void_p # MLM: double *
#lib.polygon_centroid_2d_2.restype = c_void_p # MLM: double *
#lib.polygon_centroid_3d.restype = c_void_p # MLM: double *
lib.polygon_contains_point_2d.restype = c_bool
lib.polygon_contains_point_2d_2.restype = c_bool
#lib.polygon_diameter_2d.restype = c_double
#lib.polygon_expand_2d.restype = c_void_p # MLM: double *
#lib.polygon_inrad_data_2d.restype = None
#lib.polygon_is_convex.restype = c_int
#lib.polygon_lattice_area_2d.restype = c_double
#lib.polygon_normal_3d.restype = c_void_p # MLM: double *
#lib.polygon_outrad_data_2d.restype = None
#lib.polygon_side_data_2d.restype = None
#lib.polygon_solid_angle_3d.restype = c_double
#lib.polygon_x_2d.restype = c_double
#lib.polygon_y_2d.restype = c_double
#lib.polygon_xx_2d.restype = c_double
#lib.polygon_xy_2d.restype = c_double
#lib.polygon_yy_2d.restype = c_double
#lib.polyhedron_area_3d.restype = c_double
#lib.polyhedron_centroid_3d.restype = c_void_p # MLM: double *
#lib.polyhedron_contains_point_3d.restype = c_int
#lib.polyhedron_volume_3d.restype = c_double
#lib.polyhedron_volume_3d_2.restype = c_double
#lib.polyline_arclength_nd.restype = c_void_p # MLM: double *
#lib.polyline_index_point_nd.restype = c_void_p # MLM: double *
#lib.polyline_length_nd.restype = c_double
#lib.polyline_points_nd.restype = c_void_p #MLM: double *
#lib.polyloop_arclength_nd.restype = c_void_p # MLM: double *
#lib.polyloop_points_nd.restype = c_void_p # MLM: double *
#lib.provec.restype = None
#lib.pyramid_volume_3d.restype = c_double
#lib.quad_area_2d.restype = c_double
#lib.quad_area2_2d.restype = c_double
#lib.quad_area_3d.restype = c_double
#lib.quad_contains_point_2d.restype = c_int
#lib.quad_convex_random.restype = None
#lib.quad_point_dist_2d.restype = c_double
#lib.quad_point_dist_signed_2d.restype = c_double
#lib.quad_point_near_2d.restype = None
#lib.quat_conj.restype = c_void_p # MLM: double *
#lib.quat_inv.restype = c_void_p # MLM: double *
#lib.quat_mul.restype = c_void_p # MLM: double *
#lib.quat_norm.restype = c_double
#lib.r4_abs.restype = c_float
#lib.r4_nint.restype = c_int
#lib.r8_abs.restype = c_double
#lib.r8_acos.restype = c_double
#lib.r8_asin.restype = c_double
#lib.r8_atan.restype = c_double
#lib.r8_epsilon.restype = c_double
#lib.r8_huge.restype = c_double;
#lib.r8_max.restype = c_double
#lib.r8_min.restype = c_double
#lib.r8_modp.restype = c_double
#lib.r8_nint.restype = c_int
#lib.r8_normal_01.restype = c_double
#lib.r8_pi.restype = c_double
#lib.r8_sign.restype = c_double
#lib.r8_sign_opposite_strict.restype = c_int
#lib.r8_swap.restype = None 
#lib.r8_uniform.restype =  c_double
#lib.r8_uniform_01.restype = c_double
#lib.r82vec_part_quick_a.restype = None 
#lib.r82vec_permute.restype = None
#lib.r82vec_print.restype = None
#lib.r82vec_sort_heap_index_a.restype = c_int_p 
#lib.r82vec_sort_quick_a.restype = None
#lib.r8mat_copy.restype = None
#lib.r8mat_det_2d.restype = c_double
#lib.r8mat_det_3d.restype =  c_double
#lib.r8mat_det_4d.restype = c_double
#lib.r8mat_det_5d.restype = c_double
#lib.r8mat_inverse_2d.restype = c_void_p # MLM: double *
#lib.r8mat_inverse_3d.restype = c_void_p # MLM: double *
#lib.r8mat_mv.restype = c_void_p # MLM: double *
#lib.r8mat_print.restype = None
#lib.r8mat_print_some.restype = None
#lib.r8mat_solve.restype = c_int
#lib.r8mat_solve_2d.restype = c_void_p # MLM: double *
#lib.r8mat_transpose_print.restype = None
#lib.r8mat_transpose_print_some.restype = None
#lib.r8mat_uniform_new.restype = c_void_p # MLM: double *
#lib.r8mat_uniform_01.restype = None
#lib.r8mat_uniform_01_new.restype = c_void_p # MLM: double *
#lib.r8vec_angle_3d.restype = c_double
#lib.r8vec_any_normal.restype = c_void_p # MLM; double *
#lib.r8vec_bracket.restype = None
#lib.r8vec_copy.restype = None
#lib.r8vec_cross_product_2d.restype = c_double
#lib.r8vec_cross_product_2d_affine.restype = c_double
#lib.r8vec_cross_product_3d.restype = c_void_p # MLM: double *
#lib.r8vec_cross_product_3d_affine.restype = c_void_p # MLM: double *
#lib.r8vec_distance.restype = c_double
#lib.r8vec_dot_product.restype = c_double
#lib.r8vec_dot_product_affine.restype = c_double
#lib.r8vec_eq.restype = c_int
#lib.r8vec_gt.restype = c_int
#lib.r8vec_lt.restype = c_int
#lib.r8vec_max.restype = c_double
#lib.r8vec_mean.restype = c_double
#lib.r8vec_min.restype = c_double
#lib.r8vec_negative_strict.restype = c_int
#lib.r8vec_norm.restype = c_double
#lib.r8vec_norm_affine.restype = c_double
#lib.r8vec_normal_01_new.restype = c_void_p # MLM; double *
#lib.r8vec_normsq.restype = c_double
#lib.r8vec_normsq_affine.restype = c_double
#lib.r8vec_positive_strict.restype = c_int
#lib.r8vec_print.restype = None
#lib.r8vec_print_2d.restype = None
#lib.r8vec_print_3d.restype = None
#lib.r8vec_scalar_triple_product.restype = c_double
#lib.r8vec_swap.restype = None
#lib.r8vec_uniform_new.restype = _void_p # MLM: double *;
#lib.r8vec_uniform_01_new.restype = c_void_p # MLM; double *
#lib.r8vec_uniform_unit_new.restype = c_void_p # MLM: double *
#lib.r8vec_variance.restype = c_double
#lib.r8vec_zero.restype = None
#lib.radec_distance_3d.restype = c_double
#lib.radec_to_xyz.restype = c_void_p # MLM: double *
#lib.radians_to_degrees.restype = c_double
#lib.radians_to_dms.restype = None
#lib.random_initialize.restype = c_ulong
#lib.rotation_axis_vector_3d.restype = None
#lib.rotation_axis2mat_3d.restype = None
#lib.rotation_axis2quat_3d.restype = None
#lib.rotation_mat_vector_3d.restype = None
#lib.rotation_mat2axis_3d.restype = None
#lib.rotation_mat2quat_3d.restype = None
#lib.rotation_quat_vector_3d.restype = None
#lib.rotation_quat2axis_3d.restype = None
#lib.rotation_quat2mat_3d.restype = None
#lib.rtp_to_xyz.restype = None
#lib.s_len_trim.restype = c_int
#lib.sec_deg.restype = c_double
#lib.segment_contains_point_1d.restype = None
#lib.segment_contains_point_2d.restype = None
#lib.segment_point_coords_2d.restype = None
#lib.segment_point_coords_3d.restype = None
#lib.segment_point_dist_2d.restype = c_double
#lib.segment_point_dist_3d.restype = c_double
#lib.segment_point_near_2d.restype = None
#lib.segment_point_near_3d.restype = None
#lib.segments_curvature_2d.restype = c_double
#lib.segments_dist_2d.restype = c_double
#lib.segments_dist_3d.restype = c_double 
#lib.segments_dist_3d_old.restype = c_double 
#lib.segments_int_1d.restype = c_double 
#lib.segments_int_2d.restype = None 
#lib.shape_point_dist_2d.restype = c_double 
#lib.shape_point_near_2d.restype = None 
#lib.shape_print_3d.restype = None 
#lib.shape_ray_int_2d.restype = None 
#lib.simplex_lattice_layer_point_next.restype = None 
#lib.simplex_lattice_point_next.restype = None 
#lib.simplex_unit_lattice_point_num_nd.restype = c_int 
#lib.simplex_unit_volume_nd.restype = c_double 
#lib.simplex_volume_nd.restype = c_double 
#lib.sin_deg.restype = c_double 
#lib.sin_power_int.restype = c_double 
#lib.soccer_shape_3d.restype = None 
#lib.soccer_size_3d.restype = None 
#lib.sort_heap_external.restype = None 
#lib.sphere_cap_area_2d.restype = c_double 
#lib.sphere_cap_area_3d.restype = c_double 
#lib.sphere_cap_area_nd.restype = c_double 
#lib.sphere_cap_volume_2d.restype = c_double 
#lib.sphere_cap_volume_3d.restype = c_double 
#lib.sphere_cap_volume_nd.restype = c_double 
#lib.sphere_dia2imp_3d.restype = None 
lib.sphere_distance1.restype = c_double 
lib.sphere_distance2.restype = c_double 
lib.sphere_distance3.restype = c_double 
#lib.sphere_exp_contains_point_3d.restype = c_int 
#lib.sphere_exp_point_near_3d.restype = None 
#lib.sphere_exp2imp_3d.restype = None 
#lib.sphere_exp2imp_nd.restype = None 
#lib.sphere_imp_contains_point_3d.restype = c_int 
#lib.sphere_imp_grid_icos_size.restype = None 
#lib.sphere_imp_gridfaces_3d.restype = None 
#lib.sphere_imp_gridlines_3d.restype = None 
#lib.sphere_imp_gridpoints_3d.restype = None 
#lib.sphere_imp_gridpoints_icos1.restype = None 
#lib.sphere_imp_gridpoints_icos2.restype = None 
#lib.sphere_imp_line_project_3d.restype = c_int 
#lib.sphere_imp_local2xyz_3d.restype = None 
#lib.sphere_imp_point_near_3d.restype = None 
#lib.sphere_imp_point_project_3d.restype = None 
#lib.sphere_imp_spiralpoints_3d.restype = None 
#lib.sphere_imp_volume_3d.restype = c_double 
#lib.sphere_imp_volume_nd.restype = c_double 
#lib.sphere_imp_zone_area_3d.restype = c_double 
#lib.sphere_imp_zone_volume_3d.restype = c_double 
#lib.sphere_imp2exp_3d.restype = None 
#lib.sphere_k.restype = c_double 
#lib.sphere_triangle_angles_to_area.restype = c_double 
#lib.sphere_triangle_contains_point.restype = c_double 
#lib.sphere_triangle_sides_to_angles.restype = None 
#lib.sphere_triangle_vertices_to_angles.restype = None 
#lib.sphere_triangle_vertices_to_area.restype  = c_double 
#lib.sphere_triangle_vertices_to_centroid.restype = None 
#lib.sphere_triangle_vertices_to_orientation.restype = c_int 
#lib.sphere_triangle_vertices_to_sides.restype = None 
#lib.sphere_unit_area_nd.restype = c_double 
#lib.sphere_unit_area_values.restype = None 
#lib.sphere_unit_sample_2d.restype = c_void_p # MLM: double *
#lib.sphere_unit_sample_3d.restype = c_void_p # MLM: double *
#lib.sphere_unit_sample_3d_2.restype = c_void_p # MLM; double *
#lib.sphere_unit_sample_nd.restype = c_void_p # MLM: double *
#lib.sphere_unit_sample_nd_2.restype = c_void_p # MLM: double *
#lib.sphere_unit_sample_nd_3.restype = c_void_p # MLM: double *
#lib.sphere_unit_volume_nd.restype = c_double 
#lib.sphere_unit_volume_values.restype = None 
#lib.sphere01_distance_xyz.restype = c_double 
#lib.sphere01_polygon_area.restype = c_double 
#lib.sphere01_triangle_angles_to_area.restype = c_double 
#lib.sphere01_triangle_sides_to_angles.restype = None 
#lib.sphere01_triangle_vertices_to_angles.restype = None 
#lib.sphere01_triangle_vertices_to_area.restype = c_double 
#lib.sphere01_triangle_vertices_to_midpoints.restype = None 
#lib.sphere01_triangle_vertices_to_centroid.restype = c_void_p # MLM; double *
#lib.sphere01_triangle_vertices_to_sides.restype = None 
#lib.string_2d.restype = None 
#lib.super_ellipse_points_2d.restype = None 
#lib.tan_deg.restype = c_double 
#lib.tetrahedron_barycentric_3d.restype = c_void_p # MLM; double *
#lib.tetrahedron_centroid_3d.restype = c_void_p # MLM; double *
#lib.tetrahedron_circumsphere_3d.restype = None 
#lib.tetrahedron_contains_point_3d.restype = c_int 
#lib.tetrahedron_dihedral_angles_3d.restype = c_void_p # MLM; double *
#lib.tetrahedron_edge_length_3d.restype = c_void_p # MLM; double *
#lib.tetrahedron_face_angles_3d.restype = None 
#lib.tetrahedron_face_areas_3d.restype = None 
#lib.tetrahedron_insphere_3d.restype = None 
#lib.tetrahedron_lattice_layer_point_next.restype = None 
#lib.tetrahedron_lattice_point_next.restype = None 
#lib.tetrahedron_quality1_3d.restype = c_double 
#lib.tetrahedron_quality2_3d.restype = c_double 
#lib.tetrahedron_quality3_3d.restype = c_double 
#lib.tetrahedron_quality4_3d.restype = c_double 
#lib.tetrahedron_rhombic_shape_3d.restype = None 
#lib.tetrahedron_rhombic_size_3d.restype = None 
#lib.tetrahedron_sample_3d.restype = None 
#lib.tetrahedron_shape_3d.restype = None 
#lib.tetrahedron_size_3d.restype = None 
#lib.tetrahedron_solid_angles_3d.restype = c_void_p # MLM; double *
#lib.tetrahedron_unit_lattice_point_num_3d.restype = c_int 
#lib.tetrahedron_volume_3d.restype = c_double 
#lib.theta2_adjust.restype = None 
#lib.theta3_adjust.restype = None 
#lib.timestamp.restype = None 
#lib.tmat_init.restype = None 
#lib.tmat_mxm.restype = None 
#lib.tmat_mxp.restype = None 
#lib.tmat_mxp2.restype = None 
#lib.tmat_mxv.restype = None 
#lib.tmat_rot_axis.restype = None 
#lib.tmat_rot_vector.restype = None 
#lib.tmat_scale.restype = None 
#lib.tmat_shear.restype = None 
#lib.tmat_trans.restype = None 
#lib.torus_volume_3d.restype = c_double 
#lib.tp_to_xyz.restype = c_void_p # MLM; double *;
#lib.triangle_angles_2d.restype = None 
#lib.triangle_angles_2d_new.restype = c_void_p # MLM: double *
#lib.triangle_angles_3d.restype = None 
#lib.triangle_angles_3d_new.restype = c_void_p # MLM: double *
#lib.triangle_area_2d.restype = c_double 
#lib.triangle_area_3d.restype = c_double 
#lib.triangle_area_3d_2.restype = c_double 
#lib.triangle_area_3d_3.restype = c_double 
#lib.triangle_area_heron.restype = c_double 
#lib.triangle_area_vector_3d.restype = c_void_p # MLM: double *
#lib.triangle_barycentric_2d.restype = c_void_p # MLM; double *
#lib.triangle_centroid_2d.restype = c_void_p # MLM: double *
#lib.triangle_centroid_3d.restype = c_void_p # MLM: double *
#lib.triangle_circumcenter_2d.restype = c_void_p # MLM: double *
#lib.triangle_circumcenter_2d_2.restype = c_void_p # MLM; double *
#lib.triangle_circumcenter.restype = c_void_p # MLM: double *
#lib.triangle_circumcircle_2d.restype = None
#lib.triangle_circumcircle_2d_2.restype = None 
#lib.triangle_circumradius_2d.restype = c_double 
#lib.triangle_contains_line_exp_3d.restype = None 
#lib.triangle_contains_line_par_3d.restype = None 
#lib.triangle_contains_point_2d_1.restype = c_int 
#lib.triangle_contains_point_2d_2.restype = c_int 
#lib.triangle_contains_point_2d_3.restype = c_int 
#lib.triangle_diameter_2d.restype = c_double 
#lib.triangle_edge_length_2d.restype = c_void_p # MLM: double *
#lib.triangle_gridpoints_2d.restype = None 
##lib.triangle_incenter_2d.restype = None 
#lib.triangle_incircle_2d.restype = None 
#lib.triangle_inradius_2d.restype = c_double 
#lib.triangle_is_degenerate_nd.restype = c_int 
#lib.triangle_lattice_layer_point_next.restype = None 
#lib.triangle_lattice_point_next.restype = None 
#lib.triangle_line_imp_int_2d.restype = None 
#lib.triangle_orientation_2d.restype = c_int 
#lib.triangle_orthocenter_2d.restype = None
#lib.triangle_point_dist_2d.restype = c_double 
#lib.triangle_point_dist_3d.restype = _double ;
#lib.triangle_point_dist_signed_2d.restype = c_double 
#lib.triangle_point_near_2d.restype = None 
#lib.triangle_quality_2d.restype = c_double 
#lib.triangle_right_lattice_point_num_2d.restype = c_int 
#lib.triangle_sample.restype = None 
#lib.triangle_unit_lattice_point_num_2d.restype = c_int 
#lib.triangle_xsi_to_xy_2d.restype = None 
#lib.triangle_xy_to_xsi_2d.restype = None
#lib.truncated_octahedron_shape_3d.restype = None 
#lib.truncated_octahedron_size_3d.restype = None 
#lib.tube_2d.restype = None 
#lib.tuple_next2.restype = None 
#lib.vector_directions_nd.restype = None 
#lib.vector_rotate_2d.restype = None 
#lib.vector_rotate_3d.restype = None 
#lib.vector_rotate_base_2d.restype = None 
#lib.vector_separation_2d.restype = c_double 
#lib.vector_separation_3d.restype = c_double 
#lib.vector_separation_nd.restype = c_double 
#lib.vector_unit_nd.restype = None 
#lib.voxels_dist_l1_3d.restype = c_int 
#lib.voxels_dist_l1_nd.restype = c_int 
#lib.voxels_line_3d.restype = None 
#lib.voxels_region_3d.restype = None 
#lib.voxels_step_3d.restype = None 
#lib.xy_to_polar.restype = None 
#lib.xyz_to_radec.restype = None 
#lib.xyz_to_rtp.restype = None 
#lib.xyz_to_tp.restype = None 

# Geometry

class Geometry():

    # void angle_box_2d ( double dist, double p1[2], double p2[2], double p3[2], double p4[2], double p5[2] );
    def angle_box_2d ( self, dist, p1, p2, p3, p4, p5 ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        cp4 = (c_double * 2)(*p4)
        cp5 = (c_double * 2)(*p5)
        pcp4 = pointer(cp4);
        pcp5 = pointer(cp5);
        # MLM: cp4 and cp5 are returned
        lib.angle_box_2d(c_double(dist), cp1, cp2, cp3, pcp4, pcp5)
        p4.append(pcp4.contents[0]);
        p4.append(pcp4.contents[1]);
        p5.append(pcp5.contents[0]);
        p5.append(pcp5.contents[1]);

    # int angle_contains_ray_2d ( double p1[2], double p2[2], double p3[2], double p[2] );
    def angle_contains_ray_2d ( self, p1, p2, p3, p):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        cp = (c_double * 2)(*p)
        return lib.angle_contains_ray_2d(cp1, cp2, cp3, cp)

    # double angle_deg_2d ( double p1[2], double p2[2], double p3[2] );
    def angle_deg_2d (self,  p1, p2, p3 ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        return lib.angle_deg_2d ( cp1, cp2, cp3 )

    # double *angle_half_2d ( double p1[2], double p2[2], double p3[2] );
    def angle_half_2d (self,  p1, p2, p3 ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        pcp4 = pointer((c_double * 2)(lib.angle_half_2d( cp1, cp2, cp3 )))
        cp4 = []
        cp4.append(pcp4.contents[0])
        cp4.append(pcp4.contents[1])
        return cp4 

    # double angle_rad_2d ( double p1[2], double p2[2], double p3[2] );
    def angle_rad_2d (self,  p1, p2, p3 ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        return lib.angle_rad_2d( cp1, cp2, cp3 )

    # double angle_rad_3d ( double p1[3], double p2[3], double p3[3] );
    def angle_rad_3d (self,  p1, p2, p3 ):
        cp1 = (c_double * 3)(*p1)
        cp2 = (c_double * 3)(*p2)
        cp3 = (c_double * 3)(*p3)
        return lib.angle_rad_3d( cp1, cp2, cp3 )

    # double angle_rad_nd ( int dim_num, double vec1[], double vec2[] );
    def angle_rad_nd (self,  dim_num, vec1, vec2 ):
        cvec1 = (c_double * len(vec1))(*vec1)
        cvec2 = (c_double * len(vec2))(*vec2)
        return lib.angle_rad_nd(c_int(dim_num), cvec1, cvec2)

    # double angle_turn_2d ( double p1[2], double p2[2], double p3[2] );
    def angle_turn_2d (self,  p1, p2, p3 ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        return lib.angle_turn_2d( cp1, cp2, cp3 )

    # double anglei_deg_2d ( double p1[2], double p2[2], double p3[2] );
    def anglei_deg_2d (self,  p1, p2, p3 ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        return lib.anglei_deg_2d( cp1, cp2, cp3 )

    # double anglei_rad_2d ( double p1[2], double p2[2], double p3[2] );
    def anglei_rad_2d (self,  p1, p2, p3 ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        return lib.anglei_rad_2d( cp1, cp2, cp3 )

    #double annulus_area_2d ( double r1, double r2 );
    def annulus_area_2d (self,  r1, r2 ):
        return lib.annulus_area_2d(c_double(r1),c_double(r2))

    # double annulus_sector_area_2d ( double r1, double r2, double theta1, double theta2 );
    def annulus_sector_area_2d (self,  r1, r2, theta1, theta2 ):
        return lib.annulus_sector_area_2d(c_double(r1),c_double(r2),c_double(theta1),c_double(theta2))

    # double *annulus_sector_centroid_2d ( double pc[2], double r1, double r2, double theta1, double theta2 );
    def annulus_sector_centroid_2d (self,  pc, r1, r2, theta1, theta2 ):
        cpc = (c_double * 2)(*pc)
        pcp4 = pointer((c_double * 2)(lib.annulus_sector_centroid_2d(cpc, c_double(r1), c_double(r2), c_double(theta1), c_double(theta2))))
        cp4 = []
        cp4.append(pcp4.contents[0])
        cp4.append(pcp4.contents[1])
        return cp4

    # double *ball_unit_sample_2d ( int *seed );
    def ball_unit_sample_2d (self, seed ):
        pseed = pointer(c_int(seed))
        pcp4 = pointer((c_double * 2)(lib.ball_unit_sample_2d(pseed)))
        cp4 = []
        cp4.append(pcp4.contents[0])
        cp4.append(pcp4.contents[1])
        return cp4

    # double *ball_unit_sample_3d ( int *seed );
    def ball_unit_sample_3d (self, seed ):
        pseed = pointer(c_int(seed))
        pcp4 = pointer((c_double * 3)(lib.ball_unit_sample_3d(pseed)))
        cp4 = []
        cp4.append(pcp4.contents[0])
        cp4.append(pcp4.contents[1])
        cp4.append(pcp4.contents[2])
        return cp4

    # double *ball_unit_sample_nd ( int dim_num, int *seed );
    def ball_unit_sample_nd (self, dim_num, seed ):
        pseed = pointer(c_int(seed))
        tpcp4 = lib.ball_unit_sample_nd(dim_num,pseed)
        cp4 = [dim_num]
        if tpcp4 != None:
            pcp4 = pointer((c_double * dim_num)(tpcp4))
            for i in xrange(dim_num):
                cp4.append(pcp4.contents[i])
        return cp4

    # double *basis_map_3d ( double u[3*3], double v[3*3] );
    def basis_map_3d (self, u, v ):
        dim_num = 9
        cu = (c_double * dim_num)(*u)
        cv = (c_double * dim_num)(*v)
        tpcp4 = lib.basis_map_3d(cu,cv)
        cp4 = [dim_num]
        if tpcp4 != None:
            pcp4 = pointer((c_double * dim_num)(tpcp4))
            for i in xrange(dim_num):
                cp4.append(pcp4.contents[i])
        return cp4

    # int box_01_contains_point_2d ( double p[2] );
    def box_01_contains_point_2d (self, p ):
        cp = (c_double * 2)(*p)
        return lib.box_01_contains_point_2d ( cp )


    # int box_01_contains_point_nd ( int dim_num, double p[] );
    def box_01_contains_point_nd (self, dim_num, p ):
        cp = (c_double * len(p))(*p)
        return lib.box_01_contains_point_nd(c_int(dim_num), cp)

    # int box_contains_point_2d ( double p1[2], double p2[2], double p[2] );
    def box_contains_point_2d (self, p1, p2, p ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp = (c_double * 2)(*p)
        return lib.box_contains_point_2d( p1, p2, p )

    # int box_contains_point_nd ( int dim_num, double p1[], double p2[], double p[] );
    def box_contains_point_nd (self, dim_num, p1, p2, p ):
        cp1 = (c_double * len(p1))(*p1)
        cp2 = (c_double * len(p2))(*p2)
        cp = (c_double * len(p))(*p)
        return lib.box_contains_point_nd(c_int(dim_num), cp1, cp2, cp) 


    # void box_ray_int_2d ( double p1[2], double p2[2], double pa[2], double pb[2], double pint[2] );
    def box_ray_int_2d (self, p1, p2, pa, pb, pint ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cpa = (c_double * 2)(*pa)
        cpb = (c_double * 2)(*pb)
        cpint = (c_double * 2)(*pint)
        return lib.box_ray_int_2d(cp1, cp2, cpa, cpb, cpint) 

    # int box_segment_clip_2d ( double p1[2], double p2[2], double pa[2], double pb[2] );
    def box_segment_clip_2d (self, p1, p2, pa, pb ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cpa = (c_double * 2)(*pa)
        cpb = (c_double * 2)(*pb)
        return lib.box_segment_clip_2d(cp1, cp2, cpa, cpb)


    # void circle_arc_point_near_2d ( double r, double pc[2], double theta1, double theta2, double p[2], double pn[2], double *dist );
    def circle_arc_point_near_2d (self, r, pc, theta1, theta2, p, pn, dist ):
        cpc = (c_double * 2)(*pc)
        cp = (c_double * 2)(*p)
        cpn = (c_double * 2)(*pn)

        cp2 = (c_double * 2)(*p2)
        cpa = (c_double * 2)(*pa)
        pdist = pointer(c_double(dist))

        lib.circle_arc_point_near_2d(c_double(r), cpc, c_double(theta1), c_double(theta2), cp, cpn, pdist)


    #double circle_area_2d ( double r );
    def circle_area_2d (self, r ):
        return lib.circle_area_2d(c_double(r))


    #void circle_dia2imp_2d ( double p1[2], double p2[2], double *r, double pc[2] );
    def circle_dia2imp_2d (self, p1, p2, r, pc ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cpc = (c_double * 2)(*pc)
        pr = pointer(c_double(r))
        lib.circle_dia2imp_2d(cp1, cp2, pr, cpc)


    #int circle_exp_contains_point_2d ( double p1[2], double p2[2], double p3[2], double p[2] );
    def circle_exp_contains_point_2d (self, p1, p2, p3, p ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        cp = (c_double * 2)(*p)
        return lib.circle_exp_contains_point_2d(cp1, cp2, cp3, cp)

    # void circle_exp2imp_2d ( double p1[2], double p2[2], double p3[2], double *r, double pc[2] );
    def circle_exp2imp_2d (self, p1, p2, p3, r, pc ):
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        pr = pointer(c_double(r))
        cpc = (c_double * 2)(*pc)
        return lib.circle_exp2imp_2d(cp1, cp2, cp3, pr, cpc)
        


    # int circle_imp_contains_point_2d ( double r, double pc[2], double p[2] );
    def circle_imp_contains_point_2d (self, r, pc, p ):
        cpc = (c_double * 2)(*pc)
        cp = (c_double * 2)(*p)
        return lib.circle_imp_contains_point_2d(c_double(r), cpc, cp)


    # void circle_imp_line_par_int_2d ( double r, double pc[2], double x0, double y0, double f, double g, int *int_num, double p[] );
    def circle_imp_line_par_int_2d (self, r, pc, x0, y0, f, g, int_num, p ):
        cpc = (c_double * 2)(*pc)
        pint_num = pointer(c_int(int_num))
        cp = (c_double * len(p))(*p)
        lib.circle_imp_line_par_int_2d(c_double(r), cpc, c_double(x0), c_double(y0), c_double(f), c_double(g), pint_num, cp)


    # double circle_imp_point_dist_2d ( double r, double pc[2], double p[2]  );
    def circle_imp_point_dist_2d (self, r, pc, p ):
        cpc = (c_double * 2)(*pc)
        cp = (c_double * 2)(*p)
        return lib.circle_imp_point_dist_2d(c_double(r), cpc, cp)


    # double circle_imp_point_dist_signed_2d ( double r, double pc[2], double p[2] );
    def circle_imp_point_dist_signed_2d (self, r, pc, p ):
        cpc = (c_double * 2)(*pc)
        cp = (c_double * 2)(*p)
        return lib.circle_imp_point_dist_signed_2d(c_double(r), cpc, cp)


    # double circle_imp_point_near_2d ( double r, double pc[2], double p[2], double pn[2] );
    def circle_imp_point_near_2d (self, r, pc, p, pn ):
        cpc = (c_double * 2)(*pc)
        cp = (c_double * 2)(*p)
        cpn = (c_double * 2)(*pn)
        lib.circle_imp_point_near_2d(c_double(r), cpc, cp, cpn)


    # double *circle_imp_points_2d ( double r, double pc[2], int n );
    def circle_imp_points_2d (self, r, pc, n ):
        cpc = (c_double * 2)(*pc)
        return lib.circle_imp_points_2d(c_double(r), cpc, c_int(n)) 


    # double *circle_imp_points_3d ( double r, double pc[3], double nc[3], int n );
    def circle_imp_points_3d (self, r, pc, nc, n ):
        cpc = (c_double * 3)(*pc)
        cnc = (c_double * 3)(*nc)
        return lib.circle_imp_points_3d( c_double(r), cpc, cnc, c_int(n))



    # void circle_imp_points_arc_2d ( double r, double pc[2], double theta1, double theta2, int n, double p[] );
    def circle_imp_points_arc_2d (self, r, pc, theta1, theta2, n, p ):
        cpc = (c_double * 2)(*pc)
        cp = (c_double * len(p))(*p)
        lib.circle_imp_points_arc_2d(c_double(r), cpc, c_double(theta1), c_double(theta2), c_int(n), cp)
     

    # void circle_imp_print_2d ( double r, double pc[2], char *title );
    def circle_imp_print_2d (self, r, pc, title ):
        cpc = (c_double * 2)(*pc)
        ptitle = pointer(c_char(title))
        lib.circle_imp_print_2d(c_double(r), cpc, ptitle)


    # void circle_imp_print_3d ( double r, double pc[3], double nc[3], char *title );
    def circle_imp_print_3d (self, r, pc, nc, title ):
        cpc = (c_double * 3)(*pc)
        cnc = (c_double * 3)(*nc)
        ptitle = pointer(c_char(title))
        lib.circle_imp_print_2d(c_double(r), cpc, cnc, ptitle)

    # void circle_imp2exp_2d ( double r, double pc[2], double p1[2], double p2[2], double p3[2] );
    def circle_imp2exp_2d (self, r, pc, p1, p2, p3 ):
        cpc = (c_double * 2)(*pc)
        cp1 = (c_double * 2)(*p1)
        cp2 = (c_double * 2)(*p2)
        cp3 = (c_double * 2)(*p3)
        lib.circle_imp2exp_2d(c_double(r), cpc, cp1, cp2, cp3)

    # double *circle_llr2imp_2d ( double p1[], double p2[], double q1[], double q2[], double r );
    def circle_llr2imp_2d (self, p1, p2, q1, q2, r ):
        cp1 = (c_double * len(p1))(*p1)
        cp2 = (c_double * len(p2))(*p2)
        cq1 = (c_double * len(q1))(*q1)
        cq2 = (c_double * len(q2))(*q2)
        return lib.circle_llr2imp_2d(cp1, cp2, cq1, cq2)

    # double circle_lune_area_2d ( double r, double pc[2], double theta1, double theta2 );
    def circle_lune_area_2d (self, r, pc, theta1, theta2 ):
        cpc = (c_double * 2)(*pc)
        return lib.circle_lune_area_2d(c_double(r), cpc, c_double(theta1), c_double(theta2))
        

    # double *circle_lune_centroid_2d ( double r, double pc[2], double theta1, double theta2 );
    def circle_lune_centroid_2d (self, r, pc, theta1, theta2 ):
        cpc = (c_double * 2)(*pc)
        return lib.circle_lune_centroid_2d(c_double(r), cpc, c_double(theta1), c_double(theta2) )
        

    # void circle_pppr2imp_3d ( double p1[], double p2[], double p3[], double r, double pc[], double normal[] );
    def circle_pppr2imp_3d (self, p1, p2, p3, r, pc, normal ):
        cp1 = (c_double * len(p1))(*p1)
        cp2 = (c_double * len(p2))(*p2)
        cp3 = (c_double * len(p3))(*p3)
        cpc = (c_double * len(pc))(*pc)
        cnormal = (c_double * len(normal))(*normal)
        return lib.circle_pppr2imp_3d(cp1, cp2, cp3, c_double(r), cpc, cnormal)
        
    # double *circle_ppr2imp_2d ( double p1[], double p2[], double r );
    def circle_ppr2imp_2d ( self, p1, p2, r ):
        cp1 = (c_double * len(p1))(*p1)
        cp2 = (c_double * len(p2))(*p2)
        pcp4 = pointer((c_double * 4)(lib.circle_ppr2imp_2d(cp1, cp2, c_double(r))))
        cp4 = []
        cp4.append(pcp4.contents[0])
        cp4.append(pcp4.contents[1])
        cp4.append(pcp4.contents[2])
        cp4.append(pcp4.contents[3])
        return cp4



    # double circle_sector_area_2d ( double r, double pc[2], double theta1, double theta2 );

    # double *circle_sector_centroid_2d ( double r, double pc[2], double theta1, double theta2 );

    # int circle_sector_contains_point_2d ( double r, double pc[2], double theta1, double theta2, double p[2] );

    # void circle_sector_print_2d ( double r, double pc[2], double theta1, double theta2 );

    # double circle_triangle_area_2d ( double r, double pc[2], double theta1, double theta2 );

    # void circle_triple_angles_2d ( double r1, double r2, double r3, double *angle1, double *angle2, double *angle3 );
    
    # void circles_imp_int_2d ( double r1, double pc1[2], double r2, double pc2[2], int *int_num, double p[] );

    # double cone_area_3d ( double h, double r );

    # double *cone_centroid_3d ( double r, double pc[3], double pt[3] );

    # double cone_volume_3d ( double h, double r );

    # void conv3d ( char axis, double theta, int n, double cor3[], double cor2[] );

    # double cos_deg ( double angle );

    # double cot_deg ( double angle );

    # double cot_rad ( double angle );

    # double csc_deg ( double angle );

    # void cube_shape_3d ( int point_num, int face_num, int face_order_max, double point_coord[], int face_order[], int face_point[] );

    # void cube_size_3d ( int *point_num, int *edge_num, int *face_num, int *face_order_max );

    #double cylinder_point_dist_3d ( double p1[3], double p2[3], double r, double p[3] );

    # double cylinder_point_dist_signed_3d ( double p1[3], double p2[3], double r, double p[3] );

    # int cylinder_point_inside_3d ( double p1[3], double p2[3], double r, double p[3] );

    # double *cylinder_point_near_3d ( double p1[3], double p2[3], double r, double p[3] );

    # double *cylinder_sample_3d ( double p1[3], double p2[3], double r, int n, int *seed );

    # double cylinder_volume_3d ( double p1[3], double p2[3], double r );

    # double degrees_to_radians ( double degrees );

    # double dge_det ( int n, double a[], int pivot[] );

    # int dge_fa ( int n, double a[], int pivot[] );
    # void dge_sl ( int n, double a[], int pivot[], double b[], int job );
    # double *direction_pert_3d ( double sigma, double vbase[3], int *seed );
    # double *direction_uniform_2d ( int *seed );
    # double *direction_uniform_3d ( int *seed );
    # double *direction_uniform_nd ( int dim_num, int *seed );
    # double disk_point_dist_3d ( double pc[3], double r, double axis[3], double p[3] );

    # double dms_to_radians ( int degrees, int minutes, int seconds );
    # void dodec_shape_3d ( int point_num, int face_num, int face_order_max, double point_coord[], int face_order[], int face_point[] );
    # void dodec_size_3d ( int *point_num, int *edge_num, int *face_num, int *face_order_max );
    # void dual_shape_3d ( int point_num, int face_num, int face_order_max, double point_coord[], int face_order[], int face_point[], int point_num2, int face_num2, int face_order_max2, double point_coord2[], int face_order2[], int face_point2[] );
    # void dual_size_3d ( int point_num, int edge_num, int face_num, int face_order_max, double point_coord[], int face_order[], int face_point[], int *point_num2, int *edge_num2, int *face_num2, int *face_order_max2 );

    # double ellipse_area_2d ( double r1, double r2 );
    # double ellipse_point_dist_2d ( double r1, double r2, double p[2] );
    # double *ellipse_point_near_2d ( double r1, double r2, double p[2] );
    # void ellipse_points_2d ( double pc[2], double r1, double r2, double psi, int n, double p[] );

    # void ellipse_points_arc_2d ( double pc[2], double r1, double r2, double psi, double theta1, double theta2, int n, double p[] );

    # double enorm0_nd ( int n, double x[], double y[] );
    # int get_seed ( void );
    # void glob2loc_3d ( double cospitch, double cosroll, double cosyaw, double sinpitch, double sinroll, double sinyaw, double globas[3], double glopts[3], double locpts[3] );

    # int halfplane_contains_point_2d ( double pa[2], double pb[2], double p[2] );
    # int halfspace_imp_triangle_int_3d ( double a, double b, double c, double d, double t[3*3], double p[3*4] );
    # int halfspace_norm_triangle_int_3d ( double pp[3], double pn[3], double t[3*3], double p[3*4] );
    # int halfspace_triangle_int_3d ( double dist1, double dist2, double dist3, double t[3*3], double p[3*4] );
    # double haversine ( double a );
    # void helix_shape_3d ( double a, int n, double r, double theta1, double theta2, double p[] );
    # double hexagon_area_2d ( double r );
    # int hexagon_contains_point_2d ( double v[2*6], double p[2] );
    # void hexagon_shape_2d ( double angle, double p[2] );
    # double hexagon_unit_area_2d ( );
    # void hexagon_vertices_2d ( double h[2*6] );
    # double i4_dedekind_factor ( int p, int q );
    # double i4_dedekind_sum ( int p, int q );
    # int i4_factorial2 ( int n );
    # int i4_gcd ( int i, int j );
    # int i4_lcm ( int i, int j );
    # int i4_max ( int i1, int i2 );
    # int i4_min ( int i1, int i2 );
    # int i4_modp ( int i, int j );
    # int i4_sign ( int i );
    # void i4_swap ( int *i, int *j );
    # int i4_uniform ( int a, int b, int *seed );
    # int i4_wrap ( int ival, int ilo, int ihi );
    # int i4col_compare ( int m, int n, int a[], int i, int j );
    # void i4col_find_item ( int m, int n, int a[], int item, int *row, int *col );
    # void i4col_find_pair_wrap ( int m, int n, int a[], int item1, int item2, int *row, int *col );
    # void i4col_sort_a ( int m, int n, int a[] );
    # int i4col_sorted_unique_count ( int m, int n, int a[] );
    # void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
    # void i4mat_print ( int m, int n, int a[], char *title );
    # void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi, int jhi, char *title );
    # void i4mat_transpose_print ( int m, int n, int a[], char *title );
    # void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi, int jhi, char *title );
    # int i4row_compare ( int m, int n, int a[], int i, int j );
    # void i4row_sort_a ( int m, int n, int a[] );
    # void i4row_swap ( int m, int n, int a[], int irow1, int irow2 );
    # void i4vec_copy ( int n, int a1[], int a2[] );
    # void i4vec_heap_d ( int n, int a[] );
    # int *i4vec_indicator_new ( int n );
    # int i4vec_lcm ( int n, int v[] );
    # void i4vec_print ( int n, int a[], char *title );
    # int i4vec_product ( int n, int a[] );
    # void i4vec_reverse ( int n, int a[] );
    # void i4vec_sort_heap_a ( int n, int a[] );
    # int i4vec_sorted_unique ( int n, int a[] );
    # int *i4vec_uniform_new ( int n, int a, int b, int *seed );
    # void i4vec_zero ( int n, int a[] );
    # int i4vec2_compare ( int n, int a1[], int a2[], int i, int j );
    # void i4vec2_sort_a ( int n, int a1[], int a2[] );
    # void i4vec2_sorted_unique ( int n, int a1[], int a2[], int *nuniq );
    # void icos_shape ( int point_num, int edge_num, int face_num, int face_order_max, double point_coord[], int edge_point[], int face_order[], int face_point[] );
    # void icos_size ( int *point_num, int *edge_num, int *face_num, int *face_order_max );
    # int line_exp_is_degenerate_nd ( int dim_num, double p1[], double p2[] );
    # double *line_exp_normal_2d ( double p1[2], double p2[2] );
    # double *line_exp_perp_2d ( double p1[2], double p2[2], double p3[2], int *flag );
    # double line_exp_point_dist_2d ( double p1[2], double p2[2], double p[2] );
    # double line_exp_point_dist_3d ( double p1[3], double p2[3], double p[3] );
    # double line_exp_point_dist_signed_2d ( double p1[2], double p2[2], double p[2] );
    # void line_exp_point_near_2d ( double p1[2], double p2[2], double p[2], double pn[2], double *dist, double *t );
    # void line_exp_point_near_3d ( double p1[3], double p2[3], double p[3], double pn[3], double *dist, double *t );
    # void line_exp2imp_2d ( double p1[2], double p2[2], double *a, double *b, double *c );
    # void line_exp2par_2d ( double p1[2], double p2[2], double *f, double *g, double *x0, double *y0 );
    # void line_exp2par_3d ( double p1[3], double p2[3], double *f, double *g, double *h, double *x0, double *y0, double *z0 );
    # int line_imp_is_degenerate_2d ( double a, double b, double c );
    # double line_imp_point_dist_2d ( double a, double b, double c, double p[2] );
    # double line_imp_point_dist_signed_2d ( double a, double b, double c, double p[2] );
    # void line_imp2exp_2d ( double a, double b, double c, double p1[2], double p2[2] );
    # void line_imp2par_2d ( double a, double b, double c, double *f, double *g, double *x0, double *y0 );
    # double line_par_point_dist_2d ( double f, double g, double x0, double y0, double p[2] );
    # double line_par_point_dist_3d ( double f, double g, double h, double x0, double y0, double z0, double p[3] );
    # double *line_par_point_near_2d ( double f, double g, double x0, double y0, double p[2] );
    # double *line_par_point_near_3d ( double f, double g, double h, double x0, double y0, double z0, double p[3] );
    # void line_par2exp_2d ( double f, double g, double x0, double y0, double p1[2], double p2[2] );
    # void line_par2exp_3d ( double f, double g, double h, double x0, double y0, double z0, double p1[3], double p2[3] );
    # void line_par2imp_2d ( double f, double g, double x0, double y0, double *a, double *b, double *c );
    # double lines_exp_angle_3d ( double p1[3], double p2[3], double p3[3], double p4[3] );
    # double lines_exp_angle_nd ( double p1[], double p2[], double q1[], double q2[], int n );
    # double lines_exp_dist_3d ( double p1[3], double p2[3], double q1[3], double q2[3] );
    # double lines_exp_dist_3d_2 ( double p1[3], double p2[3], double q1[3], double q2[3] );
    # int lines_exp_equal_2d ( double p1[2], double p2[2], double q1[2], double q2[2] );
    # void lines_exp_int_2d ( double p1[2], double p2[2], double p3[2], double p4[2], int *ival, double p[2] );
    # void lines_exp_near_3d ( double p1[3], double p2[3], double q1[3], double q2[3], double pn[3], double qn[3] );
    # int lines_exp_parallel_2d ( double p1[2], double p2[2], double q1[2], double q2[2] );
    # int lines_exp_parallel_3d ( double p1[3], double p2[3], double q1[3], double q2[3] );
    # double lines_imp_angle_2d ( double a1, double b1, double c1, double a2, double b2, double c2 );
    # double lines_imp_dist_2d ( double a1, double b1, double c1, double a2, double b2, double c2 );
    # void lines_imp_int_2d ( double a1, double b1, double c1, double a2, double b2, double c2, int *ival, double p[2] );
    # double lines_par_angle_2d ( double f1, double g1, double x01, double y01, double f2, double g2, double x02, double y02 );
    # double lines_par_angle_3d ( double f1, double g1, double h1, double x01, double y01, double z01, double f2, double g2, double h2, double x02, double y02, double z02 );
    # double lines_par_dist_3d ( double f1, double g1, double h1, double x01, double y01, double z01, double f2, double g2, double h2, double x02, double y02, double z02 );
    # void lines_par_int_2d ( double f1, double g1, double x1, double y1, double f2, double g2, double x2, double y2, double *t1, double *t2, double pint[2] );
    # void loc2glob_3d ( double cospitch, double cosroll, double cosyaw, double sinpitch, double sinroll, double sinyaw, double locpts[3], double globas[3], double glopts[3] );
    # int lrline ( double xu, double yu, double xv1, double yv1, double xv2, double yv2, double dv );

    # void lvec_print ( int n, int a[], char *title );
    # void minabs ( double x1, double y1, double x2, double y2, double x3, double y3, double *xmin, double *ymin );
    # int minquad ( double x1, double y1, double x2, double y2, double x3, double y3, double *xmin, double *ymin );
    # void octahedron_shape_3d ( int point_num, int face_num, int face_order_max, double point_coord[], int face_order[], int face_point[] );
    # void octahedron_size_3d ( int *point_num, int *edge_num, int *face_num, int *face_order_max );
    # int parabola_ex ( double x1, double y1, double x2, double y2, double x3, double y3, double *x, double *y );
    # int parabola_ex2 ( double x1, double y1, double x2, double y2, double x3, double y3, double *x, double *y, double *a, double *b, double *c );
    # double parallelogram_area_2d ( double p[] );
    # double parallelogram_area_3d ( double p[] );
    # int parallelogram_contains_point_2d ( double p1[2], double p2[2], double p3[2], double p[2] );
    # int parallelogram_contains_point_3d ( double p1[3], double p2[3], double p3[3], double p[3] );
    # double parallelogram_point_dist_3d ( double p1[3], double p2[3], double p3[3], double p[3] );
    # int parallelepiped_contains_point_3d ( double p1[3], double p2[3], double p3[3], double p4[3], double p[3] );
    # double parallelepiped_point_dist_3d ( double p1[3], double p2[3], double p3[3], double p4[3], double p[3] );
    # int perm_check ( int n, int p[] );
    # void perm_inv ( int n, int p[] );
    # void plane_exp_grid_3d ( double p1[3], double p2[3], double p3[3], int *ncor3, int *line_num, double cor3[], int lines[], int maxcor3, int line_max, int *ierror );
    # double plane_exp_point_dist_3d ( double p1[3], double p2[3], double p3[3], double p[3] );
    # void plane_exp_normal_3d ( double p1[3], double p2[3], double p3[3], double pn[3] );
    # void plane_exp_pro2 ( double p1[3], double p2[3], double p3[3], int npnt, double pp[], double alpha[], double beta[] );
    # void plane_exp_pro3 ( double p1[3], double p2[3], double p3[3], int n, double po[], double pp[] );
    # void plane_exp_project_3d ( double p1[3], double p2[3], double p3[3], double pf[3], int n, double po[], double pp[], int ivis[] );
    # void plane_exp2imp_3d ( double p1[3], double p2[3], double p3[3], double *a, double *b, double *c, double *d );
    # void plane_exp2normal_3d ( double p1[3], double p2[3], double p3[3], double pp[3], double pn[3] );
    # int plane_imp_is_degenerate_3d ( double a, double b, double c );
    # int plane_imp_line_par_int_3d ( double a, double b, double c, double d, double x0, double y0, double z0, double f, double g, double h, double p[3] );
    # double plane_imp_point_dist_3d ( double a, double b, double c, double d, double p[3] );
    # double plane_imp_point_dist_signed_3d ( double a, double b, double c, double d, double p[3] );
    # void plane_imp_point_near_3d ( double a, double b, double c, double d, double p[3], double pn[3] );
    # void plane_imp_segment_near_3d ( double p1[3], double p2[3], double a, double b, double c, double d, double *dist, double pnp[3], double pnl[3] );
    # void plane_imp_triangle_int_3d ( double a, double b, double c, double d, double t[3*3], int *int_num, double p[] );
    # void plane_imp_triangle_int_add_3d ( double p1[3], double p2[3], double dist1, double dist2, int *int_num, double p[] );
    # int plane_imp_triangle_near_3d ( double t[3*3], double a, double b, double c, double d, double *dist, double pn[] );
    # void plane_imp2exp_3d ( double a, double b, double c, double d, double p1[3], double p2[3], double p3[3] );
    # void plane_imp2normal_3d ( double a, double b, double c, double d, double pp[3], double pn[3] );
    # void plane_normal_basis_3d ( double pp[3], double pn[3], double pq[3], double pr[3] );
    # int plane_normal_line_exp_int_3d ( double pp[3], double normal[3], double p1[3], double p2[3], double pint[3] );
    # double *plane_normal_qr_to_xyz ( double pp[], double normal[], double pq[], double pr[], int n, double qr[] );
    # void plane_normal_tetrahedron_intersect ( double pp[3], double normal[3], double t[3*4], int *int_num, double pint[3*4] );
    # int plane_normal_triangle_int_3d ( double pp[3], double pn[3], double t[3*3], double p[3*3] );


    # void plane_normal_uniform_3d ( int *seed, double pp[3], double normal[3] );
    # void plane_normal_uniform_nd ( int dim_num, int *seed, double pp[], double normal[] );
    # double *plane_normal_xyz_to_qr ( double pp[], double normal[], double pq[], double pr[], int n, double xyz[] );
    # void plane_normal2exp_3d ( double pp[3], double pn[3], double p1[3], double p2[3], double p3[3] );
    # void plane_normal2imp_3d ( double pp[3], double pn[3], double *a, double *b, double *c, double *d );
    # double planes_imp_angle_3d ( double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2 );
    # int points_avoid_point_naive_2d ( int n, double pset[], double p[2] );
    # void points_bisect_line_imp_2d ( double p1[2], double p2[2], double *a, double *b, double *c );
    # void points_bisect_line_par_2d ( double p1[2], double p2[2], double *f, double *g, double *x, double *y );
    # int points_centroid_2d ( int n, double p[] );
    # double points_colin_2d ( double p1[2], double p2[2], double p3[2] );
    # double points_colin_3d ( double p1[3], double p2[3], double p3[3] );
    # double points_dist_2d ( double p1[2], double p2[2] );
    # double points_dist_3d ( double p1[3], double p2[3] );
    # double points_dist_nd ( int dim_num, double p1[], double p2[] );
    # void points_hull_2d ( int node_num, double node_xy[], int *hull_num, int hull[] );
    # void points_plot ( char *file_name, int node_num, double node_xy[], int node_label );
    # int points_point_near_naive_2d ( int nset, double pset[], double ptest[], double *d_min );
    # int points_point_near_naive_3d ( int nset, double pset[], double ptest[], double *d_min );
    # int points_point_near_naive_nd ( int dim_num, int nset, double pset[], double ptest[], double *d_min );
    # int *points_points_near_naive_2d ( int nset, double pset[], int ntest, double ptest[] );
    # int *points_points_near_naive_3d ( int nset, double pset[], int ntest, double ptest[] );
    # void polar_to_xy ( double r, double t, double xy[2] );
    # double polygon_1_2d ( int n, double v[] );
    # double *polygon_angles_2d ( int n, double v[] );
    # double polygon_area_2d ( int n, double v[] );
    # double polygon_area_2d_2 ( int n, double v[] );
    # double polygon_area_3d ( int n, double v[], double normal[] );
    # double polygon_area_3d_2 ( int n, double v[] );
    # double *polygon_centroid_2d ( int n, double v[] );
    # double *polygon_centroid_2d_2 ( int n, double v[] );
    # double *polygon_centroid_3d ( int n, double v[] );


    # int polygon_contains_point_2d ( int n, double v[], double p[2] );
    def polygon_contains_point_2d(self, v, p):
        n = len(v)
        p1 = (c_double * 2)(*p)
        pn = (c_double * len(v))(*v)
        return lib.polygon_contains_point_2d (c_int(n), pn, p1)

    # int polygon_contains_point_2d_2 ( int n, double v[], double p[2] );
    def polygon_contains_point_2d_2(self, v, p):
        n = len(v)
        p1 = (c_double * 2)(*p)
        pn = (c_double * len(v))(*v)
        return lib.polygon_contains_point_2d_2 (c_int(n), pn, p1)

    # double polygon_diameter_2d ( int n, double v[] );
    # double *polygon_expand_2d ( int n, double v[], double h );
    # void polygon_inrad_data_2d ( int n, double radin, double *area, double *radout, double *side );
    # int polygon_is_convex ( int n, double v[] );
    # double polygon_lattice_area_2d ( int i, int b );
    # double *polygon_normal_3d ( int n, double v[] );
    # void polygon_outrad_data_2d ( int n, double radout, double *area, double *radin, double *side );
    # void polygon_side_data_2d ( int n, double side, double *area, double *radin, double *radout );
    # double polygon_solid_angle_3d ( int n, double v[], double p[3] );
    # double polygon_x_2d ( int n, double v[] );
    # double polygon_y_2d ( int n, double v[] );
    # double polygon_xx_2d ( int n, double v[] );
    # double polygon_xy_2d ( int n, double v[] );
    # double polygon_yy_2d ( int n, double v[] );
    # double polyhedron_area_3d ( double coord[], int order_max, int face_num, int node[], int node_num, int order[] );
    # double *polyhedron_centroid_3d ( double coord[], int order_max, int face_num, int node[], int node_num, int order[] );
    # int polyhedron_contains_point_3d ( int node_num, int face_num, int face_order_max, double v[], int face_order[], int face_point[], double p[3] );
    # double polyhedron_volume_3d ( double coord[], int order_max, int face_num, int node[], int node_num, int order[] );
    # double polyhedron_volume_3d_2 ( double coord[], int order_max, int face_num, int node[], int node_num, int order[] );
    # double *polyline_arclength_nd ( int dim_num, int n, double p[] );
    # double *polyline_index_point_nd ( int dim_num, int n, double p[], double t );
    # double polyline_length_nd ( int dim_num, int n, double p[] );
    # double *polyline_points_nd ( int dim_num, int n, double p[], int nt );
    # double *polyloop_arclength_nd ( int dim_num, int nk, double pk[] );
    # double *polyloop_points_nd ( int dim_num, int nk, double pk[], int nt );
    # void provec ( int m, int n, double base[], double vecm[], double vecn[], double vecnm[] );
    # double pyramid_volume_3d ( double h, double s );
    # double quad_area_2d ( double q[2*4] );
    # double quad_area2_2d ( double q[] );
    # double quad_area_3d ( double q[] );
    # int quad_contains_point_2d ( double q[2*4], double p[2] );
    # void quad_convex_random ( int *seed, double xy[] );
    # double quad_point_dist_2d ( double q[2*4], double p[2] );
    # double quad_point_dist_signed_2d ( double q[2*4], double p[2] );
    # void quad_point_near_2d ( double q[2*4], double p[2], double pn[2], double *dist );
    # double *quat_conj ( double q[] );
    # double *quat_inv ( double q[] );
    # double *quat_mul ( double q1[], double q2[] );
    # double quat_norm ( double q[] );
    # float r4_abs ( float x );
    # int r4_nint ( float x );
    # double r8_abs ( double x );
    # double r8_acos ( double c );
    # double r8_asin ( double s );
    # double r8_atan ( double y, double x );
    # double r8_epsilon ( void );
    # double r8_huge ( void );
    # double r8_max ( double x, double y );
    # double r8_min ( double x, double y );
    # double r8_modp ( double x, double y );
    # int r8_nint ( double x );
    # double r8_normal_01 ( int *seed );
    # double r8_pi ( void );
    # double r8_sign ( double x );
    # int r8_sign_opposite_strict ( double r1, double r2 );
    # void r8_swap ( double *x, double *y );
    # double r8_uniform ( double b, double c, int *seed );
    # double r8_uniform_01 ( int *seed );
    # void r82vec_part_quick_a ( int n, double a[], int *l, int *r );
    # void r82vec_permute ( int n, double a[], int p[] );
    # void r82vec_print ( int n, double a[], char *title );
    # int *r82vec_sort_heap_index_a ( int n, double a[] );
    # void r82vec_sort_quick_a ( int n, double a[] );
    # void r8mat_copy ( int m, int n, double a1[], double a2[] );
    # double r8mat_det_2d ( double a[] );
    # double r8mat_det_3d ( double a[] );
    # double r8mat_det_4d ( double a[] );
    # double r8mat_det_5d ( double a[] );
    # double *r8mat_inverse_2d ( double a[] );
    # double *r8mat_inverse_3d ( double a[] );
    # double *r8mat_mv ( int m, int n, double a[], double x[] );
    # void r8mat_print ( int m, int n, double a[], char *title );
    # void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, int jhi, char *title );
    # int r8mat_solve ( int n, int rhs_num, double a[] );
    # double *r8mat_solve_2d ( double a[], double b[], double *det );
    # void r8mat_transpose_print ( int m, int n, double a[], char *title );
    # void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, int jhi, char *title );
    # double *r8mat_uniform_new ( int m, int n, double alo, double ahi, int *seed );
    # void r8mat_uniform_01 ( int m, int n, int *seed, double r[] );
    # double *r8mat_uniform_01_new ( int m, int n, int *seed );
    # double r8vec_angle_3d ( double u[], double v[] );
    # double *r8vec_any_normal ( int dim_num, double v1[] );
    # void r8vec_bracket ( int n, double x[], double xval, int *left, int *right );
    # void r8vec_copy ( int n, double a1[], double a2[] );
    # double r8vec_cross_product_2d ( double v1[2], double v2[2] );
    # double r8vec_cross_product_2d_affine ( double v0[2], double v1[2], double v2[2] );
    # double *r8vec_cross_product_3d ( double v1[3], double v2[3] );
    # double *r8vec_cross_product_3d_affine ( double v0[3], double v1[3], double v2[3] );
    # double r8vec_distance ( int dim_num, double v1[], double v2[] );
    # double r8vec_dot_product ( int n, double a1[], double a2[] );
    # double r8vec_dot_product_affine ( int n, double v0[], double v1[], double v2[] );
    # int r8vec_eq ( int n, double a1[], double a2[] );
    # int r8vec_gt ( int n, double a1[], double a2[] );
    # int r8vec_lt ( int n, double a1[], double a2[] );
    # double r8vec_max ( int n, double *rvec );
    # double r8vec_mean ( int n, double x[] );
    # double r8vec_min ( int n, double *rvec );
    # int r8vec_negative_strict ( int n, double a[] );
    # double r8vec_norm ( int n, double x[] );
    # double r8vec_norm_affine ( int n, double v0[], double v1[] );
    # double *r8vec_normal_01_new ( int n, int *seed );
    # double r8vec_normsq ( int n, double a[] );
    # double r8vec_normsq_affine ( int n, double v0[], double v1[] );
    # int r8vec_positive_strict ( int n, double a[] );
    # void r8vec_print ( int n, double a[], char *title );
    # void r8vec_print_2d ( double x, double y, char *title );
    # void r8vec_print_3d ( double x, double y, double z, char *title );
    # double r8vec_scalar_triple_product ( double v1[3], double v2[3], double v3[3] );
    # void r8vec_swap ( int n, double a1[], double a2[] );
    # double *r8vec_uniform_new ( int n, double b, double c, int *seed );
    # double *r8vec_uniform_01_new ( int n, int *seed );
    # double *r8vec_uniform_unit_new ( int m, int *seed );
    # double r8vec_variance ( int n, double x[] );
    # void r8vec_zero ( int n, double a[] );
    # double radec_distance_3d ( double ra1, double dec1, double ra2, double dec2 );
    # double *radec_to_xyz ( double ra, double dec );
    # double radians_to_degrees ( double angle );
    # void radians_to_dms ( double radians, int *degrees, int *minutes, int *seconds );
    # unsigned long random_initialize ( unsigned long seed );
    # void rotation_axis_vector_3d ( double axis[3], double angle, double v[3], double w[3] );
    # void rotation_axis2mat_3d ( double axis[3], double angle, double a[3*3] );
    # void rotation_axis2quat_3d ( double axis[3], double angle, double q[4] );
    # void rotation_mat_vector_3d ( double a[3*3], double v[3], double w[3] );
    # void rotation_mat2axis_3d ( double a[3*3], double axis[3], double *angle );
    # void rotation_mat2quat_3d ( double a[3*3], double q[4] );
    # void rotation_quat_vector_3d ( double q[4], double v[3], double w[3] );
    # void rotation_quat2axis_3d ( double q[4], double axis[3], double *angle );
    # void rotation_quat2mat_3d ( double q[4], double a[3*3] );
    # void rtp_to_xyz ( double r, double theta, double phi, double xyz[3] );
    # int s_len_trim ( char *s );
    # double sec_deg ( double angle );
    # void segment_contains_point_1d ( double p1, double p2, double p3, double *u );
    # void segment_contains_point_2d ( double p1[2], double p2[2], double p3[2], double u[2] );
    # void segment_point_coords_2d ( double p1[2], double p2[2], double p[2], double *s, double *t );
    # void segment_point_coords_3d ( double p1[3], double p2[3], double p[3], double *s, double *t );
    # double segment_point_dist_2d ( double p1[2], double p2[2], double p[2] );
    # double segment_point_dist_3d ( double p1[3], double p2[3], double p[3] );
    # void segment_point_near_2d ( double p1[2], double p2[2], double p[2], double pn[2], double *dist, double *t );
    # void segment_point_near_3d ( double p1[3], double p2[3], double p[3], double pn[3], double *dist, double *t );
    # double segments_curvature_2d ( double p1[2], double p2[2], double p3[2] );
    # double segments_dist_2d ( double p1[2], double p2[2], double q1[2], double q2[2] );
    # double segments_dist_3d ( double p1[3], double p2[3], double q1[3], double q2[3] );
    # double segments_dist_3d_old ( double p1[3], double p2[3], double q1[3], double q2[3] );
    # double segments_int_1d ( double p1, double p2, double q1, double q2, double *r1, double *r2 );
    # void segments_int_2d ( double p1[2], double p2[2], double p3[2], double p4[2], int *flag, double p5[2] );
    # double shape_point_dist_2d ( double pc[2], double p1[2], int nside, double p[2] );
    # void shape_point_near_2d ( double pc[2], double p1[2], int nside, double p[2], double pn[2], double *dist );
    # void shape_print_3d ( int point_num, int face_num, int face_order_max, double point_coord[], int face_order[], int face_point[] );
    # void shape_ray_int_2d ( double pc[2], double p1[2], int nside, double pa[2], double pb[2], double pi[2] );
    # void simplex_lattice_layer_point_next ( int n, int c[], int v[], int *more );
    # void simplex_lattice_point_next ( int n, int c[], int v[], int *more );
    # int simplex_unit_lattice_point_num_nd ( int d, int s );
    # double simplex_unit_volume_nd ( int ndim );
    # double simplex_volume_nd ( int ndim, double a[] );
    # double sin_deg ( double angle );
    # double sin_power_int ( double a, double b, int n );
    # void soccer_shape_3d ( int point_num, int face_num, int face_order_max, double point_coord[], int face_order[], int face_point[] );
    # void soccer_size_3d ( int *point_num, int *edge_num, int *face_num, int *face_order_max );
    # void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
    # double sphere_cap_area_2d ( double r, double h );
    # double sphere_cap_area_3d ( double r, double h );
    # double sphere_cap_area_nd ( int ndim, double r, double h );
    # double sphere_cap_volume_2d ( double r, double h );
    # double sphere_cap_volume_3d ( double r, double h );
    # double sphere_cap_volume_nd ( int ndim, double r, double h );
    # void sphere_dia2imp_3d ( double p1[3], double p2[3], double *r, double pc[3] );
#################################################

    # double sphere_distance1 ( double lat1, double long1, double lat2, double long2, double radius );
    def sphere_distance1(self, lat1, lon1, lat2, lon2, r):
        return lib.sphere_distance1(c_double(lat1), c_double(lon1), c_double(lat2), c_double(lon2), c_double(r)) 

    # double sphere_distance2 ( double lat1, double long1, double lat2, double long2, double radius );
    def sphere_distance2(self, lat1, lon1, lat2, lon2, r):
        return lib.sphere_distance2(c_double(lat1), c_double(lon1), c_double(lat2), c_double(lon2), c_double(r))

    #double sphere_distance3 ( double lat1, double long1, double lat2, double long2, double radius );
    def sphere_distance3(self, lat1, lon1, lat2, lon2, r):
        return lib.sphere_distance3(c_double(lat1), c_double(lon1), c_double(lat2), c_double(lon2), c_double(r))

#################################################

    # int sphere_exp_contains_point_3d ( double p1[3], double p2[3], double p3[3], double p4[3], double p[3] );
    # void sphere_exp_point_near_3d ( double p1[3], double p2[3], double p3[3], double p4[3], double p[3], double pn[3] );
    # void sphere_exp2imp_3d ( double p1[3], double p2[3], double p3[3], double p4[3], double *r, double pc[3] );
    # void sphere_exp2imp_nd ( int n, double p[], double *r, double pc[] );
    # double sphere_imp_area_nd ( int n, double r );
    # int sphere_imp_contains_point_3d ( double r, double pc[3], double p[3] );
    # void sphere_imp_grid_icos_size ( int factor, int *node_num,  int *edge_num, int *triangle_num );
    # void sphere_imp_gridfaces_3d ( int maxtri, int nlat, int nlong, int *ntri, int tri[] );
    # void sphere_imp_gridlines_3d ( int line_max, int nlat, int nlong, int *line_num, int line[] );
    # void sphere_imp_gridpoints_3d ( double r, double pc[3], int maxpoint, int nlat, int nlong, int *npoint, double p[] );
    # void sphere_imp_gridpoints_icos1 ( int factor, int node_num, double node_xyz[] );
    # void sphere_imp_gridpoints_icos2 ( int factor, int node_num, double node_xyz[] );
    # int sphere_imp_line_project_3d ( double r, double pc[3], int n, double p[], int maxpnt2, double pp[], double thetamin, double thetamax );
    # void sphere_imp_local2xyz_3d ( double r, double pc[3], double theta, double phi, double p[3] );
    # void sphere_imp_point_near_3d ( double r, double pc[3], double p[3], double pn[3] );
    # void sphere_imp_point_project_3d ( double r, double pc[3], double p[3], double pp[3] );
    # void sphere_imp_spiralpoints_3d ( double r, double pc[3], int n, double p[] );
    # double sphere_imp_volume_3d ( double r );
    # double sphere_imp_volume_nd ( int n, double r );
    # double sphere_imp_zone_area_3d ( double r, double h1, double h2 );
    # double sphere_imp_zone_volume_3d ( double r, double h1, double h2 );
    # void sphere_imp2exp_3d ( double r, double pc[3], double p1[3], double p2[3], double p3[3], double p4[3] );
    # double sphere_k ( int n );
    # double sphere_triangle_angles_to_area ( double r, double a, double b, double c );
    # double sphere_triangle_contains_point ( double v1[3], double v2[3], double v3[3], double p[3] );
    # void sphere_triangle_sides_to_angles ( double r, double as, double bs, double cs, double *a, double *b, double *c );
    # void sphere_triangle_vertices_to_angles ( double r, double v1[3], double v2[3], double v3[3], double *a, double *b, double *c );
    # double sphere_triangle_vertices_to_area ( double r, double v1[3], double v2[3], double v3[3] );
    # void sphere_triangle_vertices_to_centroid ( double r, double v1[3], double v2[3], double v3[3], double vs[] );
    # int sphere_triangle_vertices_to_orientation ( double v1[3], double v2[3], double v3[3] );
    # void sphere_triangle_vertices_to_sides ( double r, double v1[3], double v2[3], double v3[3], double *as, double *bs, double *cs );
    # double sphere_unit_area_nd ( int n );
    # void sphere_unit_area_values ( int *n_data, int *n, double *area );
    # double *sphere_unit_sample_2d ( int *seed );
    # double *sphere_unit_sample_3d ( int *seed );
    # double *sphere_unit_sample_3d_2 ( int *seed );
    # double *sphere_unit_sample_nd ( int n, int *seed );
    # double *sphere_unit_sample_nd_2 ( int n, int *seed );
    # double *sphere_unit_sample_nd_3 ( int n, int *seed );
    # double sphere_unit_volume_nd ( int n );
    # void sphere_unit_volume_values ( int *n_data, int *n, double *volume );
    # double sphere01_distance_xyz ( double xyz1[3], double xyz2[3] );
    # double sphere01_polygon_area ( int n, double lat[], double lon[] );
    # double sphere01_triangle_angles_to_area ( double a, double b, double c );
    # void sphere01_triangle_sides_to_angles ( double as, double bs, double cs, double *a, double *b, double *c );
    # void sphere01_triangle_vertices_to_angles ( double v1[3], double v2[3], double v3[3], double *a, double *b, double *c );
    # double sphere01_triangle_vertices_to_area ( double v1[3], double v2[3], double v3[3] );
    # void sphere01_triangle_vertices_to_midpoints ( double v1[3], double v2[3], double v3[3], double m1[3], double m2[3], double m3[3] );
    # double *sphere01_triangle_vertices_to_centroid ( double v1[3], double v2[3], double v3[3] );
    # void sphere01_triangle_vertices_to_sides ( double v1[3], double v2[3], double v3[3], double *as, double *bs, double *cs );
    # void string_2d ( int vec_num, double p1[], double p2[], int *string_num, int order[], int string[] );
    # void super_ellipse_points_2d ( double pc[2], double r1, double r2, double expo, double psi, int n, double p[] );
    # double tan_deg ( double angle );
    # double *tetrahedron_barycentric_3d ( double tetra[3*4], double p[3] );
    # double *tetrahedron_centroid_3d ( double tetra[3*4] );
    # void tetrahedron_circumsphere_3d ( double tetra[3*4], double *r, double pc[3] );
    # int tetrahedron_contains_point_3d ( double tetra[3*4], double p[3] );
    # double *tetrahedron_dihedral_angles_3d ( double tetra[] );
    # double *tetrahedron_edge_length_3d ( double tetra[3*4] );
    # void tetrahedron_face_angles_3d ( double tetra[], double angles[] );
    # void tetrahedron_face_areas_3d ( double tetra[], double areas[] );
    # void tetrahedron_insphere_3d ( double tetra[3*4], double *r, double pc[3] );
    # void tetrahedron_lattice_layer_point_next ( int c[], int v[], int *more );
    # void tetrahedron_lattice_point_next ( int c[], int v[], int *more );
    # double tetrahedron_quality1_3d ( double tetra[3*4] );
    # double tetrahedron_quality2_3d ( double tetra[3*4] );
    # double tetrahedron_quality3_3d ( double tetra[3*4] );
    # double tetrahedron_quality4_3d ( double tetra[3*4] );
    # void tetrahedron_rhombic_shape_3d ( int point_num, int face_num, int face_order_max, double point_coord[], int face_order[], int face_point[] );
    # void tetrahedron_rhombic_size_3d ( int *point_num, int *edge_num, int *face_num, int *face_order_max );
    # void tetrahedron_sample_3d ( double tetra[3*4], int n, int *seed, double p[] );
    # void tetrahedron_shape_3d ( int point_num, int face_num, int face_order_max, double point_coord[], int face_order[], int face_point[] );
    # void tetrahedron_size_3d ( int *point_num, int *edge_num, int *face_num, int *face_order_max );
    # double *tetrahedron_solid_angles_3d ( double tetra[] );
    # int tetrahedron_unit_lattice_point_num_3d ( int s );
    # double tetrahedron_volume_3d ( double tetra[3*4] );
    # void theta2_adjust ( double *theta1, double *theta2 );
    # void theta3_adjust ( double *theta1, double *theta2, double *theta3 );
    # void timestamp ( void );
    # void tmat_init ( double a[4*4] );
    # void tmat_mxm ( double a[4*4], double b[4*4], double c[4*4] );
    # void tmat_mxp ( double a[4*4], double x[4], double y[4] );
    # void tmat_mxp2 ( double a[4*4], double p1[], double p2[], int n );
    # void tmat_mxv ( double a[4*4], double x[4], double y[4] );
    # void tmat_rot_axis ( double a[4*4], double b[4*4], double angle, char axis );
    # void tmat_rot_vector ( double a[4*4], double b[4*4], double angle, double v[3] );
    # void tmat_scale ( double a[4*4], double b[4*4], double s[3] );
    # void tmat_shear ( double a[4*4], double b[4*4], char *axis, double s );
    # void tmat_trans ( double a[4*4], double b[4*4], double v[3] );
    # double torus_volume_3d ( double r1, double r2 );
    # double *tp_to_xyz ( double theta, double phi );
    # void triangle_angles_2d ( double t[2*3], double angle[3] );
    # double *triangle_angles_2d_new ( double t[2*3] );
    # void triangle_angles_3d ( double t[3*3], double angle[3] );
    # double *triangle_angles_3d_new ( double t[3*3] );
    # double triangle_area_2d ( double t[2*3] );
    # double triangle_area_3d ( double t[3*3] );
    # double triangle_area_3d_2 ( double t[3*3] );
    # double triangle_area_3d_3 ( double t[3*3] );
    # double triangle_area_heron ( double s[3] );
    # double *triangle_area_vector_3d ( double t[3*3] );
    # double *triangle_barycentric_2d ( double t[2*3], double p[2] );
    # double *triangle_centroid_2d ( double t[2*3] );
    # double *triangle_centroid_3d ( double t[3*3] );
    # double *triangle_circumcenter_2d ( double t[2*3] );
    # double *triangle_circumcenter_2d_2 ( double t[2*3] );
    # double *triangle_circumcenter ( int n, double t[] );
    # void triangle_circumcircle_2d ( double t[2*3], double *r, double pc[2] );
    # void triangle_circumcircle_2d_2 ( double t[2*3], double *r, double pc[2] );
    # double triangle_circumradius_2d ( double t[2*3] );
    # void triangle_contains_line_exp_3d ( double t[3*3], double p1[3], double p2[3], int *inside, double pint[3] );
    # void triangle_contains_line_par_3d ( double t[], double p0[], double pd[], int *inside, double p[] );
    # int triangle_contains_point_2d_1 ( double t[2*3], double p[2] );
    # int triangle_contains_point_2d_2 ( double t[2*3], double p[2] );
    # int triangle_contains_point_2d_3 ( double t[2*3], double p[2] );
    # double triangle_diameter_2d ( double t[2*3] );
    # double *triangle_edge_length_2d ( double t[2*3] );
    # void triangle_gridpoints_2d ( double t[2*3], int sub_num, int grid_max, int *grid_num, double p[] );
    # void triangle_incenter_2d ( double t[2*3], double pc[2] );
    # void triangle_incircle_2d ( double t[2*3], double pc[2], double *r );
    # double triangle_inradius_2d ( double t[2*3] );
    # int triangle_is_degenerate_nd ( int dim_num, double t[] );
    # void triangle_lattice_layer_point_next ( int c[], int v[], int *more );
    # void triangle_lattice_point_next ( int c[], int v[], int *more );
    # void triangle_line_imp_int_2d ( double t[2*3], double a, double b, double c, int *int_num, double pint[] );
    # int triangle_orientation_2d ( double t[2*3] );
    # void triangle_orthocenter_2d ( double t[2*3], double p[2], int *flag );
    # double triangle_point_dist_2d ( double t[2*3], double p[2] );
    # double triangle_point_dist_3d ( double t[3*3], double p[3] );
    # double triangle_point_dist_signed_2d ( double t[2*3], double p[2] );
    # void triangle_point_near_2d ( double t[2*3], double p[2], double pn[2], double *dist );
    # double triangle_quality_2d ( double t[2*3] );
    # int triangle_right_lattice_point_num_2d ( int a, int b );
    # void triangle_sample ( double t[2*3], int n, int *seed, double p[] );
    # int triangle_unit_lattice_point_num_2d ( int s );
    # void triangle_xsi_to_xy_2d ( double t[2*3], double xsi[3], double p[2] );
    # void triangle_xy_to_xsi_2d ( double t[2*3], double p[2], double xsi[3] );
    # void truncated_octahedron_shape_3d ( int point_num, int face_num, int face_order_max, double point_coord[], int face_order[], int face_point[] );
    # void truncated_octahedron_size_3d ( int *point_num, int *edge_num, int *face_num, int *face_order_max );
    # void tube_2d ( double dist, int n, double p[], double p1[], double p2[] );
    # void tuple_next2 ( int n, int xmin[], int xmax[], int x[], int *rank );
    # void vector_directions_nd ( int dim_num, double v[], double angle[] );
    # void vector_rotate_2d ( double v1[2], double angle, double v2[2] );
    # void vector_rotate_3d ( double p1[3], double pa[3], double angle, double p2[3] );
    # void vector_rotate_base_2d ( double p1[2], double pb[2], double angle, double p2[2] );
    # double vector_separation_2d ( double v1[], double v2[] );
    # double vector_separation_3d ( double v1[], double v2[] );
    # double vector_separation_nd ( int n, double v1[], double v2[] );
    # void vector_unit_nd ( int n, double p[] );
    # int voxels_dist_l1_3d ( int v1[3], int v2[3] );
    # int voxels_dist_l1_nd ( int dim_num, int v1[], int v2[] );
    # void voxels_line_3d ( int p1[3], int p2[3], int n, int p[] );
    # void voxels_region_3d ( int maxlist, int nx, int ny, int nz, int ishow[], int *list_num, int list[], int *nregion );
    # void voxels_step_3d ( int v1[3], int v2[3], int inc, int jnc, int knc, int v3[3] );
    # void xy_to_polar ( double xy[2], double *r, double *t );
    # void xyz_to_radec ( double p[3], double *ra, double *dec );
    # void xyz_to_rtp ( double xyz[3], double *r, double *theta, double *phi );
    # void xyz_to_tp ( double xyz[3], double *theta, double *phi );

    ########## TESTS ###########
    def test_angle_box_2d(self):
        dist = 1
        p1 = [0,0]
        p2 =[3,0]
        p3 =[4,2]
        p4 = []
        p5 = []
        self.angle_box_2d ( dist, p1, p2, p3, p4, p5 )
        assert round(p4[0],6) == 2.381966
        assert round(p4[1],6) == 1.000000
        assert round(p5[0],6) == 3.618034
        assert round(p5[1],6) == -1.000000
        dist = 1
        p1 = [0,0]
        p2 =[3,0]
        p3 =[2,-2]
        p4 = []
        p5 = []
        self.angle_box_2d ( dist, p1, p2, p3, p4, p5 )
        assert round(p4[0],6) == 3.618034
        assert round(p4[1],6) == -1.000000
        assert round(p5[0],6) == 2.381966
        assert round(p5[1],6) == 1.000000
        dist = 1
        p1 = [3,0]
        p2 =[3,0]
        p3 =[2,-2]
        p4 = []
        p5 = []
        self.angle_box_2d ( dist, p1, p2, p3, p4, p5 )
        assert round(p4[0],6) == 2.105573
        assert round(p4[1],6) == 0.447214
        assert round(p5[0],6) == 3.894427
        assert round(p5[1],6) == -0.447214


    def test_angle_contains_ray_2d(self):
        p1 = [1,0]
        p2 =[0,0]
        p3 =[1,1]
        #( self, p1, p2, p3, p):
        assert self.angle_contains_ray_2d(p1, p2, p3, [1,0]) == True
        assert self.angle_contains_ray_2d(p1, p2, p3, [0.866025,0.5]) == True
        assert self.angle_contains_ray_2d(p1, p2, p3, [0.5,0.866025]) == False
        assert self.angle_contains_ray_2d(p1, p2, p3, [6.12323e-17,1]) == False
        assert self.angle_contains_ray_2d(p1, p2, p3, [-0.5,0.866025]) == False
        assert self.angle_contains_ray_2d(p1, p2, p3, [-0.866025,0.5]) == False
        assert self.angle_contains_ray_2d(p1, p2, p3, [-1,1.22465e-16]) == False
        assert self.angle_contains_ray_2d(p1, p2, p3, [-0.866025,-0.5]) == False
        assert self.angle_contains_ray_2d(p1, p2, p3, [-0.5,-0.866025]) == False
        assert self.angle_contains_ray_2d(p1, p2, p3, [-1.83697e-16,-1]) == False
        assert self.angle_contains_ray_2d(p1, p2, p3, [0.5,-0.866025]) == False
        assert self.angle_contains_ray_2d(p1, p2, p3, [0.866025,-0.5]) == False
        assert self.angle_contains_ray_2d(p1, p2, p3, [1,-2.44929e-16]) == False

    def test_angle_deg_2d (self):
        print "Warning: angle_deg_2d is untested"
        p1 = [-1,1]
        p2 =[0,0]
        p3 =[3,2]
        print self.angle_deg_2d(p1, p2, p3)

    def test_angle_half_2d (self):
        print "Warning: angle_half_2d is untested"
        p1 = [-1,1]
        p2 =[0,0]
        p3 =[3,2]
        print self.angle_half_2d(p1,p2,p3)

    def test_angle_rad_2d (self):
        print "Warning: angle_rad_2d is untested"
        p1 = [-1,1]
        p2 =[0,0]
        p3 =[3,2]
        print self.angle_rad_2d(p1, p2, p3)

    # double angle_rad_3d ( double p1[3], double p2[3], double p3[3] );
    def test_angle_rad_3d (self):
        print "Warning: angle_rad_3d is untested"
        p1 = [-1,5,10]
        p2 =[0,0,0]
        p3 =[3,2,1]
        print self.angle_rad_3d(p1,p2,p3)

    def test_angle_rad_nd(self):
        print "Warning: angle_rad_nd is untested"
        p1 = [-1,1,0]
        p2 =[0,0,0]
        print self.angle_rad_nd(3,p1,p2)

    # double angle_turn_2d ( double p1[2], double p2[2], double p3[2] );
    def test_angle_turn_2d (self):
        print "Warning: angle_turn_2d is untested"
        p1 = [-1,1]
        p2 =[0,0]
        p3 =[3,2]
        print self.angle_turn_2d(p1,p2,p3)

    # double anglei_deg_2d ( double p1[2], double p2[2], double p3[2] );
    def test_anglei_deg_2d (self):
        print "Warning: anglei_deg_2d is untested"
        p1 = [-1,1]
        p2 =[0,0]
        p3 =[3,2]
        print self.anglei_deg_2d(p1,p2,p3)

    # double anglei_rad_2d ( double p1[2], double p2[2], double p3[2] );
    def test_anglei_rad_2d (self):
        print "Warning: anglei_rad_2d is untested"
        p1 = [-1,1]
        p2 =[0,0]
        p3 =[3,2]
        print self.anglei_rad_2d(p1,p2,p3)

    def test_annulus_area_2d (self):
        print "Warning: annulus_area_2d is untested"
        print self.annulus_area_2d(0.1,0.5)

    def test_annulus_sector_area_2d(self ):
        print "Warning: annulus_sector_area_2d is untested"
        print self.annulus_sector_area_2d(0.1,0.5,0.3,0.4)


    def test_annulus_sector_centroid_2d (self):
        print "Warning: annulus_sector_centroid_2d is untested"
        pc = [1,1]
        print self.annulus_sector_centroid_2d(pc, 0.3, 0.5, 0.8, 0.9);

    def test_ball_unit_sample_2d (self):
        print "Warning: ball_unit_sample_2d is untested"
        pseed = 38
        print self.ball_unit_sample_2d(pseed)

    def test_ball_unit_sample_3d (self):
        print "Warning: ball_unit_sample_3d is untested"
        pseed = 38
        print self.ball_unit_sample_3d(pseed)

    def test_ball_unit_sample_nd (self):
        print "Warning: ball_unit_sample_nd is untested"
        seed = 38
        print self.ball_unit_sample_nd(4,seed)

    def test_basis_map_3d (self):
        print "Warning: basis_map_3d is untested"
        u = [1,2,3,4,5,6]
        v = [6,5,4,3,2,1]
        print self.basis_map_3d(u,v)

    def test_box_01_contains_point_2d (self):
        print "Warning: box_01_contains_point_2d is untested"
        p2 =[4,3]
        print self.box_01_contains_point_2d (p2)

    def test_box_01_contains_point_nd (self):
        print "Warning: box_01_contains_point_nd is untested"
        p2 = [3,4,5,6]
        print self.box_01_contains_point_nd (4,p2)

            

    def test_circle_ppr2imp_2d(self):
        print "Warning: circle_ppr2imp_2d is untested"
        p1 = [3,3]
        p2 = [5,5]
        r = 5
        print self.circle_ppr2imp_2d(p1, p2, r)

def testPerf(numtests):
    g = Geometry()

    toradians = 0.0174532925;

    trange = range(0,numtests)

    for t in trange:
        g.sphere_distance1(39.7391500 * toradians, -104.9847000 * toradians, 34.0522300 * toradians, -118.2436800 * toradians, 6371)
        g.polygon_contains_point_2d([0,0,0,4,4,4,4,0], [1,1])



## 2000 calls to a sphere_distance1
## and 2000 calls to polygon_contains_point_2d
## = 4000 calls total
#cProfile.run('testPerf(2000)')      # ~0.062 seconds on 64bit fedora

## 20000 calls to a sphere_distance1
## and 20000 calls to polygon_contains_point_2d
## = 40000 calls total
#cProfile.run('testPerf(20000)')     # ~0.603 seconds on 64bit fedora

## 30000 calls to a sphere_distance1
## and 30000 calls to polygon_contains_point_2d
## = 60000 calls total
#cProfile.run('testPerf(30000)')     # ~0.905 seconds on 64bit fedora

## 40000 calls to a sphere_distance1
## and 40000 calls to polygon_contains_point_2d
## = 80000 calls total
#cProfile.run('testPerf(40000)')     # ~1.198 seconds on 64bit fedora








def test():
    g = Geometry()
    g.test_angle_box_2d()
    g.test_angle_contains_ray_2d()
    g.test_angle_deg_2d()
    g.test_angle_half_2d()
    g.test_angle_rad_2d()
    g.test_angle_rad_3d()
    g.test_angle_rad_nd()
    g.test_angle_turn_2d()
    g.test_anglei_deg_2d()
    g.test_anglei_rad_2d()
    g.test_annulus_area_2d()
    g.test_annulus_sector_area_2d()
    g.test_annulus_sector_centroid_2d()
    g.test_ball_unit_sample_2d()
    g.test_ball_unit_sample_3d()
    g.test_ball_unit_sample_nd ()
    g.test_basis_map_3d ()
    g.test_box_01_contains_point_2d ()
    g.test_box_01_contains_point_nd ()

    g.test_circle_ppr2imp_2d()

test();

