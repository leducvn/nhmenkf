USRDIR 	= .
NAME = libropp.a
WKD = ${USRDIR}
MKF = ${USRDIR}/Makefile

include  ../Make3.inc
FFLAGS_ARCH = 
CPPFLAGS = 

OPT1_f90 = ${FFLAGS_ARCH}   
INC1_f90 = 
SRC1_f90 = ropp_pp_version.f90
OBJ1_f90 = ${SRC1_f90:.f90=.o}
VPT1_f90 = ${USRDIR}/Pp

OPT2_f90 = ${FFLAGS_ARCH}   
INC2_f90 = 
SRC2_f90 = ropp_1dvar_version.f90
OBJ2_f90 = ${SRC2_f90:.f90=.o}
VPT2_f90 = ${USRDIR}/Odvar

OPT3_f90 = ${FFLAGS_ARCH}   
INC3_f90 = 
SRC3_f90 = ropp_io_version.f90 gtshdrs.f90 bufrutils.f90
OBJ3_f90 = ${SRC3_f90:.f90=.o}
VPT3_f90 = ${USRDIR}/Io

OPT4_f90 = ${FFLAGS_ARCH}   
INC4_f90 = 
SRC4_f90 = ropp_fm_version.f90
OBJ4_f90 = ${SRC4_f90:.f90=.o}
VPT4_f90 = ${USRDIR}/Fm

OPT5_f90 = ${FFLAGS_ARCH}   
INC5_f90 = 
SRC5_f90 = outer_and.f90 datetimetypes.f90 messages.f90 typesizes.f90
OBJ5_f90 = ${SRC5_f90:.f90=.o}
VPT5_f90 = ${USRDIR}/Utils

OPT6_f90 = ${FFLAGS_ARCH}   
INC6_f90 = -I${USRDIR}/Utils
SRC6_f90 = ropp_pp_types.f90 ropp_pp_constants.f90 ropp_pp_cov_transform.f90 ropp_pp_spline.f90 ropp_pp_msis_types.f90 ropp_pp_fft.f90 ropp_pp_utils.f90 ropp_pp_interpol_log.f90 ropp_pp_monotonous.f90 ropp_pp_satellite_velocities.f90
OBJ6_f90 = ${SRC6_f90:.f90=.o}
VPT6_f90 = ${USRDIR}/Pp:${USRDIR}/Utils

OPT7_f90 = ${FFLAGS_ARCH}   
INC7_f90 = -I${USRDIR}/Utils
SRC7_f90 = ropp_1dvar_minropp.f90 matrix_full2bm_alloc.f90 matrix_full2bm.f90 matrix_bm2full_alloc.f90 matrix_bm2full.f90 matrix_svd.f90 matrix_types.f90 matrix_toast.f90
OBJ7_f90 = ${SRC7_f90:.f90=.o}
VPT7_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils

OPT8_f90 = ${FFLAGS_ARCH}   
INC8_f90 = -I${USRDIR}/Utils
SRC8_f90 = imaxloc.f90 message_get_addinfo.f90 monthofyear.f90 message_set_routine.f90 coordinates.f90 great_circle_distance.f90 message_get_routine.f90 message_get_program.f90 caltojul.f90 copy_alloc.f90 callocate.f90 earthmod.f90 message_set_addinfo.f90 swap.f90 unitconvert.f90 assert.f90 outer_product.f90 reverse.f90 message.f90 iminloc.f90 isinrange.f90 geodesy.f90 blend.f90 datetimeprogs.f90 arrays.f90 preallocate.f90 reallocate.f90 cross_product.f90 message_set_program.f90
OBJ8_f90 = ${SRC8_f90:.f90=.o}
VPT8_f90 = ${USRDIR}/Utils

OPT9_f90 = ${FFLAGS_ARCH}   
INC9_f90 = -I${USRDIR}/Utils
SRC9_f90 = ropp_fm_gpspderivs.f90 ropp_fm_interpol.f90 ropp_fm_asinh.f90 ropp_fm_gpspderivs_ad.f90 ropp_fm_interpol_log_tl.f90 ropp_fm_interpol_ad.f90 ropp_fm_interpol_log.f90 ropp_fm_interpol_tl.f90 ropp_fm_constants.f90 ropp_fm_interpol_log_ad.f90 ropp_fm_gpspderivs_tl.f90 ropp_fm_levels.f90
OBJ9_f90 = ${SRC9_f90:.f90=.o}
VPT9_f90 = ${USRDIR}/Fm:${USRDIR}/Utils

OPT10_f90 = ${FFLAGS_ARCH}   
INC10_f90 = -I${USRDIR}/Utils
SRC10_f90 = ropp_io_ecmwf.f90 ropp_io_thin_skip.f90
OBJ10_f90 = ${SRC10_f90:.f90=.o}
VPT10_f90 = ${USRDIR}/Io:${USRDIR}/Utils

OPT11_f90 = ${FFLAGS_ARCH}   
INC11_f90 = -I${USRDIR}/Utils -I${USRDIR}/Odvar
SRC11_f90 = ropp_fm_types.f90
OBJ11_f90 = ${SRC11_f90:.f90=.o}
VPT11_f90 = ${USRDIR}/Fm:${USRDIR}/Utils:${USRDIR}/Odvar

OPT12_f90 = ${FFLAGS_ARCH}   
INC12_f90 = -I${USRDIR}/Utils -I${USRDIR}/Fm
SRC12_f90 = ropp_fm_2d_plane.f90
OBJ12_f90 = ${SRC12_f90:.f90=.o}
VPT12_f90 = ${USRDIR}/Fm:${USRDIR}/Utils

OPT13_f90 = ${FFLAGS_ARCH}   
INC13_f90 = -I${USRDIR}/Utils -I${USRDIR}/Odvar
SRC13_f90 = ropp_control2state_ad.f90 ropp_control2state.f90 matrix_assign.f90 ropp_state2control.f90
OBJ13_f90 = ${SRC13_f90:.f90=.o}
VPT13_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils

OPT14_f90 = ${FFLAGS_ARCH}   
INC14_f90 = -I${USRDIR}/Utils
SRC14_f90 = gravity.f90 geod2cart.f90 vectors.f90 ropp_unit_conversion.f90 timesince.f90 plane_coordinates.f90 ecf2eci.f90 convertdt.f90 ut_convert.f90 r_eff.f90 nruns.f90 where.f90 datetimeoffset.f90 eci2ecf.f90 unit_vector.f90 date_and_time_utc.f90 ropp_utils.f90 l2norm.f90 setminus.f90 cart2geod.f90 impact_parameter.f90 sort.f90 occ_point.f90 tangent_point.f90 geopotential2geometric.f90 curvature.f90 getrun.f90 geometric2geopotential.f90 uniq.f90 locate.f90
OBJ14_f90 = ${SRC14_f90:.f90=.o}
VPT14_f90 = ${USRDIR}/Utils

OPT15_f90 = ${FFLAGS_ARCH}   
INC15_f90 = -I${USRDIR}/Odvar -I${USRDIR}/Utils
SRC15_f90 = matrix.f90 matrix_delete.f90
OBJ15_f90 = ${SRC15_f90:.f90=.o}
VPT15_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils

OPT16_f90 = ${FFLAGS_ARCH}   
INC16_f90 = -I${USRDIR}/Utils -I${USRDIR}/Pp
SRC16_f90 = ropp_pp_sliding_polynomial.f90 ropp_pp_filter.f90 ropp_pp_bending_angle_wo.f90 ropp_pp_invert_refraction.f90 ropp_pp_bg_refraction.f90 ropp_pp_refrac_msis.f90 ropp_pp_bangle_msis.f90 ropp_pp_smooth_profile.f90 ropp_pp_interpol.f90
OBJ16_f90 = ${SRC16_f90:.f90=.o}
VPT16_f90 = ${USRDIR}/Pp:${USRDIR}/Utils

OPT17_f90 = ${FFLAGS_ARCH}   
INC17_f90 = -I${USRDIR}/Utils -I${USRDIR}/Odvar -I${USRDIR}/Fm
SRC17_f90 = ropp_1dvar_types.f90
OBJ17_f90 = ${SRC17_f90:.f90=.o}
VPT17_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils:${USRDIR}/Fm

OPT18_f90 = ${FFLAGS_ARCH}   
INC18_f90 = -I${USRDIR}/Utils
SRC18_f90 = ropp_pp_amplitude_go.f90
OBJ18_f90 = ${SRC18_f90:.f90=.o}
VPT18_f90 = ${USRDIR}/Pp:${USRDIR}/Utils

OPT19_f90 = ${FFLAGS_ARCH}   
INC19_f90 = -I${USRDIR}/Utils -I${USRDIR}/Fm
SRC19_f90 = ropp_fm_alpha2drk_ad.f90 ropp_fm_abel_ad.f90 ropp_fm_alpha2drk_tl.f90 ropp_fm_alpha2drk.f90 ropp_fm_abel_tl.f90 ropp_fm_abel.f90
OBJ19_f90 = ${SRC19_f90:.f90=.o}
VPT19_f90 = ${USRDIR}/Fm:${USRDIR}/Utils

OPT20_f90 = ${FFLAGS_ARCH}   
INC20_f90 = -I${USRDIR}/Utils -I${USRDIR}/Odvar
SRC20_f90 = matrix_solve.f90 matrix_pp2full_subset.f90 matrix_sqrt.f90 matrix_pp2full_alloc.f90 matrix_full2pp.f90 matrix_pp2full.f90 matrix_full2pp_alloc.f90
OBJ20_f90 = ${SRC20_f90:.f90=.o}
VPT20_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils

OPT21_f90 = ${FFLAGS_ARCH}   
INC21_f90 = -I${USRDIR}/Utils -I${USRDIR}/Fm -I${USRDIR}/Odvar
SRC21_f90 = ropp_fm_state2state_meto_tl.f90 ropp_fm_state2state_meto.f90 ropp_fm_state2state_ecmwf.f90 ropp_fm_compress_2d_tl.f90 ropp_fm_obs2obs.f90 ropp_fm_state2state_ecmwf_ad.f90 ropp_fm_state2state.f90 ropp_fm_compress_single_ad.f90 ropp_fm_compress_single.f90 ropp_fm_free.f90 ropp_fm_state2state_meto_ad.f90 ropp_fm_state2state_ecmwf_tl.f90 ropp_fm_compress_2d.f90 ropp_fm_compress_2d_ad.f90 ropp_fm_compress_single_tl.f90
OBJ21_f90 = ${SRC21_f90:.f90=.o}
VPT21_f90 = ${USRDIR}/Fm:${USRDIR}/Utils:${USRDIR}/Odvar

OPT22_f90 = ${FFLAGS_ARCH}   
INC22_f90 = -I${USRDIR}/Fm -I${USRDIR}/Utils -I${USRDIR}/Odvar
SRC22_f90 = ropp_fm.f90
OBJ22_f90 = ${SRC22_f90:.f90=.o}
VPT22_f90 = ${USRDIR}/Fm:${USRDIR}/Utils:${USRDIR}/Odvar

OPT23_f90 = ${FFLAGS_ARCH}   
INC23_f90 = -I${USRDIR}/Utils
SRC23_f90 = FileDelete.f90 ToLower.f90 GetIOUnit.f90 ropp_utils_version.f90 ToUpper.f90
OBJ23_f90 = ${SRC23_f90:.f90=.o}
VPT23_f90 = ${USRDIR}/Utils

OPT24_f90 = ${FFLAGS_ARCH}   
INC24_f90 = -I${USRDIR}/Odvar -I${USRDIR}/Utils
SRC24_f90 = matrix_invert.f90 matrix_operators.f90
OBJ24_f90 = ${SRC24_f90:.f90=.o}
VPT24_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils

OPT25_f90 = ${FFLAGS_ARCH}   
INC25_f90 = -I${USRDIR}/Utils -I${USRDIR}/Pp
SRC25_f90 = ropp_pp_openloop.f90 ropp_pp_abel_exp.f90 ropp_pp_abel_lin.f90 ropp_pp_linear_combination.f90 ropp_pp_modelphase.f90 ropp_pp_invert_exp.f90 ropp_pp_invert_lin.f90 ropp_pp_fit_model_refraction.f90 ropp_pp_impact2doppler.f90
OBJ25_f90 = ${SRC25_f90:.f90=.o}
VPT25_f90 = ${USRDIR}/Pp:${USRDIR}/Utils

OPT26_f90 = ${FFLAGS_ARCH}   
INC26_f90 = -I${USRDIR}/Utils
SRC26_f90 = ropp_io_thin_sg.f90 ropp_io_types.f90 ropp_io_success.f90
OBJ26_f90 = ${SRC26_f90:.f90=.o}
VPT26_f90 = ${USRDIR}/Io:${USRDIR}/Utils

OPT27_f90 = ${FFLAGS_ARCH}   
INC27_f90 = -I${USRDIR}/Utils -I${USRDIR}/Pp -I${USRDIR}/Io
SRC27_f90 = ropp_pp.f90
OBJ27_f90 = ${SRC27_f90:.f90=.o}
VPT27_f90 = ${USRDIR}/Pp:${USRDIR}/Utils:${USRDIR}/Io

OPT28_f90 = ${FFLAGS_ARCH}   
INC28_f90 = -I${USRDIR}/Utils -I${USRDIR}/Odvar -I${USRDIR}/Fm
SRC28_f90 = ropp_qc_ion_bg.f90 ropp_qc_ion_bangle.f90 ropp_1dvar_diagnostics.f90
OBJ28_f90 = ${SRC28_f90:.f90=.o}
VPT28_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils:${USRDIR}/Fm

OPT29_f90 = ${FFLAGS_ARCH}   
INC29_f90 = -I${USRDIR}/Utils -I${USRDIR}/Fm -I${USRDIR}/Odvar
SRC29_f90 = ropp_qc_pge.f90 ropp_qc_bgqc.f90
OBJ29_f90 = ${SRC29_f90:.f90=.o}
VPT29_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils:${USRDIR}/Fm

OPT30_f90 = ${FFLAGS_ARCH}   
INC30_f90 = -I${USRDIR}/Utils -I${USRDIR}/Fm -I${USRDIR}/Odvar
SRC30_f90 = ropp_fm_compress.f90 ropp_fm_refrac_1d_new_tl.f90 ropp_fm_bangle_1d_grad.f90 ropp_fm_compress_tl.f90 ropp_fm_bangle_2d_ad.f90 ropp_fm_zorro.f90 ropp_fm_bangle_2d_tl.f90 ropp_fm_bangle_1d_tl.f90 ropp_fm_refrac_1d_ad.f90 ropp_fm_iono_bangle_ad.f90 ropp_fm_dzorro_dlg.f90 ropp_fm_compress_ad.f90 ropp_fm_refrac_1d.f90 ropp_fm_refrac_1d_grad.f90 ropp_fm_refrac_1d_new_ad.f90 ropp_fm_iono_bangle_tl.f90 ropp_fm_bangle_2d.f90 ropp_fm_iono_bangle.f90 ropp_fm_refrac_1d_tl.f90 ropp_fm_refrac_1d_new.f90
OBJ30_f90 = ${SRC30_f90:.f90=.o}
VPT30_f90 = ${USRDIR}/Fm:${USRDIR}/Utils:${USRDIR}/Odvar

OPT31_f90 = ${FFLAGS_ARCH}   
INC31_f90 = -I${USRDIR}/Fm -I${USRDIR}/Utils -I${USRDIR}/Odvar
SRC31_f90 = ropp_qc.f90
OBJ31_f90 = ${SRC31_f90:.f90=.o}
VPT31_f90 = ${USRDIR}/Odvar:${USRDIR}/Fm:${USRDIR}/Utils

OPT32_f90 = ${FFLAGS_ARCH}   
INC32_f90 = -I${USRDIR}/Utils -I${USRDIR}/Io
SRC32_f90 = ropp_io_shrink.f90 ropp_io_rangecheck.f90 ropp_io_isinrange.f90 ropp_io_thin.f90 ropp2bufr.f90 ropp_io_init.f90 ropp_io_assign.f90 ropp_io_fields_compare.f90 ropp_io_occid.f90 ropp_io_ascend.f90 ropp_io_free.f90 ropp_io_addvar.f90
OBJ32_f90 = ${SRC32_f90:.f90=.o}
VPT32_f90 = ${USRDIR}/Io:${USRDIR}/Utils

OPT33_f90 = ${FFLAGS_ARCH}   
INC33_f90 = -I${USRDIR}/Io -I${USRDIR}/Utils -I${USRDIR}/Odvar -I${USRDIR}/Fm
SRC33_f90 = ropp_1dvar_iono.f90
OBJ33_f90 = ${SRC33_f90:.f90=.o}
VPT33_f90 = ${USRDIR}/Odvar:${USRDIR}/Io:${USRDIR}/Utils:${USRDIR}/Fm

OPT34_f90 = ${FFLAGS_ARCH}   
INC34_f90 = -I${USRDIR}/Io -I${USRDIR}/Utils
SRC34_f90 = ropp_fm_iono.f90
OBJ34_f90 = ${SRC34_f90:.f90=.o}
VPT34_f90 = ${USRDIR}/Fm:${USRDIR}/Io:${USRDIR}/Utils

OPT35_f90 = ${FFLAGS_ARCH}   
INC35_f90 = -I${USRDIR}/Io -I${USRDIR}/Utils -I${USRDIR}/Fm -I${USRDIR}/Odvar
SRC35_f90 = ropp_fm_copy.f90
OBJ35_f90 = ${SRC35_f90:.f90=.o}
VPT35_f90 = ${USRDIR}/Fm:${USRDIR}/Io:${USRDIR}/Utils:${USRDIR}/Odvar

OPT36_f90 = ${FFLAGS_ARCH}   
INC36_f90 = -I${USRDIR}/Io -I${USRDIR}/Utils
SRC36_f90 = convertcodes.f90 ropp_io_vlist_size.f90 ropp_io.f90 ropp_io_test.f90 ropp_io_vlist_print.f90
OBJ36_f90 = ${SRC36_f90:.f90=.o}
VPT36_f90 = ${USRDIR}/Io:${USRDIR}/Utils

OPT37_f90 = ${FFLAGS_ARCH}   
INC37_f90 = -I${USRDIR}/Io -I${USRDIR}/Utils -I${USRDIR}/Pp
SRC37_f90 = ropp_pp_copy.f90 ropp_pp_preproc.f90
OBJ37_f90 = ${SRC37_f90:.f90=.o}
VPT37_f90 = ${USRDIR}/Pp:${USRDIR}/Io:${USRDIR}/Utils

OPT38_f90 = ${FFLAGS_ARCH}   
INC38_f90 = -I${USRDIR}/Odvar -I${USRDIR}/Utils -I${USRDIR}/Fm
SRC38_f90 = ropp_1dvar.f90
OBJ38_f90 = ${SRC38_f90:.f90=.o}
VPT38_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils:${USRDIR}/Fm

OPT39_f90 = ${FFLAGS_ARCH}   
INC39_f90 = -I${USRDIR}/Io -I${USRDIR}/Utils -I${USRDIR}/Fm -I${USRDIR}/Odvar
SRC39_f90 = ropp_1dvar_copy.f90
OBJ39_f90 = ${SRC39_f90:.f90=.o}
VPT39_f90 = ${USRDIR}/Odvar:${USRDIR}/Io:${USRDIR}/Utils:${USRDIR}/Fm

OPT40_f90 = ${FFLAGS_ARCH}   
INC40_f90 = -I${USRDIR}/Utils -I${USRDIR}/Io -I${USRDIR}/Fm -I${USRDIR}/Odvar
SRC40_f90 = ropp_1dvar_diag2roprof.f90
OBJ40_f90 = ${SRC40_f90:.f90=.o}
VPT40_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils:${USRDIR}/Io:${USRDIR}/Fm

OPT41_f90 = ${FFLAGS_ARCH}   
INC41_f90 = -I${USRDIR}/Utils -I${USRDIR}/Odvar -I${USRDIR}/Fm
SRC41_f90 = ropp_qc_omb.f90 ropp_qc_genqc.f90
OBJ41_f90 = ${SRC41_f90:.f90=.o}
VPT41_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils:${USRDIR}/Fm

OPT42_f90 = ${FFLAGS_ARCH}   
INC42_f90 = -I${USRDIR}/Utils -I${USRDIR}/Pp -I${USRDIR}/Io
SRC42_f90 = ropp_pp_dct.f90 ropp_pp_interpolate_trajectory.f90 ropp_pp_bangle2phase.f90 ropp_pp_correct_L2.f90 ropp_pp_fourier_filter.f90 ropp_pp_model_refraction.f90 ropp_pp_tdry.f90 ropp_pp_geometric_optics.f90 ropp_pp_merge_profile.f90 ropp_pp_geometric_optics_adj.f90 ropp_pp_radioholographic_filter.f90 ropp_pp_ionospheric_correction.f90 ropp_pp_search_model_refraction.f90 ropp_pp_bending_angle_go.f90
OBJ42_f90 = ${SRC42_f90:.f90=.o}
VPT42_f90 = ${USRDIR}/Pp:${USRDIR}/Utils:${USRDIR}/Io

OPT43_f90 = ${FFLAGS_ARCH}   
INC43_f90 = -I${USRDIR}/Utils -I${USRDIR}/Fm -I${USRDIR}/Odvar
SRC43_f90 = ropp_1dvar_cost.f90 ropp_1dvar_solve.f90 ropp_1dvar_levmarq.f90 ropp_qc_cutoff.f90
OBJ43_f90 = ${SRC43_f90:.f90=.o}
VPT43_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils:${USRDIR}/Fm

OPT44_f90 = ${FFLAGS_ARCH}   
INC44_f90 = -I${USRDIR}/Pp -I${USRDIR}/Utils -I${USRDIR}/Io
SRC44_f90 = ropp_pp_read_config.f90
OBJ44_f90 = ${SRC44_f90:.f90=.o}
VPT44_f90 = ${USRDIR}/Pp:${USRDIR}/Utils:${USRDIR}/Io

OPT45_f90 = ${FFLAGS_ARCH}   
INC45_f90 = -I${USRDIR}/Utils -I${USRDIR}/Io -I${USRDIR}/Fm
SRC45_f90 = ropp_fm_iono_unpack_bangle.f90
OBJ45_f90 = ${SRC45_f90:.f90=.o}
VPT45_f90 = ${USRDIR}/Fm:${USRDIR}/Utils:${USRDIR}/Io

OPT46_f90 = ${FFLAGS_ARCH}   
INC46_f90 = -I${USRDIR}/Utils -I${USRDIR}/Io
SRC46_f90 = ropp_fm_set_units.f90
OBJ46_f90 = ${SRC46_f90:.f90=.o}
VPT46_f90 = ${USRDIR}/Fm:${USRDIR}/Utils:${USRDIR}/Io

OPT47_f90 = ${FFLAGS_ARCH}   
INC47_f90 = -I${USRDIR}/Utils -I${USRDIR}/Io -I${USRDIR}/Pp
SRC47_f90 = ropp_pp_preprocess_grasrs.f90 ropp_pp_tph_tdry.f90 ropp_pp_spectra.f90 ropp_pp_preprocess_cosmic.f90 ropp_pp_preprocess.f90 ropp_pp_cutoff.f90 ropp_pp_tph_bangle.f90 ropp_pp_tph_refrac.f90 ropp_pp_radiooptic_analysis.f90 ropp_pp_tph_temp.f90 ropp_pp_diag2roprof.f90 ropp_pp_cutoff_amplitude.f90
OBJ47_f90 = ${SRC47_f90:.f90=.o}
VPT47_f90 = ${USRDIR}/Pp:${USRDIR}/Utils:${USRDIR}/Io

OPT48_f90 = ${FFLAGS_ARCH}   
INC48_f90 = -I${USRDIR}/Utils -I${USRDIR}/Io
SRC48_f90 = ropp_io_thin_fixed.f90 bufr2ropp.f90 ropp_io_thin_select.f90
OBJ48_f90 = ${SRC48_f90:.f90=.o}
VPT48_f90 = ${USRDIR}/Io:${USRDIR}/Utils

OPT49_f90 = ${FFLAGS_ARCH}   
INC49_f90 = -I${USRDIR}/Utils -I${USRDIR}/Io -I${USRDIR}/Fm -I${USRDIR}/Odvar
SRC49_f90 = ropp_fm_roprof2obs.f90 ropp_fm_iono_set_default.f90 ropp_fm_obs2roprof.f90 ropp_fm_state2roprof.f90 ropp_fm_roprof2state.f90
OBJ49_f90 = ${SRC49_f90:.f90=.o}
VPT49_f90 = ${USRDIR}/Fm:${USRDIR}/Utils:${USRDIR}/Io:${USRDIR}/Odvar

OPT50_f90 = ${FFLAGS_ARCH}   
INC50_f90 = -I${USRDIR}/Utils -I${USRDIR}/Fm -I${USRDIR}/Odvar -I${USRDIR}/Io
SRC50_f90 = ropp_fm_bangle_1d_ad.f90 ropp_fm_bangle_1d.f90
OBJ50_f90 = ${SRC50_f90:.f90=.o}
VPT50_f90 = ${USRDIR}/Fm:${USRDIR}/Utils:${USRDIR}/Odvar:${USRDIR}/Io

OPT51_f90 = ${FFLAGS_ARCH}   
INC51_f90 = -I${USRDIR}/Utils -I${USRDIR}/Io -I${USRDIR}/Fm
SRC51_f90 = ropp_1dvar_iono_unpack_bangle.f90
OBJ51_f90 = ${SRC51_f90:.f90=.o}
VPT51_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils:${USRDIR}/Io:${USRDIR}/Fm

OPT52_f90 = ${FFLAGS_ARCH}   
INC52_f90 = -I${USRDIR}/Utils -I${USRDIR}/Odvar -I${USRDIR}/Io -I${USRDIR}/Fm
SRC52_f90 = ropp_1dvar_iono_repack_bg.f90 ropp_1dvar_iono_repack_bangle.f90 ropp_1dvar_covar_bg.f90 ropp_1dvar_covar_bangle.f90 ropp_1dvar_covar_refrac.f90
OBJ52_f90 = ${SRC52_f90:.f90=.o}
VPT52_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils:${USRDIR}/Io:${USRDIR}/Fm

OPT53_f90 = ${FFLAGS_ARCH}   
INC53_f90 = -I${USRDIR}/Utils -I${USRDIR}/Io
SRC53_f90 = ropp_pp_set_coordinates.f90
OBJ53_f90 = ${SRC53_f90:.f90=.o}
VPT53_f90 = ${USRDIR}/Pp:${USRDIR}/Utils:${USRDIR}/Io

OPT54_f90 = ${FFLAGS_ARCH}   
INC54_f90 = -I${USRDIR}/Odvar -I${USRDIR}/Utils -I${USRDIR}/Fm
SRC54_f90 = ropp_1dvar_read_config.f90
OBJ54_f90 = ${SRC54_f90:.f90=.o}
VPT54_f90 = ${USRDIR}/Odvar:${USRDIR}/Utils:${USRDIR}/Fm

OBJ_ALL = ${OBJ1_f90} ${OBJ2_f90} ${OBJ3_f90} ${OBJ4_f90} ${OBJ5_f90} ${OBJ6_f90} ${OBJ7_f90} ${OBJ8_f90} ${OBJ9_f90} ${OBJ10_f90} ${OBJ11_f90} ${OBJ12_f90} ${OBJ13_f90} ${OBJ14_f90} ${OBJ15_f90} ${OBJ16_f90} ${OBJ17_f90} ${OBJ18_f90} ${OBJ19_f90} ${OBJ20_f90} ${OBJ21_f90} ${OBJ22_f90} ${OBJ23_f90} ${OBJ24_f90} ${OBJ25_f90} ${OBJ26_f90} ${OBJ27_f90} ${OBJ28_f90} ${OBJ29_f90} ${OBJ30_f90} ${OBJ31_f90} ${OBJ32_f90} ${OBJ33_f90} ${OBJ34_f90} ${OBJ35_f90} ${OBJ36_f90} ${OBJ37_f90} ${OBJ38_f90} ${OBJ39_f90} ${OBJ40_f90} ${OBJ41_f90} ${OBJ42_f90} ${OBJ43_f90} ${OBJ44_f90} ${OBJ45_f90} ${OBJ46_f90} ${OBJ47_f90} ${OBJ48_f90} ${OBJ49_f90} ${OBJ50_f90} ${OBJ51_f90} ${OBJ52_f90} ${OBJ53_f90} ${OBJ54_f90}

ropp_1dvar_minropp.o : typesizes.o
matrix_full2bm_alloc.o : typesizes.o messages.o
ropp_fm_gpspderivs.o : typesizes.o
imaxloc.o : typesizes.o
message_get_addinfo.o : messages.o
ropp_fm_interpol.o : typesizes.o
ropp_fm_asinh.o : typesizes.o
monthofyear.o : datetimetypes.o
message_set_routine.o : messages.o
ropp_pp_types.o : typesizes.o
coordinates.o : typesizes.o
great_circle_distance.o : typesizes.o
ropp_io_ecmwf.o : typesizes.o
message_get_routine.o : messages.o
message_get_program.o : messages.o
matrix_full2bm.o : typesizes.o messages.o
caltojul.o : datetimetypes.o
ropp_pp_constants.o : typesizes.o
copy_alloc.o : typesizes.o
ropp_pp_cov_transform.o : typesizes.o
ropp_fm_gpspderivs_ad.o : typesizes.o
callocate.o : typesizes.o
matrix_bm2full_alloc.o : typesizes.o messages.o
earthmod.o : typesizes.o messages.o
ropp_fm_interpol_log_tl.o : typesizes.o
ropp_fm_interpol_ad.o : typesizes.o
message_set_addinfo.o : messages.o
swap.o : typesizes.o
unitconvert.o : typesizes.o
assert.o : messages.o
outer_product.o : typesizes.o
ropp_fm_interpol_log.o : typesizes.o
ropp_io_thin_skip.o : messages.o
reverse.o : typesizes.o
ropp_fm_interpol_tl.o : typesizes.o
ropp_pp_spline.o : typesizes.o
message.o : messages.o
iminloc.o : typesizes.o
isinrange.o : typesizes.o
ropp_pp_msis_types.o : typesizes.o
ropp_fm_constants.o : typesizes.o
geodesy.o : typesizes.o
blend.o : typesizes.o
ropp_pp_fft.o : typesizes.o
matrix_bm2full.o : typesizes.o messages.o
ropp_pp_utils.o : typesizes.o
datetimeprogs.o : datetimetypes.o
ropp_pp_interpol_log.o : typesizes.o
ropp_fm_interpol_log_ad.o : typesizes.o
matrix_svd.o : typesizes.o
ropp_pp_monotonous.o : typesizes.o
arrays.o : typesizes.o
preallocate.o : typesizes.o
reallocate.o : typesizes.o
cross_product.o : typesizes.o
matrix_types.o : typesizes.o
message_set_program.o : messages.o
ropp_fm_gpspderivs_tl.o : typesizes.o
ropp_pp_satellite_velocities.o : typesizes.o
ropp_fm_levels.o : typesizes.o
matrix_toast.o : typesizes.o
gravity.o : typesizes.o geodesy.o
geod2cart.o : typesizes.o coordinates.o
vectors.o : typesizes.o coordinates.o
ropp_unit_conversion.o : typesizes.o unitconvert.o
ropp_fm_2d_plane.o : typesizes.o ropp_fm_constants.o
timesince.o : datetimetypes.o datetimeprogs.o
plane_coordinates.o : typesizes.o coordinates.o
ecf2eci.o : typesizes.o datetimeprogs.o coordinates.o
ropp_pp_sliding_polynomial.o : typesizes.o ropp_pp_utils.o
convertdt.o : datetimeprogs.o
ut_convert.o : typesizes.o unitconvert.o
r_eff.o : typesizes.o geodesy.o
nruns.o : arrays.o
where.o : arrays.o
datetimeoffset.o : datetimetypes.o datetimeprogs.o
ropp_pp_filter.o : typesizes.o ropp_pp_utils.o
ropp_control2state_ad.o : typesizes.o matrix_types.o
ropp_control2state.o : typesizes.o matrix_types.o
eci2ecf.o : typesizes.o datetimeprogs.o coordinates.o
matrix.o : matrix_types.o typesizes.o
unit_vector.o : typesizes.o arrays.o
date_and_time_utc.o : datetimetypes.o datetimeprogs.o
ropp_utils.o : typesizes.o arrays.o geodesy.o coordinates.o datetimeprogs.o messages.o earthmod.o unitconvert.o
l2norm.o : typesizes.o arrays.o
ropp_fm_types.o : typesizes.o matrix_types.o
matrix_assign.o : typesizes.o matrix_types.o arrays.o messages.o
ropp_pp_bending_angle_wo.o : typesizes.o ropp_pp_types.o messages.o
ropp_pp_invert_refraction.o : typesizes.o geodesy.o ropp_pp_types.o
ropp_pp_bg_refraction.o : typesizes.o ropp_pp_types.o
setminus.o : arrays.o typesizes.o
ropp_pp_refrac_msis.o : typesizes.o ropp_pp_msis_types.o
ropp_pp_bangle_msis.o : typesizes.o ropp_pp_msis_types.o
cart2geod.o : typesizes.o coordinates.o
impact_parameter.o : typesizes.o coordinates.o
sort.o : typesizes.o arrays.o
matrix_delete.o : matrix_types.o
occ_point.o : typesizes.o coordinates.o earthmod.o
tangent_point.o : typesizes.o coordinates.o earthmod.o
ropp_pp_smooth_profile.o : typesizes.o ropp_pp_utils.o ropp_pp_types.o
ropp_state2control.o : typesizes.o matrix_types.o
geopotential2geometric.o : typesizes.o geodesy.o
curvature.o : typesizes.o coordinates.o
getrun.o : arrays.o
geometric2geopotential.o : typesizes.o geodesy.o
uniq.o : typesizes.o arrays.o
locate.o : typesizes.o arrays.o
ropp_pp_interpol.o : typesizes.o ropp_pp_spline.o
ropp_fm_alpha2drk_ad.o : typesizes.o ropp_utils.o ropp_fm_constants.o
ropp_fm_abel_ad.o : typesizes.o ropp_fm_constants.o ropp_utils.o
FileDelete.o : ropp_utils.o
ropp_fm_state2state_meto_tl.o : typesizes.o ropp_fm_constants.o ropp_fm_types.o
ropp_io_thin_sg.o : typesizes.o messages.o ropp_utils.o
ropp_fm_state2state_meto.o : typesizes.o ropp_fm_constants.o ropp_fm_types.o
ropp_pp_openloop.o : typesizes.o ropp_pp_constants.o ropp_utils.o messages.o
ropp_fm_alpha2drk_tl.o : typesizes.o ropp_utils.o ropp_fm_constants.o
ropp_fm.o : ropp_fm_types.o ropp_fm_constants.o typesizes.o
matrix_solve.o : typesizes.o messages.o matrix_types.o ropp_utils.o
ropp_pp_abel_exp.o : typesizes.o ropp_utils.o ropp_pp_constants.o
ropp_fm_state2state_ecmwf.o : typesizes.o ropp_fm_types.o ropp_fm_constants.o
ropp_pp_abel_lin.o : typesizes.o ropp_utils.o ropp_pp_constants.o
ropp_pp_linear_combination.o : typesizes.o ropp_utils.o ropp_pp_constants.o
ropp_pp_modelphase.o : typesizes.o ropp_pp_types.o ropp_pp_constants.o ropp_utils.o
ropp_fm_compress_2d_tl.o : typesizes.o ropp_fm_types.o ropp_fm_constants.o
ToLower.o : ropp_utils.o
ropp_pp_invert_exp.o : typesizes.o ropp_utils.o ropp_pp_constants.o
ropp_pp_amplitude_go.o : typesizes.o ropp_utils.o
ropp_fm_alpha2drk.o : ropp_utils.o typesizes.o ropp_fm_constants.o
ropp_fm_obs2obs.o : arrays.o ropp_fm_types.o
GetIOUnit.o : ropp_utils.o
ropp_fm_abel_tl.o : typesizes.o ropp_utils.o ropp_fm_constants.o
matrix_invert.o : matrix.o typesizes.o messages.o
ropp_fm_state2state_ecmwf_ad.o : typesizes.o ropp_fm_constants.o ropp_fm_types.o
matrix_operators.o : matrix_types.o messages.o typesizes.o matrix.o
matrix_pp2full_subset.o : typesizes.o messages.o ropp_utils.o matrix_types.o
ropp_utils_version.o : ropp_utils.o
ropp_fm_state2state.o : arrays.o ropp_fm_types.o
ropp_pp_invert_lin.o : typesizes.o ropp_utils.o ropp_pp_constants.o
matrix_sqrt.o : typesizes.o messages.o matrix.o
ropp_fm_compress_single_ad.o : typesizes.o ropp_fm_types.o ropp_fm_constants.o
ropp_1dvar_types.o : typesizes.o matrix_types.o ropp_fm_types.o
matrix_pp2full_alloc.o : typesizes.o messages.o ropp_utils.o matrix_types.o
ropp_fm_compress_single.o : typesizes.o ropp_fm_types.o ropp_fm_constants.o
ropp_fm_free.o : ropp_utils.o ropp_fm_types.o
ropp_fm_state2state_meto_ad.o : typesizes.o ropp_fm_constants.o ropp_fm_types.o
ropp_pp_fit_model_refraction.o : typesizes.o messages.o ropp_utils.o ropp_pp_utils.o ropp_pp_types.o
ropp_io_types.o : typesizes.o ropp_utils.o
ropp_fm_state2state_ecmwf_tl.o : typesizes.o ropp_fm_types.o ropp_fm_constants.o
ropp_pp_impact2doppler.o : typesizes.o ropp_pp_constants.o ropp_utils.o
matrix_full2pp.o : typesizes.o messages.o ropp_utils.o matrix_types.o
ToUpper.o : ropp_utils.o
matrix_pp2full.o : typesizes.o messages.o ropp_utils.o matrix_types.o
ropp_io_success.o : ropp_utils.o
ropp_fm_compress_2d.o : typesizes.o ropp_fm_types.o ropp_fm_constants.o
ropp_fm_abel.o : typesizes.o ropp_utils.o ropp_fm_constants.o
matrix_full2pp_alloc.o : typesizes.o messages.o ropp_utils.o matrix_types.o matrix.o
ropp_fm_compress_2d_ad.o : typesizes.o ropp_fm_types.o ropp_fm_constants.o
ropp_fm_compress_single_tl.o : typesizes.o ropp_fm_types.o ropp_fm_constants.o
ropp_fm_compress.o : typesizes.o ropp_fm.o ropp_fm_types.o ropp_fm_constants.o
ropp_1dvar.o : ropp_1dvar_types.o ropp_fm_types.o typesizes.o matrix_types.o
ropp_fm_refrac_1d_new_tl.o : typesizes.o ropp_utils.o ropp_fm.o ropp_fm_types.o ropp_fm_constants.o
ropp_fm_bangle_1d_grad.o : typesizes.o ropp_fm.o ropp_fm_types.o messages.o
ropp_fm_compress_tl.o : typesizes.o ropp_fm.o ropp_fm_types.o ropp_fm_constants.o
ropp_fm_iono.o : ropp_io_types.o
ropp_io_shrink.o : ropp_utils.o ropp_io_types.o
ropp_fm_bangle_2d_ad.o : typesizes.o ropp_fm.o geodesy.o ropp_fm_types.o ropp_fm_constants.o
ropp_qc_ion_bg.o : typesizes.o ropp_utils.o ropp_1dvar_types.o
ropp_fm_zorro.o : typesizes.o ropp_fm.o ropp_fm_constants.o
ropp_qc_pge.o : typesizes.o arrays.o messages.o ropp_utils.o ropp_fm_types.o ropp_fm_constants.o ropp_1dvar_types.o
ropp_fm_bangle_2d_tl.o : typesizes.o ropp_fm.o geodesy.o ropp_fm_types.o ropp_fm_constants.o
ropp_io_rangecheck.o : typesizes.o ropp_io_types.o ropp_utils.o geodesy.o
ropp_fm_bangle_1d_tl.o : typesizes.o ropp_fm.o ropp_utils.o ropp_fm_types.o ropp_fm_constants.o
ropp_io_isinrange.o : ropp_utils.o ropp_io_types.o
ropp_fm_refrac_1d_ad.o : typesizes.o ropp_fm.o ropp_fm_types.o ropp_fm_constants.o
ropp_io_thin.o : messages.o typesizes.o ropp_io_types.o geodesy.o ropp_utils.o
ropp_fm_iono_bangle_ad.o : typesizes.o ropp_fm.o ropp_utils.o ropp_fm_constants.o
ropp_pp_copy.o : ropp_io_types.o ropp_pp_types.o
ropp_qc.o : ropp_fm_types.o ropp_1dvar_types.o
ropp_qc_ion_bangle.o : typesizes.o ropp_utils.o ropp_1dvar_types.o
ropp2bufr.o : typesizes.o messages.o datetimetypes.o gtshdrs.o ropp_io_types.o
ropp_fm_dzorro_dlg.o : typesizes.o ropp_fm.o ropp_fm_constants.o
ropp_fm_compress_ad.o : typesizes.o ropp_fm.o ropp_fm_types.o ropp_fm_constants.o
convertcodes.o : ropp_io_types.o messages.o ropp_utils.o
ropp_io_vlist_size.o : ropp_io_types.o
ropp_1dvar_diagnostics.o : typesizes.o arrays.o messages.o ropp_utils.o matrix.o ropp_fm.o ropp_1dvar_types.o
ropp_fm_refrac_1d.o : typesizes.o ropp_utils.o ropp_fm.o ropp_fm_types.o ropp_fm_constants.o
ropp_fm_refrac_1d_grad.o : typesizes.o ropp_fm.o ropp_fm_types.o messages.o
ropp_io.o : ropp_io_types.o typesizes.o
ropp_io_init.o : ropp_utils.o ropp_io_types.o datetimeprogs.o datetimetypes.o
ropp_pp.o : typesizes.o ropp_pp_types.o ropp_pp_constants.o ropp_pp_utils.o ropp_pp_spline.o ropp_pp_msis_types.o ropp_io_types.o
ropp_io_assign.o : ropp_utils.o ropp_io_types.o
ropp_1dvar_copy.o : ropp_io_types.o ropp_fm_types.o ropp_1dvar_types.o
ropp_io_test.o : ropp_io_types.o
ropp_1dvar_iono.o : ropp_io_types.o ropp_1dvar_types.o
ropp_io_fields_compare.o : typesizes.o ropp_io_types.o messages.o
ropp_fm_refrac_1d_new_ad.o : typesizes.o ropp_utils.o ropp_fm.o ropp_fm_types.o ropp_fm_constants.o
ropp_io_vlist_print.o : ropp_io_types.o
ropp_fm_iono_bangle_tl.o : typesizes.o ropp_fm.o ropp_utils.o ropp_fm_constants.o
ropp_fm_copy.o : ropp_io_types.o ropp_fm_types.o
ropp_fm_bangle_2d.o : typesizes.o ropp_fm.o geodesy.o ropp_fm_types.o ropp_fm_constants.o
ropp_io_occid.o : ropp_utils.o ropp_io_types.o
ropp_io_ascend.o : messages.o arrays.o ropp_utils.o ropp_io_types.o
ropp_io_free.o : typesizes.o ropp_utils.o ropp_io_types.o
ropp_qc_bgqc.o : typesizes.o arrays.o messages.o ropp_utils.o ropp_fm_types.o ropp_1dvar_types.o
ropp_fm_iono_bangle.o : typesizes.o ropp_fm.o ropp_utils.o ropp_fm_constants.o
ropp_fm_refrac_1d_tl.o : typesizes.o ropp_fm.o ropp_fm_types.o ropp_fm_constants.o
ropp_io_addvar.o : typesizes.o ropp_io_types.o
ropp_pp_preproc.o : ropp_io_types.o ropp_pp_types.o typesizes.o
ropp_fm_refrac_1d_new.o : typesizes.o ropp_utils.o ropp_fm.o ropp_fm_types.o ropp_fm_constants.o
ropp_fm_roprof2obs.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_fm.o ropp_fm_types.o geodesy.o
ropp_pp_preprocess_grasrs.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_pp.o ropp_pp_spline.o ropp_pp_constants.o ropp_pp_types.o
ropp_fm_iono_set_default.o : typesizes.o ropp_utils.o ropp_io_types.o ropp_fm.o ropp_fm_copy.o
ropp_pp_dct.o : typesizes.o ropp_pp.o ropp_pp_constants.o ropp_pp_spline.o ropp_pp_utils.o ropp_pp_types.o ropp_utils.o messages.o
ropp_pp_interpolate_trajectory.o : typesizes.o ropp_pp.o ropp_utils.o
ropp_pp_set_coordinates.o : ropp_utils.o ropp_io.o ropp_io_types.o
ropp_pp_tph_tdry.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_pp.o ropp_pp_constants.o
ropp_pp_bangle2phase.o : typesizes.o ropp_pp.o ropp_utils.o
ropp_pp_spectra.o : typesizes.o ropp_utils.o ropp_io_types.o ropp_io.o ropp_pp_types.o ropp_pp_constants.o
ropp_pp_correct_L2.o : typesizes.o messages.o ropp_pp.o ropp_pp_constants.o
ropp_pp_fourier_filter.o : typesizes.o ropp_pp.o ropp_pp_constants.o
ropp_pp_model_refraction.o : typesizes.o ropp_pp.o ropp_pp_types.o
ropp_pp_preprocess_cosmic.o : typesizes.o ropp_utils.o ropp_io_types.o ropp_pp.o ropp_pp_types.o
ropp_qc_omb.o : typesizes.o matrix.o messages.o ropp_utils.o ropp_fm.o ropp_1dvar.o
ropp_1dvar_cost.o : typesizes.o messages.o arrays.o ropp_fm.o ropp_1dvar.o matrix.o ropp_fm_types.o
ropp_1dvar_iono_repack_bg.o : typesizes.o matrix.o ropp_utils.o ropp_io_types.o ropp_fm.o ropp_fm_copy.o ropp_fm_iono.o ropp_1dvar_types.o
ropp_pp_tdry.o : typesizes.o ropp_utils.o ropp_pp.o ropp_pp_constants.o
ropp_io_thin_fixed.o : typesizes.o messages.o ropp_io.o ropp_utils.o
ropp_pp_preprocess.o : typesizes.o ropp_io_types.o ropp_pp.o ropp_pp_types.o ropp_pp_constants.o ropp_utils.o
ropp_pp_cutoff.o : typesizes.o ropp_utils.o messages.o ropp_io_types.o ropp_io.o ropp_pp.o ropp_pp_types.o
ropp_1dvar_solve.o : typesizes.o ropp_utils.o ropp_fm.o ropp_1dvar.o matrix.o
ropp_pp_tph_bangle.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_pp.o
ropp_pp_tph_refrac.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_pp.o
ropp_pp_geometric_optics.o : typesizes.o ropp_pp.o ropp_utils.o
ropp_pp_merge_profile.o : typesizes.o ropp_pp.o
ropp_pp_geometric_optics_adj.o : typesizes.o ropp_pp.o ropp_utils.o
ropp_pp_radiooptic_analysis.o : typesizes.o ropp_utils.o ropp_io_types.o ropp_io.o messages.o ropp_pp.o
bufr2ropp.o : messages.o typesizes.o ropp_io_types.o ropp_io.o geodesy.o
ropp_pp_tph_temp.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_pp.o ropp_pp_constants.o
ropp_1dvar_iono_unpack_bangle.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_fm_constants.o
ropp_1dvar_diag2roprof.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_fm.o ropp_1dvar.o
ropp_pp_radioholographic_filter.o : typesizes.o ropp_pp.o ropp_utils.o
ropp_1dvar_iono_repack_bangle.o : typesizes.o matrix.o ropp_utils.o ropp_io_types.o ropp_fm.o ropp_fm_copy.o ropp_1dvar_types.o
ropp_pp_diag2roprof.o : typesizes.o ropp_io.o ropp_io_types.o ropp_pp_types.o
ropp_1dvar_levmarq.o : typesizes.o ropp_utils.o ropp_fm.o ropp_1dvar.o matrix.o
ropp_fm_set_units.o : ropp_utils.o ropp_io.o ropp_io_types.o typesizes.o unitconvert.o
ropp_fm_obs2roprof.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_fm.o ropp_fm_types.o
ropp_fm_bangle_1d_ad.o : typesizes.o ropp_utils.o ropp_fm.o ropp_fm_types.o ropp_fm_iono.o ropp_fm_constants.o geodesy.o
ropp_fm_iono_unpack_bangle.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_fm_constants.o
ropp_pp_ionospheric_correction.o : typesizes.o ropp_pp.o ropp_pp_types.o
ropp_io_thin_select.o : typesizes.o messages.o ropp_io.o ropp_utils.o
ropp_fm_bangle_1d.o : typesizes.o ropp_utils.o ropp_fm.o ropp_fm_types.o ropp_fm_iono.o ropp_fm_constants.o geodesy.o
ropp_1dvar_covar_bg.o : typesizes.o matrix.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_fm.o ropp_1dvar.o
ropp_1dvar_read_config.o : ropp_1dvar.o messages.o ropp_utils.o arrays.o
ropp_pp_cutoff_amplitude.o : typesizes.o ropp_utils.o messages.o ropp_io_types.o ropp_io.o ropp_pp.o ropp_pp_types.o
ropp_pp_search_model_refraction.o : typesizes.o messages.o ropp_utils.o ropp_pp_utils.o ropp_pp.o ropp_pp_types.o
ropp_1dvar_covar_bangle.o : typesizes.o matrix.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_fm.o ropp_1dvar.o
ropp_fm_state2roprof.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_fm.o ropp_fm_types.o ropp_fm_copy.o
ropp_qc_genqc.o : typesizes.o matrix.o ropp_utils.o ropp_fm.o ropp_1dvar.o
ropp_pp_read_config.o : ropp_pp.o ropp_pp_types.o messages.o ropp_utils.o arrays.o
ropp_1dvar_covar_refrac.o : typesizes.o matrix.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_fm.o ropp_1dvar.o
ropp_pp_bending_angle_go.o : typesizes.o ropp_pp.o messages.o
ropp_fm_roprof2state.o : typesizes.o ropp_utils.o ropp_io.o ropp_io_types.o ropp_fm.o ropp_fm_types.o
ropp_qc_cutoff.o : typesizes.o ropp_utils.o ropp_fm.o ropp_1dvar.o

.SUFFIXES: .o .c .f .f90 .F90

.c.o:
	$(CC) $(CFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.F90.o:
	$(F90) $(F90FLAGS) -c $<

clean:
	rm -f *.o *.mod

all:
	@echo OBJ_1 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ1_f90} -f ${MKF} VPATH=${VPT1_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC1_f90}" "FFLAGS=${FFLAGS}  ${OPT1_f90} ${INC1_f90}" \
#
	@echo OBJ_2 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ2_f90} -f ${MKF} VPATH=${VPT2_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC2_f90}" "FFLAGS=${FFLAGS}  ${OPT2_f90} ${INC2_f90}" \
#
	@echo OBJ_3 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ3_f90} -f ${MKF} VPATH=${VPT3_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC3_f90}" "FFLAGS=${FFLAGS}  ${OPT3_f90} ${INC3_f90}" \
#
	@echo OBJ_4 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ4_f90} -f ${MKF} VPATH=${VPT4_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC4_f90}" "FFLAGS=${FFLAGS}  ${OPT4_f90} ${INC4_f90}" \
#
	@echo OBJ_5 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ5_f90} -f ${MKF} VPATH=${VPT5_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC5_f90}" "FFLAGS=${FFLAGS}  ${OPT5_f90} ${INC5_f90}" \
#
	@echo OBJ_6 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ6_f90} -f ${MKF} VPATH=${VPT6_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC6_f90}" "FFLAGS=${FFLAGS}  ${OPT6_f90} ${INC6_f90}" \
#
	@echo OBJ_7 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ7_f90} -f ${MKF} VPATH=${VPT7_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC7_f90}" "FFLAGS=${FFLAGS}  ${OPT7_f90} ${INC7_f90}" \
#
	@echo OBJ_8 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ8_f90} -f ${MKF} VPATH=${VPT8_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC8_f90}" "FFLAGS=${FFLAGS}  ${OPT8_f90} ${INC8_f90}" \
#
	@echo OBJ_9 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ9_f90} -f ${MKF} VPATH=${VPT9_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC9_f90}" "FFLAGS=${FFLAGS}  ${OPT9_f90} ${INC9_f90}" \
#
	@echo OBJ_10 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ10_f90} -f ${MKF} VPATH=${VPT10_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC10_f90}" "FFLAGS=${FFLAGS}  ${OPT10_f90} ${INC10_f90}" \
#
	@echo OBJ_11 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ11_f90} -f ${MKF} VPATH=${VPT11_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC11_f90}" "FFLAGS=${FFLAGS}  ${OPT11_f90} ${INC11_f90}" \
#
	@echo OBJ_12 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ12_f90} -f ${MKF} VPATH=${VPT12_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC12_f90}" "FFLAGS=${FFLAGS}  ${OPT12_f90} ${INC12_f90}" \
#
	@echo OBJ_13 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ13_f90} -f ${MKF} VPATH=${VPT13_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC13_f90}" "FFLAGS=${FFLAGS}  ${OPT13_f90} ${INC13_f90}" \
#
	@echo OBJ_14 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ14_f90} -f ${MKF} VPATH=${VPT14_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC14_f90}" "FFLAGS=${FFLAGS}  ${OPT14_f90} ${INC14_f90}" \
#
	@echo OBJ_15 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ15_f90} -f ${MKF} VPATH=${VPT15_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC15_f90}" "FFLAGS=${FFLAGS}  ${OPT15_f90} ${INC15_f90}" \
#
	@echo OBJ_16 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ16_f90} -f ${MKF} VPATH=${VPT16_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC16_f90}" "FFLAGS=${FFLAGS}  ${OPT16_f90} ${INC16_f90}" \
#
	@echo OBJ_17 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ17_f90} -f ${MKF} VPATH=${VPT17_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC17_f90}" "FFLAGS=${FFLAGS}  ${OPT17_f90} ${INC17_f90}" \
#
	@echo OBJ_18 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ18_f90} -f ${MKF} VPATH=${VPT18_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC18_f90}" "FFLAGS=${FFLAGS}  ${OPT18_f90} ${INC18_f90}" \
#
	@echo OBJ_19 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ19_f90} -f ${MKF} VPATH=${VPT19_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC19_f90}" "FFLAGS=${FFLAGS}  ${OPT19_f90} ${INC19_f90}" \
#
	@echo OBJ_20 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ20_f90} -f ${MKF} VPATH=${VPT20_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC20_f90}" "FFLAGS=${FFLAGS}  ${OPT20_f90} ${INC20_f90}" \
#
	@echo OBJ_21 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ21_f90} -f ${MKF} VPATH=${VPT21_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC21_f90}" "FFLAGS=${FFLAGS}  ${OPT21_f90} ${INC21_f90}" \
#
	@echo OBJ_22 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ22_f90} -f ${MKF} VPATH=${VPT22_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC22_f90}" "FFLAGS=${FFLAGS}  ${OPT22_f90} ${INC22_f90}" \
#
	@echo OBJ_23 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ23_f90} -f ${MKF} VPATH=${VPT23_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC23_f90}" "FFLAGS=${FFLAGS}  ${OPT23_f90} ${INC23_f90}" \
#
	@echo OBJ_24 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ24_f90} -f ${MKF} VPATH=${VPT24_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC24_f90}" "FFLAGS=${FFLAGS}  ${OPT24_f90} ${INC24_f90}" \
#
	@echo OBJ_25 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ25_f90} -f ${MKF} VPATH=${VPT25_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC25_f90}" "FFLAGS=${FFLAGS}  ${OPT25_f90} ${INC25_f90}" \
#
	@echo OBJ_26 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ26_f90} -f ${MKF} VPATH=${VPT26_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC26_f90}" "FFLAGS=${FFLAGS}  ${OPT26_f90} ${INC26_f90}" \
#
	@echo OBJ_27 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ27_f90} -f ${MKF} VPATH=${VPT27_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC27_f90}" "FFLAGS=${FFLAGS}  ${OPT27_f90} ${INC27_f90}" \
#
	@echo OBJ_28 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ28_f90} -f ${MKF} VPATH=${VPT28_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC28_f90}" "FFLAGS=${FFLAGS}  ${OPT28_f90} ${INC28_f90}" \
#
	@echo OBJ_29 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ29_f90} -f ${MKF} VPATH=${VPT29_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC29_f90}" "FFLAGS=${FFLAGS}  ${OPT29_f90} ${INC29_f90}" \
#
	@echo OBJ_30 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ30_f90} -f ${MKF} VPATH=${VPT30_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC30_f90}" "FFLAGS=${FFLAGS}  ${OPT30_f90} ${INC30_f90}" \
#
	@echo OBJ_31 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ31_f90} -f ${MKF} VPATH=${VPT31_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC31_f90}" "FFLAGS=${FFLAGS}  ${OPT31_f90} ${INC31_f90}" \
#
	@echo OBJ_32 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ32_f90} -f ${MKF} VPATH=${VPT32_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC32_f90}" "FFLAGS=${FFLAGS}  ${OPT32_f90} ${INC32_f90}" \
#
	@echo OBJ_33 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ33_f90} -f ${MKF} VPATH=${VPT33_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC33_f90}" "FFLAGS=${FFLAGS}  ${OPT33_f90} ${INC33_f90}" \
#
	@echo OBJ_34 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ34_f90} -f ${MKF} VPATH=${VPT34_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC34_f90}" "FFLAGS=${FFLAGS}  ${OPT34_f90} ${INC34_f90}" \
#
	@echo OBJ_35 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ35_f90} -f ${MKF} VPATH=${VPT35_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC35_f90}" "FFLAGS=${FFLAGS}  ${OPT35_f90} ${INC35_f90}" \
#
	@echo OBJ_36 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ36_f90} -f ${MKF} VPATH=${VPT36_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC36_f90}" "FFLAGS=${FFLAGS}  ${OPT36_f90} ${INC36_f90}" \
#
	@echo OBJ_37 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ37_f90} -f ${MKF} VPATH=${VPT37_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC37_f90}" "FFLAGS=${FFLAGS}  ${OPT37_f90} ${INC37_f90}" \
#
	@echo OBJ_38 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ38_f90} -f ${MKF} VPATH=${VPT38_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC38_f90}" "FFLAGS=${FFLAGS}  ${OPT38_f90} ${INC38_f90}" \
#
	@echo OBJ_39 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ39_f90} -f ${MKF} VPATH=${VPT39_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC39_f90}" "FFLAGS=${FFLAGS}  ${OPT39_f90} ${INC39_f90}" \
#
	@echo OBJ_40 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ40_f90} -f ${MKF} VPATH=${VPT40_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC40_f90}" "FFLAGS=${FFLAGS}  ${OPT40_f90} ${INC40_f90}" \
#
	@echo OBJ_41 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ41_f90} -f ${MKF} VPATH=${VPT41_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC41_f90}" "FFLAGS=${FFLAGS}  ${OPT41_f90} ${INC41_f90}" \
#
	@echo OBJ_42 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ42_f90} -f ${MKF} VPATH=${VPT42_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC42_f90}" "FFLAGS=${FFLAGS}  ${OPT42_f90} ${INC42_f90}" \
#
	@echo OBJ_43 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ43_f90} -f ${MKF} VPATH=${VPT43_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC43_f90}" "FFLAGS=${FFLAGS}  ${OPT43_f90} ${INC43_f90}" \
#
	@echo OBJ_44 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ44_f90} -f ${MKF} VPATH=${VPT44_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC44_f90}" "FFLAGS=${FFLAGS}  ${OPT44_f90} ${INC44_f90}" \
#
	@echo OBJ_45 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ45_f90} -f ${MKF} VPATH=${VPT45_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC45_f90}" "FFLAGS=${FFLAGS}  ${OPT45_f90} ${INC45_f90}" \
#
	@echo OBJ_46 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ46_f90} -f ${MKF} VPATH=${VPT46_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC46_f90}" "FFLAGS=${FFLAGS}  ${OPT46_f90} ${INC46_f90}" \
#
	@echo OBJ_47 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ47_f90} -f ${MKF} VPATH=${VPT47_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC47_f90}" "FFLAGS=${FFLAGS}  ${OPT47_f90} ${INC47_f90}" \
#
	@echo OBJ_48 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ48_f90} -f ${MKF} VPATH=${VPT48_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC48_f90}" "FFLAGS=${FFLAGS}  ${OPT48_f90} ${INC48_f90}" \
#
	@echo OBJ_49 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ49_f90} -f ${MKF} VPATH=${VPT49_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC49_f90}" "FFLAGS=${FFLAGS}  ${OPT49_f90} ${INC49_f90}" \
#
	@echo OBJ_50 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ50_f90} -f ${MKF} VPATH=${VPT50_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC50_f90}" "FFLAGS=${FFLAGS}  ${OPT50_f90} ${INC50_f90}" \
#
	@echo OBJ_51 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ51_f90} -f ${MKF} VPATH=${VPT51_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC51_f90}" "FFLAGS=${FFLAGS}  ${OPT51_f90} ${INC51_f90}" \
#
	@echo OBJ_52 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ52_f90} -f ${MKF} VPATH=${VPT52_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC52_f90}" "FFLAGS=${FFLAGS}  ${OPT52_f90} ${INC52_f90}" \
#
	@echo OBJ_53 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ53_f90} -f ${MKF} VPATH=${VPT53_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC53_f90}" "FFLAGS=${FFLAGS}  ${OPT53_f90} ${INC53_f90}" \
#
	@echo OBJ_54 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ54_f90} -f ${MKF} VPATH=${VPT54_f90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS}  ${INC54_f90}" "FFLAGS=${FFLAGS}  ${OPT54_f90} ${INC54_f90}" \
#
	@echo LINK --------------------------------------------
	rm -f $(LIBPATH)/$(NAME)
	ar q  $(LIBPATH)/$(NAME) $(OBJ_ALL)
	@echo END =============================================

