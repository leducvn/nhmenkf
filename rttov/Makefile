USRDIR 	= .
NAME = librttov.a
WKD = ${USRDIR}
MKF = ${USRDIR}/Makefile

include ../Make3.inc
MAKE = make
CPPFLAGS = 
NETCDFINC = /work/opt/local/aarch64/apps/fj/1.2.31/netcdf-fortran/4.5.3/include
OPT_X = 
OPT_RT = ${OPT_X}
OPT_NETCDF = -DNF_REAL_IS_C_FLOAT -DNF_INT1_IS_C_SIGNED_CHAR -DNF_INT2_IS_C_SHORT
OPT_RT_ATLAS = ${OPT_RT}

OPT1_F90 = ${OPT_RT}   
INC1_F90 = 
SRC1_F90 = parkind1.F90
OBJ1_F90 = ${SRC1_F90:.F90=.o}
VPT1_F90 = ${USRDIR}/Source

OPT2_f = ${OPT_RT}   
INC2_f = 
SRC2_f = lapack.f
OBJ2_f = ${SRC2_f:.f=.o}
VPT2_f = ${USRDIR}/Source

OPT3_c = ${OPT_NETCDF}   
INC3_c = 
SRC3_c = npd_wrap_netcdf.c
OBJ3_c = ${SRC3_c:.c=.o}
VPT3_c = ${USRDIR}/Source

OPT4_F90 = ${OPT_RT}   
INC4_F90 = 
SRC4_F90 = rttov_unix_env.F90 mod_rttov_emis_atlas.F90 rttov_sublayer_tl.F90 mod_rttov_fastem5_coef.F90 rttov_deletecomment.F90 yomhook.F90 rttov_skipcommentline.F90 rttov_const.F90 mod_rttov_fastem3_coef.F90 rttov_erfcx.F90 rttov_cmpuc.F90 rttov_sublayer_ad.F90 rttov_sublayer.F90 rttov_sublayer_k.F90 rttov_intex.F90
OBJ4_F90 = ${SRC4_F90:.F90=.o}
VPT4_F90 = ${USRDIR}/Source

OPT5_F90 = ${OPT_RT}   
INC5_F90 = -I${USRDIR}/Inc
SRC5_F90 = rttov_layeravg_k.F90 rttov_layeravg.F90 rttov_layeravg_tl.F90 rttov_layeravg_ad.F90
OBJ5_F90 = ${SRC5_F90:.F90=.o}
VPT5_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT6_F90 = ${OPT_RT}   
INC6_F90 = 
SRC6_F90 = rttov_types.F90 rttov_global.F90 rttov_fastem5.F90 rttov_fastem5_ad.F90 rttov_setpressure.F90 rttov_fastem5_tl.F90 rttov_q2v.F90 rttov_fastem5_k.F90 rttov_findnextsection.F90 rttov_v2q.F90
OBJ6_F90 = ${SRC6_F90:.F90=.o}
VPT6_F90 = ${USRDIR}/Source

OPT7_F90 = ${OPT_RT}   
INC7_F90 = -I${USRDIR}/Source
SRC7_F90 = rttov_boundaryconditions.F90 rttov_boundaryconditions_tl.F90 rttov_boundaryconditions_ad.F90
OBJ7_F90 = ${SRC7_F90:.F90=.o}
VPT7_F90 = ${USRDIR}/Source

OPT8_F90 = ${OPT_RT}   
INC8_F90 = -I${USRDIR}/Inc
SRC8_F90 = rttov_intavg_chan.F90 rttov_intavg_chan_k.F90 rttov_intavg_prof_k.F90 rttov_intavg_chan_tl.F90 rttov_intavg_prof_ad.F90 rttov_intavg_prof.F90 rttov_calcemis_mw_ad.F90 rttov_calcemis_mw_k.F90 rttov_intavg_chan_ad.F90 rttov_calcemis_mw_tl.F90 rttov_refsun.F90 rttov_intavg_prof_tl.F90
OBJ8_F90 = ${SRC8_F90:.F90=.o}
VPT8_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT9_F90 = ${OPT_RT}   
INC9_F90 = 
SRC9_F90 = rttov_setpredictors_8.F90 rttov_init_auxrad_stream.F90 rttov_nullify_optpar_ir.F90 rttov_fresnel.F90 rttov_setpredictors_9_ad.F90 rttov_transmit_k.F90 rttov_getoptions.F90 rttov_integratesource_tl.F90 rttov_opdep_9_ad.F90 rttov_pcscores_rec_k.F90 mod_rttovscatt_test.F90 rttov_add_raytracing.F90 rttov_locpat_tl.F90 rttov_calcrad_tl.F90 rttov_calcbt_pc_tl.F90 rttov_setpredictors_9_solar_k.F90 rttov_scatt_setupindex.F90 rttov_pcscores.F90 rttov_copy_opdp_path.F90 rttov_copy_pccomp.F90 rttov_calcbt_ad.F90 rttov_opdep_k.F90 rttov_opdpscattir_ad.F90 rttov_opdep_9_solar_ad.F90 rttov_transmit_9_solar_tl.F90 rttov_refsun_tl.F90 rttov_reconstruct_tl.F90 rttov_opdep_ad.F90 rttov_opdep_9_solar.F90 rttov_fresnel_tl.F90 rttov_cldstr_tl.F90 rttov_calcrad.F90 rttov_copy_prof.F90 rttov_reconstruct.F90 rttov_pcscores_ad.F90 rttov_nullify_coef_scatt_ir.F90 rttov_transmit_ad.F90 rttov_setpredictors_7_k.F90 rttov_setpredictors_8_k.F90 rttov_setpredictors_9_k.F90 rttov_profaux_tl.F90 rttov_copy_raytracing.F90 rttov_mieproc_tl.F90 rttov_iniedd_ad.F90 rttov_checkinput_tl.F90 rttov_init_pccomp.F90 rttov_cldstr_k.F90 rttov_nullify_coef.F90 rttov_profaux_k.F90 rttov_setpredictors_7_tl.F90 rttov_init_transmission.F90 rttov_transmit_9_solar_k.F90 rttov_setpredictors_9_solar.F90 rttov_calcemis_ir_tl.F90 rttov_setpredictors_9.F90 rttov_locpat.F90 rttov_mieproc.F90 rttov_add_opdp_path.F90 rttov_init_raytracing.F90 rttov_init_ircld.F90 rttov_integratesource_ad.F90 rttov_profaux.F90 rttov_init_prof.F90 rttov_setpredictors_7.F90 rttov_opts_eq.F90 rttov_init_sunglint.F90 rttov_setpredictors_9_solar_ad.F90 rttov_locpat_k.F90 rttov_copy_rad.F90 rttov_locpat_ad.F90 rttov_setpredictors_9_tl.F90 rttov_init_rad.F90 rttov_calcbt_pc_ad.F90 rttov_transmit.F90 rttov_init_predictor.F90 rttov_setpredictors_8_tl.F90 rttov_transmit_9_solar_ad.F90 rttov_integratesource.F90 rttov_refsun_ad.F90 rttov_fresnel_ad.F90 rttov_opdep_9_k.F90 rttov_cldstr.F90 rttov_opdpscattir_tl.F90 rttov_opdep_9_solar_tl.F90 rttov_init_aux_prof.F90 rttov_opdep_9_tl.F90 rttov_opdep_tl.F90 rttov_transmit_9_solar.F90 rttov_nullify_coef_pccomp.F90 rttov_calcrad_ad.F90 rttov_profaux_ad.F90 rttov_calcbt.F90 rttov_pcscores_tl.F90 rttov_opdpscattir.F90 rttov_fresnel_k.F90 rttov_reconstruct_ad.F90 rttov_opdep_9.F90 rttov_alloc_scatt_prof.F90 rttov_opdpscattir_k.F90 rttov_init_trans_scatt_ir.F90 rttov_cldstr_ad.F90 rttov_setpredictors_7_ad.F90 rttov_calcbt_tl.F90 rttov_reconstruct_k.F90 rttov_calcemis_ir_ad.F90 rttov_add_prof.F90 rttov_errorhandling.F90 rttov_opdep_9_solar_k.F90 rttov_calcemis_ir_k.F90 rttov_refsun_k.F90 rttov_checkinput_ad.F90 rttov_mieproc_ad.F90 rttov_opdep.F90 rttov_pcscores_k.F90 rttov_add_aux_prof.F90 rttov_copy_aux_prof.F90 rttov_iniedd.F90 rttov_mult_profiles_k.F90 rttov_transmit_tl.F90 rttov_copy_transmission.F90 rttov_checkinput_k.F90 rttov_init_opdp_path.F90 rttov_calcrad_k.F90 rttov_setpredictors_9_solar_tl.F90 rttov_iniedd_tl.F90 rttov_init_transmission_aux.F90 rttov_setpredictors_8_ad.F90 rttov_calcbt_pc.F90
OBJ9_F90 = ${SRC9_F90:.F90=.o}
VPT9_F90 = ${USRDIR}/Source

OPT10_F90 = ${OPT_RT}   
INC10_F90 = -I${USRDIR}/Inc
SRC10_F90 = rttov_errorreport.F90 rttov_setgeometry_ad.F90 rttov_integrate_k.F90 rttov_nullify_coefs.F90 rttov_setgeometry.F90 rttov_setgeometry_k.F90 rttov_integrate_ad.F90 rttov_setgeometry_tl.F90 rttov_integrate.F90 rttov_integrate_tl.F90
OBJ10_F90 = ${SRC10_F90:.F90=.o}
VPT10_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT11_F90 = ${OPT_RT}   
INC11_F90 = -I${USRDIR}/Inc -I${USRDIR}/Source
SRC11_F90 = rttov_eddington_ad.F90 rttov_eddington.F90 rttov_eddington_tl.F90
OBJ11_F90 = ${SRC11_F90:.F90=.o}
VPT11_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT12_F90 = ${OPT_RT}   
INC12_F90 = -I${USRDIR}/Source -I$(NETCDFINC) -I${USRDIR}/Inc
SRC12_F90 = mod_iratlas.F90
OBJ12_F90 = ${SRC12_F90:.F90=.o}
VPT12_F90 = ${USRDIR}/Source:$(NETCDFINC):${USRDIR}/Inc

OPT13_F90 = ${OPT_RT}   
INC13_F90 = -I${USRDIR}/Inc
SRC13_F90 = rttov_coeffname.F90 rttov_init_coef_pccomp.F90 rttov_calcemis_ir.F90 rttov_write_binary_pccoef.F90 rttov_read_ascii_pccoef.F90 rttov_read_binary_scaercoef.F90 rttov_alloc_auxrad.F90 rttov_alloc_auxrad_stream.F90 mod_mwatlas.F90 rttov_dealloc_coef_scatt_ir.F90 rttov_opencoeff.F90 rttov_init_coef.F90 rttov_checkinput.F90 rttov_alloc_transmission.F90 rttov_alloc_pc_dimensions.F90 rttov_dealloc_coef_pccomp.F90 rttov_checkpcchan.F90 rttov_write_ascii_scaercoef.F90 rttov_alloc_transmission_aux.F90 rttov_write_binary_scaercoef.F90 rttov_write_ascii_sccldcoef.F90 rttov_coef_io_mod.F90 rttov_alloc_trans_scatt_ir.F90 rttov_write_binary_coef.F90 rttov_alloc_raytracing.F90 rttov_alloc_rad.F90 rttov_zutility.F90 rttov_alloc_prof.F90 rttov_dealloc_scattcoeffs.F90 rttov_alloc_predictor.F90 rttov_write_ascii_coef.F90 rttov_alloc_ircld.F90 rttov_read_binary_coef.F90 rttov_write_ascii_pccoef.F90 rttov_read_binary_pccoef.F90 rttov_user_profile_checkinput.F90 rttov_read_binary_sccldcoef.F90 rttov_alloc_aux_prof.F90 rttov_alloc_opdp_path.F90 rttov_write_binary_sccldcoef.F90 rttov_calcemis_mw.F90 rttov_dealloc_coef.F90 rttov_alloc_pccomp.F90
OBJ13_F90 = ${SRC13_F90:.F90=.o}
VPT13_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT14_F90 = ${OPT_RT}   
INC14_F90 = -I${USRDIR}/Inc
SRC14_F90 = rttov_init_coefs.F90 rttov_iniscatt.F90 rttov_get_pc_predictindex.F90 rttov_readscattcoeffs.F90 rttov_read_ascii_sccldcoef.F90 rttov_alloc_traj_sta.F90 rttov_iniscatt_tl.F90 rttov_dealloc_optpar_ir.F90 rttov_write_coefs.F90 rttov_read_ascii_coef.F90 rttov_alloc_traj_dyn.F90 rttov_read_ascii_scaercoef.F90 rttov_iniscatt_ad.F90 rttov_alloc_sunglint.F90 mod_cnrm_mw_atlas.F90
OBJ14_F90 = ${SRC14_F90:.F90=.o}
VPT14_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT15_F90 = ${OPT_RT}   
INC15_F90 = -I${USRDIR}/Source -I$(NETCDFINC) -I${USRDIR}/Inc
SRC15_F90 = rttov_atlas_setup.F90 rttov_get_emis.F90 rttov_deallocate_atlas.F90
OBJ15_F90 = ${SRC15_F90:.F90=.o}
VPT15_F90 = ${USRDIR}/Source:$(NETCDFINC):${USRDIR}/Inc

OPT16_F90 = ${OPT_RT}   
INC16_F90 = -I${USRDIR}/Inc
SRC16_F90 = rttov_read_coefs.F90 rttov_dealloc_coefs.F90 rttov_alloc_traj.F90
OBJ16_F90 = ${SRC16_F90:.F90=.o}
VPT16_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT17_F90 = ${OPT_RT}   
INC17_F90 = -I${USRDIR}/Inc
SRC17_F90 = rttov_setup.F90 rttov_check_traj.F90
OBJ17_F90 = ${SRC17_F90:.F90=.o}
VPT17_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT18_F90 = ${OPT_RT}   
INC18_F90 = -I${USRDIR}/Inc
SRC18_F90 = rttov_direct.F90
OBJ18_F90 = ${SRC18_F90:.F90=.o}
VPT18_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT19_F90 = ${OPT_RT}   
INC19_F90 = -I${USRDIR}/Inc
SRC19_F90 = rttov_parallel_direct.F90 rttov_tl.F90 rttov_ad.F90 rttov_k.F90
OBJ19_F90 = ${SRC19_F90:.F90=.o}
VPT19_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT20_F90 = ${OPT_RT}   
INC20_F90 = -I${USRDIR}/Inc -I${USRDIR}/Source
SRC20_F90 = rttov_scatt.F90
OBJ20_F90 = ${SRC20_F90:.F90=.o}
VPT20_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT21_F90 = ${OPT_RT}   
INC21_F90 = -I${USRDIR}/Inc
SRC21_F90 = rttov_parallel_tl.F90 rttov_parallel_k.F90 rttov_parallel_ad.F90
OBJ21_F90 = ${SRC21_F90:.F90=.o}
VPT21_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT22_F90 = ${OPT_RT}   
INC22_F90 = -I${USRDIR}/Inc -I${USRDIR}/Source
SRC22_F90 = rttov_scatt_ad.F90 rttov_scatt_tl.F90
OBJ22_F90 = ${SRC22_F90:.F90=.o}
VPT22_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OPT23_F90 = ${OPT_RT}   
INC23_F90 = -I${USRDIR}/Inc -I${USRDIR}/Source
SRC23_F90 = rttovscatt_test_one.F90
OBJ23_F90 = ${SRC23_F90:.F90=.o}
VPT23_F90 = ${USRDIR}/Source:${USRDIR}/Inc

OBJ_ALL = ${OBJ1_F90} ${OBJ2_f} ${OBJ3_c} ${OBJ4_F90} ${OBJ5_F90} ${OBJ6_F90} ${OBJ7_F90} ${OBJ8_F90} ${OBJ9_F90} ${OBJ10_F90} ${OBJ11_F90} ${OBJ12_F90} ${OBJ13_F90} ${OBJ14_F90} ${OBJ15_F90} ${OBJ16_F90} ${OBJ17_F90} ${OBJ18_F90} ${OBJ19_F90} ${OBJ20_F90} ${OBJ21_F90} ${OBJ22_F90} ${OBJ23_F90}

rttov_unix_env.o : parkind1.o
mod_rttov_emis_atlas.o : parkind1.o
rttov_sublayer_tl.o : parkind1.o
mod_rttov_fastem5_coef.o : parkind1.o
rttov_deletecomment.o : parkind1.o
yomhook.o : parkind1.o
rttov_skipcommentline.o : parkind1.o
rttov_const.o : parkind1.o
mod_rttov_fastem3_coef.o : parkind1.o
rttov_erfcx.o : parkind1.o
rttov_cmpuc.o : parkind1.o
rttov_sublayer_ad.o : parkind1.o
rttov_sublayer.o : parkind1.o
rttov_sublayer_k.o : parkind1.o
rttov_intex.o : parkind1.o
rttov_types.o : parkind1.o rttov_const.o
rttov_global.o : parkind1.o rttov_const.o
rttov_fastem5.o : mod_rttov_fastem5_coef.o parkind1.o
rttov_fastem5_ad.o : mod_rttov_fastem5_coef.o parkind1.o
rttov_layeravg_k.o : parkind1.o rttov_sublayer.h rttov_sublayer.o rttov_sublayer_k.h rttov_sublayer_k.o
rttov_layeravg.o : parkind1.o rttov_sublayer.h rttov_sublayer.o
rttov_layeravg_tl.o : parkind1.o rttov_sublayer_tl.h rttov_sublayer_tl.o
rttov_setpressure.o : parkind1.o yomhook.o
rttov_fastem5_tl.o : mod_rttov_fastem5_coef.o parkind1.o
rttov_layeravg_ad.o : parkind1.o rttov_sublayer.h rttov_sublayer.o rttov_sublayer_ad.h rttov_sublayer_ad.o
rttov_q2v.o : parkind1.o rttov_const.o
rttov_fastem5_k.o : mod_rttov_fastem5_coef.o parkind1.o
rttov_findnextsection.o : parkind1.o rttov_const.o
rttov_v2q.o : parkind1.o rttov_const.o
rttov_setpredictors_8.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_init_auxrad_stream.o : parkind1.o rttov_types.o
rttov_intavg_chan.o : parkind1.o rttov_types.o yomhook.o rttov_layeravg.h rttov_layeravg.o
rttov_nullify_optpar_ir.o : rttov_types.o
rttov_fresnel.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_setpredictors_9_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_transmit_k.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_getoptions.o : parkind1.o rttov_const.o rttov_types.o rttov_unix_env.o
rttov_integratesource_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_intavg_chan_k.o : parkind1.o rttov_types.o yomhook.o rttov_layeravg.h rttov_layeravg.o rttov_layeravg_k.h rttov_layeravg_k.o
rttov_opdep_9_ad.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_pcscores_rec_k.o : parkind1.o rttov_types.o yomhook.o
mod_rttovscatt_test.o : parkind1.o rttov_types.o
rttov_add_raytracing.o : rttov_types.o
rttov_locpat_tl.o : parkind1.o rttov_const.o rttov_types.o
rttov_calcrad_tl.o : parkind1.o rttov_types.o
rttov_calcbt_pc_tl.o : parkind1.o rttov_types.o yomhook.o
rttov_setpredictors_9_solar_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_scatt_setupindex.o : parkind1.o rttov_types.o rttov_const.o
rttov_pcscores.o : parkind1.o rttov_types.o yomhook.o
rttov_copy_opdp_path.o : rttov_types.o
rttov_intavg_prof_k.o : parkind1.o rttov_types.o yomhook.o rttov_layeravg.h rttov_layeravg.o rttov_layeravg_k.h rttov_layeravg_k.o
rttov_intavg_chan_tl.o : parkind1.o rttov_types.o yomhook.o rttov_layeravg.h rttov_layeravg.o rttov_layeravg_tl.h rttov_layeravg_tl.o
rttov_copy_pccomp.o : rttov_types.o
rttov_calcbt_ad.o : parkind1.o rttov_const.o rttov_types.o
rttov_opdep_k.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_opdpscattir_ad.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_boundaryconditions.o : parkind1.o rttov_const.o rttov_types.o yomhook.o lapack.o
rttov_opdep_9_solar_ad.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_transmit_9_solar_tl.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_refsun_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_reconstruct_tl.o : parkind1.o rttov_types.o yomhook.o
rttov_opdep_ad.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_boundaryconditions_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o lapack.o
rttov_opdep_9_solar.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_fresnel_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_intavg_prof_ad.o : parkind1.o rttov_types.o yomhook.o rttov_layeravg.h rttov_layeravg.o rttov_layeravg_ad.h rttov_layeravg_ad.o
rttov_cldstr_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_intavg_prof.o : parkind1.o rttov_types.o yomhook.o rttov_layeravg.h rttov_layeravg.o
rttov_calcrad.o : parkind1.o rttov_const.o rttov_types.o
rttov_calcemis_mw_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o mod_rttov_fastem3_coef.o rttov_fastem5_ad.h rttov_fastem5_ad.o
rttov_copy_prof.o : parkind1.o rttov_types.o
rttov_reconstruct.o : parkind1.o rttov_types.o yomhook.o
rttov_pcscores_ad.o : parkind1.o rttov_types.o yomhook.o
rttov_nullify_coef_scatt_ir.o : parkind1.o rttov_types.o
rttov_transmit_ad.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_setpredictors_7_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_setpredictors_8_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_calcemis_mw_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o mod_rttov_fastem3_coef.o rttov_fastem5_k.h rttov_fastem5_k.o
rttov_setpredictors_9_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_profaux_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_copy_raytracing.o : rttov_types.o
rttov_mieproc_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_iniedd_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_checkinput_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_init_pccomp.o : parkind1.o rttov_types.o
rttov_cldstr_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_nullify_coef.o : parkind1.o rttov_types.o
rttov_profaux_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_setpredictors_7_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_init_transmission.o : parkind1.o rttov_types.o
rttov_transmit_9_solar_k.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_setpredictors_9_solar.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_calcemis_ir_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_setpredictors_9.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_locpat.o : parkind1.o rttov_const.o rttov_types.o
rttov_mieproc.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_add_opdp_path.o : rttov_types.o
rttov_init_raytracing.o : parkind1.o rttov_types.o
rttov_init_ircld.o : parkind1.o rttov_types.o
rttov_integratesource_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_profaux.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_init_prof.o : parkind1.o rttov_types.o
rttov_setpredictors_7.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_opts_eq.o : parkind1.o rttov_types.o
rttov_init_sunglint.o : parkind1.o rttov_types.o
rttov_setpredictors_9_solar_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_locpat_k.o : parkind1.o rttov_const.o rttov_types.o
rttov_copy_rad.o : rttov_types.o
rttov_locpat_ad.o : parkind1.o rttov_const.o rttov_types.o
rttov_setpredictors_9_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_init_rad.o : parkind1.o rttov_types.o
rttov_calcbt_pc_ad.o : parkind1.o rttov_types.o yomhook.o
rttov_transmit.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_init_predictor.o : parkind1.o rttov_types.o
rttov_intavg_chan_ad.o : parkind1.o rttov_types.o yomhook.o rttov_layeravg.h rttov_layeravg.o rttov_layeravg_ad.h rttov_layeravg_ad.o
rttov_setpredictors_8_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_transmit_9_solar_ad.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_integratesource.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_refsun_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_fresnel_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_opdep_9_k.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_cldstr.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_opdpscattir_tl.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_boundaryconditions_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o lapack.o
rttov_opdep_9_solar_tl.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_init_aux_prof.o : parkind1.o rttov_types.o
rttov_opdep_9_tl.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_opdep_tl.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_transmit_9_solar.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_nullify_coef_pccomp.o : parkind1.o rttov_types.o
rttov_calcrad_ad.o : parkind1.o rttov_types.o
rttov_profaux_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_calcbt.o : parkind1.o rttov_const.o rttov_types.o
rttov_calcemis_mw_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o mod_rttov_fastem3_coef.o rttov_fastem5_tl.h rttov_fastem5_tl.o
rttov_pcscores_tl.o : parkind1.o rttov_types.o yomhook.o
rttov_opdpscattir.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_fresnel_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_reconstruct_ad.o : parkind1.o rttov_types.o yomhook.o
rttov_opdep_9.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_alloc_scatt_prof.o : parkind1.o rttov_types.o yomhook.o
rttov_opdpscattir_k.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_init_trans_scatt_ir.o : parkind1.o rttov_types.o
rttov_cldstr_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_setpredictors_7_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_calcbt_tl.o : parkind1.o rttov_const.o rttov_types.o
rttov_refsun.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_erfcx.h rttov_erfcx.o
rttov_reconstruct_k.o : parkind1.o rttov_types.o yomhook.o
rttov_calcemis_ir_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_add_prof.o : parkind1.o rttov_types.o
rttov_errorhandling.o : parkind1.o rttov_const.o rttov_global.o
rttov_opdep_9_solar_k.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_intavg_prof_tl.o : parkind1.o rttov_types.o yomhook.o rttov_layeravg.h rttov_layeravg.o rttov_layeravg_tl.h rttov_layeravg_tl.o
rttov_calcemis_ir_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_refsun_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_checkinput_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_mieproc_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_opdep.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_pcscores_k.o : parkind1.o rttov_types.o yomhook.o
rttov_add_aux_prof.o : rttov_types.o
rttov_copy_aux_prof.o : rttov_types.o
rttov_iniedd.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_mult_profiles_k.o : parkind1.o rttov_types.o
rttov_transmit_tl.o : rttov_types.o parkind1.o rttov_const.o yomhook.o
rttov_copy_transmission.o : rttov_types.o
rttov_checkinput_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_init_opdp_path.o : parkind1.o rttov_types.o
rttov_calcrad_k.o : parkind1.o rttov_types.o
rttov_setpredictors_9_solar_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_iniedd_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_init_transmission_aux.o : parkind1.o rttov_types.o
rttov_setpredictors_8_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o
rttov_calcbt_pc.o : parkind1.o rttov_types.o yomhook.o
rttov_errorreport.o : parkind1.o rttov_const.o rttov_global.o rttov_types.o rttov_errorhandling.h rttov_errorhandling.o rttov_errorreport.h
rttov_eddington_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_boundaryconditions.h rttov_boundaryconditions.o rttov_integratesource.h rttov_integratesource.o rttov_boundaryconditions_ad.h rttov_boundaryconditions_ad.o rttov_integratesource_ad.h rttov_integratesource_ad.o
rttov_setgeometry_ad.o : parkind1.o rttov_types.o yomhook.o rttov_locpat_ad.h rttov_locpat_ad.o
rttov_integrate_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_calcbt_ad.h rttov_calcbt_ad.o rttov_calcrad_k.h rttov_calcrad_k.o
rttov_nullify_coefs.o : rttov_types.o rttov_nullify_coef.h rttov_nullify_coef.o rttov_nullify_coef_pccomp.h rttov_nullify_coef_pccomp.o rttov_nullify_coef_scatt_ir.h rttov_nullify_coef_scatt_ir.o rttov_nullify_optpar_ir.h rttov_nullify_optpar_ir.o
rttov_setgeometry.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_locpat.h rttov_locpat.o
rttov_setgeometry_k.o : parkind1.o rttov_types.o yomhook.o rttov_locpat_k.h rttov_locpat_k.o
rttov_integrate_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_calcbt_ad.h rttov_calcbt_ad.o rttov_calcrad_ad.h rttov_calcrad_ad.o
rttov_eddington.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_boundaryconditions.h rttov_boundaryconditions.o rttov_integratesource.h rttov_integratesource.o
rttov_eddington_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_boundaryconditions_tl.h rttov_boundaryconditions_tl.o rttov_integratesource_tl.h rttov_integratesource_tl.o
rttov_setgeometry_tl.o : parkind1.o rttov_types.o yomhook.o rttov_locpat_tl.h rttov_locpat_tl.o
rttov_integrate.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_calcbt.h rttov_calcbt.o rttov_calcrad.h rttov_calcrad.o
rttov_integrate_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_calcbt_tl.h rttov_calcbt_tl.o rttov_calcrad_tl.h rttov_calcrad_tl.o
rttov_coeffname.o : parkind1.o rttov_const.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_init_coef_pccomp.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_q2v.h rttov_q2v.o
rttov_calcemis_ir.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o throw.h
mod_iratlas.o : lapack.o netcdf.inc parkind1.o rttov_const.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_write_binary_pccoef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_read_ascii_pccoef.o : parkind1.o rttov_const.o rttov_types.o rttov_cmpuc.h rttov_cmpuc.o rttov_deletecomment.h rttov_deletecomment.o rttov_errorreport.h rttov_errorreport.o rttov_findnextsection.h rttov_findnextsection.o rttov_nullify_coef_pccomp.h rttov_nullify_coef_pccomp.o rttov_skipcommentline.h rttov_skipcommentline.o throw.h
rttov_read_binary_scaercoef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_nullify_coef_scatt_ir.h rttov_nullify_coef_scatt_ir.o
rttov_alloc_auxrad.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_alloc_auxrad_stream.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_auxrad_stream.h rttov_init_auxrad_stream.o
mod_mwatlas.o : parkind1.o rttov_const.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_dealloc_coef_scatt_ir.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_nullify_coef_scatt_ir.h rttov_nullify_coef_scatt_ir.o
rttov_opencoeff.o : parkind1.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_init_coef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_q2v.h rttov_q2v.o
rttov_checkinput.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_alloc_transmission.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_transmission.h rttov_init_transmission.o
rttov_alloc_pc_dimensions.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_dealloc_coef_pccomp.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_nullify_coef_pccomp.h rttov_nullify_coef_pccomp.o
rttov_checkpcchan.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_write_ascii_scaercoef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_alloc_transmission_aux.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_transmission_aux.h rttov_init_transmission_aux.o
rttov_write_binary_scaercoef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_write_ascii_sccldcoef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_coef_io_mod.o : parkind1.o rttov_const.o rttov_cmpuc.h rttov_coeffname.h rttov_errorreport.h rttov_errorreport.o throw.h
rttov_alloc_trans_scatt_ir.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_trans_scatt_ir.h rttov_init_trans_scatt_ir.o
rttov_write_binary_coef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_alloc_raytracing.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_raytracing.h rttov_init_raytracing.o
rttov_alloc_rad.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_rad.h rttov_init_rad.o
rttov_zutility.o : parkind1.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_alloc_prof.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_prof.h rttov_init_prof.o
rttov_dealloc_scattcoeffs.o : parkind1.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_alloc_predictor.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_predictor.h rttov_init_predictor.o
rttov_write_ascii_coef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_alloc_ircld.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_ircld.h rttov_init_ircld.o
rttov_read_binary_coef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_write_ascii_pccoef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_read_binary_pccoef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_user_profile_checkinput.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_read_binary_sccldcoef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_nullify_coef_scatt_ir.h rttov_nullify_coef_scatt_ir.o
rttov_alloc_aux_prof.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_aux_prof.h rttov_init_aux_prof.o
rttov_alloc_opdp_path.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_opdp_path.h rttov_init_opdp_path.o
rttov_write_binary_sccldcoef.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_calcemis_mw.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o throw.h mod_rttov_fastem3_coef.o rttov_fastem5.h rttov_fastem5.o
rttov_dealloc_coef.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_nullify_coef.h rttov_nullify_coef.o
rttov_alloc_pccomp.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_pccomp.h rttov_init_pccomp.o
rttov_init_coefs.o : parkind1.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_init_coef.h rttov_init_coef.o rttov_init_coef_pccomp.h rttov_init_coef_pccomp.o
rttov_iniscatt.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o rttov_calcemis_mw.h rttov_calcemis_mw.o rttov_iniedd.h rttov_iniedd.o rttov_mieproc.h rttov_mieproc.o rttov_setgeometry.h rttov_setgeometry.o
rttov_get_pc_predictindex.o : parkind1.o rttov_const.o rttov_types.o rttov_coef_io_mod.o rttov_cmpuc.h rttov_cmpuc.o rttov_deletecomment.h rttov_deletecomment.o rttov_errorreport.h rttov_errorreport.o rttov_findnextsection.h rttov_findnextsection.o rttov_skipcommentline.h rttov_skipcommentline.o throw.h
rttov_readscattcoeffs.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o rttov_findnextsection.h rttov_findnextsection.o rttov_opencoeff.h rttov_opencoeff.o rttov_skipcommentline.h rttov_skipcommentline.o throw.h
rttov_read_ascii_sccldcoef.o : parkind1.o rttov_const.o rttov_types.o rttov_cmpuc.h rttov_cmpuc.o rttov_deletecomment.h rttov_deletecomment.o rttov_errorreport.h rttov_errorreport.o rttov_findnextsection.h rttov_findnextsection.o rttov_nullify_coef_scatt_ir.h rttov_nullify_coef_scatt_ir.o rttov_opencoeff.h rttov_opencoeff.o rttov_skipcommentline.h rttov_skipcommentline.o throw.h
rttov_alloc_traj_sta.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o rttov_alloc_auxrad.h rttov_alloc_auxrad.o rttov_alloc_pc_dimensions.h rttov_alloc_pc_dimensions.o rttov_alloc_prof.h rttov_alloc_prof.o throw.h
rttov_iniscatt_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o rttov_calcemis_mw.h rttov_calcemis_mw.o rttov_calcemis_mw_tl.h rttov_calcemis_mw_tl.o rttov_iniedd_tl.h rttov_iniedd_tl.o rttov_mieproc_tl.h rttov_mieproc_tl.o rttov_setgeometry.h rttov_setgeometry.o
rttov_dealloc_optpar_ir.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_dealloc_coef_scatt_ir.h rttov_dealloc_coef_scatt_ir.o rttov_nullify_optpar_ir.h rttov_nullify_optpar_ir.o
rttov_write_coefs.o : parkind1.o rttov_coef_io_mod.o rttov_types.o yomhook.o rttov_cmpuc.h rttov_cmpuc.o rttov_errorreport.h rttov_errorreport.o rttov_write_ascii_coef.h rttov_write_ascii_coef.o rttov_write_ascii_pccoef.h rttov_write_ascii_pccoef.o rttov_write_ascii_scaercoef.h rttov_write_ascii_scaercoef.o rttov_write_ascii_sccldcoef.h rttov_write_ascii_sccldcoef.o rttov_write_binary_coef.h rttov_write_binary_coef.o rttov_write_binary_pccoef.h rttov_write_binary_pccoef.o rttov_write_binary_scaercoef.h rttov_write_binary_scaercoef.o rttov_write_binary_sccldcoef.h rttov_write_binary_sccldcoef.o throw.h
rttov_read_ascii_coef.o : parkind1.o rttov_const.o rttov_types.o rttov_cmpuc.h rttov_cmpuc.o rttov_deletecomment.h rttov_deletecomment.o rttov_errorreport.h rttov_errorreport.o rttov_findnextsection.h rttov_findnextsection.o rttov_nullify_coef.h rttov_nullify_coef.o rttov_nullify_coef_pccomp.h rttov_nullify_coef_pccomp.o rttov_nullify_coef_scatt_ir.h rttov_nullify_coef_scatt_ir.o rttov_nullify_optpar_ir.h rttov_nullify_optpar_ir.o rttov_opencoeff.h rttov_opencoeff.o rttov_skipcommentline.h rttov_skipcommentline.o throw.h
rttov_alloc_traj_dyn.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o rttov_alloc_auxrad_stream.h rttov_alloc_auxrad_stream.o rttov_alloc_trans_scatt_ir.h rttov_alloc_trans_scatt_ir.o rttov_alloc_transmission_aux.h rttov_alloc_transmission_aux.o throw.h
rttov_read_ascii_scaercoef.o : parkind1.o rttov_const.o rttov_types.o rttov_cmpuc.h rttov_cmpuc.o rttov_deletecomment.h rttov_deletecomment.o rttov_errorreport.h rttov_errorreport.o rttov_findnextsection.h rttov_findnextsection.o rttov_nullify_coef_scatt_ir.h rttov_nullify_coef_scatt_ir.o rttov_opencoeff.h rttov_opencoeff.o rttov_skipcommentline.h rttov_skipcommentline.o throw.h
rttov_iniscatt_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o rttov_calcemis_mw.h rttov_calcemis_mw.o rttov_calcemis_mw_k.h rttov_calcemis_mw_k.o rttov_calcemis_mw_ad.h rttov_calcemis_mw_ad.o rttov_iniedd.h rttov_iniedd.o rttov_iniedd_ad.h rttov_iniedd_ad.o rttov_mieproc_ad.h rttov_mieproc_ad.o rttov_mieproc.h rttov_mieproc.o rttov_setgeometry.h rttov_setgeometry.o
rttov_alloc_sunglint.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_alloc_aux_prof.h rttov_alloc_aux_prof.o rttov_init_sunglint.h rttov_init_sunglint.o
mod_cnrm_mw_atlas.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o rttov_opencoeff.h rttov_opencoeff.o throw.h
rttov_read_coefs.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_coef_io_mod.o rttov_cmpuc.h rttov_cmpuc.o rttov_errorreport.h rttov_errorreport.o rttov_nullify_coef.h rttov_nullify_coef.o rttov_nullify_coef_pccomp.h rttov_nullify_coef_pccomp.o rttov_nullify_coef_scatt_ir.h rttov_nullify_coef_scatt_ir.o rttov_nullify_optpar_ir.h rttov_nullify_optpar_ir.o rttov_read_ascii_coef.h rttov_read_ascii_coef.o rttov_read_ascii_pccoef.h rttov_read_ascii_pccoef.o rttov_read_ascii_scaercoef.h rttov_read_ascii_scaercoef.o rttov_read_ascii_sccldcoef.h rttov_read_ascii_sccldcoef.o rttov_read_binary_coef.h rttov_read_binary_coef.o rttov_read_binary_pccoef.h rttov_read_binary_pccoef.o rttov_read_binary_scaercoef.h rttov_read_binary_scaercoef.o rttov_read_binary_sccldcoef.h rttov_read_binary_sccldcoef.o throw.h
rttov_dealloc_coefs.o : parkind1.o rttov_types.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_dealloc_coef.h rttov_dealloc_coef.o rttov_dealloc_coef_pccomp.h rttov_dealloc_coef_pccomp.o rttov_dealloc_coef_scatt_ir.h rttov_dealloc_coef_scatt_ir.o rttov_dealloc_optpar_ir.h rttov_dealloc_optpar_ir.o
rttov_atlas_setup.o : parkind1.o rttov_const.o rttov_types.o mod_iratlas.o mod_mwatlas.o mod_cnrm_mw_atlas.o mod_rttov_emis_atlas.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_alloc_traj.o : parkind1.o rttov_const.o rttov_types.o rttov_errorreport.h rttov_errorreport.o rttov_alloc_aux_prof.h rttov_alloc_aux_prof.o rttov_alloc_ircld.h rttov_alloc_ircld.o rttov_alloc_opdp_path.h rttov_alloc_opdp_path.o rttov_alloc_predictor.h rttov_alloc_predictor.o rttov_alloc_prof.h rttov_alloc_prof.o rttov_alloc_raytracing.h rttov_alloc_raytracing.o rttov_alloc_sunglint.h rttov_alloc_sunglint.o rttov_alloc_trans_scatt_ir.h rttov_alloc_trans_scatt_ir.o throw.h
rttov_get_emis.o : parkind1.o rttov_const.o rttov_types.o mod_iratlas.o mod_mwatlas.o mod_cnrm_mw_atlas.o mod_rttov_emis_atlas.o yomhook.o rttov_errorreport.h rttov_errorreport.o throw.h
rttov_deallocate_atlas.o : parkind1.o rttov_const.o rttov_types.o mod_iratlas.o mod_mwatlas.o mod_cnrm_mw_atlas.o mod_rttov_emis_atlas.o
rttov_setup.o : parkind1.o rttov_types.o yomhook.o rttov_errorhandling.h rttov_errorhandling.o rttov_errorreport.h rttov_errorreport.o rttov_init_coefs.h rttov_init_coefs.o rttov_read_coefs.h rttov_read_coefs.o throw.h
rttov_check_traj.o : parkind1.o rttov_types.o yomhook.o rttov_errorreport.h rttov_errorreport.o throw.h rttov_alloc_traj.h rttov_alloc_traj.o rttov_opts_eq.h rttov_opts_eq.o
rttov_direct.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_integrate.h rttov_integrate.o rttov_alloc_traj_dyn.h rttov_alloc_traj_dyn.o rttov_alloc_traj_sta.h rttov_alloc_traj_sta.o rttov_calcbt_pc.h rttov_calcbt_pc.o rttov_calcemis_ir.h rttov_calcemis_ir.o rttov_calcemis_mw.h rttov_calcemis_mw.o rttov_check_traj.h rttov_check_traj.o rttov_checkinput.h rttov_checkinput.o rttov_checkpcchan.h rttov_checkpcchan.o rttov_cldstr.h rttov_cldstr.o rttov_copy_aux_prof.h rttov_copy_aux_prof.o rttov_copy_opdp_path.h rttov_copy_opdp_path.o rttov_copy_prof.h rttov_copy_prof.o rttov_copy_raytracing.h rttov_copy_raytracing.o rttov_errorreport.h rttov_errorreport.o rttov_fresnel.h rttov_fresnel.o rttov_init_aux_prof.h rttov_init_aux_prof.o rttov_init_ircld.h rttov_init_ircld.o rttov_init_prof.h rttov_init_prof.o rttov_init_raytracing.h rttov_init_raytracing.o rttov_intavg_chan.h rttov_intavg_chan.o rttov_intavg_prof.h rttov_intavg_prof.o rttov_opdep.h rttov_opdep.o rttov_opdep_9.h rttov_opdep_9.o rttov_opdep_9_solar.h rttov_opdep_9_solar.o rttov_opdpscattir.h rttov_opdpscattir.o rttov_pcscores.h rttov_pcscores.o rttov_profaux.h rttov_profaux.o rttov_reconstruct.h rttov_reconstruct.o rttov_refsun.h rttov_refsun.o rttov_setgeometry.h rttov_setgeometry.o rttov_setpredictors_7.h rttov_setpredictors_7.o rttov_setpredictors_8.h rttov_setpredictors_8.o rttov_setpredictors_9.h rttov_setpredictors_9.o rttov_setpredictors_9_solar.h rttov_setpredictors_9_solar.o rttov_transmit.h rttov_transmit.o rttov_transmit_9_solar.h rttov_transmit_9_solar.o throw.h
rttov_parallel_direct.o : parkind1.o rttov_const.o rttov_types.o rttov_integrate.h rttov_integrate.h rttov_add_prof.h rttov_add_prof.o rttov_alloc_prof.h rttov_alloc_prof.o rttov_alloc_rad.h rttov_alloc_rad.o rttov_direct.h rttov_direct.o rttov_errorreport.h rttov_errorreport.o rttov_init_prof.h rttov_init_prof.o
rttov_tl.o : rttov_types.o parkind1.o rttov_const.o yomhook.o rttov_integrate.h rttov_integrate.h rttov_alloc_traj_dyn.h rttov_alloc_traj_dyn.o rttov_alloc_traj_sta.h rttov_alloc_traj_sta.o rttov_calcbt_pc_tl.h rttov_calcbt_pc_tl.o rttov_calcemis_ir_tl.h rttov_calcemis_ir_tl.o rttov_calcemis_mw_tl.h rttov_calcemis_mw_tl.o rttov_check_traj.h rttov_check_traj.o rttov_checkinput_tl.h rttov_checkinput_tl.o rttov_cldstr_tl.h rttov_cldstr_tl.o rttov_copy_aux_prof.h rttov_copy_aux_prof.o rttov_copy_opdp_path.h rttov_copy_opdp_path.o rttov_copy_prof.h rttov_copy_prof.o rttov_copy_raytracing.h rttov_copy_raytracing.o rttov_direct.h rttov_direct.o rttov_errorreport.h rttov_errorreport.o rttov_fresnel_tl.h rttov_fresnel_tl.o rttov_init_aux_prof.h rttov_init_aux_prof.o rttov_init_opdp_path.h rttov_init_opdp_path.o rttov_init_prof.h rttov_init_prof.o rttov_init_raytracing.h rttov_init_raytracing.o rttov_intavg_chan_tl.h rttov_intavg_chan_tl.o rttov_intavg_prof_tl.h rttov_intavg_prof_tl.o rttov_integrate_tl.h rttov_integrate_tl.o rttov_opdep_9_solar_tl.h rttov_opdep_9_solar_tl.o rttov_opdep_9_tl.h rttov_opdep_9_tl.o rttov_opdep_tl.h rttov_opdep_tl.o rttov_opdpscattir_tl.h rttov_opdpscattir_tl.o rttov_pcscores_tl.h rttov_pcscores_tl.o rttov_profaux_tl.h rttov_profaux_tl.o rttov_reconstruct_tl.h rttov_reconstruct_tl.o rttov_refsun_tl.h rttov_refsun_tl.o rttov_setgeometry_tl.h rttov_setgeometry_tl.o rttov_setpredictors_7_tl.h rttov_setpredictors_7_tl.o rttov_setpredictors_8_tl.h rttov_setpredictors_8_tl.o rttov_setpredictors_9_solar_tl.h rttov_setpredictors_9_solar_tl.o rttov_setpredictors_9_tl.h rttov_setpredictors_9_tl.o rttov_transmit_9_solar_tl.h rttov_transmit_9_solar_tl.o rttov_transmit_tl.h rttov_transmit_tl.o throw.h
rttov_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_integrate.h rttov_integrate.h rttov_add_aux_prof.h rttov_add_aux_prof.o rttov_add_opdp_path.h rttov_add_opdp_path.o rttov_add_prof.h rttov_add_prof.o rttov_add_raytracing.h rttov_add_raytracing.o rttov_alloc_traj_dyn.h rttov_alloc_traj_dyn.o rttov_alloc_traj_sta.h rttov_alloc_traj_sta.o rttov_calcbt_pc_ad.h rttov_calcbt_pc_ad.o rttov_calcemis_ir_ad.h rttov_calcemis_ir_ad.o rttov_calcemis_mw_ad.h rttov_calcemis_mw_ad.o rttov_check_traj.h rttov_check_traj.o rttov_checkinput_ad.h rttov_checkinput_ad.o rttov_cldstr_ad.h rttov_cldstr_ad.o rttov_direct.h rttov_direct.o rttov_errorreport.h rttov_errorreport.o rttov_fresnel_ad.h rttov_fresnel_ad.o rttov_init_aux_prof.h rttov_init_aux_prof.o rttov_init_auxrad_stream.h rttov_init_auxrad_stream.o rttov_init_ircld.h rttov_init_ircld.o rttov_init_opdp_path.h rttov_init_opdp_path.o rttov_init_predictor.h rttov_init_predictor.o rttov_init_prof.h rttov_init_prof.o rttov_init_raytracing.h rttov_init_raytracing.o rttov_init_sunglint.h rttov_init_sunglint.o rttov_init_trans_scatt_ir.h rttov_init_trans_scatt_ir.o rttov_init_transmission_aux.h rttov_init_transmission_aux.o rttov_intavg_chan_ad.h rttov_intavg_chan_ad.o rttov_intavg_prof_ad.h rttov_intavg_prof_ad.o rttov_integrate_ad.h rttov_integrate_ad.o rttov_opdep_9_ad.h rttov_opdep_9_ad.o rttov_opdep_9_solar_ad.h rttov_opdep_9_solar_ad.o rttov_opdep_ad.h rttov_opdep_ad.o rttov_opdpscattir_ad.h rttov_opdpscattir_ad.o rttov_pcscores_ad.h rttov_pcscores_ad.o rttov_profaux_ad.h rttov_profaux_ad.o rttov_reconstruct_ad.h rttov_reconstruct_ad.o rttov_refsun_ad.h rttov_refsun_ad.o rttov_setgeometry_ad.h rttov_setgeometry_ad.o rttov_setpredictors_7_ad.h rttov_setpredictors_7_ad.o rttov_setpredictors_8_ad.h rttov_setpredictors_8_ad.o rttov_setpredictors_9_ad.h rttov_setpredictors_9_ad.o rttov_setpredictors_9_solar_ad.h rttov_setpredictors_9_solar_ad.o rttov_transmit_9_solar_ad.h rttov_transmit_9_solar_ad.o rttov_transmit_ad.h rttov_transmit_ad.o throw.h
rttov_scatt.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_integrate.h rttov_integrate.o rttov_direct.h rttov_direct.o rttov_eddington.h rttov_eddington.o rttov_errorreport.h rttov_errorreport.o rttov_iniscatt.h rttov_iniscatt.o
rttov_k.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_integrate.h rttov_integrate.o rttov_add_aux_prof.h rttov_add_aux_prof.o rttov_add_opdp_path.h rttov_add_opdp_path.o rttov_add_prof.h rttov_add_prof.o rttov_add_raytracing.h rttov_add_raytracing.o rttov_alloc_traj_dyn.h rttov_alloc_traj_dyn.o rttov_alloc_traj_sta.h rttov_alloc_traj_sta.o rttov_calcbt_pc.h rttov_calcbt_pc.o rttov_calcbt_pc_ad.h rttov_calcbt_pc_ad.o rttov_calcemis_ir_k.h rttov_calcemis_ir_k.o rttov_calcemis_mw_k.h rttov_calcemis_mw_k.o rttov_check_traj.h rttov_check_traj.o rttov_checkinput_k.h rttov_checkinput_k.o rttov_cldstr_k.h rttov_cldstr_k.o rttov_direct.h rttov_direct.o rttov_errorreport.h rttov_errorreport.o rttov_fresnel_k.h rttov_fresnel_k.o rttov_init_aux_prof.h rttov_init_aux_prof.o rttov_init_auxrad_stream.h rttov_init_auxrad_stream.o rttov_init_ircld.h rttov_init_ircld.o rttov_init_opdp_path.h rttov_init_opdp_path.o rttov_init_predictor.h rttov_init_predictor.o rttov_init_prof.h rttov_init_prof.o rttov_init_rad.h rttov_init_rad.o rttov_init_raytracing.h rttov_init_raytracing.o rttov_init_sunglint.h rttov_init_sunglint.o rttov_init_trans_scatt_ir.h rttov_init_trans_scatt_ir.o rttov_init_transmission_aux.h rttov_init_transmission_aux.o rttov_intavg_chan_k.h rttov_intavg_chan_k.o rttov_intavg_prof_k.h rttov_intavg_prof_k.o rttov_integrate_k.h rttov_integrate_k.o rttov_mult_profiles_k.h rttov_mult_profiles_k.o rttov_opdep_9_k.h rttov_opdep_9_k.o rttov_opdep_9_solar_k.h rttov_opdep_9_solar_k.o rttov_opdep_k.h rttov_opdep_k.o rttov_opdpscattir_k.h rttov_opdpscattir_k.o rttov_pcscores.h rttov_pcscores.o rttov_pcscores_k.h rttov_pcscores_k.o rttov_pcscores_rec_k.h rttov_pcscores_rec_k.o rttov_profaux_k.h rttov_profaux_k.o rttov_reconstruct.h rttov_reconstruct.o rttov_reconstruct_k.h rttov_reconstruct_k.o rttov_refsun_k.h rttov_refsun_k.o rttov_setgeometry_k.h rttov_setgeometry_k.o rttov_setpredictors_7_k.h rttov_setpredictors_7_k.o rttov_setpredictors_8_k.h rttov_setpredictors_8_k.o rttov_setpredictors_9_k.h rttov_setpredictors_9_k.o rttov_setpredictors_9_solar_k.h rttov_setpredictors_9_solar_k.o rttov_transmit_9_solar_k.h rttov_transmit_9_solar_k.o rttov_transmit_k.h rttov_transmit_k.o throw.h
rttov_scatt_ad.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_integrate.h rttov_integrate.o rttov_ad.h rttov_ad.o rttov_direct.h rttov_direct.o rttov_eddington.h rttov_eddington.o rttov_eddington_ad.h rttov_eddington_ad.o rttov_errorreport.h rttov_errorreport.o rttov_iniscatt.h rttov_iniscatt.o rttov_iniscatt_ad.h rttov_iniscatt_ad.o rttov_k.h rttov_k.o
rttov_parallel_tl.o : parkind1.o rttov_const.o rttov_types.o rttov_add_prof.h rttov_add_prof.o rttov_alloc_prof.h rttov_alloc_prof.o rttov_alloc_rad.h rttov_alloc_rad.o rttov_errorreport.h rttov_errorreport.o rttov_init_prof.h rttov_init_prof.o rttov_tl.h rttov_tl.o
rttov_parallel_k.o : parkind1.o rttov_const.o rttov_types.o rttov_add_prof.h rttov_add_prof.o rttov_alloc_prof.h rttov_alloc_prof.o rttov_alloc_rad.h rttov_alloc_rad.o rttov_errorreport.h rttov_errorreport.o rttov_init_prof.h rttov_init_prof.o rttov_k.h rttov_k.o
rttov_scatt_tl.o : parkind1.o rttov_const.o rttov_types.o yomhook.o rttov_eddington_tl.h rttov_eddington_tl.o rttov_errorreport.h rttov_errorreport.o rttov_iniscatt_tl.h rttov_iniscatt_tl.o rttov_tl.h rttov_tl.o
rttov_parallel_ad.o : parkind1.o rttov_const.o rttov_types.o rttov_ad.h rttov_ad.o rttov_add_prof.h rttov_add_prof.o rttov_alloc_prof.h rttov_alloc_prof.o rttov_alloc_rad.h rttov_alloc_rad.o rttov_errorreport.h rttov_errorreport.o rttov_init_prof.h rttov_init_prof.o
rttovscatt_test_one.o : parkind1.o rttov_const.o rttov_types.o mod_rttovscatt_test.o yomhook.o rttov_alloc_rad.h rttov_alloc_rad.o rttov_alloc_prof.h rttov_alloc_prof.o rttov_alloc_scatt_prof.h rttov_alloc_scatt_prof.o rttov_scatt.h rttov_scatt.o rttov_scatt_tl.h rttov_scatt_tl.o rttov_scatt_ad.h rttov_scatt_ad.o rttov_setpressure.h rttov_setpressure.o

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
	${MAKE} ${OBJ1_F90} -f ${MKF} VPATH=${VPT1_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC1_F90}" "F90FLAGS=${FFLAGS}  ${OPT1_F90} ${INC1_F90}" \
#
	@echo OBJ_2 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ2_f} -f ${MKF} VPATH=${VPT2_f} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC2_f}" "FFLAGS=${FFLAGS}  ${OPT2_f} ${INC2_f}" \
#
	@echo OBJ_3 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ3_c} -f ${MKF} VPATH=${VPT3_c} \
		"CC=${CC}" "CPPFLAGS=${CPPFLAGS}  ${INC3_c}" "CFLAGS=${CFLAGS}  ${OPT3_c} ${INC3_c}" \
#
	@echo OBJ_4 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ4_F90} -f ${MKF} VPATH=${VPT4_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC4_F90}" "F90FLAGS=${FFLAGS}  ${OPT4_F90} ${INC4_F90}" \
#
	@echo OBJ_5 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ5_F90} -f ${MKF} VPATH=${VPT5_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC5_F90}" "F90FLAGS=${FFLAGS}  ${OPT5_F90} ${INC5_F90}" \
#
	@echo OBJ_6 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ6_F90} -f ${MKF} VPATH=${VPT6_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC6_F90}" "F90FLAGS=${FFLAGS}  ${OPT6_F90} ${INC6_F90}" \
#
	@echo OBJ_7 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ7_F90} -f ${MKF} VPATH=${VPT7_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC7_F90}" "F90FLAGS=${FFLAGS}  ${OPT7_F90} ${INC7_F90}" \
#
	@echo OBJ_8 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ8_F90} -f ${MKF} VPATH=${VPT8_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC8_F90}" "F90FLAGS=${FFLAGS}  ${OPT8_F90} ${INC8_F90}" \
#
	@echo OBJ_9 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ9_F90} -f ${MKF} VPATH=${VPT9_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC9_F90}" "F90FLAGS=${FFLAGS}  ${OPT9_F90} ${INC9_F90}" \
#
	@echo OBJ_10 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ10_F90} -f ${MKF} VPATH=${VPT10_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC10_F90}" "F90FLAGS=${FFLAGS}  ${OPT10_F90} ${INC10_F90}" \
#
	@echo OBJ_11 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ11_F90} -f ${MKF} VPATH=${VPT11_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC11_F90}" "F90FLAGS=${FFLAGS}  ${OPT11_F90} ${INC11_F90}" \
#
	@echo OBJ_12 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ12_F90} -f ${MKF} VPATH=${VPT12_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC12_F90}" "F90FLAGS=${FFLAGS}  ${OPT12_F90} ${INC12_F90}" \
#
	@echo OBJ_13 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ13_F90} -f ${MKF} VPATH=${VPT13_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC13_F90}" "F90FLAGS=${FFLAGS}  ${OPT13_F90} ${INC13_F90}" \
#
	@echo OBJ_14 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ14_F90} -f ${MKF} VPATH=${VPT14_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC14_F90}" "F90FLAGS=${FFLAGS}  ${OPT14_F90} ${INC14_F90}" \
#
	@echo OBJ_15 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ15_F90} -f ${MKF} VPATH=${VPT15_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC15_F90}" "F90FLAGS=${FFLAGS}  ${OPT15_F90} ${INC15_F90}" \
#
	@echo OBJ_16 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ16_F90} -f ${MKF} VPATH=${VPT16_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC16_F90}" "F90FLAGS=${FFLAGS}  ${OPT16_F90} ${INC16_F90}" \
#
	@echo OBJ_17 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ17_F90} -f ${MKF} VPATH=${VPT17_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC17_F90}" "F90FLAGS=${FFLAGS}  ${OPT17_F90} ${INC17_F90}" \
#
	@echo OBJ_18 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ18_F90} -f ${MKF} VPATH=${VPT18_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC18_F90}" "F90FLAGS=${FFLAGS}  ${OPT18_F90} ${INC18_F90}" \
#
	@echo OBJ_19 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ19_F90} -f ${MKF} VPATH=${VPT19_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC19_F90}" "F90FLAGS=${FFLAGS}  ${OPT19_F90} ${INC19_F90}" \
#
	@echo OBJ_20 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ20_F90} -f ${MKF} VPATH=${VPT20_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC20_F90}" "F90FLAGS=${FFLAGS}  ${OPT20_F90} ${INC20_F90}" \
#
	@echo OBJ_21 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ21_F90} -f ${MKF} VPATH=${VPT21_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC21_F90}" "F90FLAGS=${FFLAGS}  ${OPT21_F90} ${INC21_F90}" \
#
	@echo OBJ_22 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ22_F90} -f ${MKF} VPATH=${VPT22_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC22_F90}" "F90FLAGS=${FFLAGS}  ${OPT22_F90} ${INC22_F90}" \
#
	@echo OBJ_23 ------------------------------
	cd ${WKD}; \
	${MAKE} ${OBJ23_F90} -f ${MKF} VPATH=${VPT23_F90} \
		"FC=${FC}" "CPPFLAGS=${CPPFLAGS} -DRTTOV_ARCH_VECTOR ${INC23_F90}" "F90FLAGS=${FFLAGS}  ${OPT23_F90} ${INC23_F90}" \
#
	@echo LINK --------------------------------------------
	rm -f $(LIBPATH)/$(NAME)
	ar q  $(LIBPATH)/$(NAME) $(OBJ_ALL)
	@echo END =============================================

