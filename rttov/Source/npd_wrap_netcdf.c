/*
 * Wrapper for netcdff
 *  NF_INT1_IS_C_SIGNED_CHAR etc. are defined in "nfconfig.inc"
 *                                     -- T. Egawa  2012/2/16
 */

int NF_OPEN(char *path, int *mode, int *ncid) {
  return nf_open(path, mode, ncid);
}
int NF_INQ(int *ncid, int *ndims, int *nvars, int *natts, int *unlimdimid) {
  return nf_inq(ncid, ndims, nvars, natts, unlimdimid);
}
int NF_INQ_DIM(int *ncid, int *dimid, char *name, int *len) {
  return nf_inq_dim(ncid, dimid, name, len);
}
int NF_INQ_VARID(int *ncid, char *name, int *varid) {
  return nf_inq_varid(ncid, name, varid);
}
#if NF_REAL_IS_C_DOUBLE
int NF_GET_ATT_REAL(int *ncid, int *varid, char *name, double *rvals) {
#else
int NF_GET_ATT_REAL(int *ncid, int *varid, char *name, float *rvals) {
#endif
  return nf_get_att_real(ncid, varid, name, rvals);
}
#if NF_INT1_IS_C_SIGNED_CHAR
int NF_GET_VAR_INT1(int *ncid, int *varid, signed char *i1vals) {
#elif NF_INT1_IS_C_SHORT
int NF_GET_VAR_INT1(int *ncid, int *varid, short *i1vals) {
#elif NF_INT1_IS_C_INT
int NF_GET_VAR_INT1(int *ncid, int *varid, int *i1vals) {
#elif NF_INT1_IS_C_LONG
int NF_GET_VAR_INT1(int *ncid, int *varid, long *i1vals) {
#endif
  return nf_get_var_int1(ncid, varid, i1vals);
}
#if NF_INT2_IS_C_SHORT
int NF_GET_VAR_INT2(int *ncid, int *varid, short *i2vals) {
#elif NF_INT2_IS_C_INT
int NF_GET_VAR_INT2(int *ncid, int *varid, int *i2vals) {
#elif NF_INT2_IS_C_LONG
int NF_GET_VAR_INT2(int *ncid, int *varid, long *i2vals) {
#endif
  return nf_get_var_int2(ncid, varid, i2vals);
}
#if NF_REAL_IS_C_DOUBLE
int NF_GET_VAR_REAL(int *ncid, int *varid, double *rvals) {
#else
int NF_GET_VAR_REAL(int *ncid, int *varid, float *rvals) {
#endif
  return nf_get_var_real(ncid, varid, rvals);
}
int NF_CLOSE(int *ncid) {
  return nf_close(ncid);
}
