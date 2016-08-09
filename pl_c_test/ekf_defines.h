/*
EKF_CAMERA_CALC_NIS(__P, __R, __TBN, __SUBX, __X, __Z, __RET_NIS)
EKF_CAMERA_CALC_STATE(__P, __R, __TBN, __SUBX, __X, __Z, __RET_STATE)
EKF_CAMERA_CALC_COV(__P, __R, __TBN, __SUBX, __X, __Z, __RET_COV)
EKF_CAMERA_CALC_SUBX(__P, __R, __TBN, __X, __Z, __RET_SUBX)
EKF_CAMERA_CALC_INNOV(__P, __R, __TBN, __SUBX, __X, __Z, __RET_INNOV)
EKF_HEIGHT_CALC_NIS(__P, __R, __SUBX, __X, __Z, __RET_NIS)
EKF_HEIGHT_CALC_STATE(__P, __R, __SUBX, __X, __Z, __RET_STATE)
EKF_HEIGHT_CALC_COV(__P, __R, __SUBX, __X, __Z, __RET_COV)
EKF_HEIGHT_CALC_SUBX(__P, __R, __X, __Z, __RET_SUBX)
EKF_HEIGHT_CALC_INNOV(__P, __R, __SUBX, __X, __Z, __RET_INNOV)
EKF_INITIALIZATION_CALC_STATE(__TBN, __CAM_POS, __CAM_POS_R, __HGT, __HGT_R, __INIT_FOC_LEN, __INIT_FOC_LEN_R, __SUBX, __VEL, __VEL_R, __RET_STATE)
EKF_INITIALIZATION_CALC_SUBX(__TBN, __CAM_POS, __CAM_POS_R, __HGT, __HGT_R, __INIT_FOC_LEN, __INIT_FOC_LEN_R, __VEL, __VEL_R, __RET_SUBX)
EKF_INITIALIZATION_CALC_COV(__TBN, __CAM_POS, __CAM_POS_R, __HGT, __HGT_R, __INIT_FOC_LEN, __INIT_FOC_LEN_R, __SUBX, __VEL, __VEL_R, __RET_COV)
EKF_PREDICTION_CALC_STATE(__P, __DT, __SUBX, __U, __W_U_SIGMA, __X, __RET_STATE)
EKF_PREDICTION_CALC_SUBX(__P, __DT, __U, __W_U_SIGMA, __X, __RET_SUBX)
EKF_PREDICTION_CALC_COV(__P, __DT, __SUBX, __U, __W_U_SIGMA, __X, __RET_COV)
EKF_VELD_CALC_NIS(__P, __R, __SUBX, __X, __Z, __RET_NIS)
EKF_VELD_CALC_STATE(__P, __R, __SUBX, __X, __Z, __RET_STATE)
EKF_VELD_CALC_COV(__P, __R, __SUBX, __X, __Z, __RET_COV)
EKF_VELD_CALC_SUBX(__P, __R, __X, __Z, __RET_SUBX)
EKF_VELD_CALC_INNOV(__P, __R, __SUBX, __X, __Z, __RET_INNOV)
EKF_VELNE_CALC_NIS(__P, __R, __SUBX, __X, __Z, __RET_NIS)
EKF_VELNE_CALC_STATE(__P, __R, __SUBX, __X, __Z, __RET_STATE)
EKF_VELNE_CALC_COV(__P, __R, __SUBX, __X, __Z, __RET_COV)
EKF_VELNE_CALC_SUBX(__P, __R, __X, __Z, __RET_SUBX)
EKF_VELNE_CALC_INNOV(__P, __R, __SUBX, __X, __Z, __RET_INNOV)
*/

#define EKF_NUM_STATES 7
#define EKF_NUM_CONTROL_INPUTS 3
#define EKF_MAX_NUM_SUBX 88
#define EKF_STATE_IDX_PT_N 0
#define EKF_STATE_IDX_PT_E 1
#define EKF_STATE_IDX_PT_D 2
#define EKF_STATE_IDX_VT_N 3
#define EKF_STATE_IDX_VT_E 4
#define EKF_STATE_IDX_VT_D 5
#define EKF_STATE_IDX_FOC_LEN 6
#define EKF_U_IDX_DVV_N 0
#define EKF_U_IDX_DVV_E 1
#define EKF_U_IDX_DVV_D 2

#define EKF_CAMERA_CALC_NIS(__P, __R, __TBN, __SUBX, __X, __Z, __RET_NIS) \
__RET_NIS = __SUBX[24]*(__SUBX[21]*__SUBX[22]*__SUBX[24] - __SUBX[23]*__SUBX[25]) + \
__SUBX[25]*(__SUBX[20]*__SUBX[22]*__SUBX[25] - __SUBX[23]*__SUBX[24]); 

#define EKF_CAMERA_CALC_STATE(__P, __R, __TBN, __SUBX, __X, __Z, __RET_STATE) \
__RET_STATE[0] = __SUBX[24]*__SUBX[26] + __SUBX[25]*__SUBX[27] + __X[0]; __RET_STATE[1] = \
__SUBX[24]*__SUBX[28] + __SUBX[25]*__SUBX[29] + __X[1]; __RET_STATE[2] = __SUBX[24]*__SUBX[30] + \
__SUBX[25]*__SUBX[31] + __X[2]; __RET_STATE[3] = __SUBX[24]*__SUBX[32] + __SUBX[25]*__SUBX[33] + \
__X[3]; __RET_STATE[4] = __SUBX[24]*__SUBX[34] + __SUBX[25]*__SUBX[35] + __X[4]; __RET_STATE[5] = \
__SUBX[24]*__SUBX[36] + __SUBX[25]*__SUBX[37] + __X[5]; __RET_STATE[6] = __SUBX[24]*__SUBX[38] + \
__SUBX[25]*__SUBX[39] + __X[6]; 

#define EKF_CAMERA_CALC_COV(__P, __R, __TBN, __SUBX, __X, __Z, __RET_COV) \
__RET_COV[0] = __R*__SUBX[26]*__SUBX[26] + __R*__SUBX[27]*__SUBX[27] + __SUBX[40]*__SUBX[44] + \
__SUBX[41]*__SUBX[47] + __SUBX[42]*__SUBX[45] + __SUBX[43]*__SUBX[46]; __RET_COV[1] = \
__R*__SUBX[26]*__SUBX[28] + __R*__SUBX[27]*__SUBX[29] + __SUBX[44]*__SUBX[51] + __SUBX[45]*__SUBX[48] \
+ __SUBX[46]*__SUBX[50] + __SUBX[47]*__SUBX[49]; __RET_COV[2] = __R*__SUBX[26]*__SUBX[30] + \
__R*__SUBX[27]*__SUBX[31] + __SUBX[44]*__SUBX[52] + __SUBX[45]*__SUBX[55] + __SUBX[46]*__SUBX[54] + \
__SUBX[47]*__SUBX[53]; __RET_COV[3] = __P[14]*__SUBX[42] + __P[21]*__SUBX[43] + __P[3]*__SUBX[41] + \
__P[9]*__SUBX[40] + __R*__SUBX[26]*__SUBX[32] + __R*__SUBX[27]*__SUBX[33] + __SUBX[44]*__SUBX[56] + \
__SUBX[45]*__SUBX[57] + __SUBX[46]*__SUBX[59] + __SUBX[47]*__SUBX[58]; __RET_COV[4] = \
__P[10]*__SUBX[40] + __P[15]*__SUBX[42] + __P[24]*__SUBX[43] + __P[4]*__SUBX[41] + \
__R*__SUBX[26]*__SUBX[34] + __R*__SUBX[27]*__SUBX[35] + __SUBX[44]*__SUBX[60] + __SUBX[45]*__SUBX[61] \
+ __SUBX[46]*__SUBX[63] + __SUBX[47]*__SUBX[62]; __RET_COV[5] = __P[11]*__SUBX[40] + \
__P[16]*__SUBX[42] + __P[26]*__SUBX[43] + __P[5]*__SUBX[41] + __R*__SUBX[26]*__SUBX[36] + \
__R*__SUBX[27]*__SUBX[37] + __SUBX[44]*__SUBX[64] + __SUBX[45]*__SUBX[65] + __SUBX[46]*__SUBX[67] + \
__SUBX[47]*__SUBX[66]; __RET_COV[6] = __R*__SUBX[26]*__SUBX[38] + __R*__SUBX[27]*__SUBX[39] + \
__SUBX[44]*__SUBX[68] + __SUBX[45]*__SUBX[69] + __SUBX[46]*__SUBX[71] + __SUBX[47]*__SUBX[70]; \
__RET_COV[7] = __R*__SUBX[28]*__SUBX[28] + __R*__SUBX[29]*__SUBX[29] + __SUBX[48]*__SUBX[72] + \
__SUBX[49]*__SUBX[73] + __SUBX[50]*__SUBX[74] + __SUBX[51]*__SUBX[75]; __RET_COV[8] = \
__R*__SUBX[28]*__SUBX[30] + __R*__SUBX[29]*__SUBX[31] + __SUBX[52]*__SUBX[75] + __SUBX[53]*__SUBX[73] \
+ __SUBX[54]*__SUBX[74] + __SUBX[55]*__SUBX[72]; __RET_COV[9] = __P[14]*__SUBX[48] + \
__P[21]*__SUBX[50] + __P[3]*__SUBX[49] + __P[9]*__SUBX[51] + __R*__SUBX[28]*__SUBX[32] + \
__R*__SUBX[29]*__SUBX[33] + __SUBX[56]*__SUBX[75] + __SUBX[57]*__SUBX[72] + __SUBX[58]*__SUBX[73] + \
__SUBX[59]*__SUBX[74]; __RET_COV[10] = __P[10]*__SUBX[51] + __P[15]*__SUBX[48] + __P[24]*__SUBX[50] + \
__P[4]*__SUBX[49] + __R*__SUBX[28]*__SUBX[34] + __R*__SUBX[29]*__SUBX[35] + __SUBX[60]*__SUBX[75] + \
__SUBX[61]*__SUBX[72] + __SUBX[62]*__SUBX[73] + __SUBX[63]*__SUBX[74]; __RET_COV[11] = \
__P[11]*__SUBX[51] + __P[16]*__SUBX[48] + __P[26]*__SUBX[50] + __P[5]*__SUBX[49] + \
__R*__SUBX[28]*__SUBX[36] + __R*__SUBX[29]*__SUBX[37] + __SUBX[64]*__SUBX[75] + __SUBX[65]*__SUBX[72] \
+ __SUBX[66]*__SUBX[73] + __SUBX[67]*__SUBX[74]; __RET_COV[12] = __R*__SUBX[28]*__SUBX[38] + \
__R*__SUBX[29]*__SUBX[39] + __SUBX[68]*__SUBX[75] + __SUBX[69]*__SUBX[72] + __SUBX[70]*__SUBX[73] + \
__SUBX[71]*__SUBX[74]; __RET_COV[13] = __R*__SUBX[30]*__SUBX[30] + __R*__SUBX[31]*__SUBX[31] + \
__SUBX[52]*__SUBX[76] + __SUBX[53]*__SUBX[77] + __SUBX[54]*__SUBX[78] + __SUBX[55]*__SUBX[79]; \
__RET_COV[14] = __P[14]*__SUBX[55] + __P[21]*__SUBX[54] + __P[3]*__SUBX[53] + __P[9]*__SUBX[52] + \
__R*__SUBX[30]*__SUBX[32] + __R*__SUBX[31]*__SUBX[33] + __SUBX[56]*__SUBX[76] + __SUBX[57]*__SUBX[79] \
+ __SUBX[58]*__SUBX[77] + __SUBX[59]*__SUBX[78]; __RET_COV[15] = __P[10]*__SUBX[52] + \
__P[15]*__SUBX[55] + __P[24]*__SUBX[54] + __P[4]*__SUBX[53] + __R*__SUBX[30]*__SUBX[34] + \
__R*__SUBX[31]*__SUBX[35] + __SUBX[60]*__SUBX[76] + __SUBX[61]*__SUBX[79] + __SUBX[62]*__SUBX[77] + \
__SUBX[63]*__SUBX[78]; __RET_COV[16] = __P[11]*__SUBX[52] + __P[16]*__SUBX[55] + __P[26]*__SUBX[54] + \
__P[5]*__SUBX[53] + __R*__SUBX[30]*__SUBX[36] + __R*__SUBX[31]*__SUBX[37] + __SUBX[64]*__SUBX[76] + \
__SUBX[65]*__SUBX[79] + __SUBX[66]*__SUBX[77] + __SUBX[67]*__SUBX[78]; __RET_COV[17] = \
__R*__SUBX[30]*__SUBX[38] + __R*__SUBX[31]*__SUBX[39] + __SUBX[68]*__SUBX[76] + __SUBX[69]*__SUBX[79] \
+ __SUBX[70]*__SUBX[77] + __SUBX[71]*__SUBX[78]; __RET_COV[18] = __P[14]*__SUBX[57] + __P[18] + \
__P[21]*__SUBX[59] + __P[3]*__SUBX[58] + __P[9]*__SUBX[56] + __R*__SUBX[32]*__SUBX[32] + \
__R*__SUBX[33]*__SUBX[33] + __SUBX[56]*__SUBX[80] + __SUBX[57]*__SUBX[81] + __SUBX[58]*__SUBX[82] + \
__SUBX[59]*__SUBX[83]; __RET_COV[19] = __P[10]*__SUBX[56] + __P[15]*__SUBX[57] + __P[19] + \
__P[24]*__SUBX[59] + __P[4]*__SUBX[58] + __R*__SUBX[32]*__SUBX[34] + __R*__SUBX[33]*__SUBX[35] + \
__SUBX[60]*__SUBX[80] + __SUBX[61]*__SUBX[81] + __SUBX[62]*__SUBX[82] + __SUBX[63]*__SUBX[83]; \
__RET_COV[20] = __P[11]*__SUBX[56] + __P[16]*__SUBX[57] + __P[20] + __P[26]*__SUBX[59] + \
__P[5]*__SUBX[58] + __R*__SUBX[32]*__SUBX[36] + __R*__SUBX[33]*__SUBX[37] + __SUBX[64]*__SUBX[80] + \
__SUBX[65]*__SUBX[81] + __SUBX[66]*__SUBX[82] + __SUBX[67]*__SUBX[83]; __RET_COV[21] = \
__R*__SUBX[32]*__SUBX[38] + __R*__SUBX[33]*__SUBX[39] + __SUBX[68]*__SUBX[80] + __SUBX[69]*__SUBX[81] \
+ __SUBX[70]*__SUBX[82] + __SUBX[71]*__SUBX[83]; __RET_COV[22] = __P[10]*__SUBX[60] + \
__P[15]*__SUBX[61] + __P[22] + __P[24]*__SUBX[63] + __P[4]*__SUBX[62] + __R*__SUBX[34]*__SUBX[34] + \
__R*__SUBX[35]*__SUBX[35] + __SUBX[60]*__SUBX[84] + __SUBX[61]*__SUBX[85] + __SUBX[62]*__SUBX[86] + \
__SUBX[63]*__SUBX[87]; __RET_COV[23] = __P[11]*__SUBX[60] + __P[16]*__SUBX[61] + __P[23] + \
__P[26]*__SUBX[63] + __P[5]*__SUBX[62] + __R*__SUBX[34]*__SUBX[36] + __R*__SUBX[35]*__SUBX[37] + \
__SUBX[64]*__SUBX[84] + __SUBX[65]*__SUBX[85] + __SUBX[66]*__SUBX[86] + __SUBX[67]*__SUBX[87]; \
__RET_COV[24] = __R*__SUBX[34]*__SUBX[38] + __R*__SUBX[35]*__SUBX[39] + __SUBX[68]*__SUBX[84] + \
__SUBX[69]*__SUBX[85] + __SUBX[70]*__SUBX[86] + __SUBX[71]*__SUBX[87]; __RET_COV[25] = \
__P[11]*__SUBX[64] + __P[16]*__SUBX[65] + __P[25] + __P[26]*__SUBX[67] + __P[5]*__SUBX[66] + \
__R*__SUBX[36]*__SUBX[36] + __R*__SUBX[37]*__SUBX[37] + __SUBX[64]*(__P[11] + __P[12]*__SUBX[67] + \
__P[1]*__SUBX[66] + __P[7]*__SUBX[64] + __P[8]*__SUBX[65]) + __SUBX[65]*(__P[13]*__SUBX[65] + __P[16] \
+ __P[17]*__SUBX[67] + __P[2]*__SUBX[66] + __P[8]*__SUBX[64]) + __SUBX[66]*(__P[0]*__SUBX[66] + \
__P[1]*__SUBX[64] + __P[2]*__SUBX[65] + __P[5] + __P[6]*__SUBX[67]) + __SUBX[67]*(__P[12]*__SUBX[64] \
+ __P[17]*__SUBX[65] + __P[26] + __P[27]*__SUBX[67] + __P[6]*__SUBX[66]); __RET_COV[26] = \
__R*__SUBX[36]*__SUBX[38] + __R*__SUBX[37]*__SUBX[39] + __SUBX[68]*(__P[11] + __P[12]*__SUBX[67] + \
__P[1]*__SUBX[66] + __P[7]*__SUBX[64] + __P[8]*__SUBX[65]) + __SUBX[69]*(__P[13]*__SUBX[65] + __P[16] \
+ __P[17]*__SUBX[67] + __P[2]*__SUBX[66] + __P[8]*__SUBX[64]) + __SUBX[70]*(__P[0]*__SUBX[66] + \
__P[1]*__SUBX[64] + __P[2]*__SUBX[65] + __P[5] + __P[6]*__SUBX[67]) + __SUBX[71]*(__P[12]*__SUBX[64] \
+ __P[17]*__SUBX[65] + __P[26] + __P[27]*__SUBX[67] + __P[6]*__SUBX[66]); __RET_COV[27] = \
__R*__SUBX[38]*__SUBX[38] + __R*__SUBX[39]*__SUBX[39] + __SUBX[68]*(__P[12]*__SUBX[71] + \
__P[1]*__SUBX[70] + __P[7]*__SUBX[68] + __P[8]*__SUBX[69]) + __SUBX[69]*(__P[13]*__SUBX[69] + \
__P[17]*__SUBX[71] + __P[2]*__SUBX[70] + __P[8]*__SUBX[68]) + __SUBX[70]*(__P[0]*__SUBX[70] + \
__P[1]*__SUBX[68] + __P[2]*__SUBX[69] + __P[6]*__SUBX[71]) + __SUBX[71]*(__P[12]*__SUBX[68] + \
__P[17]*__SUBX[69] + __P[27]*__SUBX[71] + __P[6]*__SUBX[70]); 

#define EKF_CAMERA_CALC_SUBX(__P, __R, __TBN, __X, __Z, __RET_SUBX) \
__RET_SUBX[0] = __TBN[0][2]*__X[0] + __TBN[1][2]*__X[1] + __TBN[2][2]*__X[2]; __RET_SUBX[1] = \
__TBN[0][0]*__X[0] + __TBN[1][0]*__X[1] + __TBN[2][0]*__X[2]; __RET_SUBX[2] = \
__RET_SUBX[1]/__RET_SUBX[0]; __RET_SUBX[3] = __TBN[0][1]*__X[0] + __TBN[1][1]*__X[1] + \
__TBN[2][1]*__X[2]; __RET_SUBX[4] = __RET_SUBX[3]/__RET_SUBX[0]; __RET_SUBX[5] = \
__RET_SUBX[0]*__RET_SUBX[0]; __RET_SUBX[6] = -__RET_SUBX[1]*__RET_SUBX[5]*__TBN[0][2]*__X[6] + \
__TBN[0][0]*__X[6]/__RET_SUBX[0]; __RET_SUBX[7] = -__RET_SUBX[1]*__RET_SUBX[5]*__TBN[1][2]*__X[6] + \
__TBN[1][0]*__X[6]/__RET_SUBX[0]; __RET_SUBX[8] = -__RET_SUBX[1]*__RET_SUBX[5]*__TBN[2][2]*__X[6] + \
__TBN[2][0]*__X[6]/__RET_SUBX[0]; __RET_SUBX[9] = -__RET_SUBX[3]*__RET_SUBX[5]*__TBN[0][2]*__X[6] + \
__TBN[0][1]*__X[6]/__RET_SUBX[0]; __RET_SUBX[10] = -__RET_SUBX[3]*__RET_SUBX[5]*__TBN[1][2]*__X[6] + \
__TBN[1][1]*__X[6]/__RET_SUBX[0]; __RET_SUBX[11] = -__RET_SUBX[3]*__RET_SUBX[5]*__TBN[2][2]*__X[6] + \
__TBN[2][1]*__X[6]/__RET_SUBX[0]; __RET_SUBX[12] = __P[0]*__RET_SUBX[9] + __P[1]*__RET_SUBX[10] + \
__P[2]*__RET_SUBX[11] + __P[6]*__RET_SUBX[4]; __RET_SUBX[13] = __P[12]*__RET_SUBX[4] + \
__P[1]*__RET_SUBX[9] + __P[7]*__RET_SUBX[10] + __P[8]*__RET_SUBX[11]; __RET_SUBX[14] = \
__P[13]*__RET_SUBX[11] + __P[17]*__RET_SUBX[4] + __P[2]*__RET_SUBX[9] + __P[8]*__RET_SUBX[10]; \
__RET_SUBX[15] = __P[12]*__RET_SUBX[10] + __P[17]*__RET_SUBX[11] + __P[27]*__RET_SUBX[4] + \
__P[6]*__RET_SUBX[9]; __RET_SUBX[16] = __P[0]*__RET_SUBX[6] + __P[1]*__RET_SUBX[7] + \
__P[2]*__RET_SUBX[8] + __P[6]*__RET_SUBX[2]; __RET_SUBX[17] = __P[12]*__RET_SUBX[2] + \
__P[1]*__RET_SUBX[6] + __P[7]*__RET_SUBX[7] + __P[8]*__RET_SUBX[8]; __RET_SUBX[18] = \
__P[13]*__RET_SUBX[8] + __P[17]*__RET_SUBX[2] + __P[2]*__RET_SUBX[6] + __P[8]*__RET_SUBX[7]; \
__RET_SUBX[19] = __P[12]*__RET_SUBX[7] + __P[17]*__RET_SUBX[8] + __P[27]*__RET_SUBX[2] + \
__P[6]*__RET_SUBX[6]; __RET_SUBX[20] = __R + __RET_SUBX[16]*__RET_SUBX[6] + \
__RET_SUBX[17]*__RET_SUBX[7] + __RET_SUBX[18]*__RET_SUBX[8] + __RET_SUBX[19]*__RET_SUBX[2]; \
__RET_SUBX[21] = __R + __RET_SUBX[10]*__RET_SUBX[13] + __RET_SUBX[11]*__RET_SUBX[14] + \
__RET_SUBX[12]*__RET_SUBX[9] + __RET_SUBX[15]*__RET_SUBX[4]; __RET_SUBX[22] = \
1.0f/(__RET_SUBX[20]*__RET_SUBX[21] - __RET_SUBX[12]*__RET_SUBX[6] + __RET_SUBX[13]*__RET_SUBX[7] + \
__RET_SUBX[14]*__RET_SUBX[8] + __RET_SUBX[15]*__RET_SUBX[2]*__RET_SUBX[12]*__RET_SUBX[6] + \
__RET_SUBX[13]*__RET_SUBX[7] + __RET_SUBX[14]*__RET_SUBX[8] + __RET_SUBX[15]*__RET_SUBX[2]); \
__RET_SUBX[23] = __RET_SUBX[22]*(__RET_SUBX[12]*__RET_SUBX[6] + __RET_SUBX[13]*__RET_SUBX[7] + \
__RET_SUBX[14]*__RET_SUBX[8] + __RET_SUBX[15]*__RET_SUBX[2]); __RET_SUBX[24] = -__RET_SUBX[2]*__X[6] \
+ __Z[0]; __RET_SUBX[25] = -__RET_SUBX[4]*__X[6] + __Z[1]; __RET_SUBX[26] = \
-__RET_SUBX[12]*__RET_SUBX[23] + __RET_SUBX[16]*__RET_SUBX[21]*__RET_SUBX[22]; __RET_SUBX[27] = \
__RET_SUBX[12]*__RET_SUBX[20]*__RET_SUBX[22] - __RET_SUBX[16]*__RET_SUBX[23]; __RET_SUBX[28] = \
-__RET_SUBX[13]*__RET_SUBX[23] + __RET_SUBX[17]*__RET_SUBX[21]*__RET_SUBX[22]; __RET_SUBX[29] = \
__RET_SUBX[13]*__RET_SUBX[20]*__RET_SUBX[22] - __RET_SUBX[17]*__RET_SUBX[23]; __RET_SUBX[30] = \
-__RET_SUBX[14]*__RET_SUBX[23] + __RET_SUBX[18]*__RET_SUBX[21]*__RET_SUBX[22]; __RET_SUBX[31] = \
__RET_SUBX[14]*__RET_SUBX[20]*__RET_SUBX[22] - __RET_SUBX[18]*__RET_SUBX[23]; __RET_SUBX[32] = \
__RET_SUBX[21]*__RET_SUBX[22]*(__P[14]*__RET_SUBX[8] + __P[21]*__RET_SUBX[2] + __P[3]*__RET_SUBX[6] + \
__P[9]*__RET_SUBX[7]) - __RET_SUBX[23]*(__P[14]*__RET_SUBX[11] + __P[21]*__RET_SUBX[4] + \
__P[3]*__RET_SUBX[9] + __P[9]*__RET_SUBX[10]); __RET_SUBX[33] = \
__RET_SUBX[20]*__RET_SUBX[22]*(__P[14]*__RET_SUBX[11] + __P[21]*__RET_SUBX[4] + __P[3]*__RET_SUBX[9] \
+ __P[9]*__RET_SUBX[10]) - __RET_SUBX[23]*(__P[14]*__RET_SUBX[8] + __P[21]*__RET_SUBX[2] + \
__P[3]*__RET_SUBX[6] + __P[9]*__RET_SUBX[7]); __RET_SUBX[34] = \
__RET_SUBX[21]*__RET_SUBX[22]*(__P[10]*__RET_SUBX[7] + __P[15]*__RET_SUBX[8] + __P[24]*__RET_SUBX[2] \
+ __P[4]*__RET_SUBX[6]) - __RET_SUBX[23]*(__P[10]*__RET_SUBX[10] + __P[15]*__RET_SUBX[11] + \
__P[24]*__RET_SUBX[4] + __P[4]*__RET_SUBX[9]); __RET_SUBX[35] = \
__RET_SUBX[20]*__RET_SUBX[22]*(__P[10]*__RET_SUBX[10] + __P[15]*__RET_SUBX[11] + \
__P[24]*__RET_SUBX[4] + __P[4]*__RET_SUBX[9]) - __RET_SUBX[23]*(__P[10]*__RET_SUBX[7] + \
__P[15]*__RET_SUBX[8] + __P[24]*__RET_SUBX[2] + __P[4]*__RET_SUBX[6]); __RET_SUBX[36] = \
__RET_SUBX[21]*__RET_SUBX[22]*(__P[11]*__RET_SUBX[7] + __P[16]*__RET_SUBX[8] + __P[26]*__RET_SUBX[2] \
+ __P[5]*__RET_SUBX[6]) - __RET_SUBX[23]*(__P[11]*__RET_SUBX[10] + __P[16]*__RET_SUBX[11] + \
__P[26]*__RET_SUBX[4] + __P[5]*__RET_SUBX[9]); __RET_SUBX[37] = \
__RET_SUBX[20]*__RET_SUBX[22]*(__P[11]*__RET_SUBX[10] + __P[16]*__RET_SUBX[11] + \
__P[26]*__RET_SUBX[4] + __P[5]*__RET_SUBX[9]) - __RET_SUBX[23]*(__P[11]*__RET_SUBX[7] + \
__P[16]*__RET_SUBX[8] + __P[26]*__RET_SUBX[2] + __P[5]*__RET_SUBX[6]); __RET_SUBX[38] = \
-__RET_SUBX[15]*__RET_SUBX[23] + __RET_SUBX[19]*__RET_SUBX[21]*__RET_SUBX[22]; __RET_SUBX[39] = \
__RET_SUBX[15]*__RET_SUBX[20]*__RET_SUBX[22] - __RET_SUBX[19]*__RET_SUBX[23]; __RET_SUBX[40] = \
-__RET_SUBX[10]*__RET_SUBX[27] - __RET_SUBX[26]*__RET_SUBX[7]; __RET_SUBX[41] = \
-__RET_SUBX[26]*__RET_SUBX[6] - __RET_SUBX[27]*__RET_SUBX[9] + 1; __RET_SUBX[42] = \
-__RET_SUBX[11]*__RET_SUBX[27] - __RET_SUBX[26]*__RET_SUBX[8]; __RET_SUBX[43] = \
-__RET_SUBX[26]*__RET_SUBX[2] - __RET_SUBX[27]*__RET_SUBX[4]; __RET_SUBX[44] = __P[12]*__RET_SUBX[43] \
+ __P[1]*__RET_SUBX[41] + __P[7]*__RET_SUBX[40] + __P[8]*__RET_SUBX[42]; __RET_SUBX[45] = \
__P[13]*__RET_SUBX[42] + __P[17]*__RET_SUBX[43] + __P[2]*__RET_SUBX[41] + __P[8]*__RET_SUBX[40]; \
__RET_SUBX[46] = __P[12]*__RET_SUBX[40] + __P[17]*__RET_SUBX[42] + __P[27]*__RET_SUBX[43] + \
__P[6]*__RET_SUBX[41]; __RET_SUBX[47] = __P[0]*__RET_SUBX[41] + __P[1]*__RET_SUBX[40] + \
__P[2]*__RET_SUBX[42] + __P[6]*__RET_SUBX[43]; __RET_SUBX[48] = -__RET_SUBX[11]*__RET_SUBX[29] - \
__RET_SUBX[28]*__RET_SUBX[8]; __RET_SUBX[49] = -__RET_SUBX[28]*__RET_SUBX[6] - \
__RET_SUBX[29]*__RET_SUBX[9]; __RET_SUBX[50] = -__RET_SUBX[28]*__RET_SUBX[2] - \
__RET_SUBX[29]*__RET_SUBX[4]; __RET_SUBX[51] = -__RET_SUBX[10]*__RET_SUBX[29] - \
__RET_SUBX[28]*__RET_SUBX[7] + 1; __RET_SUBX[52] = -__RET_SUBX[10]*__RET_SUBX[31] - \
__RET_SUBX[30]*__RET_SUBX[7]; __RET_SUBX[53] = -__RET_SUBX[30]*__RET_SUBX[6] - \
__RET_SUBX[31]*__RET_SUBX[9]; __RET_SUBX[54] = -__RET_SUBX[2]*__RET_SUBX[30] - \
__RET_SUBX[31]*__RET_SUBX[4]; __RET_SUBX[55] = -__RET_SUBX[11]*__RET_SUBX[31] - \
__RET_SUBX[30]*__RET_SUBX[8] + 1; __RET_SUBX[56] = -__RET_SUBX[10]*__RET_SUBX[33] - \
__RET_SUBX[32]*__RET_SUBX[7]; __RET_SUBX[57] = -__RET_SUBX[11]*__RET_SUBX[33] - \
__RET_SUBX[32]*__RET_SUBX[8]; __RET_SUBX[58] = -__RET_SUBX[32]*__RET_SUBX[6] - \
__RET_SUBX[33]*__RET_SUBX[9]; __RET_SUBX[59] = -__RET_SUBX[2]*__RET_SUBX[32] - \
__RET_SUBX[33]*__RET_SUBX[4]; __RET_SUBX[60] = -__RET_SUBX[10]*__RET_SUBX[35] - \
__RET_SUBX[34]*__RET_SUBX[7]; __RET_SUBX[61] = -__RET_SUBX[11]*__RET_SUBX[35] - \
__RET_SUBX[34]*__RET_SUBX[8]; __RET_SUBX[62] = -__RET_SUBX[34]*__RET_SUBX[6] - \
__RET_SUBX[35]*__RET_SUBX[9]; __RET_SUBX[63] = -__RET_SUBX[2]*__RET_SUBX[34] - \
__RET_SUBX[35]*__RET_SUBX[4]; __RET_SUBX[64] = -__RET_SUBX[10]*__RET_SUBX[37] - \
__RET_SUBX[36]*__RET_SUBX[7]; __RET_SUBX[65] = -__RET_SUBX[11]*__RET_SUBX[37] - \
__RET_SUBX[36]*__RET_SUBX[8]; __RET_SUBX[66] = -__RET_SUBX[36]*__RET_SUBX[6] - \
__RET_SUBX[37]*__RET_SUBX[9]; __RET_SUBX[67] = -__RET_SUBX[2]*__RET_SUBX[36] - \
__RET_SUBX[37]*__RET_SUBX[4]; __RET_SUBX[68] = -__RET_SUBX[10]*__RET_SUBX[39] - \
__RET_SUBX[38]*__RET_SUBX[7]; __RET_SUBX[69] = -__RET_SUBX[11]*__RET_SUBX[39] - \
__RET_SUBX[38]*__RET_SUBX[8]; __RET_SUBX[70] = -__RET_SUBX[38]*__RET_SUBX[6] - \
__RET_SUBX[39]*__RET_SUBX[9]; __RET_SUBX[71] = -__RET_SUBX[2]*__RET_SUBX[38] - \
__RET_SUBX[39]*__RET_SUBX[4] + 1; __RET_SUBX[72] = __P[13]*__RET_SUBX[48] + __P[17]*__RET_SUBX[50] + \
__P[2]*__RET_SUBX[49] + __P[8]*__RET_SUBX[51]; __RET_SUBX[73] = __P[0]*__RET_SUBX[49] + \
__P[1]*__RET_SUBX[51] + __P[2]*__RET_SUBX[48] + __P[6]*__RET_SUBX[50]; __RET_SUBX[74] = \
__P[12]*__RET_SUBX[51] + __P[17]*__RET_SUBX[48] + __P[27]*__RET_SUBX[50] + __P[6]*__RET_SUBX[49]; \
__RET_SUBX[75] = __P[12]*__RET_SUBX[50] + __P[1]*__RET_SUBX[49] + __P[7]*__RET_SUBX[51] + \
__P[8]*__RET_SUBX[48]; __RET_SUBX[76] = __P[12]*__RET_SUBX[54] + __P[1]*__RET_SUBX[53] + \
__P[7]*__RET_SUBX[52] + __P[8]*__RET_SUBX[55]; __RET_SUBX[77] = __P[0]*__RET_SUBX[53] + \
__P[1]*__RET_SUBX[52] + __P[2]*__RET_SUBX[55] + __P[6]*__RET_SUBX[54]; __RET_SUBX[78] = \
__P[12]*__RET_SUBX[52] + __P[17]*__RET_SUBX[55] + __P[27]*__RET_SUBX[54] + __P[6]*__RET_SUBX[53]; \
__RET_SUBX[79] = __P[13]*__RET_SUBX[55] + __P[17]*__RET_SUBX[54] + __P[2]*__RET_SUBX[53] + \
__P[8]*__RET_SUBX[52]; __RET_SUBX[80] = __P[12]*__RET_SUBX[59] + __P[1]*__RET_SUBX[58] + \
__P[7]*__RET_SUBX[56] + __P[8]*__RET_SUBX[57] + __P[9]; __RET_SUBX[81] = __P[13]*__RET_SUBX[57] + \
__P[14] + __P[17]*__RET_SUBX[59] + __P[2]*__RET_SUBX[58] + __P[8]*__RET_SUBX[56]; __RET_SUBX[82] = \
__P[0]*__RET_SUBX[58] + __P[1]*__RET_SUBX[56] + __P[2]*__RET_SUBX[57] + __P[3] + \
__P[6]*__RET_SUBX[59]; __RET_SUBX[83] = __P[12]*__RET_SUBX[56] + __P[17]*__RET_SUBX[57] + __P[21] + \
__P[27]*__RET_SUBX[59] + __P[6]*__RET_SUBX[58]; __RET_SUBX[84] = __P[10] + __P[12]*__RET_SUBX[63] + \
__P[1]*__RET_SUBX[62] + __P[7]*__RET_SUBX[60] + __P[8]*__RET_SUBX[61]; __RET_SUBX[85] = \
__P[13]*__RET_SUBX[61] + __P[15] + __P[17]*__RET_SUBX[63] + __P[2]*__RET_SUBX[62] + \
__P[8]*__RET_SUBX[60]; __RET_SUBX[86] = __P[0]*__RET_SUBX[62] + __P[1]*__RET_SUBX[60] + \
__P[2]*__RET_SUBX[61] + __P[4] + __P[6]*__RET_SUBX[63]; __RET_SUBX[87] = __P[12]*__RET_SUBX[60] + \
__P[17]*__RET_SUBX[61] + __P[24] + __P[27]*__RET_SUBX[63] + __P[6]*__RET_SUBX[62]; 

#define EKF_CAMERA_CALC_INNOV(__P, __R, __TBN, __SUBX, __X, __Z, __RET_INNOV) \
__RET_INNOV[0] = __SUBX[24]; __RET_INNOV[1] = __SUBX[25]; 

#define EKF_HEIGHT_CALC_NIS(__P, __R, __SUBX, __X, __Z, __RET_NIS) \
__RET_NIS = __SUBX[0]*-__X[2] + __Z*-__X[2] + __Z; 

#define EKF_HEIGHT_CALC_STATE(__P, __R, __SUBX, __X, __Z, __RET_STATE) \
__RET_STATE[0] = __P[2]*__SUBX[0]*(-__X[2] + __Z) + __X[0]; __RET_STATE[1] = \
__P[8]*__SUBX[0]*(-__X[2] + __Z) + __X[1]; __RET_STATE[2] = __SUBX[1]*(-__X[2] + __Z) + __X[2]; \
__RET_STATE[3] = __P[14]*__SUBX[0]*(-__X[2] + __Z) + __X[3]; __RET_STATE[4] = \
__P[15]*__SUBX[0]*(-__X[2] + __Z) + __X[4]; __RET_STATE[5] = __P[16]*__SUBX[0]*(-__X[2] + __Z) + \
__X[5]; __RET_STATE[6] = __P[17]*__SUBX[0]*(-__X[2] + __Z) + __X[6]; 

#define EKF_HEIGHT_CALC_COV(__P, __R, __SUBX, __X, __Z, __RET_COV) \
__RET_COV[0] = __P[0] + __P[2]*__P[2]*__R*__SUBX[2] - __P[2]*__P[2]*__SUBX[0] - \
__P[2]*__SUBX[0]*__SUBX[3]; __RET_COV[1] = __P[1] + __P[2]*__P[8]*__R*__SUBX[2] - \
__P[2]*__P[8]*__SUBX[0] - __P[8]*__SUBX[0]*__SUBX[3]; __RET_COV[2] = __P[2]*__SUBX[5] + \
__SUBX[3]*__SUBX[4]; __RET_COV[3] = __P[14]*__P[2]*__R*__SUBX[2] - __P[14]*__P[2]*__SUBX[0] - \
__P[14]*__SUBX[0]*__SUBX[3] + __P[3]; __RET_COV[4] = __P[15]*__P[2]*__R*__SUBX[2] - \
__P[15]*__P[2]*__SUBX[0] - __P[15]*__SUBX[0]*__SUBX[3] + __P[4]; __RET_COV[5] = \
__P[16]*__P[2]*__R*__SUBX[2] - __P[16]*__P[2]*__SUBX[0] - __P[16]*__SUBX[0]*__SUBX[3] + __P[5]; \
__RET_COV[6] = __P[17]*__P[2]*__R*__SUBX[2] - __P[17]*__P[2]*__SUBX[0] - __P[17]*__SUBX[0]*__SUBX[3] \
+ __P[6]; __RET_COV[7] = __P[7] + __P[8]*__P[8]*__R*__SUBX[2] - __P[8]*__P[8]*__SUBX[0] - \
__P[8]*__SUBX[0]*__SUBX[6]; __RET_COV[8] = __P[8]*__SUBX[5] + __SUBX[4]*__SUBX[6]; __RET_COV[9] = \
__P[14]*__P[8]*__R*__SUBX[2] - __P[14]*__P[8]*__SUBX[0] - __P[14]*__SUBX[0]*__SUBX[6] + __P[9]; \
__RET_COV[10] = __P[10] + __P[15]*__P[8]*__R*__SUBX[2] - __P[15]*__P[8]*__SUBX[0] - \
__P[15]*__SUBX[0]*__SUBX[6]; __RET_COV[11] = __P[11] + __P[16]*__P[8]*__R*__SUBX[2] - \
__P[16]*__P[8]*__SUBX[0] - __P[16]*__SUBX[0]*__SUBX[6]; __RET_COV[12] = __P[12] + \
__P[17]*__P[8]*__R*__SUBX[2] - __P[17]*__P[8]*__SUBX[0] - __P[17]*__SUBX[0]*__SUBX[6]; __RET_COV[13] \
= __P[13]*__P[13]*__R*__SUBX[2] + __P[13]*__SUBX[4]*__SUBX[4]; __RET_COV[14] = \
-__P[14]*__SUBX[1]*__SUBX[4] + __P[14]*__SUBX[4] + __P[14]*__SUBX[5]; __RET_COV[15] = \
-__P[15]*__SUBX[1]*__SUBX[4] + __P[15]*__SUBX[4] + __P[15]*__SUBX[5]; __RET_COV[16] = \
-__P[16]*__SUBX[1]*__SUBX[4] + __P[16]*__SUBX[4] + __P[16]*__SUBX[5]; __RET_COV[17] = \
-__P[17]*__SUBX[1]*__SUBX[4] + __P[17]*__SUBX[4] + __P[17]*__SUBX[5]; __RET_COV[18] = \
__P[14]*__P[14]*__R*__SUBX[2] - __P[14]*__P[14]*__SUBX[0] - __P[14]*__SUBX[0]*(-__P[14]*__SUBX[1] + \
__P[14]) + __P[18]; __RET_COV[19] = __P[14]*__P[15]*__R*__SUBX[2] - __P[14]*__P[15]*__SUBX[0] - \
__P[15]*__SUBX[0]*(-__P[14]*__SUBX[1] + __P[14]) + __P[19]; __RET_COV[20] = \
__P[14]*__P[16]*__R*__SUBX[2] - __P[14]*__P[16]*__SUBX[0] - __P[16]*__SUBX[0]*(-__P[14]*__SUBX[1] + \
__P[14]) + __P[20]; __RET_COV[21] = __P[14]*__P[17]*__R*__SUBX[2] - __P[14]*__P[17]*__SUBX[0] - \
__P[17]*__SUBX[0]*(-__P[14]*__SUBX[1] + __P[14]) + __P[21]; __RET_COV[22] = \
__P[15]*__P[15]*__R*__SUBX[2] - __P[15]*__P[15]*__SUBX[0] - __P[15]*__SUBX[0]*(-__P[15]*__SUBX[1] + \
__P[15]) + __P[22]; __RET_COV[23] = __P[15]*__P[16]*__R*__SUBX[2] - __P[15]*__P[16]*__SUBX[0] - \
__P[16]*__SUBX[0]*(-__P[15]*__SUBX[1] + __P[15]) + __P[23]; __RET_COV[24] = \
__P[15]*__P[17]*__R*__SUBX[2] - __P[15]*__P[17]*__SUBX[0] - __P[17]*__SUBX[0]*(-__P[15]*__SUBX[1] + \
__P[15]) + __P[24]; __RET_COV[25] = __P[16]*__P[16]*__R*__SUBX[2] - __P[16]*__P[16]*__SUBX[0] - \
__P[16]*__SUBX[0]*(-__P[16]*__SUBX[1] + __P[16]) + __P[25]; __RET_COV[26] = \
__P[16]*__P[17]*__R*__SUBX[2] - __P[16]*__P[17]*__SUBX[0] - __P[17]*__SUBX[0]*(-__P[16]*__SUBX[1] + \
__P[16]) + __P[26]; __RET_COV[27] = __P[17]*__P[17]*__R*__SUBX[2] - __P[17]*__P[17]*__SUBX[0] - \
__P[17]*__SUBX[0]*(-__P[17]*__SUBX[1] + __P[17]) + __P[27]; 

#define EKF_HEIGHT_CALC_SUBX(__P, __R, __X, __Z, __RET_SUBX) \
__RET_SUBX[0] = 1.0f/(__P[13] + __R); __RET_SUBX[1] = __P[13]*__RET_SUBX[0]; __RET_SUBX[2] = \
__RET_SUBX[0]*__RET_SUBX[0]; __RET_SUBX[3] = -__P[2]*__RET_SUBX[1] + __P[2]; __RET_SUBX[4] = \
-__RET_SUBX[1] + 1; __RET_SUBX[5] = __P[13]*__R*__RET_SUBX[2]; __RET_SUBX[6] = -__P[8]*__RET_SUBX[1] \
+ __P[8]; 

#define EKF_HEIGHT_CALC_INNOV(__P, __R, __SUBX, __X, __Z, __RET_INNOV) \
__RET_INNOV = -__X[2] + __Z; 

#define EKF_INITIALIZATION_CALC_STATE(__TBN, __CAM_POS, __CAM_POS_R, __HGT, __HGT_R, __INIT_FOC_LEN, __INIT_FOC_LEN_R, __SUBX, __VEL, __VEL_R, __RET_STATE) \
__RET_STATE[0] = __HGT*__SUBX[3]; __RET_STATE[1] = __SUBX[4]*__SUBX[5]; __RET_STATE[2] = __HGT; \
__RET_STATE[3] = -__VEL[0]; __RET_STATE[4] = -__VEL[1]; __RET_STATE[5] = -__VEL[2]; __RET_STATE[6] = \
__INIT_FOC_LEN; 

#define EKF_INITIALIZATION_CALC_SUBX(__TBN, __CAM_POS, __CAM_POS_R, __HGT, __HGT_R, __INIT_FOC_LEN, __INIT_FOC_LEN_R, __VEL, __VEL_R, __RET_SUBX) \
__RET_SUBX[0] = __CAM_POS[0]*__TBN[2][0] + __CAM_POS[1]*__TBN[2][1] + 1.0f*__TBN[2][2]; __RET_SUBX[1] \
= 1.0f/__RET_SUBX[0]; __RET_SUBX[2] = __CAM_POS[0]*__TBN[0][0] + __CAM_POS[1]*__TBN[0][1] + \
1.0f*__TBN[0][2]; __RET_SUBX[3] = __RET_SUBX[1]*__RET_SUBX[2]; __RET_SUBX[4] = \
__CAM_POS[0]*__TBN[1][0] + __CAM_POS[1]*__TBN[1][1] + 1.0f*__TBN[1][2]; __RET_SUBX[5] = \
__HGT*__RET_SUBX[1]; __RET_SUBX[6] = __RET_SUBX[0]*__RET_SUBX[0]; __RET_SUBX[7] = \
__HGT*__RET_SUBX[2]*__RET_SUBX[6]; __RET_SUBX[8] = __RET_SUBX[5]*__TBN[0][0] - \
__RET_SUBX[7]*__TBN[2][0]; __RET_SUBX[9] = __RET_SUBX[5]*__TBN[0][1] - __RET_SUBX[7]*__TBN[2][1]; \
__RET_SUBX[10] = __HGT_R*__RET_SUBX[6]; __RET_SUBX[11] = __HGT*__RET_SUBX[4]*__RET_SUBX[6]; \
__RET_SUBX[12] = -__RET_SUBX[11]*__TBN[2][0] + __RET_SUBX[5]*__TBN[1][0]; __RET_SUBX[13] = \
-__RET_SUBX[11]*__TBN[2][1] + __RET_SUBX[5]*__TBN[1][1]; 

#define EKF_INITIALIZATION_CALC_COV(__TBN, __CAM_POS, __CAM_POS_R, __HGT, __HGT_R, __INIT_FOC_LEN, __INIT_FOC_LEN_R, __SUBX, __VEL, __VEL_R, __RET_COV) \
__RET_COV[0] = __CAM_POS_R*__SUBX[8]*__SUBX[8] + __CAM_POS_R*__SUBX[9]*__SUBX[9] + \
__SUBX[10]*__SUBX[2]*__SUBX[2]; __RET_COV[1] = __CAM_POS_R*__SUBX[12]*__SUBX[8] + \
__CAM_POS_R*__SUBX[13]*__SUBX[9] + __SUBX[10]*__SUBX[2]*__SUBX[4]; __RET_COV[2] = __HGT_R*__SUBX[3]; \
__RET_COV[3] = 0; __RET_COV[4] = 0; __RET_COV[5] = 0; __RET_COV[6] = 0; __RET_COV[7] = \
__CAM_POS_R*__SUBX[12]*__SUBX[12] + __CAM_POS_R*__SUBX[13]*__SUBX[13] + \
__SUBX[10]*__SUBX[4]*__SUBX[4]; __RET_COV[8] = __HGT_R*__SUBX[1]*__SUBX[4]; __RET_COV[9] = 0; \
__RET_COV[10] = 0; __RET_COV[11] = 0; __RET_COV[12] = 0; __RET_COV[13] = __HGT_R; __RET_COV[14] = 0; \
__RET_COV[15] = 0; __RET_COV[16] = 0; __RET_COV[17] = 0; __RET_COV[18] = __VEL_R[0]; __RET_COV[19] = \
0; __RET_COV[20] = 0; __RET_COV[21] = 0; __RET_COV[22] = __VEL_R[1]; __RET_COV[23] = 0; __RET_COV[24] \
= 0; __RET_COV[25] = __VEL_R[2]; __RET_COV[26] = 0; __RET_COV[27] = __INIT_FOC_LEN_R; 

#define EKF_PREDICTION_CALC_STATE(__P, __DT, __SUBX, __U, __W_U_SIGMA, __X, __RET_STATE) \
__RET_STATE[0] = __DT*__X[3] + __X[0]; __RET_STATE[1] = __DT*__X[4] + __X[1]; __RET_STATE[2] = \
__DT*__X[5] + __X[2]; __RET_STATE[3] = -__U[0] + __X[3]; __RET_STATE[4] = -__U[1] + __X[4]; \
__RET_STATE[5] = -__U[2] + __X[5]; __RET_STATE[6] = __X[6]; 

#define EKF_PREDICTION_CALC_SUBX(__P, __DT, __U, __W_U_SIGMA, __X, __RET_SUBX) \


#define EKF_PREDICTION_CALC_COV(__P, __DT, __SUBX, __U, __W_U_SIGMA, __X, __RET_COV) \
__RET_COV[0] = __DT*__P[3] + __DT*(__DT*__P[18] + __P[3]) + __P[0]; __RET_COV[1] = __DT*__P[9] + \
__DT*(__DT*__P[19] + __P[4]) + __P[1]; __RET_COV[2] = __DT*__P[14] + __DT*(__DT*__P[20] + __P[5]) + \
__P[2]; __RET_COV[3] = __DT*__P[18] + __P[3]; __RET_COV[4] = __DT*__P[19] + __P[4]; __RET_COV[5] = \
__DT*__P[20] + __P[5]; __RET_COV[6] = __DT*__P[21] + __P[6]; __RET_COV[7] = __DT*__P[10] + \
__DT*(__DT*__P[22] + __P[10]) + __P[7]; __RET_COV[8] = __DT*__P[15] + __DT*(__DT*__P[23] + __P[11]) + \
__P[8]; __RET_COV[9] = __DT*__P[19] + __P[9]; __RET_COV[10] = __DT*__P[22] + __P[10]; __RET_COV[11] = \
__DT*__P[23] + __P[11]; __RET_COV[12] = __DT*__P[24] + __P[12]; __RET_COV[13] = __DT*__P[16] + \
__DT*(__DT*__P[25] + __P[16]) + __P[13]; __RET_COV[14] = __DT*__P[20] + __P[14]; __RET_COV[15] = \
__DT*__P[23] + __P[15]; __RET_COV[16] = __DT*__P[25] + __P[16]; __RET_COV[17] = __DT*__P[26] + \
__P[17]; __RET_COV[18] = __P[18] + __W_U_SIGMA[0]*__W_U_SIGMA[0]; __RET_COV[19] = __P[19]; \
__RET_COV[20] = __P[20]; __RET_COV[21] = __P[21]; __RET_COV[22] = __P[22] + \
__W_U_SIGMA[1]*__W_U_SIGMA[1]; __RET_COV[23] = __P[23]; __RET_COV[24] = __P[24]; __RET_COV[25] = \
__P[25] + __W_U_SIGMA[2]*__W_U_SIGMA[2]; __RET_COV[26] = __P[26]; __RET_COV[27] = __P[27]; 

#define EKF_VELD_CALC_NIS(__P, __R, __SUBX, __X, __Z, __RET_NIS) \
__RET_NIS = __SUBX[0]*__X[5] + __Z*__X[5] + __Z; 

#define EKF_VELD_CALC_STATE(__P, __R, __SUBX, __X, __Z, __RET_STATE) \
__RET_STATE[0] = -__P[5]*__SUBX[0]*(__X[5] + __Z) + __X[0]; __RET_STATE[1] = \
-__P[11]*__SUBX[0]*(__X[5] + __Z) + __X[1]; __RET_STATE[2] = -__P[16]*__SUBX[0]*(__X[5] + __Z) + \
__X[2]; __RET_STATE[3] = -__P[20]*__SUBX[0]*(__X[5] + __Z) + __X[3]; __RET_STATE[4] = \
-__P[23]*__SUBX[0]*(__X[5] + __Z) + __X[4]; __RET_STATE[5] = -__SUBX[1]*(__X[5] + __Z) + __X[5]; \
__RET_STATE[6] = -__P[26]*__SUBX[0]*(__X[5] + __Z) + __X[6]; 

#define EKF_VELD_CALC_COV(__P, __R, __SUBX, __X, __Z, __RET_COV) \
__RET_COV[0] = __P[0] + __P[5]*__P[5]*__R*__SUBX[2] - __P[5]*__P[5]*__SUBX[0] - \
__P[5]*__SUBX[0]*__SUBX[3]; __RET_COV[1] = __P[11]*__P[5]*__R*__SUBX[2] - __P[11]*__P[5]*__SUBX[0] - \
__P[11]*__SUBX[0]*__SUBX[3] + __P[1]; __RET_COV[2] = -__P[16]*__P[5]*__SUBX[0] - \
__P[16]*__SUBX[0]*__SUBX[3] + __P[2] + __P[5]*__SUBX[4]; __RET_COV[3] = __P[20]*__P[5]*__R*__SUBX[2] \
- __P[20]*__P[5]*__SUBX[0] - __P[20]*__SUBX[0]*__SUBX[3] + __P[3]; __RET_COV[4] = \
__P[23]*__P[5]*__R*__SUBX[2] - __P[23]*__P[5]*__SUBX[0] - __P[23]*__SUBX[0]*__SUBX[3] + __P[4]; \
__RET_COV[5] = __P[25]*__P[5]*__R*__SUBX[2] + __SUBX[3]*(-__SUBX[1] + 1); __RET_COV[6] = \
__P[26]*__P[5]*__R*__SUBX[2] - __P[26]*__P[5]*__SUBX[0] - __P[26]*__SUBX[0]*__SUBX[3] + __P[6]; \
__RET_COV[7] = __P[11]*__P[11]*__R*__SUBX[2] - __P[11]*__P[11]*__SUBX[0] - \
__P[11]*__SUBX[0]*__SUBX[5] + __P[7]; __RET_COV[8] = -__P[11]*__P[16]*__SUBX[0] + __P[11]*__SUBX[4] - \
__P[16]*__SUBX[0]*__SUBX[5] + __P[8]; __RET_COV[9] = __P[11]*__P[20]*__R*__SUBX[2] - \
__P[11]*__P[20]*__SUBX[0] - __P[20]*__SUBX[0]*__SUBX[5] + __P[9]; __RET_COV[10] = __P[10] + \
__P[11]*__P[23]*__R*__SUBX[2] - __P[11]*__P[23]*__SUBX[0] - __P[23]*__SUBX[0]*__SUBX[5]; \
__RET_COV[11] = __P[11]*__P[25]*__R*__SUBX[2] + __SUBX[5]*(-__SUBX[1] + 1); __RET_COV[12] = \
__P[11]*__P[26]*__R*__SUBX[2] - __P[11]*__P[26]*__SUBX[0] + __P[12] - __P[26]*__SUBX[0]*__SUBX[5]; \
__RET_COV[13] = __P[13] + __P[16]*__P[16]*__R*__SUBX[2] - __P[16]*__P[16]*__SUBX[0] - \
__P[16]*__SUBX[0]*(-__P[16]*__SUBX[1] + __P[16]); __RET_COV[14] = __P[14] - __P[16]*__P[20]*__SUBX[0] \
- __P[20]*__SUBX[0]*(-__P[16]*__SUBX[1] + __P[16]) + __P[20]*__SUBX[4]; __RET_COV[15] = __P[15] - \
__P[16]*__P[23]*__SUBX[0] - __P[23]*__SUBX[0]*(-__P[16]*__SUBX[1] + __P[16]) + __P[23]*__SUBX[4]; \
__RET_COV[16] = __P[25]*__SUBX[4] + (-__SUBX[1] + 1)*(-__P[16]*__SUBX[1] + __P[16]); __RET_COV[17] = \
-__P[16]*__P[26]*__SUBX[0] + __P[17] - __P[26]*__SUBX[0]*(-__P[16]*__SUBX[1] + __P[16]) + \
__P[26]*__SUBX[4]; __RET_COV[18] = __P[18] + __P[20]*__P[20]*__R*__SUBX[2] - \
__P[20]*__P[20]*__SUBX[0] - __P[20]*__SUBX[0]*(-__P[20]*__SUBX[1] + __P[20]); __RET_COV[19] = __P[19] \
+ __P[20]*__P[23]*__R*__SUBX[2] - __P[20]*__P[23]*__SUBX[0] - __P[23]*__SUBX[0]*(-__P[20]*__SUBX[1] + \
__P[20]); __RET_COV[20] = __P[20]*__P[25]*__R*__SUBX[2] + (-__SUBX[1] + 1)*(-__P[20]*__SUBX[1] + \
__P[20]); __RET_COV[21] = __P[20]*__P[26]*__R*__SUBX[2] - __P[20]*__P[26]*__SUBX[0] + __P[21] - \
__P[26]*__SUBX[0]*(-__P[20]*__SUBX[1] + __P[20]); __RET_COV[22] = __P[22] + \
__P[23]*__P[23]*__R*__SUBX[2] - __P[23]*__P[23]*__SUBX[0] - __P[23]*__SUBX[0]*(-__P[23]*__SUBX[1] + \
__P[23]); __RET_COV[23] = __P[23]*__P[25]*__R*__SUBX[2] + (-__SUBX[1] + 1)*(-__P[23]*__SUBX[1] + \
__P[23]); __RET_COV[24] = __P[23]*__P[26]*__R*__SUBX[2] - __P[23]*__P[26]*__SUBX[0] + __P[24] - \
__P[26]*__SUBX[0]*(-__P[23]*__SUBX[1] + __P[23]); __RET_COV[25] = __P[25]*__P[25]*__R*__SUBX[2] + \
__P[25]*-__SUBX[1] + 1*-__SUBX[1] + 1; __RET_COV[26] = __P[25]*__P[26]*__R*__SUBX[2] - \
__P[26]*__SUBX[1]*(-__SUBX[1] + 1) + __P[26]*(-__SUBX[1] + 1); __RET_COV[27] = \
__P[26]*__P[26]*__R*__SUBX[2] - __P[26]*__P[26]*__SUBX[0] - __P[26]*__SUBX[0]*(-__P[26]*__SUBX[1] + \
__P[26]) + __P[27]; 

#define EKF_VELD_CALC_SUBX(__P, __R, __X, __Z, __RET_SUBX) \
__RET_SUBX[0] = 1.0f/(__P[25] + __R); __RET_SUBX[1] = __P[25]*__RET_SUBX[0]; __RET_SUBX[2] = \
__RET_SUBX[0]*__RET_SUBX[0]; __RET_SUBX[3] = -__P[5]*__RET_SUBX[1] + __P[5]; __RET_SUBX[4] = \
__P[16]*__R*__RET_SUBX[2]; __RET_SUBX[5] = -__P[11]*__RET_SUBX[1] + __P[11]; 

#define EKF_VELD_CALC_INNOV(__P, __R, __SUBX, __X, __Z, __RET_INNOV) \
__RET_INNOV = __X[5] + __Z; 

#define EKF_VELNE_CALC_NIS(__P, __R, __SUBX, __X, __Z, __RET_NIS) \
__RET_NIS = __SUBX[2]*(__SUBX[0]*__SUBX[2]*(__P[22] + __R) - __SUBX[1]*__SUBX[3]) + \
__SUBX[3]*(__SUBX[0]*__SUBX[3]*(__P[18] + __R) - __SUBX[1]*__SUBX[2]); 

#define EKF_VELNE_CALC_STATE(__P, __R, __SUBX, __X, __Z, __RET_STATE) \
__RET_STATE[0] = __SUBX[2]*__SUBX[5] + __SUBX[3]*__SUBX[4] + __X[0]; __RET_STATE[1] = \
__SUBX[2]*__SUBX[7] + __SUBX[3]*__SUBX[6] + __X[1]; __RET_STATE[2] = __SUBX[2]*__SUBX[9] + \
__SUBX[3]*__SUBX[8] + __X[2]; __RET_STATE[3] = __SUBX[10]*__SUBX[2] + __SUBX[11]*__SUBX[3] + __X[3]; \
__RET_STATE[4] = __SUBX[12]*__SUBX[3] + __SUBX[13]*__SUBX[2] + __X[4]; __RET_STATE[5] = \
__SUBX[14]*__SUBX[3] + __SUBX[15]*__SUBX[2] + __X[5]; __RET_STATE[6] = __SUBX[16]*__SUBX[3] + \
__SUBX[17]*__SUBX[2] + __X[6]; 

#define EKF_VELNE_CALC_COV(__P, __R, __SUBX, __X, __Z, __RET_COV) \
__RET_COV[0] = __P[0] + __P[3]*__SUBX[5] + __P[4]*__SUBX[4] + __R*__SUBX[4]*__SUBX[4] + \
__R*__SUBX[5]*__SUBX[5] + __SUBX[18]*__SUBX[4] + __SUBX[19]*__SUBX[5]; __RET_COV[1] = \
__P[10]*__SUBX[4] + __P[1] + __P[9]*__SUBX[5] + __R*__SUBX[4]*__SUBX[6] + __R*__SUBX[5]*__SUBX[7] + \
__SUBX[18]*__SUBX[6] + __SUBX[19]*__SUBX[7]; __RET_COV[2] = __P[14]*__SUBX[5] + __P[15]*__SUBX[4] + \
__P[2] + __R*__SUBX[4]*__SUBX[8] + __R*__SUBX[5]*__SUBX[9] + __SUBX[18]*__SUBX[8] + \
__SUBX[19]*__SUBX[9]; __RET_COV[3] = __R*__SUBX[10]*__SUBX[5] + __R*__SUBX[11]*__SUBX[4] + \
__SUBX[11]*__SUBX[18] + __SUBX[19]*__SUBX[20]; __RET_COV[4] = __R*__SUBX[12]*__SUBX[4] + \
__R*__SUBX[13]*__SUBX[5] + __SUBX[13]*__SUBX[19] + __SUBX[18]*__SUBX[21]; __RET_COV[5] = \
__P[20]*__SUBX[5] + __P[23]*__SUBX[4] + __P[5] + __R*__SUBX[14]*__SUBX[4] + __R*__SUBX[15]*__SUBX[5] \
+ __SUBX[14]*__SUBX[18] + __SUBX[15]*__SUBX[19]; __RET_COV[6] = __P[21]*__SUBX[5] + __P[24]*__SUBX[4] \
+ __P[6] + __R*__SUBX[16]*__SUBX[4] + __R*__SUBX[17]*__SUBX[5] + __SUBX[16]*__SUBX[18] + \
__SUBX[17]*__SUBX[19]; __RET_COV[7] = __P[10]*__SUBX[6] + __P[7] + __P[9]*__SUBX[7] + \
__R*__SUBX[6]*__SUBX[6] + __R*__SUBX[7]*__SUBX[7] + __SUBX[22]*__SUBX[6] + __SUBX[23]*__SUBX[7]; \
__RET_COV[8] = __P[14]*__SUBX[7] + __P[15]*__SUBX[6] + __P[8] + __R*__SUBX[6]*__SUBX[8] + \
__R*__SUBX[7]*__SUBX[9] + __SUBX[22]*__SUBX[8] + __SUBX[23]*__SUBX[9]; __RET_COV[9] = \
__R*__SUBX[10]*__SUBX[7] + __R*__SUBX[11]*__SUBX[6] + __SUBX[11]*__SUBX[22] + __SUBX[20]*__SUBX[23]; \
__RET_COV[10] = __R*__SUBX[12]*__SUBX[6] + __R*__SUBX[13]*__SUBX[7] + __SUBX[13]*__SUBX[23] + \
__SUBX[21]*__SUBX[22]; __RET_COV[11] = __P[11] + __P[20]*__SUBX[7] + __P[23]*__SUBX[6] + \
__R*__SUBX[14]*__SUBX[6] + __R*__SUBX[15]*__SUBX[7] + __SUBX[14]*__SUBX[22] + __SUBX[15]*__SUBX[23]; \
__RET_COV[12] = __P[12] + __P[21]*__SUBX[7] + __P[24]*__SUBX[6] + __R*__SUBX[16]*__SUBX[6] + \
__R*__SUBX[17]*__SUBX[7] + __SUBX[16]*__SUBX[22] + __SUBX[17]*__SUBX[23]; __RET_COV[13] = __P[13] + \
__P[14]*__SUBX[9] + __P[15]*__SUBX[8] + __R*__SUBX[8]*__SUBX[8] + __R*__SUBX[9]*__SUBX[9] + \
__SUBX[24]*__SUBX[8] + __SUBX[25]*__SUBX[9]; __RET_COV[14] = __R*__SUBX[10]*__SUBX[9] + \
__R*__SUBX[11]*__SUBX[8] + __SUBX[11]*__SUBX[24] + __SUBX[20]*__SUBX[25]; __RET_COV[15] = \
__R*__SUBX[12]*__SUBX[8] + __R*__SUBX[13]*__SUBX[9] + __SUBX[13]*__SUBX[25] + __SUBX[21]*__SUBX[24]; \
__RET_COV[16] = __P[16] + __P[20]*__SUBX[9] + __P[23]*__SUBX[8] + __R*__SUBX[14]*__SUBX[8] + \
__R*__SUBX[15]*__SUBX[9] + __SUBX[14]*__SUBX[24] + __SUBX[15]*__SUBX[25]; __RET_COV[17] = __P[17] + \
__P[21]*__SUBX[9] + __P[24]*__SUBX[8] + __R*__SUBX[16]*__SUBX[8] + __R*__SUBX[17]*__SUBX[9] + \
__SUBX[16]*__SUBX[24] + __SUBX[17]*__SUBX[25]; __RET_COV[18] = __R*__SUBX[10]*__SUBX[10] + \
__R*__SUBX[11]*__SUBX[11] + __SUBX[11]*(__P[19]*__SUBX[20] + __P[22]*__SUBX[11]) + \
__SUBX[20]*(__P[18]*__SUBX[20] + __P[19]*__SUBX[11]); __RET_COV[19] = __R*__SUBX[10]*__SUBX[13] + \
__R*__SUBX[11]*__SUBX[12] + __SUBX[13]*(__P[18]*__SUBX[20] + __P[19]*__SUBX[11]) + \
__SUBX[21]*(__P[19]*__SUBX[20] + __P[22]*__SUBX[11]); __RET_COV[20] = __P[20]*__SUBX[20] + \
__P[23]*__SUBX[11] + __R*__SUBX[10]*__SUBX[15] + __R*__SUBX[11]*__SUBX[14] + \
__SUBX[14]*(__P[19]*__SUBX[20] + __P[22]*__SUBX[11]) + __SUBX[15]*(__P[18]*__SUBX[20] + \
__P[19]*__SUBX[11]); __RET_COV[21] = __P[21]*__SUBX[20] + __P[24]*__SUBX[11] + \
__R*__SUBX[10]*__SUBX[17] + __R*__SUBX[11]*__SUBX[16] + __SUBX[16]*(__P[19]*__SUBX[20] + \
__P[22]*__SUBX[11]) + __SUBX[17]*(__P[18]*__SUBX[20] + __P[19]*__SUBX[11]); __RET_COV[22] = \
__R*__SUBX[12]*__SUBX[12] + __R*__SUBX[13]*__SUBX[13] + __SUBX[13]*(__P[18]*__SUBX[13] + \
__P[19]*__SUBX[21]) + __SUBX[21]*(__P[19]*__SUBX[13] + __P[22]*__SUBX[21]); __RET_COV[23] = \
__P[20]*__SUBX[13] + __P[23]*__SUBX[21] + __R*__SUBX[12]*__SUBX[14] + __R*__SUBX[13]*__SUBX[15] + \
__SUBX[14]*(__P[19]*__SUBX[13] + __P[22]*__SUBX[21]) + __SUBX[15]*(__P[18]*__SUBX[13] + \
__P[19]*__SUBX[21]); __RET_COV[24] = __P[21]*__SUBX[13] + __P[24]*__SUBX[21] + \
__R*__SUBX[12]*__SUBX[16] + __R*__SUBX[13]*__SUBX[17] + __SUBX[16]*(__P[19]*__SUBX[13] + \
__P[22]*__SUBX[21]) + __SUBX[17]*(__P[18]*__SUBX[13] + __P[19]*__SUBX[21]); __RET_COV[25] = \
__P[20]*__SUBX[15] + __P[23]*__SUBX[14] + __P[25] + __R*__SUBX[14]*__SUBX[14] + \
__R*__SUBX[15]*__SUBX[15] + __SUBX[14]*(__P[19]*__SUBX[15] + __P[22]*__SUBX[14] + __P[23]) + \
__SUBX[15]*(__P[18]*__SUBX[15] + __P[19]*__SUBX[14] + __P[20]); __RET_COV[26] = __P[21]*__SUBX[15] + \
__P[24]*__SUBX[14] + __P[26] + __R*__SUBX[14]*__SUBX[16] + __R*__SUBX[15]*__SUBX[17] + \
__SUBX[16]*(__P[19]*__SUBX[15] + __P[22]*__SUBX[14] + __P[23]) + __SUBX[17]*(__P[18]*__SUBX[15] + \
__P[19]*__SUBX[14] + __P[20]); __RET_COV[27] = __P[21]*__SUBX[17] + __P[24]*__SUBX[16] + __P[27] + \
__R*__SUBX[16]*__SUBX[16] + __R*__SUBX[17]*__SUBX[17] + __SUBX[16]*(__P[19]*__SUBX[17] + \
__P[22]*__SUBX[16] + __P[24]) + __SUBX[17]*(__P[18]*__SUBX[17] + __P[19]*__SUBX[16] + __P[21]); 

#define EKF_VELNE_CALC_SUBX(__P, __R, __X, __Z, __RET_SUBX) \
__RET_SUBX[0] = 1.0f/(-__P[19]*__P[19] + (__P[18] + __R)*(__P[22] + __R)); __RET_SUBX[1] = \
__P[19]*__RET_SUBX[0]; __RET_SUBX[2] = __X[3] + __Z[0]; __RET_SUBX[3] = __X[4] + __Z[1]; \
__RET_SUBX[4] = __P[3]*__RET_SUBX[1] - __P[4]*__RET_SUBX[0]*(__P[18] + __R); __RET_SUBX[5] = \
-__P[3]*__RET_SUBX[0]*(__P[22] + __R) + __P[4]*__RET_SUBX[1]; __RET_SUBX[6] = \
-__P[10]*__RET_SUBX[0]*(__P[18] + __R) + __P[9]*__RET_SUBX[1]; __RET_SUBX[7] = __P[10]*__RET_SUBX[1] \
- __P[9]*__RET_SUBX[0]*(__P[22] + __R); __RET_SUBX[8] = __P[14]*__RET_SUBX[1] - \
__P[15]*__RET_SUBX[0]*(__P[18] + __R); __RET_SUBX[9] = -__P[14]*__RET_SUBX[0]*(__P[22] + __R) + \
__P[15]*__RET_SUBX[1]; __RET_SUBX[10] = -__P[18]*__RET_SUBX[0]*(__P[22] + __R) + \
__P[19]*__P[19]*__RET_SUBX[0]; __RET_SUBX[11] = __P[18]*__RET_SUBX[1] - \
__P[19]*__RET_SUBX[0]*(__P[18] + __R); __RET_SUBX[12] = __P[19]*__P[19]*__RET_SUBX[0] - \
__P[22]*__RET_SUBX[0]*(__P[18] + __R); __RET_SUBX[13] = __P[22]*__RET_SUBX[1] - \
__RET_SUBX[1]*(__P[22] + __R); __RET_SUBX[14] = __P[20]*__RET_SUBX[1] - \
__P[23]*__RET_SUBX[0]*(__P[18] + __R); __RET_SUBX[15] = -__P[20]*__RET_SUBX[0]*(__P[22] + __R) + \
__P[23]*__RET_SUBX[1]; __RET_SUBX[16] = __P[21]*__RET_SUBX[1] - __P[24]*__RET_SUBX[0]*(__P[18] + \
__R); __RET_SUBX[17] = -__P[21]*__RET_SUBX[0]*(__P[22] + __R) + __P[24]*__RET_SUBX[1]; __RET_SUBX[18] \
= __P[19]*__RET_SUBX[5] + __P[22]*__RET_SUBX[4] + __P[4]; __RET_SUBX[19] = __P[18]*__RET_SUBX[5] + \
__P[19]*__RET_SUBX[4] + __P[3]; __RET_SUBX[20] = __RET_SUBX[10] + 1; __RET_SUBX[21] = __RET_SUBX[12] \
+ 1; __RET_SUBX[22] = __P[10] + __P[19]*__RET_SUBX[7] + __P[22]*__RET_SUBX[6]; __RET_SUBX[23] = \
__P[18]*__RET_SUBX[7] + __P[19]*__RET_SUBX[6] + __P[9]; __RET_SUBX[24] = __P[15] + \
__P[19]*__RET_SUBX[9] + __P[22]*__RET_SUBX[8]; __RET_SUBX[25] = __P[14] + __P[18]*__RET_SUBX[9] + \
__P[19]*__RET_SUBX[8]; 

#define EKF_VELNE_CALC_INNOV(__P, __R, __SUBX, __X, __Z, __RET_INNOV) \
__RET_INNOV[0] = __SUBX[2]; __RET_INNOV[1] = __SUBX[3]; 
