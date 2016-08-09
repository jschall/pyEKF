from sympy import *
import numpy as np
from math import *

EKF_NUM_STATES = 7
EKF_NUM_CONTROL_INPUTS = 3
EKF_MAX_NUM_SUBX = 72
EKF_STATE_IDX_PT_N = 0
EKF_STATE_IDX_PT_E = 1
EKF_STATE_IDX_PT_D = 2
EKF_STATE_IDX_VT_N = 3
EKF_STATE_IDX_VT_E = 4
EKF_STATE_IDX_VT_D = 5
EKF_STATE_IDX_ACC_D_B = 6
EKF_U_IDX_VDV_N = 0
EKF_U_IDX_VDV_E = 1
EKF_U_IDX_VDV_D = 2

def EKF_CAMERA_CALC_NIS(P, R, Tbn, subx, x, z):
    ret_NIS = np.zeros((1,1))
    ret_NIS = subx[20]*(subx[17]*subx[18]*subx[20] - subx[19]*subx[21]) + subx[21]*(subx[16]*subx[18]*subx[21] - subx[19]*subx[20])

    return ret_NIS

def EKF_CAMERA_CALC_STATE(P, R, Tbn, subx, x, z):
    ret_state = np.zeros((7,1))
    ret_state[0] = subx[20]*subx[22] + subx[21]*subx[23] + x[0]
    ret_state[1] = subx[20]*subx[24] + subx[21]*subx[25] + x[1]
    ret_state[2] = subx[20]*subx[26] + subx[21]*subx[27] + x[2]
    ret_state[3] = subx[20]*subx[28] + subx[21]*subx[29] + x[3]
    ret_state[4] = subx[20]*subx[30] + subx[21]*subx[31] + x[4]
    ret_state[5] = subx[20]*subx[32] + subx[21]*subx[33] + x[5]
    ret_state[6] = subx[20]*subx[34] + subx[21]*subx[35] + x[6]

    return ret_state

def EKF_CAMERA_CALC_COV(P, R, Tbn, subx, x, z):
    ret_cov = np.zeros((28,1))
    ret_cov[0] = R*subx[22]**2 + R*subx[23]**2 + subx[36]*subx[39] + subx[37]*subx[41] + subx[38]*subx[40]
    ret_cov[1] = R*subx[22]*subx[24] + R*subx[23]*subx[25] + subx[39]*subx[44] + subx[40]*subx[43] + subx[41]*subx[42]
    ret_cov[2] = R*subx[22]*subx[26] + R*subx[23]*subx[27] + subx[39]*subx[46] + subx[40]*subx[47] + subx[41]*subx[45]
    ret_cov[3] = P[14]*subx[38] + P[3]*subx[37] + P[9]*subx[36] + R*subx[22]*subx[28] + R*subx[23]*subx[29] + subx[39]*subx[49] + subx[40]*subx[50] + subx[41]*subx[48]
    ret_cov[4] = P[10]*subx[36] + P[15]*subx[38] + P[4]*subx[37] + R*subx[22]*subx[30] + R*subx[23]*subx[31] + subx[39]*subx[52] + subx[40]*subx[53] + subx[41]*subx[51]
    ret_cov[5] = P[11]*subx[36] + P[16]*subx[38] + P[5]*subx[37] + R*subx[22]*subx[32] + R*subx[23]*subx[33] + subx[39]*subx[55] + subx[40]*subx[56] + subx[41]*subx[54]
    ret_cov[6] = P[12]*subx[36] + P[17]*subx[38] + P[6]*subx[37] + R*subx[22]*subx[34] + R*subx[23]*subx[35] + subx[39]*subx[58] + subx[40]*subx[59] + subx[41]*subx[57]
    ret_cov[7] = R*subx[24]**2 + R*subx[25]**2 + subx[42]*subx[60] + subx[43]*subx[61] + subx[44]*subx[62]
    ret_cov[8] = R*subx[24]*subx[26] + R*subx[25]*subx[27] + subx[45]*subx[60] + subx[46]*subx[62] + subx[47]*subx[61]
    ret_cov[9] = P[14]*subx[43] + P[3]*subx[42] + P[9]*subx[44] + R*subx[24]*subx[28] + R*subx[25]*subx[29] + subx[48]*subx[60] + subx[49]*subx[62] + subx[50]*subx[61]
    ret_cov[10] = P[10]*subx[44] + P[15]*subx[43] + P[4]*subx[42] + R*subx[24]*subx[30] + R*subx[25]*subx[31] + subx[51]*subx[60] + subx[52]*subx[62] + subx[53]*subx[61]
    ret_cov[11] = P[11]*subx[44] + P[16]*subx[43] + P[5]*subx[42] + R*subx[24]*subx[32] + R*subx[25]*subx[33] + subx[54]*subx[60] + subx[55]*subx[62] + subx[56]*subx[61]
    ret_cov[12] = P[12]*subx[44] + P[17]*subx[43] + P[6]*subx[42] + R*subx[24]*subx[34] + R*subx[25]*subx[35] + subx[57]*subx[60] + subx[58]*subx[62] + subx[59]*subx[61]
    ret_cov[13] = R*subx[26]**2 + R*subx[27]**2 + subx[45]*subx[63] + subx[46]*subx[64] + subx[47]*subx[65]
    ret_cov[14] = P[14]*subx[47] + P[3]*subx[45] + P[9]*subx[46] + R*subx[26]*subx[28] + R*subx[27]*subx[29] + subx[48]*subx[63] + subx[49]*subx[64] + subx[50]*subx[65]
    ret_cov[15] = P[10]*subx[46] + P[15]*subx[47] + P[4]*subx[45] + R*subx[26]*subx[30] + R*subx[27]*subx[31] + subx[51]*subx[63] + subx[52]*subx[64] + subx[53]*subx[65]
    ret_cov[16] = P[11]*subx[46] + P[16]*subx[47] + P[5]*subx[45] + R*subx[26]*subx[32] + R*subx[27]*subx[33] + subx[54]*subx[63] + subx[55]*subx[64] + subx[56]*subx[65]
    ret_cov[17] = P[12]*subx[46] + P[17]*subx[47] + P[6]*subx[45] + R*subx[26]*subx[34] + R*subx[27]*subx[35] + subx[57]*subx[63] + subx[58]*subx[64] + subx[59]*subx[65]
    ret_cov[18] = P[14]*subx[50] + P[18] + P[3]*subx[48] + P[9]*subx[49] + R*subx[28]**2 + R*subx[29]**2 + subx[48]*subx[66] + subx[49]*subx[67] + subx[50]*subx[68]
    ret_cov[19] = P[10]*subx[49] + P[15]*subx[50] + P[19] + P[4]*subx[48] + R*subx[28]*subx[30] + R*subx[29]*subx[31] + subx[51]*subx[66] + subx[52]*subx[67] + subx[53]*subx[68]
    ret_cov[20] = P[11]*subx[49] + P[16]*subx[50] + P[20] + P[5]*subx[48] + R*subx[28]*subx[32] + R*subx[29]*subx[33] + subx[54]*subx[66] + subx[55]*subx[67] + subx[56]*subx[68]
    ret_cov[21] = P[12]*subx[49] + P[17]*subx[50] + P[21] + P[6]*subx[48] + R*subx[28]*subx[34] + R*subx[29]*subx[35] + subx[57]*subx[66] + subx[58]*subx[67] + subx[59]*subx[68]
    ret_cov[22] = P[10]*subx[52] + P[15]*subx[53] + P[22] + P[4]*subx[51] + R*subx[30]**2 + R*subx[31]**2 + subx[51]*subx[69] + subx[52]*subx[70] + subx[53]*subx[71]
    ret_cov[23] = P[11]*subx[52] + P[16]*subx[53] + P[23] + P[5]*subx[51] + R*subx[30]*subx[32] + R*subx[31]*subx[33] + subx[54]*subx[69] + subx[55]*subx[70] + subx[56]*subx[71]
    ret_cov[24] = P[12]*subx[52] + P[17]*subx[53] + P[24] + P[6]*subx[51] + R*subx[30]*subx[34] + R*subx[31]*subx[35] + subx[57]*subx[69] + subx[58]*subx[70] + subx[59]*subx[71]
    ret_cov[25] = P[11]*subx[55] + P[16]*subx[56] + P[25] + P[5]*subx[54] + R*subx[32]**2 + R*subx[33]**2 + subx[54]*(P[0]*subx[54] + P[1]*subx[55] + P[2]*subx[56] + P[5]) + subx[55]*(P[11] + P[1]*subx[54] + P[7]*subx[55] + P[8]*subx[56]) + subx[56]*(P[13]*subx[56] + P[16] + P[2]*subx[54] + P[8]*subx[55])
    ret_cov[26] = P[12]*subx[55] + P[17]*subx[56] + P[26] + P[6]*subx[54] + R*subx[32]*subx[34] + R*subx[33]*subx[35] + subx[57]*(P[0]*subx[54] + P[1]*subx[55] + P[2]*subx[56] + P[5]) + subx[58]*(P[11] + P[1]*subx[54] + P[7]*subx[55] + P[8]*subx[56]) + subx[59]*(P[13]*subx[56] + P[16] + P[2]*subx[54] + P[8]*subx[55])
    ret_cov[27] = P[12]*subx[58] + P[17]*subx[59] + P[27] + P[6]*subx[57] + R*subx[34]**2 + R*subx[35]**2 + subx[57]*(P[0]*subx[57] + P[1]*subx[58] + P[2]*subx[59] + P[6]) + subx[58]*(P[12] + P[1]*subx[57] + P[7]*subx[58] + P[8]*subx[59]) + subx[59]*(P[13]*subx[59] + P[17] + P[2]*subx[57] + P[8]*subx[58])

    return ret_cov

def EKF_CAMERA_CALC_SUBX(P, R, Tbn, x, z):
    ret_subx = np.zeros((72,1))
    ret_subx[0] = Tbn[0,2]*x[0] + Tbn[1,2]*x[1] + Tbn[2,2]*x[2]
    ret_subx[1] = Tbn[0,0]*x[0] + Tbn[1,0]*x[1] + Tbn[2,0]*x[2]
    ret_subx[2] = Tbn[0,1]*x[0] + Tbn[1,1]*x[1] + Tbn[2,1]*x[2]
    ret_subx[3] = ret_subx[0]**(-2)
    ret_subx[4] = Tbn[0,0]/ret_subx[0] - Tbn[0,2]*ret_subx[1]*ret_subx[3]
    ret_subx[5] = Tbn[1,0]/ret_subx[0] - Tbn[1,2]*ret_subx[1]*ret_subx[3]
    ret_subx[6] = Tbn[2,0]/ret_subx[0] - Tbn[2,2]*ret_subx[1]*ret_subx[3]
    ret_subx[7] = Tbn[0,1]/ret_subx[0] - Tbn[0,2]*ret_subx[2]*ret_subx[3]
    ret_subx[8] = Tbn[1,1]/ret_subx[0] - Tbn[1,2]*ret_subx[2]*ret_subx[3]
    ret_subx[9] = Tbn[2,1]/ret_subx[0] - Tbn[2,2]*ret_subx[2]*ret_subx[3]
    ret_subx[10] = P[0]*ret_subx[7] + P[1]*ret_subx[8] + P[2]*ret_subx[9]
    ret_subx[11] = P[1]*ret_subx[7] + P[7]*ret_subx[8] + P[8]*ret_subx[9]
    ret_subx[12] = P[13]*ret_subx[9] + P[2]*ret_subx[7] + P[8]*ret_subx[8]
    ret_subx[13] = P[0]*ret_subx[4] + P[1]*ret_subx[5] + P[2]*ret_subx[6]
    ret_subx[14] = P[1]*ret_subx[4] + P[7]*ret_subx[5] + P[8]*ret_subx[6]
    ret_subx[15] = P[13]*ret_subx[6] + P[2]*ret_subx[4] + P[8]*ret_subx[5]
    ret_subx[16] = R + ret_subx[13]*ret_subx[4] + ret_subx[14]*ret_subx[5] + ret_subx[15]*ret_subx[6]
    ret_subx[17] = R + ret_subx[10]*ret_subx[7] + ret_subx[11]*ret_subx[8] + ret_subx[12]*ret_subx[9]
    ret_subx[18] = 1/(ret_subx[16]*ret_subx[17] - (ret_subx[10]*ret_subx[4] + ret_subx[11]*ret_subx[5] + ret_subx[12]*ret_subx[6])**2)
    ret_subx[19] = ret_subx[18]*(ret_subx[10]*ret_subx[4] + ret_subx[11]*ret_subx[5] + ret_subx[12]*ret_subx[6])
    ret_subx[20] = z[0] - ret_subx[1]/ret_subx[0]
    ret_subx[21] = z[1] - ret_subx[2]/ret_subx[0]
    ret_subx[22] = -ret_subx[10]*ret_subx[19] + ret_subx[13]*ret_subx[17]*ret_subx[18]
    ret_subx[23] = ret_subx[10]*ret_subx[16]*ret_subx[18] - ret_subx[13]*ret_subx[19]
    ret_subx[24] = -ret_subx[11]*ret_subx[19] + ret_subx[14]*ret_subx[17]*ret_subx[18]
    ret_subx[25] = ret_subx[11]*ret_subx[16]*ret_subx[18] - ret_subx[14]*ret_subx[19]
    ret_subx[26] = -ret_subx[12]*ret_subx[19] + ret_subx[15]*ret_subx[17]*ret_subx[18]
    ret_subx[27] = ret_subx[12]*ret_subx[16]*ret_subx[18] - ret_subx[15]*ret_subx[19]
    ret_subx[28] = ret_subx[17]*ret_subx[18]*(P[14]*ret_subx[6] + P[3]*ret_subx[4] + P[9]*ret_subx[5]) - ret_subx[19]*(P[14]*ret_subx[9] + P[3]*ret_subx[7] + P[9]*ret_subx[8])
    ret_subx[29] = ret_subx[16]*ret_subx[18]*(P[14]*ret_subx[9] + P[3]*ret_subx[7] + P[9]*ret_subx[8]) - ret_subx[19]*(P[14]*ret_subx[6] + P[3]*ret_subx[4] + P[9]*ret_subx[5])
    ret_subx[30] = ret_subx[17]*ret_subx[18]*(P[10]*ret_subx[5] + P[15]*ret_subx[6] + P[4]*ret_subx[4]) - ret_subx[19]*(P[10]*ret_subx[8] + P[15]*ret_subx[9] + P[4]*ret_subx[7])
    ret_subx[31] = ret_subx[16]*ret_subx[18]*(P[10]*ret_subx[8] + P[15]*ret_subx[9] + P[4]*ret_subx[7]) - ret_subx[19]*(P[10]*ret_subx[5] + P[15]*ret_subx[6] + P[4]*ret_subx[4])
    ret_subx[32] = ret_subx[17]*ret_subx[18]*(P[11]*ret_subx[5] + P[16]*ret_subx[6] + P[5]*ret_subx[4]) - ret_subx[19]*(P[11]*ret_subx[8] + P[16]*ret_subx[9] + P[5]*ret_subx[7])
    ret_subx[33] = ret_subx[16]*ret_subx[18]*(P[11]*ret_subx[8] + P[16]*ret_subx[9] + P[5]*ret_subx[7]) - ret_subx[19]*(P[11]*ret_subx[5] + P[16]*ret_subx[6] + P[5]*ret_subx[4])
    ret_subx[34] = ret_subx[17]*ret_subx[18]*(P[12]*ret_subx[5] + P[17]*ret_subx[6] + P[6]*ret_subx[4]) - ret_subx[19]*(P[12]*ret_subx[8] + P[17]*ret_subx[9] + P[6]*ret_subx[7])
    ret_subx[35] = ret_subx[16]*ret_subx[18]*(P[12]*ret_subx[8] + P[17]*ret_subx[9] + P[6]*ret_subx[7]) - ret_subx[19]*(P[12]*ret_subx[5] + P[17]*ret_subx[6] + P[6]*ret_subx[4])
    ret_subx[36] = -ret_subx[22]*ret_subx[5] - ret_subx[23]*ret_subx[8]
    ret_subx[37] = -ret_subx[22]*ret_subx[4] - ret_subx[23]*ret_subx[7] + 1
    ret_subx[38] = -ret_subx[22]*ret_subx[6] - ret_subx[23]*ret_subx[9]
    ret_subx[39] = P[1]*ret_subx[37] + P[7]*ret_subx[36] + P[8]*ret_subx[38]
    ret_subx[40] = P[13]*ret_subx[38] + P[2]*ret_subx[37] + P[8]*ret_subx[36]
    ret_subx[41] = P[0]*ret_subx[37] + P[1]*ret_subx[36] + P[2]*ret_subx[38]
    ret_subx[42] = -ret_subx[24]*ret_subx[4] - ret_subx[25]*ret_subx[7]
    ret_subx[43] = -ret_subx[24]*ret_subx[6] - ret_subx[25]*ret_subx[9]
    ret_subx[44] = -ret_subx[24]*ret_subx[5] - ret_subx[25]*ret_subx[8] + 1
    ret_subx[45] = -ret_subx[26]*ret_subx[4] - ret_subx[27]*ret_subx[7]
    ret_subx[46] = -ret_subx[26]*ret_subx[5] - ret_subx[27]*ret_subx[8]
    ret_subx[47] = -ret_subx[26]*ret_subx[6] - ret_subx[27]*ret_subx[9] + 1
    ret_subx[48] = -ret_subx[28]*ret_subx[4] - ret_subx[29]*ret_subx[7]
    ret_subx[49] = -ret_subx[28]*ret_subx[5] - ret_subx[29]*ret_subx[8]
    ret_subx[50] = -ret_subx[28]*ret_subx[6] - ret_subx[29]*ret_subx[9]
    ret_subx[51] = -ret_subx[30]*ret_subx[4] - ret_subx[31]*ret_subx[7]
    ret_subx[52] = -ret_subx[30]*ret_subx[5] - ret_subx[31]*ret_subx[8]
    ret_subx[53] = -ret_subx[30]*ret_subx[6] - ret_subx[31]*ret_subx[9]
    ret_subx[54] = -ret_subx[32]*ret_subx[4] - ret_subx[33]*ret_subx[7]
    ret_subx[55] = -ret_subx[32]*ret_subx[5] - ret_subx[33]*ret_subx[8]
    ret_subx[56] = -ret_subx[32]*ret_subx[6] - ret_subx[33]*ret_subx[9]
    ret_subx[57] = -ret_subx[34]*ret_subx[4] - ret_subx[35]*ret_subx[7]
    ret_subx[58] = -ret_subx[34]*ret_subx[5] - ret_subx[35]*ret_subx[8]
    ret_subx[59] = -ret_subx[34]*ret_subx[6] - ret_subx[35]*ret_subx[9]
    ret_subx[60] = P[0]*ret_subx[42] + P[1]*ret_subx[44] + P[2]*ret_subx[43]
    ret_subx[61] = P[13]*ret_subx[43] + P[2]*ret_subx[42] + P[8]*ret_subx[44]
    ret_subx[62] = P[1]*ret_subx[42] + P[7]*ret_subx[44] + P[8]*ret_subx[43]
    ret_subx[63] = P[0]*ret_subx[45] + P[1]*ret_subx[46] + P[2]*ret_subx[47]
    ret_subx[64] = P[1]*ret_subx[45] + P[7]*ret_subx[46] + P[8]*ret_subx[47]
    ret_subx[65] = P[13]*ret_subx[47] + P[2]*ret_subx[45] + P[8]*ret_subx[46]
    ret_subx[66] = P[0]*ret_subx[48] + P[1]*ret_subx[49] + P[2]*ret_subx[50] + P[3]
    ret_subx[67] = P[1]*ret_subx[48] + P[7]*ret_subx[49] + P[8]*ret_subx[50] + P[9]
    ret_subx[68] = P[13]*ret_subx[50] + P[14] + P[2]*ret_subx[48] + P[8]*ret_subx[49]
    ret_subx[69] = P[0]*ret_subx[51] + P[1]*ret_subx[52] + P[2]*ret_subx[53] + P[4]
    ret_subx[70] = P[10] + P[1]*ret_subx[51] + P[7]*ret_subx[52] + P[8]*ret_subx[53]
    ret_subx[71] = P[13]*ret_subx[53] + P[15] + P[2]*ret_subx[51] + P[8]*ret_subx[52]

    return ret_subx

def EKF_CAMERA_CALC_INNOV(P, R, Tbn, subx, x, z):
    ret_innov = np.zeros((2,1))
    ret_innov[0] = subx[20]
    ret_innov[1] = subx[21]

    return ret_innov

def EKF_HEIGHT_CALC_NIS(P, R, subx, x, z):
    ret_NIS = np.zeros((1,1))
    ret_NIS = subx[0]*(-x[2] + z)**2

    return ret_NIS

def EKF_HEIGHT_CALC_STATE(P, R, subx, x, z):
    ret_state = np.zeros((7,1))
    ret_state[0] = P[2]*subx[0]*(-x[2] + z) + x[0]
    ret_state[1] = P[8]*subx[0]*(-x[2] + z) + x[1]
    ret_state[2] = subx[1]*(-x[2] + z) + x[2]
    ret_state[3] = P[14]*subx[0]*(-x[2] + z) + x[3]
    ret_state[4] = P[15]*subx[0]*(-x[2] + z) + x[4]
    ret_state[5] = P[16]*subx[0]*(-x[2] + z) + x[5]
    ret_state[6] = P[17]*subx[0]*(-x[2] + z) + x[6]

    return ret_state

def EKF_HEIGHT_CALC_COV(P, R, subx, x, z):
    ret_cov = np.zeros((28,1))
    ret_cov[0] = P[0] + P[2]**2*R*subx[2] - P[2]**2*subx[0] - P[2]*subx[0]*subx[3]
    ret_cov[1] = P[1] + P[2]*P[8]*R*subx[2] - P[2]*P[8]*subx[0] - P[8]*subx[0]*subx[3]
    ret_cov[2] = P[2]*subx[5] + subx[3]*subx[4]
    ret_cov[3] = P[14]*P[2]*R*subx[2] - P[14]*P[2]*subx[0] - P[14]*subx[0]*subx[3] + P[3]
    ret_cov[4] = P[15]*P[2]*R*subx[2] - P[15]*P[2]*subx[0] - P[15]*subx[0]*subx[3] + P[4]
    ret_cov[5] = P[16]*P[2]*R*subx[2] - P[16]*P[2]*subx[0] - P[16]*subx[0]*subx[3] + P[5]
    ret_cov[6] = P[17]*P[2]*R*subx[2] - P[17]*P[2]*subx[0] - P[17]*subx[0]*subx[3] + P[6]
    ret_cov[7] = P[7] + P[8]**2*R*subx[2] - P[8]**2*subx[0] - P[8]*subx[0]*subx[6]
    ret_cov[8] = P[8]*subx[5] + subx[4]*subx[6]
    ret_cov[9] = P[14]*P[8]*R*subx[2] - P[14]*P[8]*subx[0] - P[14]*subx[0]*subx[6] + P[9]
    ret_cov[10] = P[10] + P[15]*P[8]*R*subx[2] - P[15]*P[8]*subx[0] - P[15]*subx[0]*subx[6]
    ret_cov[11] = P[11] + P[16]*P[8]*R*subx[2] - P[16]*P[8]*subx[0] - P[16]*subx[0]*subx[6]
    ret_cov[12] = P[12] + P[17]*P[8]*R*subx[2] - P[17]*P[8]*subx[0] - P[17]*subx[0]*subx[6]
    ret_cov[13] = P[13]**2*R*subx[2] + P[13]*subx[4]**2
    ret_cov[14] = -P[14]*subx[1]*subx[4] + P[14]*subx[4] + P[14]*subx[5]
    ret_cov[15] = -P[15]*subx[1]*subx[4] + P[15]*subx[4] + P[15]*subx[5]
    ret_cov[16] = -P[16]*subx[1]*subx[4] + P[16]*subx[4] + P[16]*subx[5]
    ret_cov[17] = -P[17]*subx[1]*subx[4] + P[17]*subx[4] + P[17]*subx[5]
    ret_cov[18] = P[14]**2*R*subx[2] - P[14]**2*subx[0] - P[14]*subx[0]*(-P[14]*subx[1] + P[14]) + P[18]
    ret_cov[19] = P[14]*P[15]*R*subx[2] - P[14]*P[15]*subx[0] - P[15]*subx[0]*(-P[14]*subx[1] + P[14]) + P[19]
    ret_cov[20] = P[14]*P[16]*R*subx[2] - P[14]*P[16]*subx[0] - P[16]*subx[0]*(-P[14]*subx[1] + P[14]) + P[20]
    ret_cov[21] = P[14]*P[17]*R*subx[2] - P[14]*P[17]*subx[0] - P[17]*subx[0]*(-P[14]*subx[1] + P[14]) + P[21]
    ret_cov[22] = P[15]**2*R*subx[2] - P[15]**2*subx[0] - P[15]*subx[0]*(-P[15]*subx[1] + P[15]) + P[22]
    ret_cov[23] = P[15]*P[16]*R*subx[2] - P[15]*P[16]*subx[0] - P[16]*subx[0]*(-P[15]*subx[1] + P[15]) + P[23]
    ret_cov[24] = P[15]*P[17]*R*subx[2] - P[15]*P[17]*subx[0] - P[17]*subx[0]*(-P[15]*subx[1] + P[15]) + P[24]
    ret_cov[25] = P[16]**2*R*subx[2] - P[16]**2*subx[0] - P[16]*subx[0]*(-P[16]*subx[1] + P[16]) + P[25]
    ret_cov[26] = P[16]*P[17]*R*subx[2] - P[16]*P[17]*subx[0] - P[17]*subx[0]*(-P[16]*subx[1] + P[16]) + P[26]
    ret_cov[27] = P[17]**2*R*subx[2] - P[17]**2*subx[0] - P[17]*subx[0]*(-P[17]*subx[1] + P[17]) + P[27]

    return ret_cov

def EKF_HEIGHT_CALC_SUBX(P, R, x, z):
    ret_subx = np.zeros((7,1))
    ret_subx[0] = 1/(P[13] + R)
    ret_subx[1] = P[13]*ret_subx[0]
    ret_subx[2] = ret_subx[0]**2
    ret_subx[3] = -P[2]*ret_subx[1] + P[2]
    ret_subx[4] = -ret_subx[1] + 1
    ret_subx[5] = P[13]*R*ret_subx[2]
    ret_subx[6] = -P[8]*ret_subx[1] + P[8]

    return ret_subx

def EKF_HEIGHT_CALC_INNOV(P, R, subx, x, z):
    ret_innov = np.zeros((1,1))
    ret_innov = -x[2] + z

    return ret_innov

def EKF_INITIALIZATION_CALC_STATE(Tbn, cam_pos, cam_pos_R, hgt, hgt_R, subx, vdv_b_R, vel, vel_R):
    ret_state = np.zeros((7,1))
    ret_state[0] = hgt*subx[3]
    ret_state[1] = subx[4]*subx[5]
    ret_state[2] = hgt
    ret_state[3] = -vel[0]
    ret_state[4] = -vel[1]
    ret_state[5] = -vel[2]
    ret_state[6] = 0.0

    return ret_state

def EKF_INITIALIZATION_CALC_SUBX(Tbn, cam_pos, cam_pos_R, hgt, hgt_R, vdv_b_R, vel, vel_R):
    ret_subx = np.zeros((14,1))
    ret_subx[0] = Tbn[2,0]*cam_pos[0] + Tbn[2,1]*cam_pos[1] + 1.0*Tbn[2,2]
    ret_subx[1] = 1/ret_subx[0]
    ret_subx[2] = Tbn[0,0]*cam_pos[0] + Tbn[0,1]*cam_pos[1] + 1.0*Tbn[0,2]
    ret_subx[3] = ret_subx[1]*ret_subx[2]
    ret_subx[4] = Tbn[1,0]*cam_pos[0] + Tbn[1,1]*cam_pos[1] + 1.0*Tbn[1,2]
    ret_subx[5] = hgt*ret_subx[1]
    ret_subx[6] = ret_subx[0]**(-2)
    ret_subx[7] = hgt*ret_subx[2]*ret_subx[6]
    ret_subx[8] = Tbn[0,0]*ret_subx[5] - Tbn[2,0]*ret_subx[7]
    ret_subx[9] = Tbn[0,1]*ret_subx[5] - Tbn[2,1]*ret_subx[7]
    ret_subx[10] = hgt_R*ret_subx[6]
    ret_subx[11] = hgt*ret_subx[4]*ret_subx[6]
    ret_subx[12] = Tbn[1,0]*ret_subx[5] - Tbn[2,0]*ret_subx[11]
    ret_subx[13] = Tbn[1,1]*ret_subx[5] - Tbn[2,1]*ret_subx[11]

    return ret_subx

def EKF_INITIALIZATION_CALC_COV(Tbn, cam_pos, cam_pos_R, hgt, hgt_R, subx, vdv_b_R, vel, vel_R):
    ret_cov = np.zeros((28,1))
    ret_cov[0] = cam_pos_R*subx[8]**2 + cam_pos_R*subx[9]**2 + subx[10]*subx[2]**2
    ret_cov[1] = cam_pos_R*subx[12]*subx[8] + cam_pos_R*subx[13]*subx[9] + subx[10]*subx[2]*subx[4]
    ret_cov[2] = hgt_R*subx[3]
    ret_cov[3] = 0
    ret_cov[4] = 0
    ret_cov[5] = 0
    ret_cov[6] = 0
    ret_cov[7] = cam_pos_R*subx[12]**2 + cam_pos_R*subx[13]**2 + subx[10]*subx[4]**2
    ret_cov[8] = hgt_R*subx[1]*subx[4]
    ret_cov[9] = 0
    ret_cov[10] = 0
    ret_cov[11] = 0
    ret_cov[12] = 0
    ret_cov[13] = hgt_R
    ret_cov[14] = 0
    ret_cov[15] = 0
    ret_cov[16] = 0
    ret_cov[17] = 0
    ret_cov[18] = vel_R[0]
    ret_cov[19] = 0
    ret_cov[20] = 0
    ret_cov[21] = 0
    ret_cov[22] = vel_R[1]
    ret_cov[23] = 0
    ret_cov[24] = 0
    ret_cov[25] = vel_R[2]
    ret_cov[26] = 0
    ret_cov[27] = vdv_b_R

    return ret_cov

def EKF_PREDICTION_CALC_STATE(P, dt, subx, u, w_u_sigma, x):
    ret_state = np.zeros((7,1))
    ret_state[0] = dt*x[3] + x[0]
    ret_state[1] = dt*x[4] + x[1]
    ret_state[2] = dt*x[5] + x[2]
    ret_state[3] = -u[0] + x[3]
    ret_state[4] = -u[1] + x[4]
    ret_state[5] = dt*x[6] - u[2] + x[5]
    ret_state[6] = x[6]

    return ret_state

def EKF_PREDICTION_CALC_SUBX(P, dt, u, w_u_sigma, x):
    ret_subx = np.zeros((0,1))

    return ret_subx

def EKF_PREDICTION_CALC_COV(P, dt, subx, u, w_u_sigma, x):
    ret_cov = np.zeros((28,1))
    ret_cov[0] = P[0] + P[3]*dt + dt*(P[18]*dt + P[3])
    ret_cov[1] = P[1] + P[9]*dt + dt*(P[19]*dt + P[4])
    ret_cov[2] = P[14]*dt + P[2] + dt*(P[20]*dt + P[5])
    ret_cov[3] = P[18]*dt + P[3]
    ret_cov[4] = P[19]*dt + P[4]
    ret_cov[5] = P[20]*dt + P[5] + dt*(P[21]*dt + P[6])
    ret_cov[6] = P[21]*dt + P[6]
    ret_cov[7] = P[10]*dt + P[7] + dt*(P[10] + P[22]*dt)
    ret_cov[8] = P[15]*dt + P[8] + dt*(P[11] + P[23]*dt)
    ret_cov[9] = P[19]*dt + P[9]
    ret_cov[10] = P[10] + P[22]*dt
    ret_cov[11] = P[11] + P[23]*dt + dt*(P[12] + P[24]*dt)
    ret_cov[12] = P[12] + P[24]*dt
    ret_cov[13] = P[13] + P[16]*dt + dt*(P[16] + P[25]*dt)
    ret_cov[14] = P[14] + P[20]*dt
    ret_cov[15] = P[15] + P[23]*dt
    ret_cov[16] = P[16] + P[25]*dt + dt*(P[17] + P[26]*dt)
    ret_cov[17] = P[17] + P[26]*dt
    ret_cov[18] = P[18] + w_u_sigma[0]**2
    ret_cov[19] = P[19]
    ret_cov[20] = P[20] + P[21]*dt
    ret_cov[21] = P[21]
    ret_cov[22] = P[22] + w_u_sigma[1]**2
    ret_cov[23] = P[23] + P[24]*dt
    ret_cov[24] = P[24]
    ret_cov[25] = P[25] + P[26]*dt + dt*(P[26] + P[27]*dt) + w_u_sigma[2]**2
    ret_cov[26] = P[26] + P[27]*dt
    ret_cov[27] = P[27]

    return ret_cov

def EKF_VELD_CALC_NIS(P, R, subx, x, z):
    ret_NIS = np.zeros((1,1))
    ret_NIS = subx[0]*(x[5] + z)**2

    return ret_NIS

def EKF_VELD_CALC_STATE(P, R, subx, x, z):
    ret_state = np.zeros((7,1))
    ret_state[0] = -P[5]*subx[0]*(x[5] + z) + x[0]
    ret_state[1] = -P[11]*subx[0]*(x[5] + z) + x[1]
    ret_state[2] = -P[16]*subx[0]*(x[5] + z) + x[2]
    ret_state[3] = -P[20]*subx[0]*(x[5] + z) + x[3]
    ret_state[4] = -P[23]*subx[0]*(x[5] + z) + x[4]
    ret_state[5] = -subx[1]*(x[5] + z) + x[5]
    ret_state[6] = -P[26]*subx[0]*(x[5] + z) + x[6]

    return ret_state

def EKF_VELD_CALC_COV(P, R, subx, x, z):
    ret_cov = np.zeros((28,1))
    ret_cov[0] = P[0] + P[5]**2*R*subx[2] - P[5]**2*subx[0] - P[5]*subx[0]*subx[3]
    ret_cov[1] = P[11]*P[5]*R*subx[2] - P[11]*P[5]*subx[0] - P[11]*subx[0]*subx[3] + P[1]
    ret_cov[2] = -P[16]*P[5]*subx[0] - P[16]*subx[0]*subx[3] + P[2] + P[5]*subx[4]
    ret_cov[3] = P[20]*P[5]*R*subx[2] - P[20]*P[5]*subx[0] - P[20]*subx[0]*subx[3] + P[3]
    ret_cov[4] = P[23]*P[5]*R*subx[2] - P[23]*P[5]*subx[0] - P[23]*subx[0]*subx[3] + P[4]
    ret_cov[5] = P[25]*P[5]*R*subx[2] + subx[3]*(-subx[1] + 1)
    ret_cov[6] = P[26]*P[5]*R*subx[2] - P[26]*P[5]*subx[0] - P[26]*subx[0]*subx[3] + P[6]
    ret_cov[7] = P[11]**2*R*subx[2] - P[11]**2*subx[0] - P[11]*subx[0]*subx[5] + P[7]
    ret_cov[8] = -P[11]*P[16]*subx[0] + P[11]*subx[4] - P[16]*subx[0]*subx[5] + P[8]
    ret_cov[9] = P[11]*P[20]*R*subx[2] - P[11]*P[20]*subx[0] - P[20]*subx[0]*subx[5] + P[9]
    ret_cov[10] = P[10] + P[11]*P[23]*R*subx[2] - P[11]*P[23]*subx[0] - P[23]*subx[0]*subx[5]
    ret_cov[11] = P[11]*P[25]*R*subx[2] + subx[5]*(-subx[1] + 1)
    ret_cov[12] = P[11]*P[26]*R*subx[2] - P[11]*P[26]*subx[0] + P[12] - P[26]*subx[0]*subx[5]
    ret_cov[13] = P[13] + P[16]**2*R*subx[2] - P[16]**2*subx[0] - P[16]*subx[0]*(-P[16]*subx[1] + P[16])
    ret_cov[14] = P[14] - P[16]*P[20]*subx[0] - P[20]*subx[0]*(-P[16]*subx[1] + P[16]) + P[20]*subx[4]
    ret_cov[15] = P[15] - P[16]*P[23]*subx[0] - P[23]*subx[0]*(-P[16]*subx[1] + P[16]) + P[23]*subx[4]
    ret_cov[16] = P[25]*subx[4] + (-subx[1] + 1)*(-P[16]*subx[1] + P[16])
    ret_cov[17] = -P[16]*P[26]*subx[0] + P[17] - P[26]*subx[0]*(-P[16]*subx[1] + P[16]) + P[26]*subx[4]
    ret_cov[18] = P[18] + P[20]**2*R*subx[2] - P[20]**2*subx[0] - P[20]*subx[0]*(-P[20]*subx[1] + P[20])
    ret_cov[19] = P[19] + P[20]*P[23]*R*subx[2] - P[20]*P[23]*subx[0] - P[23]*subx[0]*(-P[20]*subx[1] + P[20])
    ret_cov[20] = P[20]*P[25]*R*subx[2] + (-subx[1] + 1)*(-P[20]*subx[1] + P[20])
    ret_cov[21] = P[20]*P[26]*R*subx[2] - P[20]*P[26]*subx[0] + P[21] - P[26]*subx[0]*(-P[20]*subx[1] + P[20])
    ret_cov[22] = P[22] + P[23]**2*R*subx[2] - P[23]**2*subx[0] - P[23]*subx[0]*(-P[23]*subx[1] + P[23])
    ret_cov[23] = P[23]*P[25]*R*subx[2] + (-subx[1] + 1)*(-P[23]*subx[1] + P[23])
    ret_cov[24] = P[23]*P[26]*R*subx[2] - P[23]*P[26]*subx[0] + P[24] - P[26]*subx[0]*(-P[23]*subx[1] + P[23])
    ret_cov[25] = P[25]**2*R*subx[2] + P[25]*(-subx[1] + 1)**2
    ret_cov[26] = P[25]*P[26]*R*subx[2] - P[26]*subx[1]*(-subx[1] + 1) + P[26]*(-subx[1] + 1)
    ret_cov[27] = P[26]**2*R*subx[2] - P[26]**2*subx[0] - P[26]*subx[0]*(-P[26]*subx[1] + P[26]) + P[27]

    return ret_cov

def EKF_VELD_CALC_SUBX(P, R, x, z):
    ret_subx = np.zeros((6,1))
    ret_subx[0] = 1/(P[25] + R)
    ret_subx[1] = P[25]*ret_subx[0]
    ret_subx[2] = ret_subx[0]**2
    ret_subx[3] = -P[5]*ret_subx[1] + P[5]
    ret_subx[4] = P[16]*R*ret_subx[2]
    ret_subx[5] = -P[11]*ret_subx[1] + P[11]

    return ret_subx

def EKF_VELD_CALC_INNOV(P, R, subx, x, z):
    ret_innov = np.zeros((1,1))
    ret_innov = x[5] + z

    return ret_innov

def EKF_VELNE_CALC_NIS(P, R, subx, x, z):
    ret_NIS = np.zeros((1,1))
    ret_NIS = subx[2]*(subx[0]*subx[2]*(P[22] + R) - subx[1]*subx[3]) + subx[3]*(subx[0]*subx[3]*(P[18] + R) - subx[1]*subx[2])

    return ret_NIS

def EKF_VELNE_CALC_STATE(P, R, subx, x, z):
    ret_state = np.zeros((7,1))
    ret_state[0] = subx[2]*subx[5] + subx[3]*subx[4] + x[0]
    ret_state[1] = subx[2]*subx[7] + subx[3]*subx[6] + x[1]
    ret_state[2] = subx[2]*subx[9] + subx[3]*subx[8] + x[2]
    ret_state[3] = subx[10]*subx[2] + subx[11]*subx[3] + x[3]
    ret_state[4] = subx[12]*subx[3] + subx[13]*subx[2] + x[4]
    ret_state[5] = subx[14]*subx[3] + subx[15]*subx[2] + x[5]
    ret_state[6] = subx[16]*subx[3] + subx[17]*subx[2] + x[6]

    return ret_state

def EKF_VELNE_CALC_COV(P, R, subx, x, z):
    ret_cov = np.zeros((28,1))
    ret_cov[0] = P[0] + P[3]*subx[5] + P[4]*subx[4] + R*subx[4]**2 + R*subx[5]**2 + subx[18]*subx[4] + subx[19]*subx[5]
    ret_cov[1] = P[10]*subx[4] + P[1] + P[9]*subx[5] + R*subx[4]*subx[6] + R*subx[5]*subx[7] + subx[18]*subx[6] + subx[19]*subx[7]
    ret_cov[2] = P[14]*subx[5] + P[15]*subx[4] + P[2] + R*subx[4]*subx[8] + R*subx[5]*subx[9] + subx[18]*subx[8] + subx[19]*subx[9]
    ret_cov[3] = R*subx[10]*subx[5] + R*subx[11]*subx[4] + subx[11]*subx[18] + subx[19]*subx[20]
    ret_cov[4] = R*subx[12]*subx[4] + R*subx[13]*subx[5] + subx[13]*subx[19] + subx[18]*subx[21]
    ret_cov[5] = P[20]*subx[5] + P[23]*subx[4] + P[5] + R*subx[14]*subx[4] + R*subx[15]*subx[5] + subx[14]*subx[18] + subx[15]*subx[19]
    ret_cov[6] = P[21]*subx[5] + P[24]*subx[4] + P[6] + R*subx[16]*subx[4] + R*subx[17]*subx[5] + subx[16]*subx[18] + subx[17]*subx[19]
    ret_cov[7] = P[10]*subx[6] + P[7] + P[9]*subx[7] + R*subx[6]**2 + R*subx[7]**2 + subx[22]*subx[6] + subx[23]*subx[7]
    ret_cov[8] = P[14]*subx[7] + P[15]*subx[6] + P[8] + R*subx[6]*subx[8] + R*subx[7]*subx[9] + subx[22]*subx[8] + subx[23]*subx[9]
    ret_cov[9] = R*subx[10]*subx[7] + R*subx[11]*subx[6] + subx[11]*subx[22] + subx[20]*subx[23]
    ret_cov[10] = R*subx[12]*subx[6] + R*subx[13]*subx[7] + subx[13]*subx[23] + subx[21]*subx[22]
    ret_cov[11] = P[11] + P[20]*subx[7] + P[23]*subx[6] + R*subx[14]*subx[6] + R*subx[15]*subx[7] + subx[14]*subx[22] + subx[15]*subx[23]
    ret_cov[12] = P[12] + P[21]*subx[7] + P[24]*subx[6] + R*subx[16]*subx[6] + R*subx[17]*subx[7] + subx[16]*subx[22] + subx[17]*subx[23]
    ret_cov[13] = P[13] + P[14]*subx[9] + P[15]*subx[8] + R*subx[8]**2 + R*subx[9]**2 + subx[24]*subx[8] + subx[25]*subx[9]
    ret_cov[14] = R*subx[10]*subx[9] + R*subx[11]*subx[8] + subx[11]*subx[24] + subx[20]*subx[25]
    ret_cov[15] = R*subx[12]*subx[8] + R*subx[13]*subx[9] + subx[13]*subx[25] + subx[21]*subx[24]
    ret_cov[16] = P[16] + P[20]*subx[9] + P[23]*subx[8] + R*subx[14]*subx[8] + R*subx[15]*subx[9] + subx[14]*subx[24] + subx[15]*subx[25]
    ret_cov[17] = P[17] + P[21]*subx[9] + P[24]*subx[8] + R*subx[16]*subx[8] + R*subx[17]*subx[9] + subx[16]*subx[24] + subx[17]*subx[25]
    ret_cov[18] = R*subx[10]**2 + R*subx[11]**2 + subx[11]*(P[19]*subx[20] + P[22]*subx[11]) + subx[20]*(P[18]*subx[20] + P[19]*subx[11])
    ret_cov[19] = R*subx[10]*subx[13] + R*subx[11]*subx[12] + subx[13]*(P[18]*subx[20] + P[19]*subx[11]) + subx[21]*(P[19]*subx[20] + P[22]*subx[11])
    ret_cov[20] = P[20]*subx[20] + P[23]*subx[11] + R*subx[10]*subx[15] + R*subx[11]*subx[14] + subx[14]*(P[19]*subx[20] + P[22]*subx[11]) + subx[15]*(P[18]*subx[20] + P[19]*subx[11])
    ret_cov[21] = P[21]*subx[20] + P[24]*subx[11] + R*subx[10]*subx[17] + R*subx[11]*subx[16] + subx[16]*(P[19]*subx[20] + P[22]*subx[11]) + subx[17]*(P[18]*subx[20] + P[19]*subx[11])
    ret_cov[22] = R*subx[12]**2 + R*subx[13]**2 + subx[13]*(P[18]*subx[13] + P[19]*subx[21]) + subx[21]*(P[19]*subx[13] + P[22]*subx[21])
    ret_cov[23] = P[20]*subx[13] + P[23]*subx[21] + R*subx[12]*subx[14] + R*subx[13]*subx[15] + subx[14]*(P[19]*subx[13] + P[22]*subx[21]) + subx[15]*(P[18]*subx[13] + P[19]*subx[21])
    ret_cov[24] = P[21]*subx[13] + P[24]*subx[21] + R*subx[12]*subx[16] + R*subx[13]*subx[17] + subx[16]*(P[19]*subx[13] + P[22]*subx[21]) + subx[17]*(P[18]*subx[13] + P[19]*subx[21])
    ret_cov[25] = P[20]*subx[15] + P[23]*subx[14] + P[25] + R*subx[14]**2 + R*subx[15]**2 + subx[14]*(P[19]*subx[15] + P[22]*subx[14] + P[23]) + subx[15]*(P[18]*subx[15] + P[19]*subx[14] + P[20])
    ret_cov[26] = P[21]*subx[15] + P[24]*subx[14] + P[26] + R*subx[14]*subx[16] + R*subx[15]*subx[17] + subx[16]*(P[19]*subx[15] + P[22]*subx[14] + P[23]) + subx[17]*(P[18]*subx[15] + P[19]*subx[14] + P[20])
    ret_cov[27] = P[21]*subx[17] + P[24]*subx[16] + P[27] + R*subx[16]**2 + R*subx[17]**2 + subx[16]*(P[19]*subx[17] + P[22]*subx[16] + P[24]) + subx[17]*(P[18]*subx[17] + P[19]*subx[16] + P[21])

    return ret_cov

def EKF_VELNE_CALC_SUBX(P, R, x, z):
    ret_subx = np.zeros((26,1))
    ret_subx[0] = 1/(-P[19]**2 + (P[18] + R)*(P[22] + R))
    ret_subx[1] = P[19]*ret_subx[0]
    ret_subx[2] = x[3] + z[0]
    ret_subx[3] = x[4] + z[1]
    ret_subx[4] = P[3]*ret_subx[1] - P[4]*ret_subx[0]*(P[18] + R)
    ret_subx[5] = -P[3]*ret_subx[0]*(P[22] + R) + P[4]*ret_subx[1]
    ret_subx[6] = -P[10]*ret_subx[0]*(P[18] + R) + P[9]*ret_subx[1]
    ret_subx[7] = P[10]*ret_subx[1] - P[9]*ret_subx[0]*(P[22] + R)
    ret_subx[8] = P[14]*ret_subx[1] - P[15]*ret_subx[0]*(P[18] + R)
    ret_subx[9] = -P[14]*ret_subx[0]*(P[22] + R) + P[15]*ret_subx[1]
    ret_subx[10] = -P[18]*ret_subx[0]*(P[22] + R) + P[19]**2*ret_subx[0]
    ret_subx[11] = P[18]*ret_subx[1] - P[19]*ret_subx[0]*(P[18] + R)
    ret_subx[12] = P[19]**2*ret_subx[0] - P[22]*ret_subx[0]*(P[18] + R)
    ret_subx[13] = P[22]*ret_subx[1] - ret_subx[1]*(P[22] + R)
    ret_subx[14] = P[20]*ret_subx[1] - P[23]*ret_subx[0]*(P[18] + R)
    ret_subx[15] = -P[20]*ret_subx[0]*(P[22] + R) + P[23]*ret_subx[1]
    ret_subx[16] = P[21]*ret_subx[1] - P[24]*ret_subx[0]*(P[18] + R)
    ret_subx[17] = -P[21]*ret_subx[0]*(P[22] + R) + P[24]*ret_subx[1]
    ret_subx[18] = P[19]*ret_subx[5] + P[22]*ret_subx[4] + P[4]
    ret_subx[19] = P[18]*ret_subx[5] + P[19]*ret_subx[4] + P[3]
    ret_subx[20] = ret_subx[10] + 1
    ret_subx[21] = ret_subx[12] + 1
    ret_subx[22] = P[10] + P[19]*ret_subx[7] + P[22]*ret_subx[6]
    ret_subx[23] = P[18]*ret_subx[7] + P[19]*ret_subx[6] + P[9]
    ret_subx[24] = P[15] + P[19]*ret_subx[9] + P[22]*ret_subx[8]
    ret_subx[25] = P[14] + P[18]*ret_subx[9] + P[19]*ret_subx[8]

    return ret_subx

def EKF_VELNE_CALC_INNOV(P, R, subx, x, z):
    ret_innov = np.zeros((2,1))
    ret_innov[0] = subx[2]
    ret_innov[1] = subx[3]

    return ret_innov

