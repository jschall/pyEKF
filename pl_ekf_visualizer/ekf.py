import numpy as np
from math import *

EKF_NUM_STATES = 6
EKF_NUM_CONTROL_INPUTS = 3
EKF_MAX_NUM_SUBX = 70
EKF_STATE_IDX_AZ = 0
EKF_STATE_IDX_EL = 1
EKF_STATE_IDX_RANGE_INV = 2
EKF_STATE_IDX_VT_N = 3
EKF_STATE_IDX_VT_E = 4
EKF_STATE_IDX_VT_D = 5
EKF_U_IDX_DVV_N = 0
EKF_U_IDX_DVV_E = 1
EKF_U_IDX_DVV_D = 2

def EKF_CAMERA_CALC_NIS(P, R, Tbn, subx, x, z):
    ret_NIS = 0
    ret_NIS = subx[26]*(subx[23]*subx[24]*subx[26] - subx[25]*subx[27]) + subx[27]*(subx[22]*subx[24]*subx[27] - subx[25]*subx[26])

    return ret_NIS

def EKF_CAMERA_CALC_STATE(P, R, Tbn, subx, x, z):
    ret_state = np.zeros(6)
    ret_state[0] = subx[26]*subx[28] + subx[27]*subx[29] + x[0]
    ret_state[1] = subx[26]*subx[30] + subx[27]*subx[31] + x[1]
    ret_state[2] = subx[26]*subx[32] + subx[27]*subx[33] + x[2]
    ret_state[3] = subx[26]*subx[34] + subx[27]*subx[35] + x[3]
    ret_state[4] = subx[26]*subx[36] + subx[27]*subx[37] + x[4]
    ret_state[5] = subx[26]*subx[38] + subx[27]*subx[39] + x[5]

    return ret_state

def EKF_CAMERA_CALC_COV(P, R, Tbn, subx, x, z):
    ret_cov = np.zeros(21)
    ret_cov[0] = R*subx[28]**2 + R*subx[29]**2 + subx[40]*subx[43] + subx[41]*subx[45] + subx[42]*subx[44]
    ret_cov[1] = R*subx[28]*subx[30] + R*subx[29]*subx[31] + subx[43]*subx[48] + subx[44]*subx[47] + subx[45]*subx[46]
    ret_cov[2] = R*subx[28]*subx[32] + R*subx[29]*subx[33] + subx[43]*subx[50] + subx[44]*subx[51] + subx[45]*subx[49]
    ret_cov[3] = P[12]*subx[42] + P[3]*subx[41] + P[8]*subx[40] + R*subx[28]*subx[34] + R*subx[29]*subx[35] + subx[43]*subx[53] + subx[44]*subx[54] + subx[45]*subx[52]
    ret_cov[4] = P[13]*subx[42] + P[4]*subx[41] + P[9]*subx[40] + R*subx[28]*subx[36] + R*subx[29]*subx[37] + subx[43]*subx[56] + subx[44]*subx[57] + subx[45]*subx[55]
    ret_cov[5] = P[10]*subx[40] + P[14]*subx[42] + P[5]*subx[41] + R*subx[28]*subx[38] + R*subx[29]*subx[39] + subx[43]*subx[59] + subx[44]*subx[60] + subx[45]*subx[58]
    ret_cov[6] = R*subx[30]**2 + R*subx[31]**2 + subx[46]*subx[61] + subx[47]*subx[62] + subx[48]*subx[63]
    ret_cov[7] = R*subx[30]*subx[32] + R*subx[31]*subx[33] + subx[49]*subx[61] + subx[50]*subx[63] + subx[51]*subx[62]
    ret_cov[8] = P[12]*subx[47] + P[3]*subx[46] + P[8]*subx[48] + R*subx[30]*subx[34] + R*subx[31]*subx[35] + subx[52]*subx[61] + subx[53]*subx[63] + subx[54]*subx[62]
    ret_cov[9] = P[13]*subx[47] + P[4]*subx[46] + P[9]*subx[48] + R*subx[30]*subx[36] + R*subx[31]*subx[37] + subx[55]*subx[61] + subx[56]*subx[63] + subx[57]*subx[62]
    ret_cov[10] = P[10]*subx[48] + P[14]*subx[47] + P[5]*subx[46] + R*subx[30]*subx[38] + R*subx[31]*subx[39] + subx[58]*subx[61] + subx[59]*subx[63] + subx[60]*subx[62]
    ret_cov[11] = R*subx[32]**2 + R*subx[33]**2 + subx[49]*subx[64] + subx[50]*subx[65] + subx[51]*subx[66]
    ret_cov[12] = P[12]*subx[51] + P[3]*subx[49] + P[8]*subx[50] + R*subx[32]*subx[34] + R*subx[33]*subx[35] + subx[52]*subx[64] + subx[53]*subx[65] + subx[54]*subx[66]
    ret_cov[13] = P[13]*subx[51] + P[4]*subx[49] + P[9]*subx[50] + R*subx[32]*subx[36] + R*subx[33]*subx[37] + subx[55]*subx[64] + subx[56]*subx[65] + subx[57]*subx[66]
    ret_cov[14] = P[10]*subx[50] + P[14]*subx[51] + P[5]*subx[49] + R*subx[32]*subx[38] + R*subx[33]*subx[39] + subx[58]*subx[64] + subx[59]*subx[65] + subx[60]*subx[66]
    ret_cov[15] = P[12]*subx[54] + P[15] + P[3]*subx[52] + P[8]*subx[53] + R*subx[34]**2 + R*subx[35]**2 + subx[52]*subx[67] + subx[53]*subx[68] + subx[54]*subx[69]
    ret_cov[16] = P[13]*subx[54] + P[16] + P[4]*subx[52] + P[9]*subx[53] + R*subx[34]*subx[36] + R*subx[35]*subx[37] + subx[55]*subx[67] + subx[56]*subx[68] + subx[57]*subx[69]
    ret_cov[17] = P[10]*subx[53] + P[14]*subx[54] + P[17] + P[5]*subx[52] + R*subx[34]*subx[38] + R*subx[35]*subx[39] + subx[58]*subx[67] + subx[59]*subx[68] + subx[60]*subx[69]
    ret_cov[18] = P[13]*subx[57] + P[18] + P[4]*subx[55] + P[9]*subx[56] + R*subx[36]**2 + R*subx[37]**2 + subx[55]*(P[0]*subx[55] + P[1]*subx[56] + P[2]*subx[57] + P[4]) + subx[56]*(P[1]*subx[55] + P[6]*subx[56] + P[7]*subx[57] + P[9]) + subx[57]*(P[11]*subx[57] + P[13] + P[2]*subx[55] + P[7]*subx[56])
    ret_cov[19] = P[10]*subx[56] + P[14]*subx[57] + P[19] + P[5]*subx[55] + R*subx[36]*subx[38] + R*subx[37]*subx[39] + subx[58]*(P[0]*subx[55] + P[1]*subx[56] + P[2]*subx[57] + P[4]) + subx[59]*(P[1]*subx[55] + P[6]*subx[56] + P[7]*subx[57] + P[9]) + subx[60]*(P[11]*subx[57] + P[13] + P[2]*subx[55] + P[7]*subx[56])
    ret_cov[20] = P[10]*subx[59] + P[14]*subx[60] + P[20] + P[5]*subx[58] + R*subx[38]**2 + R*subx[39]**2 + subx[58]*(P[0]*subx[58] + P[1]*subx[59] + P[2]*subx[60] + P[5]) + subx[59]*(P[10] + P[1]*subx[58] + P[6]*subx[59] + P[7]*subx[60]) + subx[60]*(P[11]*subx[60] + P[14] + P[2]*subx[58] + P[7]*subx[59])

    return ret_cov

def EKF_CAMERA_CALC_SUBX(P, R, Tbn, x, z):
    ret_subx = np.zeros(70)
    ret_subx[0] = 1/x[2]
    ret_subx[1] = sin(x[1])
    ret_subx[2] = cos(x[1])
    ret_subx[3] = 1.0*ret_subx[0]*ret_subx[2]*cos(x[0])
    ret_subx[4] = 1.0*ret_subx[0]*ret_subx[2]*sin(x[0])
    ret_subx[5] = Tbn[0,2]*ret_subx[3] + Tbn[1,2]*ret_subx[4] - 1.0*Tbn[2,2]*ret_subx[0]*ret_subx[1]
    ret_subx[6] = Tbn[0,0]*ret_subx[3] + Tbn[1,0]*ret_subx[4] - 1.0*Tbn[2,0]*ret_subx[0]*ret_subx[1]
    ret_subx[7] = Tbn[0,1]*ret_subx[3] + Tbn[1,1]*ret_subx[4] - 1.0*Tbn[2,1]*ret_subx[0]*ret_subx[1]
    ret_subx[8] = ret_subx[5]**(-2)
    ret_subx[9] = ret_subx[6]*ret_subx[8]*(Tbn[0,2]*ret_subx[4] - Tbn[1,2]*ret_subx[3]) + (-Tbn[0,0]*ret_subx[4] + Tbn[1,0]*ret_subx[3])/ret_subx[5]
    ret_subx[10] = ret_subx[6]*ret_subx[8]*(1.0*Tbn[0,2]*ret_subx[0]*ret_subx[1]*cos(x[0]) + 1.0*Tbn[1,2]*ret_subx[0]*ret_subx[1]*sin(x[0]) + 1.0*Tbn[2,2]*ret_subx[0]*ret_subx[2]) + (-1.0*Tbn[0,0]*ret_subx[0]*ret_subx[1]*cos(x[0]) - 1.0*Tbn[1,0]*ret_subx[0]*ret_subx[1]*sin(x[0]) - 1.0*Tbn[2,0]*ret_subx[0]*ret_subx[2])/ret_subx[5]
    ret_subx[11] = x[2]**(-2)
    ret_subx[12] = ret_subx[6]*ret_subx[8]*(1.0*Tbn[0,2]*ret_subx[11]*ret_subx[2]*cos(x[0]) + 1.0*Tbn[1,2]*ret_subx[11]*ret_subx[2]*sin(x[0]) - 1.0*Tbn[2,2]*ret_subx[11]*ret_subx[1]) + (-1.0*Tbn[0,0]*ret_subx[11]*ret_subx[2]*cos(x[0]) - 1.0*Tbn[1,0]*ret_subx[11]*ret_subx[2]*sin(x[0]) + 1.0*Tbn[2,0]*ret_subx[11]*ret_subx[1])/ret_subx[5]
    ret_subx[13] = ret_subx[7]*ret_subx[8]*(Tbn[0,2]*ret_subx[4] - Tbn[1,2]*ret_subx[3]) + (-Tbn[0,1]*ret_subx[4] + Tbn[1,1]*ret_subx[3])/ret_subx[5]
    ret_subx[14] = ret_subx[7]*ret_subx[8]*(1.0*Tbn[0,2]*ret_subx[0]*ret_subx[1]*cos(x[0]) + 1.0*Tbn[1,2]*ret_subx[0]*ret_subx[1]*sin(x[0]) + 1.0*Tbn[2,2]*ret_subx[0]*ret_subx[2]) + (-1.0*Tbn[0,1]*ret_subx[0]*ret_subx[1]*cos(x[0]) - 1.0*Tbn[1,1]*ret_subx[0]*ret_subx[1]*sin(x[0]) - 1.0*Tbn[2,1]*ret_subx[0]*ret_subx[2])/ret_subx[5]
    ret_subx[15] = ret_subx[7]*ret_subx[8]*(1.0*Tbn[0,2]*ret_subx[11]*ret_subx[2]*cos(x[0]) + 1.0*Tbn[1,2]*ret_subx[11]*ret_subx[2]*sin(x[0]) - 1.0*Tbn[2,2]*ret_subx[11]*ret_subx[1]) + (-1.0*Tbn[0,1]*ret_subx[11]*ret_subx[2]*cos(x[0]) - 1.0*Tbn[1,1]*ret_subx[11]*ret_subx[2]*sin(x[0]) + 1.0*Tbn[2,1]*ret_subx[11]*ret_subx[1])/ret_subx[5]
    ret_subx[16] = P[0]*ret_subx[13] + P[1]*ret_subx[14] + P[2]*ret_subx[15]
    ret_subx[17] = P[11]*ret_subx[15] + P[2]*ret_subx[13] + P[7]*ret_subx[14]
    ret_subx[18] = P[1]*ret_subx[13] + P[6]*ret_subx[14] + P[7]*ret_subx[15]
    ret_subx[19] = P[0]*ret_subx[9] + P[1]*ret_subx[10] + P[2]*ret_subx[12]
    ret_subx[20] = P[11]*ret_subx[12] + P[2]*ret_subx[9] + P[7]*ret_subx[10]
    ret_subx[21] = P[1]*ret_subx[9] + P[6]*ret_subx[10] + P[7]*ret_subx[12]
    ret_subx[22] = R + ret_subx[10]*ret_subx[21] + ret_subx[12]*ret_subx[20] + ret_subx[19]*ret_subx[9]
    ret_subx[23] = R + ret_subx[13]*ret_subx[16] + ret_subx[14]*ret_subx[18] + ret_subx[15]*ret_subx[17]
    ret_subx[24] = 1/(ret_subx[22]*ret_subx[23] - (ret_subx[10]*ret_subx[18] + ret_subx[12]*ret_subx[17] + ret_subx[16]*ret_subx[9])**2)
    ret_subx[25] = ret_subx[24]*(ret_subx[10]*ret_subx[18] + ret_subx[12]*ret_subx[17] + ret_subx[16]*ret_subx[9])
    ret_subx[26] = z[0] - ret_subx[6]/ret_subx[5]
    ret_subx[27] = z[1] - ret_subx[7]/ret_subx[5]
    ret_subx[28] = -ret_subx[16]*ret_subx[25] + ret_subx[19]*ret_subx[23]*ret_subx[24]
    ret_subx[29] = ret_subx[16]*ret_subx[22]*ret_subx[24] - ret_subx[19]*ret_subx[25]
    ret_subx[30] = -ret_subx[18]*ret_subx[25] + ret_subx[21]*ret_subx[23]*ret_subx[24]
    ret_subx[31] = ret_subx[18]*ret_subx[22]*ret_subx[24] - ret_subx[21]*ret_subx[25]
    ret_subx[32] = -ret_subx[17]*ret_subx[25] + ret_subx[20]*ret_subx[23]*ret_subx[24]
    ret_subx[33] = ret_subx[17]*ret_subx[22]*ret_subx[24] - ret_subx[20]*ret_subx[25]
    ret_subx[34] = ret_subx[23]*ret_subx[24]*(P[12]*ret_subx[12] + P[3]*ret_subx[9] + P[8]*ret_subx[10]) - ret_subx[25]*(P[12]*ret_subx[15] + P[3]*ret_subx[13] + P[8]*ret_subx[14])
    ret_subx[35] = ret_subx[22]*ret_subx[24]*(P[12]*ret_subx[15] + P[3]*ret_subx[13] + P[8]*ret_subx[14]) - ret_subx[25]*(P[12]*ret_subx[12] + P[3]*ret_subx[9] + P[8]*ret_subx[10])
    ret_subx[36] = ret_subx[23]*ret_subx[24]*(P[13]*ret_subx[12] + P[4]*ret_subx[9] + P[9]*ret_subx[10]) - ret_subx[25]*(P[13]*ret_subx[15] + P[4]*ret_subx[13] + P[9]*ret_subx[14])
    ret_subx[37] = ret_subx[22]*ret_subx[24]*(P[13]*ret_subx[15] + P[4]*ret_subx[13] + P[9]*ret_subx[14]) - ret_subx[25]*(P[13]*ret_subx[12] + P[4]*ret_subx[9] + P[9]*ret_subx[10])
    ret_subx[38] = ret_subx[23]*ret_subx[24]*(P[10]*ret_subx[10] + P[14]*ret_subx[12] + P[5]*ret_subx[9]) - ret_subx[25]*(P[10]*ret_subx[14] + P[14]*ret_subx[15] + P[5]*ret_subx[13])
    ret_subx[39] = ret_subx[22]*ret_subx[24]*(P[10]*ret_subx[14] + P[14]*ret_subx[15] + P[5]*ret_subx[13]) - ret_subx[25]*(P[10]*ret_subx[10] + P[14]*ret_subx[12] + P[5]*ret_subx[9])
    ret_subx[40] = -ret_subx[10]*ret_subx[28] - ret_subx[14]*ret_subx[29]
    ret_subx[41] = -ret_subx[13]*ret_subx[29] - ret_subx[28]*ret_subx[9] + 1
    ret_subx[42] = -ret_subx[12]*ret_subx[28] - ret_subx[15]*ret_subx[29]
    ret_subx[43] = P[1]*ret_subx[41] + P[6]*ret_subx[40] + P[7]*ret_subx[42]
    ret_subx[44] = P[11]*ret_subx[42] + P[2]*ret_subx[41] + P[7]*ret_subx[40]
    ret_subx[45] = P[0]*ret_subx[41] + P[1]*ret_subx[40] + P[2]*ret_subx[42]
    ret_subx[46] = -ret_subx[13]*ret_subx[31] - ret_subx[30]*ret_subx[9]
    ret_subx[47] = -ret_subx[12]*ret_subx[30] - ret_subx[15]*ret_subx[31]
    ret_subx[48] = -ret_subx[10]*ret_subx[30] - ret_subx[14]*ret_subx[31] + 1
    ret_subx[49] = -ret_subx[13]*ret_subx[33] - ret_subx[32]*ret_subx[9]
    ret_subx[50] = -ret_subx[10]*ret_subx[32] - ret_subx[14]*ret_subx[33]
    ret_subx[51] = -ret_subx[12]*ret_subx[32] - ret_subx[15]*ret_subx[33] + 1
    ret_subx[52] = -ret_subx[13]*ret_subx[35] - ret_subx[34]*ret_subx[9]
    ret_subx[53] = -ret_subx[10]*ret_subx[34] - ret_subx[14]*ret_subx[35]
    ret_subx[54] = -ret_subx[12]*ret_subx[34] - ret_subx[15]*ret_subx[35]
    ret_subx[55] = -ret_subx[13]*ret_subx[37] - ret_subx[36]*ret_subx[9]
    ret_subx[56] = -ret_subx[10]*ret_subx[36] - ret_subx[14]*ret_subx[37]
    ret_subx[57] = -ret_subx[12]*ret_subx[36] - ret_subx[15]*ret_subx[37]
    ret_subx[58] = -ret_subx[13]*ret_subx[39] - ret_subx[38]*ret_subx[9]
    ret_subx[59] = -ret_subx[10]*ret_subx[38] - ret_subx[14]*ret_subx[39]
    ret_subx[60] = -ret_subx[12]*ret_subx[38] - ret_subx[15]*ret_subx[39]
    ret_subx[61] = P[0]*ret_subx[46] + P[1]*ret_subx[48] + P[2]*ret_subx[47]
    ret_subx[62] = P[11]*ret_subx[47] + P[2]*ret_subx[46] + P[7]*ret_subx[48]
    ret_subx[63] = P[1]*ret_subx[46] + P[6]*ret_subx[48] + P[7]*ret_subx[47]
    ret_subx[64] = P[0]*ret_subx[49] + P[1]*ret_subx[50] + P[2]*ret_subx[51]
    ret_subx[65] = P[1]*ret_subx[49] + P[6]*ret_subx[50] + P[7]*ret_subx[51]
    ret_subx[66] = P[11]*ret_subx[51] + P[2]*ret_subx[49] + P[7]*ret_subx[50]
    ret_subx[67] = P[0]*ret_subx[52] + P[1]*ret_subx[53] + P[2]*ret_subx[54] + P[3]
    ret_subx[68] = P[1]*ret_subx[52] + P[6]*ret_subx[53] + P[7]*ret_subx[54] + P[8]
    ret_subx[69] = P[11]*ret_subx[54] + P[12] + P[2]*ret_subx[52] + P[7]*ret_subx[53]

    return ret_subx

def EKF_CAMERA_CALC_INNOV(P, R, Tbn, subx, x, z):
    ret_innov = np.zeros(2)
    ret_innov[0] = subx[26]
    ret_innov[1] = subx[27]

    return ret_innov

def EKF_HEIGHT_CALC_NIS(P, R, subx, x, z):
    ret_NIS = 0
    ret_NIS = subx[4]*subx[5]**2

    return ret_NIS

def EKF_HEIGHT_CALC_STATE(P, R, subx, x, z):
    ret_state = np.zeros(6)
    ret_state[0] = subx[4]*subx[5]*subx[6] + x[0]
    ret_state[1] = subx[3]*subx[4]*subx[5] + x[1]
    ret_state[2] = subx[2]*subx[4]*subx[5] + x[2]
    ret_state[3] = subx[4]*subx[5]*subx[7] + x[3]
    ret_state[4] = subx[4]*subx[5]*subx[8] + x[4]
    ret_state[5] = subx[4]*subx[5]*subx[9] + x[5]

    return ret_state

def EKF_HEIGHT_CALC_COV(P, R, subx, x, z):
    ret_cov = np.zeros(21)
    ret_cov[0] = P[0] + P[1]*subx[0]*subx[4]*subx[6] - P[2]*subx[1]*subx[4]*subx[6] + R*subx[10]*subx[6]**2 + subx[0]*subx[12]*subx[4]*subx[6] - subx[13]*subx[1]*subx[4]*subx[6]
    ret_cov[1] = R*subx[10]*subx[3]*subx[6] + subx[12]*subx[14] - subx[13]*subx[1]*subx[3]*subx[4]
    ret_cov[2] = R*subx[10]*subx[2]*subx[6] + subx[12]*subx[16] + subx[13]*subx[15]
    ret_cov[3] = -P[12]*subx[1]*subx[4]*subx[6] + P[3] + P[8]*subx[0]*subx[4]*subx[6] + R*subx[10]*subx[6]*subx[7] + subx[12]*subx[17] - subx[13]*subx[18]
    ret_cov[4] = -P[13]*subx[1]*subx[4]*subx[6] + P[4] + P[9]*subx[0]*subx[4]*subx[6] + R*subx[10]*subx[6]*subx[8] + subx[12]*subx[19] - subx[13]*subx[20]
    ret_cov[5] = P[10]*subx[0]*subx[4]*subx[6] - P[14]*subx[1]*subx[4]*subx[6] + P[5] + R*subx[10]*subx[6]*subx[9] + subx[12]*subx[21] - subx[13]*subx[22]
    ret_cov[6] = R*subx[10]*subx[3]**2 + subx[14]*subx[23] - subx[1]*subx[24]*subx[3]*subx[4]
    ret_cov[7] = R*subx[10]*subx[2]*subx[3] + subx[0]*subx[23]*subx[2]*subx[4] + subx[15]*subx[24]
    ret_cov[8] = -P[12]*subx[1]*subx[3]*subx[4] + P[8]*subx[14] + R*subx[10]*subx[3]*subx[7] + subx[17]*subx[23] - subx[18]*subx[24]
    ret_cov[9] = -P[13]*subx[1]*subx[3]*subx[4] + P[9]*subx[14] + R*subx[10]*subx[3]*subx[8] + subx[0]*subx[23]*subx[4]*subx[8] - subx[20]*subx[24]
    ret_cov[10] = P[10]*subx[14] - P[14]*subx[1]*subx[3]*subx[4] + R*subx[10]*subx[3]*subx[9] + subx[21]*subx[23] - subx[22]*subx[24]
    ret_cov[11] = R*subx[10]*subx[2]**2 + subx[15]*(P[11]*subx[15] + P[7]*subx[16]) + subx[16]*(P[7]*subx[15] + subx[11]*subx[2])
    ret_cov[12] = P[12]*subx[15] + P[8]*subx[0]*subx[2]*subx[4] + R*subx[10]*subx[2]*subx[7] + subx[17]*(P[7]*subx[15] + subx[11]*subx[2]) - subx[18]*(P[11]*subx[15] + P[7]*subx[16])
    ret_cov[13] = P[13]*subx[15] + P[9]*subx[0]*subx[2]*subx[4] + R*subx[10]*subx[2]*subx[8] + subx[19]*(P[7]*subx[15] + subx[11]*subx[2]) - subx[20]*(P[11]*subx[15] + P[7]*subx[16])
    ret_cov[14] = P[10]*subx[0]*subx[2]*subx[4] + P[14]*subx[15] + R*subx[10]*subx[2]*subx[9] + subx[21]*(P[7]*subx[15] + subx[11]*subx[2]) - subx[22]*(P[11]*subx[15] + P[7]*subx[16])
    ret_cov[15] = -P[12]*subx[1]*subx[4]*subx[7] + P[15] + P[8]*subx[0]*subx[4]*subx[7] + R*subx[10]*subx[7]**2 + subx[17]*(P[6]*subx[17] - P[7]*subx[18] + P[8]) - subx[18]*(-P[11]*subx[18] + P[12] + P[7]*subx[17])
    ret_cov[16] = -P[13]*subx[1]*subx[4]*subx[7] + P[16] + P[9]*subx[0]*subx[4]*subx[7] + R*subx[10]*subx[7]*subx[8] + subx[19]*(P[6]*subx[17] - P[7]*subx[18] + P[8]) - subx[20]*(-P[11]*subx[18] + P[12] + P[7]*subx[17])
    ret_cov[17] = P[10]*subx[0]*subx[4]*subx[7] - P[14]*subx[1]*subx[4]*subx[7] + P[17] + R*subx[10]*subx[7]*subx[9] + subx[21]*(P[6]*subx[17] - P[7]*subx[18] + P[8]) - subx[22]*(-P[11]*subx[18] + P[12] + P[7]*subx[17])
    ret_cov[18] = -P[13]*subx[1]*subx[4]*subx[8] + P[18] + P[9]*subx[0]*subx[4]*subx[8] + R*subx[10]*subx[8]**2 + subx[19]*(-P[7]*subx[20] + P[9] + subx[11]*subx[8]) - subx[20]*(-P[11]*subx[20] + P[13] + P[7]*subx[19])
    ret_cov[19] = P[10]*subx[0]*subx[4]*subx[8] - P[14]*subx[1]*subx[4]*subx[8] + P[19] + R*subx[10]*subx[8]*subx[9] + subx[21]*(-P[7]*subx[20] + P[9] + subx[11]*subx[8]) - subx[22]*(-P[11]*subx[20] + P[13] + P[7]*subx[19])
    ret_cov[20] = P[10]*subx[0]*subx[4]*subx[9] - P[14]*subx[1]*subx[4]*subx[9] + P[20] + R*subx[10]*subx[9]**2 + subx[21]*(P[10] + P[6]*subx[21] - P[7]*subx[22]) - subx[22]*(-P[11]*subx[22] + P[14] + P[7]*subx[21])

    return ret_cov

def EKF_HEIGHT_CALC_SUBX(P, R, x, z):
    ret_subx = np.zeros(25)
    ret_subx[0] = 1.0*cos(x[1])/x[2]
    ret_subx[1] = 1.0*sin(x[1])/x[2]**2
    ret_subx[2] = P[11]*ret_subx[1] - P[7]*ret_subx[0]
    ret_subx[3] = -P[6]*ret_subx[0] + P[7]*ret_subx[1]
    ret_subx[4] = 1/(R - ret_subx[0]*ret_subx[3] + ret_subx[1]*ret_subx[2])
    ret_subx[5] = z + 1.0*sin(x[1])/x[2]
    ret_subx[6] = -P[1]*ret_subx[0] + P[2]*ret_subx[1]
    ret_subx[7] = P[12]*ret_subx[1] - P[8]*ret_subx[0]
    ret_subx[8] = P[13]*ret_subx[1] - P[9]*ret_subx[0]
    ret_subx[9] = -P[10]*ret_subx[0] + P[14]*ret_subx[1]
    ret_subx[10] = ret_subx[4]**2
    ret_subx[11] = P[6]*ret_subx[0]*ret_subx[4]
    ret_subx[12] = P[1] - P[7]*ret_subx[1]*ret_subx[4]*ret_subx[6] + ret_subx[11]*ret_subx[6]
    ret_subx[13] = -P[11]*ret_subx[1]*ret_subx[4]*ret_subx[6] + P[2] + P[7]*ret_subx[0]*ret_subx[4]*ret_subx[6]
    ret_subx[14] = ret_subx[0]*ret_subx[3]*ret_subx[4] + 1
    ret_subx[15] = -ret_subx[1]*ret_subx[2]*ret_subx[4] + 1
    ret_subx[16] = ret_subx[0]*ret_subx[2]*ret_subx[4]
    ret_subx[17] = ret_subx[0]*ret_subx[4]*ret_subx[7]
    ret_subx[18] = ret_subx[1]*ret_subx[4]*ret_subx[7]
    ret_subx[19] = ret_subx[0]*ret_subx[4]*ret_subx[8]
    ret_subx[20] = ret_subx[1]*ret_subx[4]*ret_subx[8]
    ret_subx[21] = ret_subx[0]*ret_subx[4]*ret_subx[9]
    ret_subx[22] = ret_subx[1]*ret_subx[4]*ret_subx[9]
    ret_subx[23] = P[6]*ret_subx[14] - P[7]*ret_subx[1]*ret_subx[3]*ret_subx[4]
    ret_subx[24] = -P[11]*ret_subx[1]*ret_subx[3]*ret_subx[4] + P[7]*ret_subx[14]

    return ret_subx

def EKF_HEIGHT_CALC_INNOV(P, R, subx, x, z):
    ret_innov = 0
    ret_innov = subx[5]

    return ret_innov

def EKF_INITIALIZATION_CALC_STATE(Tbn, cam_pos, cam_pos_R, hgt, hgt_R, subx, vel, vel_R):
    ret_state = np.zeros(6)
    ret_state[0] = atan2(subx[2]*(Tbn[1,2] + subx[0]), subx[2]*(Tbn[0,2] + subx[3]))
    ret_state[1] = atan2(-hgt, sqrt(subx[8]))
    ret_state[2] = subx[9]
    ret_state[3] = -vel[0]
    ret_state[4] = -vel[1]
    ret_state[5] = -vel[2]

    return ret_state

def EKF_INITIALIZATION_CALC_SUBX(Tbn, cam_pos, cam_pos_R, hgt, hgt_R, vel, vel_R):
    ret_subx = np.zeros(20)
    ret_subx[0] = Tbn[1,0]*cam_pos[0] + Tbn[1,1]*cam_pos[1]
    ret_subx[1] = Tbn[2,0]*cam_pos[0] + Tbn[2,1]*cam_pos[1]
    ret_subx[2] = 1.0*hgt/(Tbn[2,2] + ret_subx[1])
    ret_subx[3] = Tbn[0,0]*cam_pos[0] + Tbn[0,1]*cam_pos[1]
    ret_subx[4] = (1.0*Tbn[0,2] + ret_subx[3])**2 + (1.0*Tbn[1,2] + ret_subx[0])**2
    ret_subx[5] = hgt**2
    ret_subx[6] = 1.0*Tbn[2,2] + ret_subx[1]
    ret_subx[7] = ret_subx[6]**2
    ret_subx[8] = ret_subx[4]*ret_subx[5]/ret_subx[7]
    ret_subx[9] = 1/sqrt(ret_subx[5]*(ret_subx[4] + ret_subx[7])/ret_subx[7])
    ret_subx[10] = (Tbn[2,2] + ret_subx[1])**(-2)
    ret_subx[11] = 1/(1.0*ret_subx[10]*ret_subx[5]*(Tbn[0,2] + ret_subx[3])**2 + 1.0*ret_subx[10]*ret_subx[5]*(Tbn[1,2] + ret_subx[0])**2)
    ret_subx[12] = 1.0*hgt*ret_subx[11]*(Tbn[0,2] + ret_subx[3])*(Tbn[1,0]*ret_subx[2] - 1.0*Tbn[2,0]*hgt*ret_subx[10]*(Tbn[1,2] + ret_subx[0]))/(Tbn[2,2] + ret_subx[1]) - 1.0*hgt*ret_subx[11]*(Tbn[1,2] + ret_subx[0])*(Tbn[0,0]*ret_subx[2] - 1.0*Tbn[2,0]*hgt*ret_subx[10]*(Tbn[0,2] + ret_subx[3]))/(Tbn[2,2] + ret_subx[1])
    ret_subx[13] = 1.0*hgt*ret_subx[11]*(Tbn[0,2] + ret_subx[3])*(Tbn[1,1]*ret_subx[2] - 1.0*Tbn[2,1]*hgt*ret_subx[10]*(Tbn[1,2] + ret_subx[0]))/(Tbn[2,2] + ret_subx[1]) - 1.0*hgt*ret_subx[11]*(Tbn[1,2] + ret_subx[0])*(Tbn[0,1]*ret_subx[2] - 1.0*Tbn[2,1]*hgt*ret_subx[10]*(Tbn[0,2] + ret_subx[3]))/(Tbn[2,2] + ret_subx[1])
    ret_subx[14] = 2*Tbn[0,0]*cam_pos[0] + 2*Tbn[0,1]*cam_pos[1] + 2.0*Tbn[0,2]
    ret_subx[15] = 2*Tbn[1,0]*cam_pos[0] + 2*Tbn[1,1]*cam_pos[1] + 2.0*Tbn[1,2]
    ret_subx[16] = ret_subx[5]/(2*ret_subx[7])
    ret_subx[17] = ret_subx[4]*ret_subx[5]/ret_subx[6]**3
    ret_subx[18] = Tbn[2,0]*ret_subx[5]*(ret_subx[4] + ret_subx[7])/ret_subx[6]**3 - ret_subx[16]*(Tbn[0,0]*ret_subx[14] + Tbn[1,0]*ret_subx[15] + Tbn[2,0]*(2*Tbn[2,0]*cam_pos[0] + 2*Tbn[2,1]*cam_pos[1] + 2.0*Tbn[2,2]))
    ret_subx[19] = Tbn[2,1]*ret_subx[5]*(ret_subx[4] + ret_subx[7])/ret_subx[6]**3 - ret_subx[16]*(Tbn[0,1]*ret_subx[14] + Tbn[1,1]*ret_subx[15] + Tbn[2,1]*(2*Tbn[2,0]*cam_pos[0] + 2*Tbn[2,1]*cam_pos[1] + 2.0*Tbn[2,2]))

    return ret_subx

def EKF_INITIALIZATION_CALC_COV(Tbn, cam_pos, cam_pos_R, hgt, hgt_R, subx, vel, vel_R):
    ret_cov = np.zeros(21)
    ret_cov[0] = cam_pos_R*subx[12]**2 + cam_pos_R*subx[13]**2
    ret_cov[1] = cam_pos_R*subx[12]*subx[7]*sqrt(subx[8])*(-Tbn[2,0]*subx[17] + subx[16]*(Tbn[0,0]*subx[14] + Tbn[1,0]*subx[15]))/(hgt*subx[4]*(subx[5] + subx[8])) + cam_pos_R*subx[13]*subx[7]*sqrt(subx[8])*(-Tbn[2,1]*subx[17] + subx[16]*(Tbn[0,1]*subx[14] + Tbn[1,1]*subx[15]))/(hgt*subx[4]*(subx[5] + subx[8]))
    ret_cov[2] = cam_pos_R*subx[12]*subx[18]*subx[7]*subx[9]/(subx[5]*(subx[4] + subx[7])) + cam_pos_R*subx[13]*subx[19]*subx[7]*subx[9]/(subx[5]*(subx[4] + subx[7]))
    ret_cov[3] = 0
    ret_cov[4] = 0
    ret_cov[5] = 0
    ret_cov[6] = cam_pos_R*subx[6]**4*subx[8]*(-Tbn[2,0]*subx[17] + subx[16]*(Tbn[0,0]*subx[14] + Tbn[1,0]*subx[15]))**2/(subx[4]**2*subx[5]*(subx[5] + subx[8])**2) + cam_pos_R*subx[6]**4*subx[8]*(-Tbn[2,1]*subx[17] + subx[16]*(Tbn[0,1]*subx[14] + Tbn[1,1]*subx[15]))**2/(subx[4]**2*subx[5]*(subx[5] + subx[8])**2)
    ret_cov[7] = cam_pos_R*subx[18]*subx[6]**4*sqrt(subx[8])*subx[9]*(-Tbn[2,0]*subx[17] + subx[16]*(Tbn[0,0]*subx[14] + Tbn[1,0]*subx[15]))/(hgt**3*subx[4]*(subx[4] + subx[7])*(subx[5] + subx[8])) + cam_pos_R*subx[19]*subx[6]**4*sqrt(subx[8])*subx[9]*(-Tbn[2,1]*subx[17] + subx[16]*(Tbn[0,1]*subx[14] + Tbn[1,1]*subx[15]))/(hgt**3*subx[4]*(subx[4] + subx[7])*(subx[5] + subx[8]))
    ret_cov[8] = 0
    ret_cov[9] = 0
    ret_cov[10] = 0
    ret_cov[11] = cam_pos_R*subx[18]**2*subx[6]**4*subx[7]/(hgt**4*subx[5]*(subx[4] + subx[7])**3) + cam_pos_R*subx[19]**2*subx[6]**4*subx[7]/(hgt**4*subx[5]*(subx[4] + subx[7])**3) + hgt_R*subx[7]/(subx[5]**2*(subx[4] + subx[7]))
    ret_cov[12] = 0
    ret_cov[13] = 0
    ret_cov[14] = 0
    ret_cov[15] = vel_R[0]
    ret_cov[16] = 0
    ret_cov[17] = 0
    ret_cov[18] = vel_R[1]
    ret_cov[19] = 0
    ret_cov[20] = vel_R[2]

    return ret_cov

def EKF_POLAR2NED_CALC_NED(polar):
    ret_ned = np.zeros(3)
    ret_ned[0] = 1.0*cos(polar[0])*cos(polar[1])/polar[2]
    ret_ned[1] = 1.0*sin(polar[0])*cos(polar[1])/polar[2]
    ret_ned[2] = -1.0*sin(polar[1])/polar[2]

    return ret_ned

def EKF_PREDICTION_CALC_STATE(P, dt, subx, u, w_u_sigma, x):
    ret_state = np.zeros(6)
    ret_state[0] = atan2(dt*x[4] + subx[0]*subx[1]*subx[3], dt*x[3] + subx[0]*subx[3]*subx[4])
    ret_state[1] = atan2(1.0*subx[0]*subx[5] - 1.0*subx[6], subx[10])
    ret_state[2] = subx[12]
    ret_state[3] = -u[0] + x[3]
    ret_state[4] = -u[1] + x[4]
    ret_state[5] = -u[2] + x[5]

    return ret_state

def EKF_PREDICTION_CALC_SUBX(P, dt, u, w_u_sigma, x):
    ret_subx = np.zeros(48)
    ret_subx[0] = 1/x[2]
    ret_subx[1] = sin(x[0])
    ret_subx[2] = cos(x[1])
    ret_subx[3] = 1.0*ret_subx[2]
    ret_subx[4] = cos(x[0])
    ret_subx[5] = sin(x[1])
    ret_subx[6] = dt*x[5]
    ret_subx[7] = dt*x[2]*x[3] + ret_subx[3]*ret_subx[4]
    ret_subx[8] = dt*x[2]*x[4] + ret_subx[1]*ret_subx[3]
    ret_subx[9] = ret_subx[7]**2 + ret_subx[8]**2
    ret_subx[10] = sqrt(ret_subx[9]/x[2]**2)
    ret_subx[11] = ret_subx[9] + (-1.0*ret_subx[5] + ret_subx[6]*x[2])**2
    ret_subx[12] = 1/sqrt(ret_subx[11]/x[2]**2)
    ret_subx[13] = 1/((dt*x[3] + ret_subx[0]*ret_subx[3]*ret_subx[4])**2 + (dt*x[4] + ret_subx[0]*ret_subx[1]*ret_subx[3])**2)
    ret_subx[14] = -ret_subx[0]*ret_subx[13]*ret_subx[1]*ret_subx[3]*(-dt*x[4] - ret_subx[0]*ret_subx[1]*ret_subx[3]) + ret_subx[0]*ret_subx[13]*ret_subx[3]*ret_subx[4]*(dt*x[3] + ret_subx[0]*ret_subx[3]*ret_subx[4])
    ret_subx[15] = -1.0*ret_subx[0]*ret_subx[13]*ret_subx[1]*ret_subx[5]*(dt*x[3] + ret_subx[0]*ret_subx[3]*ret_subx[4]) - 1.0*ret_subx[0]*ret_subx[13]*ret_subx[4]*ret_subx[5]*(-dt*x[4] - ret_subx[0]*ret_subx[1]*ret_subx[3])
    ret_subx[16] = -1.0*ret_subx[13]*ret_subx[1]*ret_subx[2]*(dt*x[3] + ret_subx[0]*ret_subx[3]*ret_subx[4])/x[2]**2 - 1.0*ret_subx[13]*ret_subx[2]*ret_subx[4]*(-dt*x[4] - ret_subx[0]*ret_subx[1]*ret_subx[3])/x[2]**2
    ret_subx[17] = dt*ret_subx[13]*(-dt*x[4] - ret_subx[0]*ret_subx[1]*ret_subx[3])
    ret_subx[18] = dt*ret_subx[13]*(dt*x[3] + ret_subx[0]*ret_subx[3]*ret_subx[4])
    ret_subx[19] = P[0]*ret_subx[14] + P[1]*ret_subx[15] + P[2]*ret_subx[16] + P[3]*ret_subx[17] + P[4]*ret_subx[18]
    ret_subx[20] = P[1]*ret_subx[14] + P[6]*ret_subx[15] + P[7]*ret_subx[16] + P[8]*ret_subx[17] + P[9]*ret_subx[18]
    ret_subx[21] = P[11]*ret_subx[16] + P[12]*ret_subx[17] + P[13]*ret_subx[18] + P[2]*ret_subx[14] + P[7]*ret_subx[15]
    ret_subx[22] = P[13]*ret_subx[16] + P[16]*ret_subx[17] + P[18]*ret_subx[18] + P[4]*ret_subx[14] + P[9]*ret_subx[15]
    ret_subx[23] = P[12]*ret_subx[16] + P[15]*ret_subx[17] + P[16]*ret_subx[18] + P[3]*ret_subx[14] + P[8]*ret_subx[15]
    ret_subx[24] = 1/(ret_subx[9]/x[2]**2 + (1.0*ret_subx[0]*ret_subx[5] - 1.0*ret_subx[6])**2)
    ret_subx[25] = 1.0*ret_subx[0]*ret_subx[10]*ret_subx[24]*ret_subx[2] + ret_subx[10]*ret_subx[24]*(-1.0*ret_subx[0]*ret_subx[5] + 1.0*ret_subx[6])*(-2.0*ret_subx[1]*ret_subx[5]*ret_subx[8] - 2.0*ret_subx[4]*ret_subx[5]*ret_subx[7])/(2*ret_subx[9])
    ret_subx[26] = dt*x[3]*(2*dt*x[2]*x[3] + 2.0*ret_subx[2]*ret_subx[4]) + dt*x[4]*(2*dt*x[2]*x[4] + 2.0*ret_subx[1]*ret_subx[2])
    ret_subx[27] = -1.0*ret_subx[10]*ret_subx[24]*ret_subx[5]/x[2]**2 + ret_subx[10]*ret_subx[24]*x[2]**2*(-1.0*ret_subx[0]*ret_subx[5] + 1.0*ret_subx[6])*(ret_subx[26]/(2*x[2]**2) - ret_subx[9]/x[2]**3)/ret_subx[9]
    ret_subx[28] = P[10]*ret_subx[15] + P[14]*ret_subx[16] + P[17]*ret_subx[17] + P[19]*ret_subx[18] + P[5]*ret_subx[14]
    ret_subx[29] = 1.0*dt*ret_subx[10]*ret_subx[24]
    ret_subx[30] = ret_subx[10]*ret_subx[24]*(-1.0*ret_subx[0]*ret_subx[5] + 1.0*ret_subx[6])*(-2.0*ret_subx[1]*ret_subx[2]*ret_subx[7] + 2.0*ret_subx[2]*ret_subx[4]*ret_subx[8])/(2*ret_subx[9])
    ret_subx[31] = dt*ret_subx[10]*ret_subx[24]*ret_subx[7]*x[2]*(-1.0*ret_subx[0]*ret_subx[5] + 1.0*ret_subx[6])/ret_subx[9]
    ret_subx[32] = dt*ret_subx[10]*ret_subx[24]*ret_subx[8]*x[2]*(-1.0*ret_subx[0]*ret_subx[5] + 1.0*ret_subx[6])/ret_subx[9]
    ret_subx[33] = ret_subx[12]*(-2.0*ret_subx[1]*ret_subx[2]*ret_subx[7] + 2.0*ret_subx[2]*ret_subx[4]*ret_subx[8])/(2*ret_subx[11])
    ret_subx[34] = ret_subx[12]*(-2.0*ret_subx[1]*ret_subx[5]*ret_subx[8] - 2.0*ret_subx[2]*(-1.0*ret_subx[5] + ret_subx[6]*x[2]) - 2.0*ret_subx[4]*ret_subx[5]*ret_subx[7])/(2*ret_subx[11])
    ret_subx[35] = ret_subx[12]*x[2]**2*(ret_subx[11]/x[2]**3 - (ret_subx[26] + ret_subx[6]*(2*dt*x[2]*x[5] - 2.0*ret_subx[5]))/(2*x[2]**2))/ret_subx[11]
    ret_subx[36] = dt*ret_subx[12]*x[2]*(-1.0*ret_subx[5] + ret_subx[6]*x[2])/ret_subx[11]
    ret_subx[37] = dt*ret_subx[12]*ret_subx[7]*x[2]/ret_subx[11]
    ret_subx[38] = dt*ret_subx[12]*ret_subx[8]*x[2]/ret_subx[11]
    ret_subx[39] = -P[10]*ret_subx[29] + P[1]*ret_subx[30] + P[6]*ret_subx[25] + P[7]*ret_subx[27] + P[8]*ret_subx[31] + P[9]*ret_subx[32]
    ret_subx[40] = P[11]*ret_subx[27] + P[12]*ret_subx[31] + P[13]*ret_subx[32] - P[14]*ret_subx[29] + P[2]*ret_subx[30] + P[7]*ret_subx[25]
    ret_subx[41] = P[10]*ret_subx[25] + P[14]*ret_subx[27] + P[17]*ret_subx[31] + P[19]*ret_subx[32] - P[20]*ret_subx[29] + P[5]*ret_subx[30]
    ret_subx[42] = P[0]*ret_subx[30] + P[1]*ret_subx[25] + P[2]*ret_subx[27] + P[3]*ret_subx[31] + P[4]*ret_subx[32] - P[5]*ret_subx[29]
    ret_subx[43] = P[12]*ret_subx[27] + P[15]*ret_subx[31] + P[16]*ret_subx[32] - P[17]*ret_subx[29] + P[3]*ret_subx[30] + P[8]*ret_subx[25]
    ret_subx[44] = P[13]*ret_subx[27] + P[16]*ret_subx[31] + P[18]*ret_subx[32] - P[19]*ret_subx[29] + P[4]*ret_subx[30] + P[9]*ret_subx[25]
    ret_subx[45] = -P[10]*ret_subx[34] + P[14]*ret_subx[35] - P[17]*ret_subx[37] - P[19]*ret_subx[38] - P[20]*ret_subx[36] - P[5]*ret_subx[33]
    ret_subx[46] = P[12]*ret_subx[35] - P[15]*ret_subx[37] - P[16]*ret_subx[38] - P[17]*ret_subx[36] - P[3]*ret_subx[33] - P[8]*ret_subx[34]
    ret_subx[47] = P[13]*ret_subx[35] - P[16]*ret_subx[37] - P[18]*ret_subx[38] - P[19]*ret_subx[36] - P[4]*ret_subx[33] - P[9]*ret_subx[34]

    return ret_subx

def EKF_PREDICTION_CALC_COV(P, dt, subx, u, w_u_sigma, x):
    ret_cov = np.zeros(21)
    ret_cov[0] = subx[14]*subx[19] + subx[15]*subx[20] + subx[16]*subx[21] + subx[17]*subx[23] + subx[18]*subx[22]
    ret_cov[1] = subx[19]*subx[30] + subx[20]*subx[25] + subx[21]*subx[27] + subx[22]*subx[32] + subx[23]*subx[31] - subx[28]*subx[29]
    ret_cov[2] = -subx[19]*subx[33] - subx[20]*subx[34] + subx[21]*subx[35] - subx[22]*subx[38] - subx[23]*subx[37] - subx[28]*subx[36]
    ret_cov[3] = subx[23]
    ret_cov[4] = subx[22]
    ret_cov[5] = subx[28]
    ret_cov[6] = subx[25]*subx[39] + subx[27]*subx[40] - subx[29]*subx[41] + subx[30]*subx[42] + subx[31]*subx[43] + subx[32]*subx[44]
    ret_cov[7] = -subx[33]*subx[42] - subx[34]*subx[39] + subx[35]*subx[40] - subx[36]*subx[41] - subx[37]*subx[43] - subx[38]*subx[44]
    ret_cov[8] = subx[43]
    ret_cov[9] = subx[44]
    ret_cov[10] = subx[41]
    ret_cov[11] = -subx[33]*(-P[0]*subx[33] - P[1]*subx[34] + P[2]*subx[35] - P[3]*subx[37] - P[4]*subx[38] - P[5]*subx[36]) - subx[34]*(-P[10]*subx[36] - P[1]*subx[33] - P[6]*subx[34] + P[7]*subx[35] - P[8]*subx[37] - P[9]*subx[38]) + subx[35]*(P[11]*subx[35] - P[12]*subx[37] - P[13]*subx[38] - P[14]*subx[36] - P[2]*subx[33] - P[7]*subx[34]) - subx[36]*subx[45] - subx[37]*subx[46] - subx[38]*subx[47]
    ret_cov[12] = subx[46]
    ret_cov[13] = subx[47]
    ret_cov[14] = subx[45]
    ret_cov[15] = P[15] + w_u_sigma[0]**2
    ret_cov[16] = P[16]
    ret_cov[17] = P[17]
    ret_cov[18] = P[18] + w_u_sigma[1]**2
    ret_cov[19] = P[19]
    ret_cov[20] = P[20] + w_u_sigma[2]**2

    return ret_cov

def EKF_VELD_CALC_NIS(P, R, subx, x, z):
    ret_NIS = 0
    ret_NIS = subx[0]*(x[5] + z)**2

    return ret_NIS

def EKF_VELD_CALC_STATE(P, R, subx, x, z):
    ret_state = np.zeros(6)
    ret_state[0] = -P[5]*subx[0]*(x[5] + z) + x[0]
    ret_state[1] = -P[10]*subx[0]*(x[5] + z) + x[1]
    ret_state[2] = -P[14]*subx[0]*(x[5] + z) + x[2]
    ret_state[3] = -P[17]*subx[0]*(x[5] + z) + x[3]
    ret_state[4] = -P[19]*subx[0]*(x[5] + z) + x[4]
    ret_state[5] = -subx[1]*(x[5] + z) + x[5]

    return ret_state

def EKF_VELD_CALC_COV(P, R, subx, x, z):
    ret_cov = np.zeros(21)
    ret_cov[0] = P[0] + P[5]**2*R*subx[2] - P[5]**2*subx[0] - P[5]*subx[0]*subx[3]
    ret_cov[1] = P[10]*P[5]*R*subx[2] - P[10]*P[5]*subx[0] - P[10]*subx[0]*subx[3] + P[1]
    ret_cov[2] = P[14]*P[5]*R*subx[2] - P[14]*P[5]*subx[0] - P[14]*subx[0]*subx[3] + P[2]
    ret_cov[3] = P[17]*P[5]*R*subx[2] - P[17]*P[5]*subx[0] - P[17]*subx[0]*subx[3] + P[3]
    ret_cov[4] = P[19]*P[5]*R*subx[2] - P[19]*P[5]*subx[0] - P[19]*subx[0]*subx[3] + P[4]
    ret_cov[5] = P[20]*P[5]*R*subx[2] + subx[3]*(-subx[1] + 1)
    ret_cov[6] = P[10]**2*R*subx[2] - P[10]**2*subx[0] - P[10]*subx[0]*(-P[10]*subx[1] + P[10]) + P[6]
    ret_cov[7] = P[10]*P[14]*R*subx[2] - P[10]*P[14]*subx[0] - P[14]*subx[0]*(-P[10]*subx[1] + P[10]) + P[7]
    ret_cov[8] = P[10]*P[17]*R*subx[2] - P[10]*P[17]*subx[0] - P[17]*subx[0]*(-P[10]*subx[1] + P[10]) + P[8]
    ret_cov[9] = P[10]*P[19]*R*subx[2] - P[10]*P[19]*subx[0] - P[19]*subx[0]*(-P[10]*subx[1] + P[10]) + P[9]
    ret_cov[10] = P[10]*P[20]*R*subx[2] + (-subx[1] + 1)*(-P[10]*subx[1] + P[10])
    ret_cov[11] = P[11] + P[14]**2*R*subx[2] - P[14]**2*subx[0] - P[14]*subx[0]*(-P[14]*subx[1] + P[14])
    ret_cov[12] = P[12] + P[14]*P[17]*R*subx[2] - P[14]*P[17]*subx[0] - P[17]*subx[0]*(-P[14]*subx[1] + P[14])
    ret_cov[13] = P[13] + P[14]*P[19]*R*subx[2] - P[14]*P[19]*subx[0] - P[19]*subx[0]*(-P[14]*subx[1] + P[14])
    ret_cov[14] = P[14]*P[20]*R*subx[2] + (-subx[1] + 1)*(-P[14]*subx[1] + P[14])
    ret_cov[15] = P[15] + P[17]**2*R*subx[2] - P[17]**2*subx[0] - P[17]*subx[0]*(-P[17]*subx[1] + P[17])
    ret_cov[16] = P[16] + P[17]*P[19]*R*subx[2] - P[17]*P[19]*subx[0] - P[19]*subx[0]*(-P[17]*subx[1] + P[17])
    ret_cov[17] = P[17]*P[20]*R*subx[2] + (-subx[1] + 1)*(-P[17]*subx[1] + P[17])
    ret_cov[18] = P[18] + P[19]**2*R*subx[2] - P[19]**2*subx[0] - P[19]*subx[0]*(-P[19]*subx[1] + P[19])
    ret_cov[19] = P[19]*P[20]*R*subx[2] + (-subx[1] + 1)*(-P[19]*subx[1] + P[19])
    ret_cov[20] = P[20]**2*R*subx[2] + P[20]*(-subx[1] + 1)**2

    return ret_cov

def EKF_VELD_CALC_SUBX(P, R, x, z):
    ret_subx = np.zeros(4)
    ret_subx[0] = 1/(P[20] + R)
    ret_subx[1] = P[20]*ret_subx[0]
    ret_subx[2] = ret_subx[0]**2
    ret_subx[3] = -P[5]*ret_subx[1] + P[5]

    return ret_subx

def EKF_VELD_CALC_INNOV(P, R, subx, x, z):
    ret_innov = 0
    ret_innov = x[5] + z

    return ret_innov

def EKF_VELNE_CALC_NIS(P, R, subx, x, z):
    ret_NIS = 0
    ret_NIS = (x[3] + z[0])*(subx[0]*(P[18] + R)*(x[3] + z[0]) - subx[1]*(x[4] + z[1])) + (x[4] + z[1])*(subx[0]*(P[15] + R)*(x[4] + z[1]) - subx[1]*(x[3] + z[0]))

    return ret_NIS

def EKF_VELNE_CALC_STATE(P, R, subx, x, z):
    ret_state = np.zeros(6)
    ret_state[0] = subx[2]*(x[4] + z[1]) + subx[3]*(x[3] + z[0]) + x[0]
    ret_state[1] = subx[4]*(x[4] + z[1]) + subx[5]*(x[3] + z[0]) + x[1]
    ret_state[2] = subx[6]*(x[4] + z[1]) + subx[7]*(x[3] + z[0]) + x[2]
    ret_state[3] = subx[8]*(x[3] + z[0]) + subx[9]*(x[4] + z[1]) + x[3]
    ret_state[4] = subx[10]*(x[4] + z[1]) + subx[11]*(x[3] + z[0]) + x[4]
    ret_state[5] = subx[12]*(x[4] + z[1]) + subx[13]*(x[3] + z[0]) + x[5]

    return ret_state

def EKF_VELNE_CALC_COV(P, R, subx, x, z):
    ret_cov = np.zeros(21)
    ret_cov[0] = P[0] + P[3]*subx[3] + P[4]*subx[2] + R*subx[2]**2 + R*subx[3]**2 + subx[14]*subx[2] + subx[15]*subx[3]
    ret_cov[1] = P[1] + P[8]*subx[3] + P[9]*subx[2] + R*subx[2]*subx[4] + R*subx[3]*subx[5] + subx[14]*subx[4] + subx[15]*subx[5]
    ret_cov[2] = P[12]*subx[3] + P[13]*subx[2] + P[2] + R*subx[2]*subx[6] + R*subx[3]*subx[7] + subx[14]*subx[6] + subx[15]*subx[7]
    ret_cov[3] = R*subx[2]*subx[9] + R*subx[3]*subx[8] + subx[14]*subx[9] + subx[15]*subx[16]
    ret_cov[4] = R*subx[10]*subx[2] + R*subx[11]*subx[3] + subx[11]*subx[15] + subx[14]*(subx[10] + 1)
    ret_cov[5] = P[17]*subx[3] + P[19]*subx[2] + P[5] + R*subx[12]*subx[2] + R*subx[13]*subx[3] + subx[12]*subx[14] + subx[13]*subx[15]
    ret_cov[6] = P[6] + P[8]*subx[5] + P[9]*subx[4] + R*subx[4]**2 + R*subx[5]**2 + subx[17]*subx[4] + subx[18]*subx[5]
    ret_cov[7] = P[12]*subx[5] + P[13]*subx[4] + P[7] + R*subx[4]*subx[6] + R*subx[5]*subx[7] + subx[17]*subx[6] + subx[18]*subx[7]
    ret_cov[8] = R*subx[4]*subx[9] + R*subx[5]*subx[8] + subx[16]*subx[18] + subx[17]*subx[9]
    ret_cov[9] = R*subx[10]*subx[4] + R*subx[11]*subx[5] + subx[11]*subx[18] + subx[17]*(subx[10] + 1)
    ret_cov[10] = P[10] + P[17]*subx[5] + P[19]*subx[4] + R*subx[12]*subx[4] + R*subx[13]*subx[5] + subx[12]*subx[17] + subx[13]*subx[18]
    ret_cov[11] = P[11] + P[12]*subx[7] + P[13]*subx[6] + R*subx[6]**2 + R*subx[7]**2 + subx[19]*subx[6] + subx[20]*subx[7]
    ret_cov[12] = R*subx[6]*subx[9] + R*subx[7]*subx[8] + subx[16]*subx[20] + subx[19]*subx[9]
    ret_cov[13] = R*subx[10]*subx[6] + R*subx[11]*subx[7] + subx[11]*subx[20] + subx[19]*(subx[10] + 1)
    ret_cov[14] = P[14] + P[17]*subx[7] + P[19]*subx[6] + R*subx[12]*subx[6] + R*subx[13]*subx[7] + subx[12]*subx[19] + subx[13]*subx[20]
    ret_cov[15] = R*subx[8]**2 + R*subx[9]**2 + subx[16]*(P[15]*subx[16] + P[16]*subx[9]) + subx[9]*(P[16]*subx[16] + P[18]*subx[9])
    ret_cov[16] = R*subx[10]*subx[9] + R*subx[11]*subx[8] + subx[11]*(P[15]*subx[16] + P[16]*subx[9]) + (subx[10] + 1)*(P[16]*subx[16] + P[18]*subx[9])
    ret_cov[17] = P[17]*subx[16] + P[19]*subx[9] + R*subx[12]*subx[9] + R*subx[13]*subx[8] + subx[12]*(P[16]*subx[16] + P[18]*subx[9]) + subx[13]*(P[15]*subx[16] + P[16]*subx[9])
    ret_cov[18] = R*subx[10]**2 + R*subx[11]**2 + subx[11]*(P[15]*subx[11] + P[16]*(subx[10] + 1)) + (subx[10] + 1)*(P[16]*subx[11] + P[18]*(subx[10] + 1))
    ret_cov[19] = P[17]*subx[11] + P[19]*(subx[10] + 1) + R*subx[10]*subx[12] + R*subx[11]*subx[13] + subx[12]*(P[16]*subx[11] + P[18]*(subx[10] + 1)) + subx[13]*(P[15]*subx[11] + P[16]*(subx[10] + 1))
    ret_cov[20] = P[17]*subx[13] + P[19]*subx[12] + P[20] + R*subx[12]**2 + R*subx[13]**2 + subx[12]*(P[16]*subx[13] + P[18]*subx[12] + P[19]) + subx[13]*(P[15]*subx[13] + P[16]*subx[12] + P[17])

    return ret_cov

def EKF_VELNE_CALC_SUBX(P, R, x, z):
    ret_subx = np.zeros(21)
    ret_subx[0] = 1/(-P[16]**2 + (P[15] + R)*(P[18] + R))
    ret_subx[1] = P[16]*ret_subx[0]
    ret_subx[2] = P[3]*ret_subx[1] - P[4]*ret_subx[0]*(P[15] + R)
    ret_subx[3] = -P[3]*ret_subx[0]*(P[18] + R) + P[4]*ret_subx[1]
    ret_subx[4] = P[8]*ret_subx[1] - P[9]*ret_subx[0]*(P[15] + R)
    ret_subx[5] = -P[8]*ret_subx[0]*(P[18] + R) + P[9]*ret_subx[1]
    ret_subx[6] = P[12]*ret_subx[1] - P[13]*ret_subx[0]*(P[15] + R)
    ret_subx[7] = -P[12]*ret_subx[0]*(P[18] + R) + P[13]*ret_subx[1]
    ret_subx[8] = -P[15]*ret_subx[0]*(P[18] + R) + P[16]**2*ret_subx[0]
    ret_subx[9] = P[15]*ret_subx[1] - P[16]*ret_subx[0]*(P[15] + R)
    ret_subx[10] = P[16]**2*ret_subx[0] - P[18]*ret_subx[0]*(P[15] + R)
    ret_subx[11] = P[18]*ret_subx[1] - ret_subx[1]*(P[18] + R)
    ret_subx[12] = P[17]*ret_subx[1] - P[19]*ret_subx[0]*(P[15] + R)
    ret_subx[13] = -P[17]*ret_subx[0]*(P[18] + R) + P[19]*ret_subx[1]
    ret_subx[14] = P[16]*ret_subx[3] + P[18]*ret_subx[2] + P[4]
    ret_subx[15] = P[15]*ret_subx[3] + P[16]*ret_subx[2] + P[3]
    ret_subx[16] = ret_subx[8] + 1
    ret_subx[17] = P[16]*ret_subx[5] + P[18]*ret_subx[4] + P[9]
    ret_subx[18] = P[15]*ret_subx[5] + P[16]*ret_subx[4] + P[8]
    ret_subx[19] = P[13] + P[16]*ret_subx[7] + P[18]*ret_subx[6]
    ret_subx[20] = P[12] + P[15]*ret_subx[7] + P[16]*ret_subx[6]

    return ret_subx

def EKF_VELNE_CALC_INNOV(P, R, subx, x, z):
    ret_innov = np.zeros(2)
    ret_innov[0] = x[3] + z[0]
    ret_innov[1] = x[4] + z[1]

    return ret_innov

