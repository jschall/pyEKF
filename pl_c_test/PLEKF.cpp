#include "PLEKF.h"
#include <string.h>

void PLEKF::initialize(const struct initialize_params& params)
{
    const float hgt_R = PLEKF_HGT_R;
    const float cam_pos_R = PLEKF_CAM_POS_R;
    const float vel_R[3] = {PLEKF_VELNE_R,PLEKF_VELNE_R,PLEKF_VELD_R};
    const float init_foc_len = 1.0f;
    const float init_foc_len_R = PLEKF_INIT_FOC_LEN_R;

    EKF_INITIALIZATION_CALC_SUBX(params.Tbn, params.cam_pos, cam_pos_R, params.height, hgt_R, init_foc_len, init_foc_len_R, params.velNED, vel_R, _subx);
    EKF_INITIALIZATION_CALC_STATE(params.Tbn, params.cam_pos, cam_pos_R, params.height, hgt_R, init_foc_len, init_foc_len_R, _subx, params.velNED, vel_R, _next_state);
    EKF_INITIALIZATION_CALC_COV(params.Tbn, params.cam_pos, cam_pos_R, params.height, hgt_R, init_foc_len, init_foc_len_R, _subx, params.velNED, vel_R, _next_cov);

    commitOperation();
}

void PLEKF::predict(const struct predict_params& params)
{
    float u[EKF_NUM_CONTROL_INPUTS];
    u[EKF_U_IDX_DVV_N] = params.delVelNED[0];
    u[EKF_U_IDX_DVV_E] = params.delVelNED[1];
    u[EKF_U_IDX_DVV_D] = params.delVelNED[2];

    float w_u_sigma[EKF_NUM_CONTROL_INPUTS];
    w_u_sigma[EKF_U_IDX_DVV_N] = PLEKF_VEHICLE_ACCEL_SIGMA*params.dt;
    w_u_sigma[EKF_U_IDX_DVV_E] = PLEKF_VEHICLE_ACCEL_SIGMA*params.dt;
    w_u_sigma[EKF_U_IDX_DVV_D] = PLEKF_VEHICLE_ACCEL_SIGMA*params.dt;

    EKF_PREDICTION_CALC_SUBX(_cov, params.dt, u, w_u_sigma, _state, _subx)
    EKF_PREDICTION_CALC_STATE(_cov, params.dt, _subx, u, w_u_sigma, _state, _next_state)
    EKF_PREDICTION_CALC_COV(_cov, params.dt, _subx, u, w_u_sigma, _state, _next_cov)
    
    commitOperation();
}

float PLEKF::prepFuseCamera(const struct fuse_camera_params& params)
{
    const float cam_pos_R = PLEKF_CAM_POS_R;
    float NIS;
    EKF_CAMERA_CALC_SUBX(_cov, cam_pos_R, params.Tbn, _state, params.cam_pos, _subx)
    EKF_CAMERA_CALC_STATE(_cov, cam_pos_R, params.Tbn, _subx, _state, params.cam_pos, _next_state)
    EKF_CAMERA_CALC_COV(_cov, cam_pos_R, params.Tbn, _subx, _state, params.cam_pos, _next_cov)
    EKF_CAMERA_CALC_NIS(_cov, cam_pos_R, params.Tbn, _subx, _state, params.cam_pos, NIS)
    
    return NIS;
}

float PLEKF::prepFuseVelNE(const struct fuse_velNE_params& params)
{
    const float velNE_R = PLEKF_VELNE_R;
    float NIS;
    EKF_VELNE_CALC_SUBX(_cov, velNE_R, _state, params.velNE, _subx)
    EKF_VELNE_CALC_STATE(_cov, velNE_R, _subx, _state, params.velNE, _next_state)
    EKF_VELNE_CALC_COV(_cov, velNE_R, _subx, _state, params.velNE, _next_cov)
    EKF_VELNE_CALC_NIS(_cov, velNE_R, _subx, _state, params.velNE, NIS);
    return NIS;
}

float PLEKF::prepFuseVelD(const struct fuse_velD_params& params)
{
    const float velD_R = PLEKF_VELD_R;
    float NIS;
    EKF_VELD_CALC_SUBX(_cov, velD_R, _state, params.velD, _subx)
    EKF_VELD_CALC_STATE(_cov, velD_R, _subx, _state, params.velD, _next_state)
    EKF_VELD_CALC_COV(_cov, velD_R, _subx, _state, params.velD, _next_cov)
    EKF_VELD_CALC_NIS(_cov, velD_R, _subx, _state, params.velD, NIS);
    return NIS;
}

float PLEKF::prepFuseHeight(const struct fuse_height_params& params)
{
    const float hgt_R = PLEKF_HGT_R;
    float NIS;
    EKF_HEIGHT_CALC_SUBX(_cov, hgt_R, _state, params.height, _subx)
    EKF_HEIGHT_CALC_STATE(_cov, hgt_R, _subx, _state, params.height, _next_state)
    EKF_HEIGHT_CALC_COV(_cov, hgt_R, _subx, _state, params.height, _next_cov)
    EKF_HEIGHT_CALC_NIS(_cov, hgt_R, _subx, _state, params.height, NIS);
    return NIS;
}

void PLEKF::commitOperation()
{
    memcpy(_state, _next_state, sizeof(_state));
    memcpy(_cov, _next_cov, sizeof(_cov));
}
