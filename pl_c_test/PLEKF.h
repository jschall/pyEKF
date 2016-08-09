#pragma once

#include <math.h>
#include "ekf_defines.h"

#define P_ARRAY_SIZE (EKF_NUM_STATES*(EKF_NUM_STATES-1)/2+EKF_NUM_STATES)
#define P_IDX(__ROW,__COL) (__ROW<=__COL) ? (__ROW*EKF_NUM_STATES-(__ROW-1)*__ROW/2+__COL-__ROW) : (__COL*EKF_NUM_STATES-(__COL-1)*__COL/2+__ROW-__COL)

#define PLEKF_CAM_POS_R ((float)M_PI/180.0f * (float)M_PI/180.0f)
#define PLEKF_HGT_R (4.0f*4.0f)
#define PLEKF_VELNE_R (1.0f*1.0f)
#define PLEKF_VELD_R (1.0f*1.0f)
#define PLEKF_INIT_FOC_LEN_R (0.02f*0.02f)
#define PLEKF_VEHICLE_ACCEL_SIGMA 0.25f

class PLEKF {
public:
    struct initialize_params {
        float Tbn[3][3];
        float cam_pos[2];
        float height;
        float velNED[3];
    };
    
    struct predict_params {
        float dt;
        float delVelNED[3];
    };
    
    struct fuse_camera_params {
        float Tbn[3][3];
        float cam_pos[2];
    };
    
    struct fuse_velNE_params {
        float velNE[2];
    };
    
    struct fuse_velD_params {
        float velD;
    };
    
    struct fuse_height_params {
        float height;
    };
    
    void initialize(const struct initialize_params& params);
    void predict(const struct predict_params& params);
    float prepFuseCamera(const struct fuse_camera_params& params);
    float prepFuseVelNE(const struct fuse_velNE_params& params);
    float prepFuseVelD(const struct fuse_velD_params& params);
    float prepFuseHeight(const struct fuse_height_params& params);

    void commitOperation();
//private:
    float _state[EKF_NUM_STATES];
    float _cov[P_ARRAY_SIZE];

    float _subx[EKF_MAX_NUM_SUBX];
    float _next_state[EKF_NUM_STATES];
    float _next_cov[P_ARRAY_SIZE];
};
