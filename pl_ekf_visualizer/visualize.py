from ekf import *
import numpy as np
import math
import sympy as sp
import visual as vpy
from time import sleep

class PLEKF:
    def __init__(self,dt):
        self.ang_R = radians(1.)**2
        self.hgt_R = 4.**2
        self.init_ang_sca = .98
        self.init_ang_sca_R = 0.02**2
        self.vel_R = np.array([0.2**2,0.2**2,0.1**2])
        self.w_u_sigma = np.array([.25*dt,.25*dt,.25*dt])

        self.state = np.zeros(EKF_NUM_STATES)
        self.cov = np.zeros(EKF_NUM_STATES*(EKF_NUM_STATES+1)/2)

    def initialize(self,Tbn,cam,hgt,vel):
        subx = EKF_INITIALIZATION_CALC_SUBX(Tbn, cam, self.ang_R, hgt, self.hgt_R, self.init_ang_sca, self.init_ang_sca_R, vel, self.vel_R)
        state = EKF_INITIALIZATION_CALC_STATE(Tbn, cam, self.ang_R, hgt, self.hgt_R, self.init_ang_sca, self.init_ang_sca_R, subx, vel, self.vel_R)
        cov = EKF_INITIALIZATION_CALC_COV(Tbn, cam, self.ang_R, hgt, self.hgt_R, self.init_ang_sca, self.init_ang_sca_R, subx, vel, self.vel_R)

        self.state = state
        self.cov = cov

    def predict(self, dt, delVel):
        subx = EKF_PREDICTION_CALC_SUBX(self.cov, dt, delVel, self.w_u_sigma, self.state)
        state = EKF_PREDICTION_CALC_STATE(self.cov, dt, subx, delVel, self.w_u_sigma, self.state)
        cov = EKF_PREDICTION_CALC_COV(self.cov, dt, subx, delVel, self.w_u_sigma, self.state)
        self.state = state
        self.cov = cov

    def fuseCamera(self, Tbn, cam):
        subx = EKF_CAMERA_CALC_SUBX(self.cov, self.ang_R, Tbn, self.state, cam)
        state = EKF_CAMERA_CALC_STATE(self.cov, self.ang_R, Tbn, subx, self.state, cam)
        cov = EKF_CAMERA_CALC_COV(self.cov, self.ang_R, Tbn, subx, self.state, cam)

        self.state = state
        self.cov = cov

    def fuseHeight(self, hgt):
        subx = EKF_HEIGHT_CALC_SUBX(self.cov, self.hgt_R, self.state, hgt)
        state = EKF_HEIGHT_CALC_STATE(self.cov, self.hgt_R, subx, self.state, hgt)
        cov = EKF_HEIGHT_CALC_COV(self.cov, self.hgt_R, subx, self.state, hgt)

        self.state = state
        self.cov = cov

    def fuseVertVel(self, vertVel):
        subx = EKF_VELD_CALC_SUBX(self.cov, self.vel_R[2], self.state, vertVel)
        state = EKF_VELD_CALC_STATE(self.cov, self.vel_R[2], subx, self.state, vertVel)
        cov = EKF_VELD_CALC_COV(self.cov, self.vel_R[2], subx, self.state, vertVel)

        self.state = state
        self.cov = cov

    def get_pos(self):
        return np.array([self.state[EKF_STATE_IDX_PT_N],self.state[EKF_STATE_IDX_PT_E],self.state[EKF_STATE_IDX_PT_D]])

    def get_vel(self):
        return np.array([self.state[EKF_STATE_IDX_VT_N],self.state[EKF_STATE_IDX_VT_E],self.state[EKF_STATE_IDX_VT_D]])

    def get_state(self):
        return self.state

    def get_cov(self):
        N = EKF_NUM_STATES
        ret = np.zeros((N,N))
        r = lambda k: int(math.floor((2*N+1-math.sqrt((2*N+1)*(2*N+1)-8*k))/2))
        c = lambda k: int(k - N*r(k) + r(k)*(r(k)-1)/2 + r(k))

        for i in range(len(self.cov)):
            ret[r(i),c(i)] = ret[c(i),r(i)] = self.cov[i]

        return ret

    def get_pos_cov(self):
        N = EKF_NUM_STATES

        stateindices = [EKF_STATE_IDX_PT_N,EKF_STATE_IDX_PT_E,EKF_STATE_IDX_PT_D]
        ret = np.zeros((3,3))
        cov = self.get_cov()
        for r_in,r_out in zip(stateindices,range(3)):
            for c_in,c_out in zip(stateindices,range(3)):
                ret[r_out,c_out] = cov[r_in,c_in]
        return ret

class Copter:
    def __init__(self, ekf):
        self._computeLambdas()
        self.pos = np.array([[0.],[10.],[-10.]])
        self.vel = np.array([[0.],[0.],[0.]])
        self.height_dem = -10.
        self.t = 0.
        self.camera_meas_noise = math.radians(1.)
        self.accel_meas_noise = .25
        self.height_meas_noise = .2
        self.vel_meas_noise = .5
        self.ekf = ekf
        self.do_camera_meas()
        self.do_height_meas()
        self.do_accel_meas()
        self.do_vel_meas()

    def _computeLambdas(self):
        t_sym = sp.Symbol('t')
        freq = 0.1
        spd = 3.
        copterPos = sp.Matrix([0.,spd*sp.cos(freq*2.*math.pi*t_sym)/(2.*math.pi*freq),-10.+0.5*t_sym])
        #copterPos = sp.Matrix([0.,0.,-10.+0.5*t])

        copterVel = sp.diff(copterPos,t_sym)
        copterAccel = sp.diff(copterVel,t_sym)

        self.copterPosLambda = lambdify(t_sym,copterPos)
        self.copterVelLambda = lambdify(t_sym,copterVel)
        self.copterAccelLambda = lambdify(t_sym,copterAccel)

    def update(self,dt):
        self.t += dt
        self.vel += self.get_accel()*dt
        self.pos += self.vel*dt

    def do_accel_meas(self):
        self.accel_meas = self.get_accel()+np.array([[random.gauss(0.,self.accel_meas_noise)] for _ in range(3)])

    def get_accel_meas(self):
        return self.accel_meas

    def do_camera_meas(self):
        relTargetPos = -self.get_pos()
        self.camera_meas = np.matrix([relTargetPos[0]/relTargetPos[2], relTargetPos[1]/relTargetPos[2]])+np.array([[random.gauss(0.,self.camera_meas_noise)] for _ in range(2)])

    def get_camera_meas(self):
        return self.camera_meas

    def do_height_meas(self):
        self.height_meas = -self.get_pos()[2]+random.gauss(0.,self.height_meas_noise)

    def get_height_meas(self):
        return self.height_meas

    def do_vel_meas(self):
        self.vel_meas = self.get_vel()+np.array([[random.gauss(0.,self.vel_meas_noise)] for _ in range(3)])

    def get_vel_meas(self):
        return self.vel_meas

    def get_pos(self):
        return self.pos
        #return self.copterPosLambda(self.t)

    def get_vel(self):
        return self.vel
        #return self.copterVelLambda(self.t)

    def get_accel(self):
        kP = 1.
        kD = 2.
        estPos = self.ekf.get_pos()
        estVel = self.ekf.get_vel()
        return np.array([estPos[0]*kP+estVel[0]*kD, estPos[1]*kP+estVel[1]*kD, 0.5-self.vel[2]])
        #return self.copterAccelLambda(self.t)

from helpers import tovpy
import random
import subprocess
class EKFVis:
    def __init__(self,ekf,copter):
        self.frame = 0
        self.ekf = ekf
        self.copter = copter
        self.scene = vpy.display(title = 'EKF visualizer', width=1920+9, height=1080+30, center = tovpy((0.,0.,-5.)))
        self.scene.autoscale = False
        self.scene.range = 16
        self.scene.forward = tovpy((cos(math.radians(0.)),0.,sin(math.radians(0.))))
        self.scene.ambient = vpy.color.gray(0.6)
        self.scene.background = vpy.color.white
        self.copterEst = vpy.sphere(radius=0.1, color=vpy.color.blue)
        self.copterTruth = vpy.sphere(radius=0.1, color=vpy.color.red)
        self.ground = vpy.cylinder(pos=tovpy((0.,0.,0.1)), axis=tovpy((0.,0.,0.1)), radius=12., material=vpy.materials.wood)
        self.landingPad = vpy.cylinder(axis=tovpy((0.,0.,0.1)), radius=0.3, color=vpy.color.red)
        self.heightLines = [vpy.cylinder(pos=tovpy((0.,-10.,z)), axis=tovpy((0.,20.,0.)), radius=0.04, opacity=0.4) for z in np.linspace(0.,-10.,10.)]

        self.cov_points = [np.matrix([[random.gauss(0.,1.)],[random.gauss(0.,1.)],[random.gauss(0.,1.)]]) for _ in range(200)]
        self.cov_cloud = [vpy.sphere(radius=0.03, color=vpy.color.black,opacity=0.1) for _ in range(len(self.cov_points))]

        self.camMeasInd = vpy.arrow(shaftwidth=0.08, color=vpy.color.green, opacity=0.4)
        self.screenshotProcess = None

    def update(self):
        if self.screenshotProcess is not None:
            self.screenshotProcess.wait()
        self.copterEst.pos = tovpy(-self.ekf.get_pos())
        self.copterTruth.pos = tovpy(self.copter.get_pos())
        cameraMeas = self.copter.get_camera_meas()
        self.camMeasInd.axis = tovpy((cameraMeas[0],cameraMeas[1],1.))
        self.camMeasInd.pos = self.copterTruth.pos

        noise_transform = np.linalg.cholesky(self.ekf.get_pos_cov())

        for i in range(len(self.cov_cloud)):
            self.cov_cloud[i].pos = tovpy(noise_transform*self.cov_points[i]-self.ekf.get_pos())

        vpy.rate(60)

        self.frame += 1
        #self.screenshotProcess = subprocess.Popen(['import', '-window', '0x5400003', 'frames/vp%04u.png' % (self.frame,)], shell=False)

t = 0.
dt = 1./60.

ekf = PLEKF(dt)
copter = Copter(ekf)
ekf.initialize(np.eye(3),copter.get_camera_meas(), copter.get_height_meas()-5., copter.get_vel())

vis = EKFVis(ekf, copter)
vis.update()

# fly 1 m/s for 10m
while -copter.get_pos()[2] > 0.:
    for _ in range(2):
        t += dt
        copter.update(dt)
        copter.do_accel_meas()
        ekf.predict(dt, (copter.get_accel_meas())*dt)
        copter.do_vel_meas()
        ekf.fuseVertVel(copter.get_vel_meas()[2])
        vis.update()
    copter.do_height_meas()
    ekf.fuseHeight(copter.get_height_meas())
    copter.do_camera_meas()
    ekf.fuseCamera(np.eye(3), copter.get_camera_meas())
