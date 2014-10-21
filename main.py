# -*- coding: utf-8 -*-
"""
Simulation program
"""
import numpy as np
import pylab as pl

class Controller:
    def __init__(self, dt):
        self.dt = dt
        #for pid
        self.error_integral = 0
        self.error_diff = 0
        self.error_diff_ = 0
        self.error_ = 0
        #for DOB
        self.tmp1_dob = 0.0
        self.tmp2_dob = 0.0            
            
    def pid_controller(self, error, Kp, Ki, Kd, dt, wc_diff):
        self.error_diff = (self.error_diff_ + wc_diff * (error - self.error_)) / (1.0 + wc_diff * self.dt) 
        self.error_integral = self.error_integral + error * self.dt
        
        input = Kp * error + Kd * self.error_diff + Ki * self.error_integral;

        self.error_diff_ = self.error_diff
        self.error_ = error        

        return input        
        
    def disturbance_observer(self, Jm, tau, dq, wc_dob):
        self.tmp1_dob = tau + Jm * wc_dob * dq;
        self.tmp2_dob = (self.tmp2_dob + self.dt * wc_dob * self.tmp1_dob) / (1 + self.dt * wc_dob);
        dist = self.tmp2_dob - Jm * wc_dob * dq; 

        return dist

class RigidRotor:

    def __init__(self,
                 inertia,
                 damping,
                 torque,
                 sampling_time,
                 control_sampling_time
                 ):
        # plant parameters and simulation conditions
        self.inertia = inertia
        self.damping = damping
        self.torque = torque
        self.sampling_time = sampling_time
        self.control_sampling_time = control_sampling_time

        self.xvec = [0.0, 0.0]

        # for controller
        self.control_sampling_time = 0.001

    # calculate derivative of state value
    def calc_derivative(self, torque_reac):
        return (self.xvec[1], - self.damping / self.inertia * self.xvec[1] + (1.0 / self.inertia) * (self.torque - torque_reac) )

    # update 
    def update(self):
        self.xvec[0] = self.xvec[0] + self.dxvec[0] * self.sampling_time
        self.xvec[1] = self.xvec[1] + self.dxvec[1] * self.sampling_time

# simulation parameters
inertia=0.1
damping=0.0
torque = 0.0

sampling_time = 0.001
control_sampling_time=0.001
simulation_time = 3

# state parameters and data for plot
xvec0_data=[]
xvec1_data=[]
cmd_data=[]

# simulation object
sim = RigidRotor(inertia, damping, torque, sampling_time, control_sampling_time)
controller = Controller(control_sampling_time)

# main loop 10[sec]
for i in range(simulation_time*(int)(1/sim.sampling_time)):
    time = i * sim.sampling_time

    control_delay = (int)(sim.control_sampling_time/sim.sampling_time) #[sample]
    if i % control_delay == 0:
        """ controller """
        # definition for control parameters
        if i == 0 :
            dist = 0.0
            torque_reac = 0.0

        # position controller
        # pole assignment (second order vibration system)
        wn = 10.0
        zeta = 0.5
        kp = wn**2 * sim.inertia
        kd = 2.0 * zeta * wn * sim.inertia - sim.damping

        # PID control
        theta_cmd = 1.0
        theta_res = sim.xvec[0]
        error = theta_cmd - theta_res
        sim.torque = controller.pid_controller(error, kp, 0, kd, sim.control_sampling_time, 500.0)        

        sim.torque = sim.torque + dist

        # DOB
        dist = controller.disturbance_observer(inertia, sim.torque, sim.xvec[1], 300.0)

        #data update
        xvec0_data.append(sim.xvec[0])
        xvec1_data.append(sim.xvec[1])
        cmd_data.append(theta_cmd)

        """ controller end """

    """ plant """
    #reaction torque
    if time > 1.5:
        torque_reac = 10.0

    # derivative calculation
    sim.dxvec = sim.calc_derivative(torque_reac)
    # euler-integration
    sim.update()

    """ plant end """

# data plot
time_data = np.arange(0, 3, sim.control_sampling_time)
pl.plot(time_data, cmd_data[:], label="theta cmd")
pl.plot(time_data, xvec0_data[:], label="theta res")
pl.legend()
pl.grid()
pl.xlabel('time [s]')
pl.ylabel('angle [rad]')
pl.show()
