# -*- coding: utf-8 -*-
"""
Simulation program
"""
import numpy as np
import pylab as pl

class RigidRotor:

    def __init__(self,
                 inertia,
                 friction,
                 torque,
                 sampling_time,
                 control_sampling_time
                 ):
        # plant parameters and simulation conditions
        self.inertia = inertia
        self.friction = friction
        self.torque = torque
        self.sampling_time = sampling_time
        self.control_sampling_time = control_sampling_time

        self.xvec = [0.0, 0.0]

        # for controller
        self.control_sampling_time = 0.001

    # calculate derivative of state value
    def calc_derivative(self, torque_reac):
        return (self.xvec[1], - self.friction / self.inertia * self.xvec[1] + (1.0 / self.inertia) * (self.torque - torque_reac) )

    # update
    def update(self):
        self.xvec[0] = self.xvec[0] + self.dxvec[0] * self.sampling_time
        self.xvec[1] = self.xvec[1] + self.dxvec[1] * self.sampling_time

# simulation parameters
inertia=0.1
friction=0.0
torque = 0.0

sampling_time = 0.001
control_sampling_time=0.001
simulation_time = 3

# state parameters and data for plot
xvec0_data=[]
xvec1_data=[]
cmd_data=[]



# simulation object
sim = RigidRotor(inertia, friction, torque, sampling_time, control_sampling_time)

# main loop 10[sec]
for i in range(simulation_time*(int)(1/sim.sampling_time)):
    time = i * sim.sampling_time

    control_delay = (int)(sim.control_sampling_time/sim.sampling_time) #[sample]
    if i%control_delay == 0:
        """ controller """
        # definition for control parameters
        if i == 0 :
            pass

        # position controller
        # pole assignment (second order vibration system)
        wn = 10.0
        zeta = 0.3
        kp = wn**2 * sim.inertia
        kd = 2.0 * zeta * wn * sim.inertia - sim.friction

        # P control
        theta_cmd = 1.0
        theta_res = sim.xvec[0]
        sim.torque = kp * (theta_cmd - theta_res)
        # D control
        dtheta_cmd = 0.0
        dtheta_res = sim.xvec[1]
        sim.torque = sim.torque + kd * (dtheta_cmd - dtheta_res)

        #data update
        xvec0_data.append(sim.xvec[0])
        xvec1_data.append(sim.xvec[1])
        cmd_data.append(theta_cmd)

        """ controller end """

    """ plant """
    #reaction torque
    torque_reac = 0.0

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
