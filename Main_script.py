# Label code, initialize all class attributes, extract solution

import os
from copy import deepcopy
import numpy as np
import Classes
from Functions import RocketDynModel
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

Sim = Classes.SimulationProperties()
Global = Classes.GlobalProperties()
Rocket = Classes.RocketProperties()

StateVecInit= np.hstack((
    np.array([0,0,Rocket.LaunchAlt]) + Rocket.InitDistAlongRail*Rocket.LaunchVec,
    np.array([0,-Rocket.LaunchAngle,0]),
    np.zeros(3),
    np.zeros(3)))

sol = solve_ivp(RocketDynModel,[0,Sim.SimTime],StateVecInit,args=(Global,Rocket),t_eval=np.arange(0,Sim.SimTime + Sim.TimeStepSize,Sim.TimeStepSize))

# Extract solution
TimeStep = []
for TimeStepNum in range(Sim.NumTimeSteps):
    Time = sol.t[TimeStepNum]
    StateVec = sol.y[:,TimeStepNum]
    RocketDynModel(Time,StateVec,Global,Rocket)
    TimeStep.append(Classes.TimeSeriesData(Time,deepcopy(Global),deepcopy(Rocket)))
    pass

for TimeStepNum in range(Sim.NumTimeSteps):
    print(TimeStep[TimeStepNum].Rocket.PosVec)
    pass

#plt.figure()
#plt.plot(Sim.TimeStepVec,Rocket.PosVecTraj[2,:])
#plt.show()