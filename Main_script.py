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

plt.figure()
plt.plot(sol.t,sol.y[2,:])
plt.show()