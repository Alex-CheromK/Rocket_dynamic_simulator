### Main script descript:
# This file simulates the full flight trajectory of a liquid propellant rocket.
# The code is written in the Python programming language.
# The software has been created by Dr. Ali C. Kheirabadi for the purpose of...
# simulating and validating the range of motion of the University of British Columbia's...
# Rocket Design Team's liquid propellant rocket.
# The simulation relies on standard 3D nonlinear dynamics modelling techniques...
# including Newtonian mechanics for generating the translational and rotational...
# equations of motion of the rocket.

### Import standard Python libraries
# "copy" library for making deep copies of class objects for resulting purposes
from copy import deepcopy
# "numpy" for linear algebra operations 
import numpy as np
# "scipy" integration tools for numerically solving ODE
from scipy.integrate import solve_ivp
# "matplotlib" library for outputting result plots
import matplotlib.pyplot as plt

### Import custom modules
# "Classes" module which contains definitions of all class objects for the simulation
import Classes
# "Functions" module which contains all function defitions for simulating the rocket's motion
from Functions import RocketDynModel

### Define class objects
# Object that stores all simulation parameters such as duration and time-step size
Sim = Classes.SimulationProperties()
# Object that stores global properties such as air and gravitational properties
Global = Classes.GlobalProperties()
# Object that stores all rocket parameters such as load mass and dimensions
Rocket = Classes.RocketProperties()

### Define the rocket's initial conditions for the simulation
StateVecInit= np.hstack((
    # The initial position and orientation are computed based on the...
    # characteristics of the launch rail
    np.array([0,0,Rocket.LaunchAlt]) + Rocket.InitDistAlongRail*Rocket.LaunchVec,
    np.array([0,-Rocket.LaunchAngle,0]),
    # The initial translational and angular velocities are assumed to be zero...
    # (i.e. the rocket is initially stationary)
    np.zeros(3),
    np.zeros(3)))

### Use SciPy initial value problem solver to run simulation
# The "RocketDynModel" function outputs the time-derivatives of the...
# rocket's states.
# "solve_ivp" then simply integrates these derivatives to compute...
# the trajectories of the states.
sol = solve_ivp(RocketDynModel,
                # Time-span of the simulation
                [0,Sim.SimTime],
                # Initial condition of the simulation
                StateVecInit,
                # Additional inputs to "RocketDynModel" that are not relevant for integration...
                # but that are required for calculations within "RocketDynModel"
                args=(Global,Rocket),
                # Define the time-steps at which to store results
                t_eval=np.arange(0,Sim.SimTime + Sim.TimeStepSize,Sim.TimeStepSize))

### Extract simulation results and store them within a class object list "Results"...
# wherein in each list entry corresponds to a unique time-step.
# For example, to access the position vector of the rocket at time-step...
# number 5, we may write Results[5].Rocket.PosVec.
Results = []
for TimeStepNum in range(Sim.NumTimeSteps):
    # From the initial value problem solution, extract the time and state vector...
    # corresponding to the current time step
    Time = sol.t[TimeStepNum]
    StateVec = sol.y[:,TimeStepNum]
    # Re-run "RocketDynModel" to compute all simulation values corresponding to...
    # the current values of time and the state vector
    RocketDynModel(Time,StateVec,Global,Rocket)
    # Since the class objects containing global and rocket properties have just been...
    # updated within "RocketDynModel", store these updated values in the "Results"...
    # class object at the list entry corresponding to the current time step number
    Results.append(Classes.TimeSeriesData(Time,deepcopy(Global),deepcopy(Rocket)))
    # Note: "deepcopy" is used here to make a copy of the global and rocket class...
    #       structures. Otherwise, future changes in these class objects will alter...
    #       the values stores in "Results"
    pass

### Generate output plots
# Create a new figure containing 6 subplots
fig, axs = plt.subplots(nrows=3,ncols=2,sharex=True)
# Plot rocket positions in 3D space
axs[0,0].plot(sol.t,sol.y[0,:]/1000) # Note: Division by 1000 to plot km instead of m
axs[1,0].plot(sol.t,sol.y[1,:]/1000)
axs[2,0].plot(sol.t,sol.y[2,:]/1000)
# Plot rocket velocities
axs[0,1].plot(sol.t,sol.y[6,:])
axs[1,1].plot(sol.t,sol.y[7,:])
axs[2,1].plot(sol.t,sol.y[8,:])
# Define shared x-axis labels
axs[2,0].set_xlabel('Time (sec)')
axs[2,1].set_xlabel('Time (sec)')
# Define y-axis labels
axs[0,0].set_ylabel('x (km)')
axs[1,0].set_ylabel('y (km)')
axs[2,0].set_ylabel('z (km)')
axs[0,1].set_ylabel('v_x (m/s)')
axs[1,1].set_ylabel('v_y (m/s)')
axs[2,1].set_ylabel('v_z (m/s)')
# Show plot
plt.show()