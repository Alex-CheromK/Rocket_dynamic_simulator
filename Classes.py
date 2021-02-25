### Module description:
# This module contains all class definitions and class object...
# attribute initilizations.

### Import standard Python libraries
# "numpy" for linear algebra operations
import numpy as np
# "pandas" library for reading values from Excel
import pandas as pd

# Define a class that stores all simulation parameters
class SimulationProperties:
    # Define (or initialize) all class attributes
    def __init__(self):
        # Simulation time step size (sec)
        self.TimeStepSize = 1
        # Simulation duration (sec)
        self.SimTime = 600
        # Number of simulation time steps
        self.NumTimeSteps = int(self.SimTime/self.TimeStepSize) + 1
    pass

# Define a class that stores all global simulation properties
class GlobalProperties:
    # Define (or initialize) all class attributes
    def __init__(self):
        # The Earth's radius (m)
        self.EarthRad = 6371.009e3
        # Gravitational acceleration at sea-level (m/s^2)
        self.GravAccelSL = 9.81
        # Wind velocity vector used in the simulation (m/s)
        self.WindVelVec = np.array([0,-10,0])
        # Air's ration of specific heat capacities
        self.SpHeatRatio = 1.4
        pass
    pass

# Define a class that stores all rocket parameters
class RocketProperties:
    # Define (or initialize) all class attributes
    def __init__(self):
        # Duration of engine thrust (sec)
        self.BurnTime = 34
        # Initial mass of the rocket with fuel (kg)
        self.LoadMass = 687.6761773
        # Fuselage diameter (m)
        self.FuseDia = 0.48
        # Fuselage length (m)
        self.FuseLength = 8.013640494
        # Flowrate of fuel out of the rocket during engine thrust (kg/s)
        self.PropFlowRate = 12
        # Engine parameters
        self.NozzleEff = 0.98
        self.C_Star = 1580.684648
        self.ExitPressure = 93914.05538 # Pa
        self.ChamberPressure = 1000000 # Pa
        self.ExpAreaRatio = 2.3
        # Reference area for aerodynamic drag calculations (m^2)
        self.DragArea = np.pi/4*self.FuseDia**2
        # Position vector from the rocket's base to its aerodynamic center (m)
        self.ACRelBasePosVec_B = np.array([0.625,0,0])
        # The rocket's direction vector in its body frame (this information...
        # is used for rotation and coordinate frame procedures)
        self.DirVec_B = np.array([1,0,0])
        # Read drag coefficient vs. Mach number data from external Excel spreadsheet...
        # and convert this data into numoy arrays
        self.DragDataTable = pd.read_excel(r'Rocket drag data.xlsx')
        self.MachNumData = self.DragDataTable['Mach number'].to_numpy()
        self.DragCoeffData = self.DragDataTable['Drag coefficient'].to_numpy()
        # Launch rail angle above the horizontal frame (deg)
        self.LaunchAngle = 80
        # The rocket's initial position along the launch rail relative to...
        # the rail's base (m)
        self.InitDistAlongRail = 0
        # The altitude of the launch rail's base above sea-level (m)
        self.LaunchAlt = 1401
        # The length of the launch rail (m)
        self.RailLength = 15
        # The friction coefficient between the launch rail and the rocket
        self.RailFricCoeff = 0
        # A unit vector denoting the orientation of the launch rail
        self.LaunchVec = np.array([
            np.cos(self.LaunchAngle*np.pi/180),
            0,
            np.sin(self.LaunchAngle*np.pi/180)])
        # Deployment altitude of the ballute (m)
        self.BalluteAlt = 75000
        # Deployment altitude of the main chute (m)
        self.MainChuteAlt = 3000
        # Drag coefficients of the ballute and main chute
        self.BalluteDragCoeff = 0.75
        self.MainChuteDragCoeff = 0.53
        # Attachment line lengths of the ballute and main chute (m)
        self.BalluteLineLength = 5
        self.MainChuteLineLength = 7.148
        # Deployed diameters of the ballute and main chute (m)
        self.BalluteDia = 1
        self.MainChuteDia = 4.13
        # Total number of active ballutes and chutes
        self.NumBallutes = 3
        self.NumMainChutes = 3
        # Drag reference areas of the ballute and main chute (m^2)
        self.BalluteArea = self.NumBallutes*np.pi/4*(self.BalluteDia**2)
        self.MainChuteArea = self.NumMainChutes*np.pi/4*(self.MainChuteDia**2)
        # Attachment locations of the ballute and main chute on the rocket fuselage...
        # relative to the rocket's base defined within the rocket's body frame (m)
        self.ChuteRelBasePosVec_B = np.array([5.66,0,0])
        pass
    pass

# Define a class that stores all time series results
class TimeSeriesData:
    # Define (or initialize) all class attributes
    def __init__(self,Time,Global,Rocket):
        # Simulation time at current time-step number
        self.Time = Time
        # Global and rocket properties at current time-step number
        self.Global = Global
        self.Rocket = Rocket
        pass
    pass