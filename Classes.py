from sys import exit
import numpy as np
import pandas as pd

class SimulationProperties:
    def __init__(self):
        self.TimeStepSize = 1
        self.SimTime = 600
        self.PlotAnimation = True
        if self.SimTime % self.TimeStepSize > 0:
            exit('In class SimulationProperties, the variable self.SimTime '
                 'must be divisble by self.TimeStepSize...')
            pass
        else:
            self.NumTimeSteps = self.SimTime // self.TimeStepSize + 1
            pass
    pass

class GlobalProperties:
    def __init__(self):
        self.EarthRad = 6371.009e3
        self.GravAccelSL = 9.81
        self.WindVelVec = np.array([0,-10,0])
        self.SpHeatRatio = 1.4
        pass
    pass

class RocketProperties:
    def __init__(self):
        self.BurnTime = 34
        self.LoadMass = 687.6761773
        self.FuseDia = 0.48
        self.FuseLength = 8.013640494
        self.PropFlowRate = 12
        self.NozzleEff = 0.98
        self.C_Star = 1580.684648
        self.ExitPressure = 93914.05538
        self.ChamberPressure = 1000000
        self.ExpAreaRatio = 2.3
        self.DragArea = np.pi/4*self.FuseDia**2
        self.ACRelBasePosVec_B = np.array([0.625,0,0])
        self.DirVec_B = np.array([1,0,0])
        self.DragDataTable = pd.read_excel(r'Rocket drag data.xlsx')
        self.MachNumData = self.DragDataTable['Mach number'].to_numpy()
        self.DragCoeffData = self.DragDataTable['Drag coefficient'].to_numpy()
        self.LaunchAngle = 80
        self.InitDistAlongRail = 0
        self.LaunchAlt = 1401
        self.RailLength = 15
        self.RailFricCoeff = 0
        self.LaunchVec = np.array([
            np.cos(self.LaunchAngle*np.pi/180),
            0,
            np.sin(self.LaunchAngle*np.pi/180)])
        self.BalluteAlt = 75000
        self.MainChuteAlt = 3000
        self.BalluteDragCoeff = 0.75
        self.MainChuteDragCoeff = 0.53
        self.BalluteLineLength = 5
        self.MainChuteLineLength = 7.148
        self.BalluteDia = 1
        self.MainChuteDia = 4.13
        self.NumBallutes = 3
        self.NumMainChutes = 3
        self.BalluteArea = self.NumBallutes*np.pi/4*(self.BalluteDia**2)
        self.MainChuteArea = self.NumMainChutes*np.pi/4*(self.MainChuteDia**2)
        self.ChuteRelBasePosVec_B = np.array([5.66,0,0])
        self.PosVec = np.zeros(3)
        self.RollAngle = 0
        self.ElevAngle = 0
        self.HeadAngle = 0
        self.VelVec = np.zeros(3)
        self.AngVelVec_B = np.zeros(3)
        pass
    pass