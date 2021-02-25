import os
import numpy as np
from pyatmos import coesa76

def atmos(z):
    rho,T,P = coesa76(z/1000)
    R = 287
    k = 1.4
    a = np.sqrt(k*R*T)
    return rho,T,P,a

def ComputeRocketMass(Time,Rocket):
    if Time <= Rocket.BurnTime:
        Rocket.Mass = Rocket.LoadMass - Rocket.PropFlowRate*Time
        pass
    else:
        Rocket.Mass = Rocket.LoadMass - Rocket.PropFlowRate*Rocket.BurnTime
        pass
    pass

def ComputeCGRelBase(Time,Rocket):
    if Time <= Rocket.BurnTime:
        Rocket.CGRelBasePosVec_B = np.array([(11313424894001152*Time**2)/18632593421892121875 - (187339189186985984*Time)/2661799060270303125 + 45613417578875388465333449617306391/16782748158355863008585176842240000,
                                             0,
                                             0])
        Rocket.CGRelBaseVelVec_B = np.array([(22626849788002304*Time)/18632593421892121875 - 187339189186985984/2661799060270303125,
                                             0,
                                             0])
        Rocket.CGRelBaseAccVec_B = np.array([22626849788002304/18632593421892121875,
                                             0,
                                             0])
        pass
    else:
        Rocket.CGRelBasePosVec_B = np.array([(11313424894001152*Rocket.BurnTime**2)/18632593421892121875 - (187339189186985984*Rocket.BurnTime)/2661799060270303125 + 45613417578875388465333449617306391/16782748158355863008585176842240000,
                                             0,
                                             0])
        Rocket.CGRelBaseVelVec_B = np.zeros(3)
        Rocket.CGRelBaseAccVec_B = np.zeros(3)
        pass
    pass

def ComputeThrustMag(Time,Global,Rocket):
    if Time <= Rocket.BurnTime:
        Rocket.SpImpulse = \
            Rocket.NozzleEff*(Rocket.C_Star*Global.SpHeatRatio/Global.GravAccelSL*np.sqrt((2/(Global.SpHeatRatio - 1)) \
            *((2/(Global.SpHeatRatio + 1))**((Global.SpHeatRatio + 1)/(Global.SpHeatRatio - 1))) \
            *(1 - ((Rocket.ExitPressure/Rocket.ChamberPressure)**((Global.SpHeatRatio - 1)/Global.SpHeatRatio)))) \
            + Rocket.C_Star*Rocket.ExpAreaRatio/(Global.GravAccelSL*Rocket.ChamberPressure) \
            *(Rocket.ExitPressure - Global.AtmPressure))
        Rocket.ThrustForce = Rocket.SpImpulse*Rocket.PropFlowRate*Global.GravAccelSL
        pass
    else:
        Rocket.ThrustForce = 0
        pass
    pass

def ComputeDragForceVec(Global,Rocket):
    Rocket.RelWindVelVec = Global.WindVelVec - Rocket.ACVelVec
    Rocket.MachNum = np.linalg.norm(Rocket.RelWindVelVec)/Global.SoundSpeed
    Rocket.DragCoeff = np.interp(Rocket.MachNum,Rocket.MachNumData,Rocket.DragCoeffData)
    Rocket.DragForceVec = \
        0.5*Rocket.DragCoeff*Global.AirDensity*Rocket.DragArea*np.linalg.norm(Rocket.RelWindVelVec)*Rocket.RelWindVelVec
    pass

def RocketDynModel(Time,StateVec,Global,Rocket):
    # Extract states
    Rocket.PosVec = StateVec[0:3]
    Rocket.RollAngle = StateVec[3]
    Rocket.ElevAngle = StateVec[4]
    Rocket.HeadAngle = StateVec[5]
    Rocket.VelVec = StateVec[6:9]
    Rocket.AngVelVec_B = StateVec[9:12]

    # Rocket inertial properties
    ComputeRocketMass(Time,Rocket)
    ComputeCGRelBase(Time,Rocket)
    Rocket.MOIMat = np.array([[(1/2)*Rocket.Mass*(Rocket.FuseDia/2)**2,0,0],
                              [0,(1/12)*Rocket.Mass*(3*(Rocket.FuseDia/2)**2 + Rocket.FuseLength**2),0],
                              [0,0,(1/12)*Rocket.Mass*(3*(Rocket.FuseDia/2)**2 + Rocket.FuseLength**2)]])

    # Global properties
    Global.AirDensity,_,Global.AtmPressure,Global.SoundSpeed = atmos(Rocket.PosVec[2])
    Global.GravAccel = Global.GravAccelSL*((Global.EarthRad/(Global.EarthRad + Rocket.PosVec[2]))**2)

    # Kinematics
    Rocket.HeadRotMat = np.array([[np.cos(Rocket.HeadAngle*np.pi/180),-np.sin(Rocket.HeadAngle*np.pi/180),0],
                                  [np.sin(Rocket.HeadAngle*np.pi/180),np.cos(Rocket.HeadAngle*np.pi/180),0],
                                  [0,0,1]])
    Rocket.ElevRotMat = np.array([[np.cos(Rocket.ElevAngle*np.pi/180),0,np.sin(Rocket.ElevAngle*np.pi/180)],
                                  [0,1,0],
                                  [-np.sin(Rocket.ElevAngle*np.pi/180),0,np.cos(Rocket.ElevAngle*np.pi/180)]])
    Rocket.RollRotMat = np.array([[1,0,0],
                                  [0,np.cos(Rocket.RollAngle*np.pi/180),-np.sin(Rocket.RollAngle*np.pi/180)],
                                  [0,np.sin(Rocket.RollAngle*np.pi/180),np.cos(Rocket.RollAngle*np.pi/180)]])
    Rocket.TotalRotMat = Rocket.HeadRotMat @ Rocket.ElevRotMat @ Rocket.RollRotMat

    Rocket.EulerRotMat = np.array([[1,0,-np.sin(Rocket.ElevAngle*np.pi/180)],
                                   [0,np.cos(Rocket.RollAngle*np.pi/180),np.cos(Rocket.ElevAngle*np.pi/180)*np.sin(Rocket.RollAngle*np.pi/180)],
                                   [0,-np.sin(Rocket.RollAngle*np.pi/180),np.cos(Rocket.ElevAngle*np.pi/180)*np.cos(Rocket.RollAngle*np.pi/180)]])

    Rocket.HeadTransMat_F = np.transpose(Rocket.HeadRotMat)
    Rocket.ElevTransMat_F = np.transpose(Rocket.ElevRotMat)
    Rocket.RollTransMat_F = np.transpose(Rocket.RollRotMat)
    Rocket.TotalTransMat_F = Rocket.RollTransMat_F @ Rocket.ElevTransMat_F @ Rocket.HeadTransMat_F

    Rocket.HeadTransMat_B = Rocket.HeadRotMat
    Rocket.ElevTransMat_B = Rocket.ElevRotMat
    Rocket.RollTransMat_B = Rocket.RollRotMat
    Rocket.TotalTransMat_B = Rocket.HeadTransMat_B @ Rocket.ElevTransMat_B @ Rocket.RollTransMat_B

    Rocket.DirVec = Rocket.TotalRotMat @ Rocket.DirVec_B

    Rocket.ACVelVec = \
        Rocket.VelVec + \
        Rocket.TotalTransMat_B @ \
        (-Rocket.CGRelBaseVelVec_B + np.cross(Rocket.AngVelVec_B,Rocket.ACRelBasePosVec_B - Rocket.CGRelBasePosVec_B))
    Rocket.ChuteVelVec = \
        Rocket.VelVec + \
        Rocket.TotalTransMat_B @ \
        (Rocket.CGRelBaseVelVec_B + np.cross(Rocket.AngVelVec_B,Rocket.ChuteRelBasePosVec_B - Rocket.CGRelBasePosVec_B))

    # Kinetics
    if Rocket.PosVec[2] < Rocket.RailLength*np.sin(Rocket.LaunchAngle*np.pi/180) + Rocket.LaunchAlt:
        ## Compute external forces
        # Thrust force
        ComputeThrustMag(Time,Global,Rocket)
 
        # Gravitational force
        Rocket.GravForce = -Rocket.Mass*Global.GravAccel*np.sin(Rocket.LaunchAngle*np.pi/180)
 
        # Rail friction force
        Rocket.FricForce = -Rocket.RailFricCoeff*Rocket.Mass*Global.GravAccel*np.cos(Rocket.LaunchAngle*np.pi/180)
 
        # Aerodynamic force
        ComputeDragForceVec(Global,Rocket)
        Rocket.AeroForce = np.dot(Rocket.DragForceVec,Rocket.LaunchVec)
 
        # Net external force
        Rocket.TotalForce = Rocket.ThrustForce + Rocket.GravForce + Rocket.FricForce + Rocket.AeroForce
 
        ## Compute state derivatives
        Rocket.AccVec = Rocket.TotalForce/Rocket.Mass*Rocket.LaunchVec
        Rocket.AngAccVec_B = np.zeros(3)
        pass
    else:
        ## Kinetics - Forces
        # Gravitational force
        Rocket.GravForceVec = np.array([0,0,-Rocket.Mass*Global.GravAccel])
 
        # Thrust force
        ComputeThrustMag(Time,Global,Rocket)
        Rocket.ThrustForceVec = Rocket.ThrustForce*Rocket.DirVec
 
        # Aerodynamic force
        ComputeDragForceVec(Global,Rocket)
        Rocket.AeroForceVec = Rocket.DragForceVec
 
        # Ballute force (currently using massless parachute model)
        Rocket.ChuteRelWindVelVec = Global.WindVelVec - Rocket.ChuteVelVec
        if (Rocket.VelVec[2] < 0) and (Rocket.PosVec[2] < Rocket.BalluteAlt) and (Rocket.PosVec[2] > Rocket.MainChuteAlt):
            Rocket.BalluteForceVec = \
                0.5*Rocket.BalluteDragCoeff*Global.AirDensity*Rocket.BalluteArea \
                *np.linalg.norm(Rocket.ChuteRelWindVelVec)*Rocket.ChuteRelWindVelVec
            pass
        else:
            Rocket.BalluteForceVec = np.zeros(3)
            pass
 
        # Main chute force (currently using massless parachute model)
        if (Rocket.VelVec[2] < 0) and (Rocket.PosVec[2] < Rocket.MainChuteAlt):
            Rocket.MainChuteForceVec = \
                0.5*Rocket.MainChuteDragCoeff*Global.AirDensity*Rocket.MainChuteArea \
                *np.linalg.norm(Rocket.ChuteRelWindVelVec)*Rocket.ChuteRelWindVelVec
            pass
        else:
            Rocket.MainChuteForceVec = np.zeros(3)
            pass
        # Total force vector
        Rocket.ForceVec = \
            Rocket.GravForceVec + Rocket.ThrustForceVec + Rocket.AeroForceVec \
            + Rocket.BalluteForceVec + Rocket.MainChuteForceVec
 
        ## Kinetics - Moments
        # Thrust damping moment
        if Time < Rocket.BurnTime:
            Rocket.ThrustDampMom = \
                -Rocket.PropFlowRate*np.cross(-Rocket.CGRelBasePosVec_B,np.cross(Rocket.AngVelVec_B,-Rocket.CGRelBasePosVec_B))
            pass
        else:
            Rocket.ThrustDampMom = np.zeros(3)
            pass

        # Fuselage aerodynamic moment
        Rocket.AeroMomVec = \
            np.cross(Rocket.ACRelBasePosVec_B - Rocket.CGRelBasePosVec_B,Rocket.TotalTransMat_F @ Rocket.AeroForceVec)
 
        # Ballute moment
        Rocket.BalluteMomVec = \
            np.cross(Rocket.ChuteRelBasePosVec_B - Rocket.CGRelBasePosVec_B,Rocket.TotalTransMat_F @ Rocket.BalluteForceVec)
 
        # Main chute moment
        Rocket.MainChuteMomVec = \
            np.cross(Rocket.ChuteRelBasePosVec_B - Rocket.CGRelBasePosVec_B,Rocket.TotalTransMat_F @ Rocket.MainChuteForceVec)
 
        # Total aerodynamic moment
        Rocket.MomVec = \
            Rocket.AeroMomVec + Rocket.BalluteMomVec + Rocket.MainChuteMomVec + Rocket.ThrustDampMom
 
        ## Compute state derivatives
        Rocket.AccVec = Rocket.ForceVec/Rocket.Mass
        Rocket.AngAccVec_B = np.linalg.inv(Rocket.MOIMat) @ \
            (Rocket.MomVec - np.cross(Rocket.AngVelVec_B,Rocket.MOIMat @ Rocket.AngVelVec_B))
        pass



    ## Compile state derivative vector
    StateDerivVec = np.hstack((
        Rocket.VelVec,
        np.linalg.inv(Rocket.EulerRotMat) @ Rocket.AngVelVec_B,
        Rocket.AccVec,
        Rocket.AngAccVec_B))
    return StateDerivVec