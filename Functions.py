### Module description:
# This module contains all custom functions necessary for simulating...
# the rocket's dynamics.

### Import standard Python libraries
# "numpy" for linear algebra operations
import numpy as np
# "pyatmos" library for estimating air properties based on...
# the standard atmosphere model
from pyatmos import coesa76

# Define a custom function that calls "coesa76" from "pyatmos"...
# but only outputs information useful for the current application
# The following function takes the rockets geomtric altitude z (m) as input
def atmos(z):
    # Call "coesda76" from "pyatmos" to compute the...
    # density (kg/m^3), temperature (K), and pressure (Pa) of air at altitude z
    rho,T,P = coesa76(z/1000)
    # Compute the speed of sound in air (m/s) at the current altitude
    a = np.sqrt(1.4*287*T)
    # Return air properties
    return rho,T,P,a

# Define a custom function that computes the rocket's mass...
# at a particular instance in time.
# The inputs to this function as the simulation time, and a class...
# object that contains all rocket properties.
def ComputeRocketMass(Time,Rocket):
    # If the engine combustion process is still active and the rocket...
    # is burning fuel...
    if Time <= Rocket.BurnTime:
        # The mass of the rocket is equal to the initial load mass...
        # minues the mass of propellant that has flowed out of the engine...
        # up until the current instant in time
        Rocket.Mass = Rocket.LoadMass - Rocket.PropFlowRate*Time
        pass
    # If the combustion process has ended...
    else:
        # The rocket mass remains at the value it possessed once the engine...
        # operation ceased
        Rocket.Mass = Rocket.LoadMass - Rocket.PropFlowRate*Rocket.BurnTime
        pass
    pass

# Define a custom function that computes the location of the...
# rocket's center-of-gravity relative to its base at a given instance in time.
# This function requires as inputs the current time and a class object conatining...
# all rocket properties.
def ComputeCGRelBase(Time,Rocket):
    # If the engine combustion process is still active and the rocket...
    # is burning fuel...
    if Time <= Rocket.BurnTime:
        # The following expressions have been derived for the specific rocket...
        # being simulated at UBC.
        # Other users will have to generate their own expressions for the rocket's...
        # center-of-mass location relative to the rocket base over time.
        Rocket.CGRelBasePosVec_B = np.array([(11313424894001152*Time**2)/18632593421892121875 - (187339189186985984*Time)/2661799060270303125 + 45613417578875388465333449617306391/16782748158355863008585176842240000,
                                             0,
                                             0])
        # Based on the above positions, the velocity and acceelration of the...
        # center-of-gravity relatvie to the rocket's based may be computed using differentiation.
        Rocket.CGRelBaseVelVec_B = np.array([(22626849788002304*Time)/18632593421892121875 - 187339189186985984/2661799060270303125,
                                             0,
                                             0])
        Rocket.CGRelBaseAccVec_B = np.array([22626849788002304/18632593421892121875,
                                             0,
                                             0])
        pass
    # If the combustion process has ended...
    else:
        # The center-of-gravity location relative to the rocket base...
        # remains fixed at the value corresponding...
        # to the instant at which combustion ended.
        Rocket.CGRelBasePosVec_B = np.array([(11313424894001152*Rocket.BurnTime**2)/18632593421892121875 - (187339189186985984*Rocket.BurnTime)/2661799060270303125 + 45613417578875388465333449617306391/16782748158355863008585176842240000,
                                             0,
                                             0])
        # The center-of-gravity no longer moves relative to the rocket base,...
        # therefore its relative velocity and acceleration are zero.
        Rocket.CGRelBaseVelVec_B = np.zeros(3)
        Rocket.CGRelBaseAccVec_B = np.zeros(3)
        pass
    pass

# Define a custom function that computes the magnitude of the...
# rocket's thrust force.
# This function depends on the simulation time and class objects...
# containing global and rocket properties.
def ComputeThrustMag(Time,Global,Rocket):
    # If the engine combustion process is still active and the rocket...
    # is burning fuel...
    if Time <= Rocket.BurnTime:
        # The following equations are based on standard rocket theory,...
        # and may be obtained from "Chapter 3" of "Rocket and Spacecraft Propulsion"...
        # by Martin J. L. Turner.
        Rocket.SpImpulse = \
            Rocket.NozzleEff*(Rocket.C_Star*Global.SpHeatRatio/Global.GravAccelSL*np.sqrt((2/(Global.SpHeatRatio - 1)) \
            *((2/(Global.SpHeatRatio + 1))**((Global.SpHeatRatio + 1)/(Global.SpHeatRatio - 1))) \
            *(1 - ((Rocket.ExitPressure/Rocket.ChamberPressure)**((Global.SpHeatRatio - 1)/Global.SpHeatRatio)))) \
            + Rocket.C_Star*Rocket.ExpAreaRatio/(Global.GravAccelSL*Rocket.ChamberPressure) \
            *(Rocket.ExitPressure - Global.AtmPressure))
        Rocket.ThrustForce = Rocket.SpImpulse*Rocket.PropFlowRate*Global.GravAccelSL
        pass
    # If the combustion process has ended, there is zero thrust force...
    else:
        Rocket.ThrustForce = 0
        pass
    pass

# Define a custom function that computes the drag force vector...
# acting on the rocket's aerodynamics center.
# This function requires only global and rocket properties.
def ComputeDragForceVec(Global,Rocket):
    # Compute the wind velocity experienced by an observed who is...
    # fixed to the rocket's aerodynamics center.
    Rocket.RelWindVelVec = Global.WindVelVec - Rocket.ACVelVec
    # Compute the Mach number experienced by this observer.
    Rocket.MachNum = np.linalg.norm(Rocket.RelWindVelVec)/Global.SoundSpeed
    # Using provided drag coefficient vs. Mach number data, compute...
    # the drag coefficient of the rocket.
    Rocket.DragCoeff = np.interp(Rocket.MachNum,Rocket.MachNumData,Rocket.DragCoeffData)
    # Note: The available drag coefficient vs. Mach number data is intended for...
    #       flow that is directly parallel to the rocket's fuselage.
    #       Throughout this simulation, the incoming relative flow will clearly alter...
    #       in direction and gain a component that is perpendicular to the rocket's fuselage.
    #       However, given the high speed of the rocket relative to any wind, this perpendicular...
    #       will remain small, which makes the current simplified approach for computig drag valid.

    # With the drag coefficient know, the drag force vector may be...
    # computed using Morison's equation.
    Rocket.DragForceVec = \
        0.5*Rocket.DragCoeff*Global.AirDensity*Rocket.DragArea*np.linalg.norm(Rocket.RelWindVelVec)*Rocket.RelWindVelVec

    pass

# Finally define a custom function which takes in the simulation time,...
# the rocket's state vector, and other global and rocket properties, and...
# outputs the time-derivative of the rocket's state vector for numerical...
# integration to solve the ODE.
def RocketDynModel(Time,StateVec,Global,Rocket):
    # Assign the state vector variables to uniqe rocket states.
    Rocket.PosVec = StateVec[0:3]
    Rocket.RollAngle = StateVec[3]
    Rocket.ElevAngle = StateVec[4]
    Rocket.HeadAngle = StateVec[5]
    Rocket.VelVec = StateVec[6:9]
    Rocket.AngVelVec_B = StateVec[9:12]
    Rocket.AngVelVec_B[0] = 0

    ### Compute the rocket's inertial properties.
    # Mass
    ComputeRocketMass(Time,Rocket)
    # Location of its CG relative to its base
    ComputeCGRelBase(Time,Rocket)
    Rocket.MOIMat = np.array([[(1/2)*Rocket.Mass*(Rocket.FuseDia/2)**2,0,0],
                              [0,(1/12)*Rocket.Mass*(3*(Rocket.FuseDia/2)**2 + Rocket.FuseLength**2),0],
                              [0,0,(1/12)*Rocket.Mass*(3*(Rocket.FuseDia/2)**2 + Rocket.FuseLength**2)]])

    # Compute global properties based on the rocket's gemoetric altitude
    Global.AirDensity,_,Global.AtmPressure,Global.SoundSpeed = atmos(Rocket.PosVec[2])
    Global.GravAccel = Global.GravAccelSL*((Global.EarthRad/(Global.EarthRad + Rocket.PosVec[2]))**2)

    ### Kinematics - Rotations, velocities, and accelerations of the rocket
    # Define rotation matrices that compute the rocket's orientation based...
    # on it current Euler angles
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

    # Define a roation matrix that converts the time-derivatives of...
    # the rocket's Euler angles to angular velocites in the body frame.
    Rocket.EulerRotMat = np.array([[1,0,-np.sin(Rocket.ElevAngle*np.pi/180)],
                                   [0,np.cos(Rocket.RollAngle*np.pi/180),np.cos(Rocket.ElevAngle*np.pi/180)*np.sin(Rocket.RollAngle*np.pi/180)],
                                   [0,-np.sin(Rocket.RollAngle*np.pi/180),np.cos(Rocket.ElevAngle*np.pi/180)*np.cos(Rocket.RollAngle*np.pi/180)]])

    # Define reference frame transformation matrices going...
    # "forward" from the global frame to the body frame.
    Rocket.HeadTransMat_F = np.transpose(Rocket.HeadRotMat)
    Rocket.ElevTransMat_F = np.transpose(Rocket.ElevRotMat)
    Rocket.RollTransMat_F = np.transpose(Rocket.RollRotMat)
    Rocket.TotalTransMat_F = Rocket.RollTransMat_F @ Rocket.ElevTransMat_F @ Rocket.HeadTransMat_F

    # Define reference frame transformation matrices going...
    # "backward" from the body frame to the global frame.
    Rocket.HeadTransMat_B = Rocket.HeadRotMat
    Rocket.ElevTransMat_B = Rocket.ElevRotMat
    Rocket.RollTransMat_B = Rocket.RollRotMat
    Rocket.TotalTransMat_B = Rocket.HeadTransMat_B @ Rocket.ElevTransMat_B @ Rocket.RollTransMat_B

    # Define a unit vector which is parallel to the rocket's fuselage...
    # at the current instant in time.
    # This vector is ncessary for projecting forces along the rocket's body.
    Rocket.DirVec = Rocket.TotalRotMat @ Rocket.DirVec_B

    # Compute the velocity of the rocket's aerodynamic center
    Rocket.ACVelVec = \
        Rocket.VelVec + \
        Rocket.TotalTransMat_B @ \
        (-Rocket.CGRelBaseVelVec_B + np.cross(Rocket.AngVelVec_B,Rocket.ACRelBasePosVec_B - Rocket.CGRelBasePosVec_B))
    # Compute the velocity of the location along the rocket where...
    # parachutes are attached.
    Rocket.ChuteVelVec = \
        Rocket.VelVec + \
        Rocket.TotalTransMat_B @ \
        (Rocket.CGRelBaseVelVec_B + np.cross(Rocket.AngVelVec_B,Rocket.ChuteRelBasePosVec_B - Rocket.CGRelBasePosVec_B))

    ### Kinetics - Forces and moments that act on the rocket
    # First, we must check to see of the rocket has detached from the launch rail.
    # If so, the equations of motion contain only one DOF aligned with the rail.

    # If the rocket is still moving along the rail...
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
 
        # Compute the rocket's translational and angular accelerations
        Rocket.AccVec = Rocket.TotalForce/Rocket.Mass*Rocket.LaunchVec
        # Note: Multiplying the acceleration by "Rocket.LaunchVec" converts the scalar...
        #       acceleration value projected along the rail to a 3D vector in the global frame.

        # While attached to the launch rail, the rocket cannot rotate
        Rocket.AngAccVec_B = np.zeros(3)
        pass

    # If the rocket has detached from the launch rail...
    else:
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
 
        # Compute the rocket's translational and angular accelerations
        Rocket.AccVec = Rocket.ForceVec/Rocket.Mass
        Rocket.AngAccVec_B = np.linalg.inv(Rocket.MOIMat) @ \
            (Rocket.MomVec - np.cross(Rocket.AngVelVec_B,Rocket.MOIMat @ Rocket.AngVelVec_B))
        pass

    # Compile and return the state derivative vector
    StateDerivVec = np.hstack((
        Rocket.VelVec,
        np.linalg.inv(Rocket.EulerRotMat) @ Rocket.AngVelVec_B*180/np.pi,
        Rocket.AccVec,
        Rocket.AngAccVec_B))
    return StateDerivVec