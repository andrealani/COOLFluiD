#!/usr/bin/env python
import datetime
import math
import numpy as np

### write here the parameters that will be needed to writey the CFcase file:
author = "N. Ozak"  # creator of the CFcase file
now = datetime.datetime.now()  # will print when the CFcase file was created
purpose = "Electron + ion fluid model of the circular polarized wave"  # write the goal of the CFcase file
working_dir = "./"  # working directory
CFCasefile = "circularPolarized.CFcase" #name of CFcase file
results_dir = "./Results_two_fluid_wave"  # name of folder for outputs
inter_file = "./TwoFluid_wave.inter"  # name of interactive oarameter file to be read
read_rate = 10  # Rate at which you want the inter file to be read (after this many time steps)
TwoDHalf = "true"  # for 2.5D case --> true for 2D --> false  for 3D???
nbSpecies = 2  # Number of species to be included in the model
gamma = 5. / 3.  # adiabatic constant of both fluids
divBclean = 1.0  # Include divB cleaning? yes --> 1.0
divEclean = 1.0  # Include divE cleaning? yes --> 1.0
Limiter = True
restart = False  # True # If you want to restart the simulation set this to "True" and set up the start_time which is the time 0 when the simulation will start AND the mesh that should be read initially
Mesh_file = "./Mesh_2x1_128x64_centered_periodic.CFmesh"  # Mesh that you want to use in your run

if Mesh_file == "./Mesh_2x1_128x64_centered_periodic.CFmesh":
    mesh_save_rate = 1000  # Save the mesh after these many time steps
    data_save_rate = 100  # Save the data to a .plt file in results_dir after these many time steps
    MaxSteps = 5000  # Maximum number of time steps to run in the simulation
elif Mesh_file == "./Mesh_2x1_64x32_centered_periodic.CFmesh":
    mesh_save_rate = 500  # Save the mesh after these many time steps
    data_save_rate = 50  # Save the data to a .plt file in results_dir after these many time steps
    MaxSteps = 2500  # Maximum number of time steps to run in the simulation
elif Mesh_file == "./Mesh_2x1_32x16_centered_periodic.CFmesh":
    mesh_save_rate = 250  # Save the mesh after these many time steps
    data_save_rate = 25  # Save the data to a .plt file in results_dir after these many time steps
    MaxSteps = 1250  # Maximum number of time steps to run in the simulation
elif Mesh_file == "./Mesh_2x1_16x8_centered_periodic.CFmesh":
    mesh_save_rate = 125  # Save the mesh after these many time steps
    data_save_rate = 25  # Save the data to a .plt file in results_dir after these many time steps
    MaxSteps = 625  # Maximum number of time steps to run in the simulation
else:
    print ("WRONG MESH NAME")

# Scale_factor = 1.0849417947911578e-06  # Scaling factor for the mesh
if restart:
    Scale_factor = 1.0
    start_time = 0.0018625737
Source_term = "DriftWaves2DHalfTwoFluid"
Collisional = "false"
Norm_param = -6  # Up to what decimal do you want the convergence to go to, this is in log scale

# Global constants:
mu_0 = 1.2566370614e-6  # pc.mu_0
lightSpeed = 2.99792485e8
epsilon_0 = 8.854187817e-12
k_B = 1.38064852e-23  # pc.Boltzmann
charge = 1.60217662e-19  # Electric charge to be used
m_p = 1.67262177774e-27  # pc.m_p

# Chosen conditions for normalization:
massRatio = 256  # mi/me
B_0 = 6.e-9  # [T] normalization magnetic field
n_0 = 1.e7  # [m^-3]
beta_0 = 0.1 # Plasma beta
eta_0  = 0.1 # perturbation amplitude
#omega_norm = 5.72
kappa_norm =  np.sqrt(5.)#np.sqrt(0.5**2 + 1)
mass2 = m_p  # 1.67262177774e-27 # Mass of second species (ions)
mass1 = m_p / massRatio  # .1093821545e-31 # Mass of first species (electrons)


def norm_vars(bfield, density):
    omega_pi = np.sqrt(density * charge ** 2. / (epsilon_0 * m_p))  # ion plasma frequency wpi ^-1 is unit time scale
    rho_norm = m_p * density  # [kg/m3] norm mass density
    v_norm = bfield / np.sqrt(mu_0 * rho_norm)  # [m/s] = V_Alfven
    c_norm = 10 * v_norm  # numerical speed of light, set to 10*V_A
    l_norm = lightSpeed / omega_pi  # c_num/omega_pi == ion skin depth is unit length scale
    v_norm = bfield / np.sqrt(rho_norm * mu_0)  # Alfven speed of ions
    p_norm = bfield ** 2. / (mu_0)  # magnetic pressure
    t_norm = l_norm / v_norm  # normalized time from norm dist and vel
    E_norm = v_norm * bfield  # normalized electric field from v*B
    T_norm = p_norm / (k_B * density)  # normalized temperature obtained from norm pressure
    return rho_norm, v_norm, l_norm, c_norm, t_norm, E_norm, T_norm


rho_0, v_0, l_0, c_num, t_0, E_0, T_0 = norm_vars(B_0, n_0)

# Reference Values:
epsilon_num = 1./(mu_0*c_num**2.)
omega_ci = charge*B_0/mass2
omega_ce = -charge*B_0/mass1
omega_pi = np.sqrt(n_0 * charge ** 2. / (epsilon_num * mass2))  # ion plasma frequency wpi ^-1 is unit time scale
omega_pe = np.sqrt(n_0 * charge ** 2. / (epsilon_num * mass1))  # ion plasma frequency wpi ^-1 is unit time scale


def dispersion_relation_4th(kappa, omega_ce, omega_ci, omega_pe, omega_pi):
    Coeff_e = (omega_pe/(kappa*c_num))**2.
    Coeff_i = (omega_pi/(kappa*c_num))**2.
    A4 = -1./(c_num**2.*kappa**2)
    A3 = -(omega_ce + omega_ci)/(c_num**2.*kappa**2)
    A2 = (Coeff_i + Coeff_e - (omega_ce*omega_ci)/(c_num**2.*kappa**2) + 1)
    A1 = (omega_ce*(1 + Coeff_i) + omega_ci*(1 + Coeff_e))
    A0 = omega_ce*omega_ci

    coeffs = np.array([A4, A3, A2, A1, A0])

    solution  = np.roots(coeffs)

    return solution

def dispersion_relation_2nd(kappa, omega_ce, omega_ci, omega_pe, omega_pi):
    Coeff_e = (omega_pe/(kappa*c_num))**2.
    Coeff_i = (omega_pi/(kappa*c_num))**2.
    A2 = (Coeff_i + Coeff_e + 1)
    A1 = (omega_ce*(1 + Coeff_i) + omega_ci*(1 + Coeff_e))
    A0 = omega_ce*omega_ci

    coeffs = np.array([A2, A1, A0])

    solution  = np.roots(coeffs)

    return solution

#omega    = omega_norm*omega_ci
kappa    = kappa_norm/l_0
solution = dispersion_relation_4th(kappa, omega_ce, omega_ci, omega_pe, omega_pi)
omega = abs(solution[2])

print ('kappa_norm = ', kappa*l_0)
print ('solution0  = ', solution[0]/omega_ci)
print ('solution1  = ', solution[1]/omega_ci)
print ('solution2  = ', solution[2]/omega_ci)
print ('solution3  = ', solution[3]/omega_ci)
print ('l_0 = ', l_0)


V_e      = abs(eta_0*omega_ce/(omega + omega_ce)*omega/kappa)
V_i      = abs(eta_0*omega_ci/(omega + omega_ci)*omega/kappa)

# Scale factor and length of the domain
L_x = np.pi * l_0  # length in x of the rectangular domain
L_y = L_x/2. #4.*np.pi * l_0  # length in y of the rectangular domain

if not restart:
    Scale_factor = 1./L_x

Bx_ref = B_0
By_ref = eta_0*B_0
Bz_ref = eta_0*B_0
Ex_ref = abs(omega/kappa*eta_0*B_0)
Ey_ref = abs(omega/kappa*eta_0*B_0)
Ez_ref = abs(omega/kappa*eta_0*B_0)
Psi_ref = E_0
Phi_ref = B_0
rhoe_ref = mass1 * n_0
rhoi_ref = mass2 * n_0
ue_ref = V_e
ve_ref = V_e
we_ref = V_e
ui_ref = V_i
vi_ref = V_i
wi_ref = V_i
Te_ref = beta_0/4. * T_0
Ti_ref = beta_0/4. * T_0
c_se   = np.sqrt(k_B*Te_ref/mass1) #Sound speed of electrons

# TimeStepSettings
delta_X = L_x/128
# Time Step settings
TimeStep = abs(2*np.pi/(1000*omega)) # Time Step defined for the case of 128
CFL = c_num*TimeStep/delta_X # CFL defined for the 128 point case

if Mesh_file == "./Mesh_2x1_128x64_centered_periodic.CFmesh":
    TimeStep = TimeStep
elif Mesh_file == "./Mesh_2x1_64x32_centered_periodic.CFmesh":
    TimeStep = TimeStep*2. # We keep the CFL constant for all simulations
elif Mesh_file == "./Mesh_2x1_32x16_centered_periodic.CFmesh":
    TimeStep = TimeStep*4. # We keep the CFL constant for all simulations
elif Mesh_file == "./Mesh_2x1_16x8_centered_periodic.CFmesh":
    TimeStep = TimeStep*8. # We keep the CFL constant for all simulations

print ("CFL = ", CFL)
print ("TimeStep = ", TimeStep)
print ("omega    = ", omega)


f = open(CFCasefile, 'w')
f.write("# COOLFluiD Startfile\n")
f.write("# Comments begin with '#' \n")
f.write("# Test file created by: " + author + ", on: " + now.strftime("%Y-%m-%d") + "\n")
f.write("# Info: " + purpose + "\n")
f.write("################################################################################ \n")
f.write("# Assertion For Debugging\n")
f.write("#\n")
f.write("#CFEnv.ExceptionLogLevel    = 1000\n")
f.write("#CFEnv.DoAssertions         = true\n")
f.write("#CFEnv.AssertionDumps       = true\n")
f.write("#CFEnv.AssertionThrows      = true\n")
f.write("#CFEnv.AssertThrows         = true\n")
f.write("#CFEnv.AssertDumps          = true\n")
f.write("#CFEnv.ExceptionAborts      = true\n")
f.write("#CFEnv.ExceptionDumps       = true\n")
f.write("#CFEnv.ExceptionOutputs     = true\n")
f.write("#CFEnv.RegistSignalHandlers = false\n")
f.write("#CFEnv.TraceToStdOut        = true\n")
f.write("#CFEnv.TraceActive          = true\n")
f.write("#\n")
f.write("# this will always fail with GAMBIT\n")
f.write("#CFEnv.ErrorOnUnusedConfig = true\n")
f.write("#\n")
f.write("################################################################################ \n")
f.write("#\n")
f.write("# SubSystem Modules\n")
f.write(
    "Simulator.Modules.Libs = libCFmeshFileReader libCFmeshFileWriter libTecplotWriter libNavierStokes libMaxwell libMultiFluidMHD libFiniteVolume libNewtonMethod libFiniteVolumeNavierStokes libFiniteVolumeMaxwell libFiniteVolumeMultiFluidMHD libGambit2CFmesh libForwardEuler libPetscI \n")
f.write("\n")
f.write("#SubSystem Parameters\n")
f.write("Simulator.Paths.WorkingDir       = " + working_dir + "\n")
f.write("Simulator.Paths.ResultsDir       = " + results_dir + "\n")
f.write("\n")
f.write("#SubSystem Parameters\n")
f.write("Simulator.SubSystem.InteractiveParamReader.FileName = " + inter_file + "\n")
f.write("Simulator.SubSystem.InteractiveParamReader.readRate = " + str(read_rate) + "\n")
f.write("#\n")
f.write("################################################################################ \n")
f.write("#\n")
f.write("# Physical Model\n")
f.write("# order of the reference values is Bx, By, Bz, Ex, Ey, Ez, Psi, Phi, rhoe, rhoi, ue, ve, we, ui, vi, wi, Te, Ti\n")
f.write("\n")
f.write("Simulator.SubSystem.Default.PhysicalModelType = MultiFluidMHD2D #You choose between 2D and 3D\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.refValues = " + str(Bx_ref) + " " + str(By_ref) + " " + str(Bz_ref) + " " + str(Ex_ref) + " " + str(
    Ey_ref) + " " + \
        str(Ez_ref) + " " + str(Psi_ref) + " " + str(Phi_ref) + " " + str(rhoe_ref) + " " + str(rhoi_ref) + " " + str(ue_ref) + " " + str(
    ve_ref) + " " + str(we_ref) + " " + \
        str(ui_ref) + " " + str(vi_ref) + " " + str(wi_ref) + " " + str(Te_ref) + " " + str(
    Ti_ref) + " " + "# Characteristic value of the solution for the computation of the numerical Jacobian of the convective Flux\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.refLength = 1.0 # Not changed\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.Is2DHalf  = " + TwoDHalf + " # true for 2.5D\n")
f.write("\n")
f.write("# Fluids Properties\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.nbSpecies = " + str(nbSpecies) + "\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.ConvTerm.molecularMass1 = " + str(mass1) + " # Electron mass\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.ConvTerm.molecularMass2 = " + str(mass2) + " # Ion Mass\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.ConvTerm.gamma          = " + str(gamma) + "\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.ConvTerm.nonInducedElectromagnetic = 0.0 0.0 0.0 0.0 0.0 0.0 # Bx By Bz Ex Ey Ez \n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.ConvTerm.divBCleaningConst  = " + str(divBclean) + "\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.ConvTerm.divECleaningConst  = " + str(divEclean) + "\n")
f.write("#Simulator.SubSystem.MultiFluidMHD2D.DiffTerm.nbSpecies         = 2\n")
f.write("#Simulator.SubSystem.MultiFluidMHD2D.DiffTerm.dynViscosity      = 0. 0.\n")
f.write("#Simulator.SubSystem.MultiFluidMHD2D.DiffTerm.thermConductivity = 0. 0. \n")
f.write("#Simulator.SubSystem.MultiFluidMHD2D.DiffTerm.BraginskiiTransport = false \n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.ConvTerm.lightSpeedMax      = " + str(
    c_num) + " # Tunning of the speed of light to speed up the numerical simulation\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.ConvTerm.lightSpeedMF       = " + str(c_num) + "\n")
f.write("Simulator.SubSystem.MultiFluidMHD2D.ConvTerm.IsLeake            = false # Just when the first fluid is ions and the second is neutrals\n")
f.write("\n")
f.write("################################################################################ \n")
f.write("# Output\n")
f.write("\n")
f.write("Simulator.SubSystem.OutputFormat        = Tecplot CFmesh \n")
f.write("Simulator.SubSystem.CFmesh.FileName     = Test-sol.CFmesh \n")
f.write("Simulator.SubSystem.CFmesh.SaveRate     = " + str(mesh_save_rate) + "\n")
f.write("Simulator.SubSystem.CFmesh.AppendTime   = true \n")
f.write("Simulator.SubSystem.CFmesh.AppendIter   = false")
f.write("\n")
f.write("Simulator.SubSystem.Tecplot.FileName       = multiFluid.plt \n")
f.write("Simulator.SubSystem.Tecplot.Data.outputVar = RhoiViTi \n")
f.write("Simulator.SubSystem.Tecplot.SaveRate       = " + str(data_save_rate) + "\n")
f.write("#Simulator.SubSystem.Tecplot.Data.printExtraValues = true \n")
f.write("#Simulator.SubSystem.Tecplot.Data.SurfaceTRS = x0 y0 \n")
f.write("Simulator.SubSystem.Tecplot.AppendTime     = true #false \n")
f.write("Simulator.SubSystem.Tecplot.AppendIter     = true \n")
f.write("#Simulator.SubSystem.Tecplot.WriteSol = ParWriteSolutionBlock \n")
f.write("#Simulator.SubSystem.Tecplot.Data.DataHandleOutput.CCSocketNames   = Qtot GradPyi GradPye \n")
f.write("#Simulator.SubSystem.Tecplot.Data.DataHandleOutput.CCVariableNames = Charge GradPyi GradPye \n")
f.write("#Simulator.SubSystem.Tecplot.Data.DataHandleOutput.CCBlockSize     = 1 1 1 \n")
f.write("\n")
f.write("################################################################################ \n")
f.write("# Time Marching\n")
f.write("\n")
if restart:
    f.write("Simulator.SubSystem.InitialTime = " + str(start_time) + "\n")
f.write("#Simulator.SubSystem.SubSystemStatus.TimeStep = 1.0e-3\n")
f.write("#Simulator.SubSystem.SubSystemStatus.ComputeDT = FunctionDT\n")
f.write("#Simulator.SubSystem.SubSystemStatus.FunctionDT.Vars = i\n")
f.write("#Simulator.SubSystem.SubSystemStatus.FunctionDT.Def = if(i<101,5.0e-3,if(i<111,1e-2,2e-2))\n")
f.write("\n")
f.write("# Stop Condition\n")
f.write("Simulator.SubSystem.StopCondition       = MaxNumberSteps\n")
f.write("Simulator.SubSystem.MaxNumberSteps.nbSteps = " + str(MaxSteps) + "\n")
f.write("\n")
f.write("#Simulator.SubSystem.StopCondition   = MaxTime\n")
f.write("#Simulator.SubSystem.MaxTime.maxTime = 16\n")
f.write("\n")
f.write("#Simulator.SubSystem.StopCondition       = Norm\n")
f.write("#Simulator.SubSystem.Norm.valueNorm      = -20.0\n")
f.write("\n")
f.write("# Linear System\n")
f.write("Simulator.SubSystem.LinearSystemSolver = PETSC\n")
f.write("Simulator.SubSystem.LSSNames = NewtonIteratorLSS\n")
f.write("Simulator.SubSystem.NewtonIteratorLSS.Data.PCType = PCASM\n")
f.write("Simulator.SubSystem.NewtonIteratorLSS.Data.KSPType = KSPGMRES\n")
f.write("Simulator.SubSystem.NewtonIteratorLSS.Data.MatOrderingType = MATORDERING_RCM\n")
f.write("Simulator.SubSystem.NewtonIteratorLSS.Data.MaxIter = 1000\n")
f.write("Simulator.SubSystem.NewtonIteratorLSS.Data.NbKrylovSpaces = 150\n")
f.write("Simulator.SubSystem.NewtonIteratorLSS.Data.RelativeTolerance = 1e-4\n")
f.write("\n")
f.write("# Explicit Solver\n")
f.write("#Simulator.SubSystem.ConvergenceMethod = FwdEuler\n")
f.write("#Simulator.SubSystem.FwdEuler.Data.CFL.Value = 0.1\n")
f.write("#Simulator.SubSystem.FwddEuler.Data.CFL.ComputeCFL =  Interactive\n")
f.write("\n")
f.write("# Implicit first Order\n")
f.write("#Simulator.SubSystem.ConvergenceMethod = NewtonIterator\n")
f.write("#Simulator.SubSystem.NewtonIterator.UpdateSol = StdUpdateSol\n")
f.write("#Simulator.SubSystem.NewtonIterator.StdUpdateSol.Relaxation= 1.\n")
f.write("\n")
f.write("# CFL definition\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.CFL.ComputeCFL = Interactive\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.CFL.ComputeCFL = Function\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.CFL.Function.Def =\ \n")
f.write("#if(i<160,1e4,if(i<250,1e5,if(i<670,1e6,if(i<2690,1e7,1e8))))\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.CFL.Value = 1e4\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.CFL.Interactive.CFL = 1.0\n")
f.write("\n")
f.write("#First order in time\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.MaxSteps = 30\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.L2.MonitoredVarID = 15\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.FilterState = Max\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.Max.maskIDs = 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.Max.minValues = 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 10. 10.\n")
f.write("#Simulator.SubSystem.NewtonIterator.Data.Norm = -15\n")
f.write("\n")
f.write("# Implicit Second Order\n")
f.write("Simulator.SubSystem.ConvergenceMethod = BDF2\n")
f.write("Simulator.SubSystem.BDF2.ShowRate = 1\n")
f.write("# CFL definition\n")
f.write("#Simulator.SubSystem.BDF2.Data.CFL.Value = 1.0\n")
f.write("#Simulator.SubSystem.BDF2.Data.CFL.ComputeCFL = Interactive\n")
f.write("#Simulator.SubSystem.BDF2.ConvergenceFile = convergence_UnsteadyMagnetosphereACAImplPrim0_85READFROMFILE.plt \n")
f.write("Simulator.SubSystem.BDF2.Data.MaxSteps = 1000\n")
f.write("Simulator.SubSystem.BDF2.Data.L2.GlobalRes = true\n")
f.write("Simulator.SubSystem.BDF2.Data.L2.NormalizedRes = true\n")
f.write("Simulator.SubSystem.BDF2.Data.L2.RefVals = " + str(Bx_ref) + " " + str(By_ref) + " " + str(Bz_ref) + " " + str(Ex_ref) + " " + str(Ey_ref) + " " + str(Ez_ref) + " " +
        str(Psi_ref) + " " + str(Phi_ref) + " " + str(rhoe_ref) + " " + str(rhoi_ref) + " " + str(ue_ref) + " " + str(ve_ref) + " " + str(we_ref) + " " + str(ui_ref) + " " +
        str(vi_ref) + " " + str(wi_ref) + " " + str(Te_ref) + " " + str(Ti_ref) + " # Characteristic value of the solution for the computation of the numerical Jacobian of the convective Flux\n")
f.write("Simulator.SubSystem.BDF2.Data.Norm = " + str(Norm_param) + "\n")
f.write("\n")
f.write("################################################################################ \n")
f.write("# Mesh Reader\n")
f.write("Simulator.SubSystem.Default.listTRS = PeriodicX PeriodicY # x0 x1\n")
f.write("\n")
f.write("Simulator.SubSystem.MeshCreator = CFmeshFileReader\n")
f.write("Simulator.SubSystem.CFmeshFileReader.Data.FileName = " + Mesh_file + "\n")
f.write("Simulator.SubSystem.CFmeshFileReader.Data.ScalingFactor = " + str(Scale_factor) + "\n")
f.write("\n")
if not restart:
    f.write("# comment this out to Restart\n")
    f.write("Simulator.SubSystem.CFmeshFileReader.Gambit2CFmesh.Discontinuous = true\n")
    f.write("Simulator.SubSystem.CFmeshFileReader.Gambit2CFmesh.SolutionOrder = P0\n")
    f.write("Simulator.SubSystem.CFmeshFileReader.convertFrom = Gambit2CFmesh \n")
    f.write("Simulator.SubSystem.CFmeshFileReader.ParReadCFmesh.ParCFmeshFileReader.ParMetis.NCommonNodes = 2\n")
    f.write("\n")
f.write("################################################################################ \n")
f.write("# Space Method\n")
f.write("Simulator.SubSystem.SpaceMethod = CellCenterFVM\n")
f.write("Simulator.SubSystem.CellCenterFVM.ComputeRHS = NumJacob\n")
f.write("#Simulator.SubSystem.CellCenterFVM.NumJacob.FreezeDiffCoeff = true\n")
f.write("\n")
f.write("# First Order Time stepping\n")
f.write("#Simulator.SubSystem.CellCenterFVM.ComputeTimeRHS = PseudoSteadyTimeRhs\n")
f.write("#Simulator.SubSystem.CellCenterFVM.PseudoSteadyTimeRhs.zeroDiagValue = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n")
f.write("#Simulator.SubSystem.CellCenterFVM.PseudoSteadyTimeRhs.useGlobalDT = false\n")
f.write("\n")
f.write("# second order Time stepping\n")
f.write("Simulator.SubSystem.CellCenterFVM.ComputeTimeRHS = BDF2TimeRhsLimited\n")
f.write("Simulator.SubSystem.CellCenterFVM.BDF2TimeRhsLimited.TimeLimiter = MinMod\n")
f.write("Simulator.SubSystem.CellCenterFVM.BDF2TimeRhsLimited.MinMod.SlopeRatio = 3.\n")
f.write("Simulator.SubSystem.CellCenterFVM.BDF2TimeRhsLimited.zeroDiagValue = 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0\n")
f.write("\n")
f.write("### second order: uncomment this\n")
f.write("Simulator.SubSystem.CellCenterFVM.SetupCom = LeastSquareP1Setup\n")
f.write("Simulator.SubSystem.CellCenterFVM.SetupNames = Setup1\n")
f.write("Simulator.SubSystem.CellCenterFVM.Setup1.stencil = FaceVertexPlusGhost\n")
f.write("Simulator.SubSystem.CellCenterFVM.UnSetupCom = LeastSquareP1UnSetup\n")
f.write("Simulator.SubSystem.CellCenterFVM.UnSetupNames = UnSetup1\n")
f.write("\n")
f.write("#Simulator.SubSystem.CellCenterFVM.Data.PolyRec = Constant\n")
f.write("\n")
f.write("## second order: uncomment this\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.PolyRec = LinearLS2DPeriodic\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.LinearLS2DPeriodic.limitRes = -4.0\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.LinearLS2DPeriodic.gradientFactor = 1.\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.LinearLS2DPeriodic.PeriodicBCNames = Jet1 Jet2\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.Limiter = Venktn2D\n")
if Limiter:
    f.write("Simulator.SubSystem.CellCenterFVM.Data.Venktn2D.isMFMHD = true\n")
    f.write("Simulator.SubSystem.CellCenterFVM.Data.Venktn2D.coeffEps = 1.\n")
    f.write("Simulator.SubSystem.CellCenterFVM.Data.Venktn2D.length =" + str(l_0) + "\n")
    f.write("Simulator.SubSystem.CellCenterFVM.Data.Venktn2D.magnitudeValues = " + str(Bx_ref) + " " + str(By_ref) + " " + str(Bz_ref) + " " + str(Ex_ref) + " " + str(Ey_ref) + " " + str(Ez_ref) + " " +
        str(Psi_ref) + " " + str(Phi_ref) + " " + str(rhoe_ref) + " " + str(rhoi_ref) + " " + str(c_se) + " " + str(c_se) + " " + str(c_se) + " " + str(ui_ref) + " " +
        str(vi_ref) + " " + str(wi_ref) + " " + str(Te_ref) + " " + str(Ti_ref) + " # \n")
f.write("# Simulator.SubSystem.CellCenterFVM.Data.NodalExtrapolation = DistanceBased\n")
f.write("# Simulator.SubSystem.CellCenterFVM.Data.DistanceBased.TrsPriorityList = \ \n")
f.write("#y0 PeriodicX y1\n")
f.write("\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.FluxSplitter = AUSMPlusUpMultiFluid2D\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.AUSMPlusUpMultiFluid2D.useMacCormackScaling = false\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.AUSMPlusUpMultiFluid2D.choiceA12 = 1\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.AUSMPlusUpMultiFluid2D.machInf = 1. 1.\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.AUSMPlusUpMultiFluid2D.Coeff = 1.e4 \n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.AUSMPlusUpMultiFluid2D.Bdiss = 1. \n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.AUSMPlusUpMultiFluid2D.Ediss = 1. \n")
f.write("\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.UpdateVar = RhoiViTi\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.SolutionVar = Cons\n")
f.write("#Simulator.SubSystem.CellCenterFVM.Data.DiffusiveVar = RhoiViTi\n")
f.write("#Simulator.SubSystem.CellCenterFVM.Data.DerivativeStrategy = Corrected2D\n")
f.write("\n")
f.write("#Simulator.SubSystem.CellCenterFVM.Data.DiffusiveFlux = NavierStokesMF\n")
f.write("#Simulator.SubSystem.CellCenterFVM.Data.isAxisymm = true\n")
f.write("\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.SourceTerm = " + Source_term + "\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.DriftWaves2DHalfTwoFluid.isCollisional = " + Collisional + "\n")
f.write("Simulator.SubSystem.CellCenterFVM.Data.DriftWaves2DHalfTwoFluid.ElectricCharge = " + str(charge) + "\n")
f.write("\n")
f.write("################################################################################ \n")
f.write("# Initial Conditions Bx, By, Bz, Ex, Ey, Ez, Psi, Phi, rhoe, rhoi, ue, ve, we, ui, vi, wi, Te, Ti\n")
f.write("#order of species is first electrons then of ions!\n")
if restart:
    f.write("Simulator.SubSystem.CellCenterFVM.Restart = true\n")
if not restart:
    f.write("Simulator.SubSystem.CellCenterFVM.InitComds = InitStateAddVar \n")
    f.write("Simulator.SubSystem.CellCenterFVM.InitNames = InField \n")
    f.write("Simulator.SubSystem.CellCenterFVM.InField.applyTRS = InnerFaces \n")
    f.write("Simulator.SubSystem.CellCenterFVM.InField.InitVars = x y \n")
    f.write("Simulator.SubSystem.CellCenterFVM.InField.Vars = x y xStar Bpar Bperp Epar Eperp Uepar Ueperp Uipar Uiperp sinalpha cosalpha \n")
    f.write("Simulator.SubSystem.CellCenterFVM.InField.InitDef = (x*sqrt(1./5.)+y*sqrt(4./5.))\ \n") #sqrt(1./5.) = cos(alpha) and sqrt(1./5.) = sin(alpha)
    f.write("                                            " + str(B_0) + "\ \n")
    f.write("                                            " + str(By_ref) + "*cos(" + str(kappa) + "*(x*sqrt(1./5.)+y*sqrt(4./5.)))\ \n")
    f.write("                                            0.\ \n")
    f.write("                                            " + str(-Ey_ref) + "*sin(" + str(kappa) + "*(x*sqrt(1./5.)+y*sqrt(4./5.)))\ \n")
    f.write("                                            0.\ \n")
    f.write("                                            " + str(ue_ref) + "*cos(" + str(kappa) + "*(x*sqrt(1./5.)+y*sqrt(4./5.)))\ \n")
    f.write("                                            0.\ \n")
    f.write("                                            " + str(ui_ref) + "*cos(" + str(kappa) + "*(x*sqrt(1./5.)+y*sqrt(4./5.)))\ \n")
    f.write("                                            sqrt(4./5.)\ \n")
    f.write("                                            sqrt(1./5.) \n")
    f.write("Simulator.SubSystem.CellCenterFVM.InField.Def =\ \n")
    f.write("                                               Bpar*cosalpha-Bperp*sinalpha\ \n") # this is all B_x
    f.write("                                               Bpar*sinalpha+Bperp*cosalpha\ \n") # this is all B_y
    f.write("                                               " + str(Bz_ref) + "*sin(" + str(kappa) + "*xStar)\ \n") # this is all B_z
    f.write("                                               Epar*cosalpha-Eperp*sinalpha\ \n") # this is the value for E_x
    f.write("                                               Epar*sinalpha+Eperp*cosalpha\ \n") # this is all E_y
    f.write("                                               " + str(Ez_ref) + "*cos(" + str(kappa) + "*xStar)\ \n") # this is all E_z
    f.write("                                               0.0\ \n") # this is the value for Psi
    f.write("                                               0.0\ \n") # this is the value for Phi
    f.write("                                               " + str(rhoe_ref) + "\ \n") # this is all rho_e
    f.write("                                               " + str(rhoi_ref) + "\ \n") # this is all rho_i
    f.write("                                               Uepar*cosalpha-Ueperp*sinalpha\ \n") # this is the value for u_e
    f.write("                                               Uepar*sinalpha+Ueperp*cosalpha\ \n") # this is all v_e
    f.write("                                               " + str(ue_ref) + "*sin(" + str(kappa) + "*xStar)\ \n") # this is all w_e
    f.write("                                               Uipar*cosalpha-Uiperp*sinalpha\ \n") # this is the value for u_i
    f.write("                                               Uipar*sinalpha+Uiperp*cosalpha\ \n") # this is all v_i
    f.write("                                               " + str(ui_ref) + "*sin(" + str(kappa) + "*xStar)\ \n") # this is all w_i
    f.write("                                               " + str(Te_ref) +  "\ \n") # this is the value for T_e
    f.write("                                               " + str(Ti_ref) +  "\n") # this is the value for T_i
    f.write("\n")
f.write("################################################################################ \n")
f.write("# Boundary Conditions\n")
f.write("Simulator.SubSystem.CellCenterFVM.BcComds = BCPeriodicFVMCC BCPeriodicFVMCC ##periodic boundary conditions\n")
f.write("Simulator.SubSystem.CellCenterFVM.BcNames = Jet1 Jet2 \n")
f.write("\n")
f.write("# Bottom Condition(periodic)\n")
f.write("Simulator.SubSystem.CellCenterFVM.Jet1.applyTRS = PeriodicY\n")
f.write("Simulator.SubSystem.CellCenterFVM.Jet1.Threshold = 1\n")
f.write("Simulator.SubSystem.CellCenterFVM.Jet1.TranslationVector = 0. " + str(2*L_y) +"  \n")
f.write("Simulator.SubSystem.CellCenterFVM.Jet1.ZeroGradientFlags = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n") # Depends on the function that we are putting in!
f.write("\n")
f.write("# Left Condition (periodic)\n")
f.write("Simulator.SubSystem.CellCenterFVM.Jet2.applyTRS = PeriodicX\n")
f.write("Simulator.SubSystem.CellCenterFVM.Jet2.Threshold = 1\n")
f.write("Simulator.SubSystem.CellCenterFVM.Jet2.TranslationVector = " + str(2*L_x) + " 0. \n")
f.write("Simulator.SubSystem.CellCenterFVM.Jet2.ZeroGradientFlags = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n") # Depends on the function that we are putting in!
f.write("\n")
f.write("################################################################################ \n")
f.write("# DataProcessing\n")
f.write("\n")
f.write("Simulator.SubSystem.DataPostProcessing = DataProcessing\n")
f.write("Simulator.SubSystem.DataProcessing.Data.updateVar = RhoiViTi\n")
f.write("Simulator.SubSystem.DataProcessing.Comds = GridConvergenceTwoFluid\n")
f.write("Simulator.SubSystem.DataProcessing.Names = GridConvergenceTwoFluid\n")
f.write("Simulator.SubSystem.DataProcessing.GridConvergenceTwoFluid.IsTransversal = true\n")
f.write("Simulator.SubSystem.DataProcessing.GridConvergenceTwoFluid.OutputFileError = ./Error.plt\n")
f.write("Simulator.SubSystem.DataProcessing.ProcessRate = 1 \n")
f.write("#Simulator.SubSystem.DataProcessing.DivMonitoring.options = ./DivMonitoring.plt\n")
### END OF FILE WRITING

f.close()
### Now write the .inter file
g = open(inter_file,'w')
g.write("Simulator.SubSystem.SubSystemStatus.TimeStep = " + str(TimeStep) + "\n")
g.write("Simulator.SubSystem.CellCenterFVM.Data.LinearLS2DPeriodic.gradientFactor = 1.\n")
g.write("Simulator.SubSystem.CellCenterFVM.Data.LinearLS2DPeriodic.limitRes = -20.0\n")
g.close()
