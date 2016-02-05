include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("GenMatrices.jl")
using PyPlot


#Mode = 1 for synthdata Mode = 2 for real data!
DATA = Import_Data(2)


vdm_path = "synth/TDF.csv"
vdm_act = readcsv(vdm_path)


#Initialize ADMM Parameters
max_delay = 20.0				#Default: 20.0
Pars = init_Params()
Pars.it = 2																													#Current Iteration (Starts at 2)
Pars.nits =80 																											#Number of Iterations Allowed
Pars.tau = 1.2
Pars.sigma=0.75
Pars.G = 10.0
Pars.num_lines = DATA.num_lines
Pars.num_tdf_times = 50																							#Number of points in reconstructed TDF
Pars.tdf_values = zeros(Pars.num_tdf_times)
Pars.tdf_times=collect(1.0:((max_delay-1.0)/(Pars.num_tdf_times-1)):max_delay)
Pars.initial_psi = 0.0																							#Initial value for TDFs
Pars.eps_abs = 0.0	
Pars.eps_rel = 0.1
Pars.alpha = 1.0																										#Muliplier's throttle
Pars.rho0 = 50.0																										#Initial Penalty Parameter
Pars.mu_spec = 5000.0																									#Spectral Regularization Weight
Pars.mu_temp = 5000.0																									#Temporal Regularization Weight
Pars.mu_smoo = 5000.0    																							#Smoothing Regularization Weight (TIKHONOV)



#Initialize Matrices
Mats = Gen_Mats(DATA,Pars)
X = ones((50,50))
res = discon(Mats.H,X)
println(res-Mats.H)

#println("H: ",size(Mats.H))
println("Done")


