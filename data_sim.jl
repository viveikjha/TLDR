include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("GenMatrices.jl")


using PyPlot

sim_wavelenghts = readcsv("simulation/new_wavelengths.csv")
vdm_simulated = readcsv("simulation/simulated_vdm.csv")
DATA = Import_Data(2)


#Initialize ADMM Parameters
max_delay = 50.0				#Default: 20.0
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
Pars.initial_psi = 0.1																							#Initial value for TDFs
Pars.eps_abs = 0.01	
Pars.eps_rel = 0.1
Pars.alpha = 1.0																										#Muliplier's throttle
Pars.rho0 = 50.0																										#Initial Penalty Parameter
Pars.mu_spec = 10000.0																									#Spectral Regularization Weight
Pars.mu_temp = 10000.0																									#Temporal Regularization Weight
Pars.mu_smoo = 10000.0    																							#Smoothing Regularization Weight (TIKHONOV)




println("size of sim_waves: ", length(sim_wavelenghts))
println("shape of sim_vdm: ", size(vdm_simulated))

#Initialize Matrices
Mats = Gen_Mats(DATA,Pars)

#Spectra = Mats.H*vdm_simulated
Spectra = Mats.H*vdm_simulated
writecsv("Spectrac_sim.csv", Spectra)
#println(Spectra[1:10])
println()



dims = size(Spectra)
println("dimensions:", dims)
#Create fake sigmas.
sigma = 5.0

n = randn((dims))*sigma #GENERATE NOISE
println("-----------")


println("total: ", sum(n.^2)/(dims[1]*dims[2]*sigma^2))
#println(size(Spectra))

Noisy_Spectra = Spectra+n #ADD NOISE
#println(Spectra[1:10])

sig_arr = ones(dims)*sigma
println("> ",sum((Noisy_Spectra-Spectra).^2)/(dims[1]*dims[2]))

writecsv("Spectran_sim.csv", Noisy_Spectra)



println("Array Sizes for Output Files:")
println("		Simulated VDM: ", size(vdm_simulated))
println("		Mapping Matrix: ", size(Mats.H))
println("		Simulated Spectr: ", size(Spectra'))



println("		Simulated Sigmas: ",size(sig_arr'))

writecsv("simulation/errspectra_simulated.csv",sig_arr')
writecsv("simulation/spectra_simulated.csv",Noisy_Spectra')
println("Simulated Spectra Written.")
figure()
imshow((vdm_simulated),cmap="Reds",aspect="auto",origin="lower",extent=[minimum(sim_wavelenghts),maximum(sim_wavelenghts),0,max_delay],interpolation="None")
title("Simulated TDF")
xlabel("Wavelength")
ylabel("Delay")
savefig("SimTDF.png",format="png")
#show()
