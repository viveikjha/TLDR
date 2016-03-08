include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("GenMatrices.jl")


using PyPlot

sim_wavelenghts = readcsv("simulation/new_wavelengths.csv")
vdm_simulated = readcsv("simulation/simulated_vdm.csv")
DATA = Import_Data(2)


#Initialize ADMM Parameters
Pars = init_Params()
																								#Initial Penalty Parameters
Pars.mu_spec = 1.0															#Spectral Regularization Weight
Pars.mu_temp = 1.0															#Temporal Regularization Weight
Pars.mu_smoo = 1.0															#Smoothing Regularization Weight (TIKHONOV)
max_delay=50



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
