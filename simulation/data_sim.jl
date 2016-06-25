include("../RMLib.jl")
include("../RMTypes.jl")
include("../DataImport.jl")
include("../GenMatrices.jl")
include("../DataImportNEW.jl")

using PyPlot

#sim_wavelenghts = readcsv("new_wavelengths.csv")
vdm_simulated = readcsv("Spiral/simulated_vdm.csv")
#DATA = Import_Data(2)
DATA=Import_DataN("../data/","rvm_wavelengths.csv","rvm_fluxes.csv","rvm_errfluxes.csv","rvm_dates.csv","arp151.b.dat")

#Initialize ADMM Parameters
Pars = init_Params()
																								#Initial Penalty Parameters
#Pars.mu_spec = 1.0															#Spectral Regularization Weight
#Pars.mu_temp = 1.0															#Temporal Regularization Weight
#Pars.mu_smoo = 1.0															#Smoothing Regularization Weight (TIKHONOV)
#max_delay=50


#Initialize Matrices
Mats = Gen_Mats(DATA,Pars)

#CALCULATE SNR ARRAY FOR ACTUAL DATA


SNR = abs(DATA.L./DATA.EL)
println("SNR size: ", size(SNR))
println("Mats.H ", size(Mats.H))
H2=readcsv("H.csv")
println(maximum(abs(Mats.H-H2)))
Spectra = Mats.H*vdm_simulated #NOISELESS SPECTRA
println("Spectra Size: ", size(Spectra))
#writecsv("Spectrac_sim.csv", Spectra)



nbase=vec(readcsv("../data/rvm_errfluxes.csv"))
n = rand(nbase,length(Spectra)).*rand([-1.0,1.0],length(Spectra))
n= reshape(n,size(Spectra))
println("-----------")
Error= n



println("Zeros in SNR: ", size(find(SNR.==0.0)))
println("Checking for Nans in SNR: ",any(isnan,SNR))
println("Checking for Nans in Error: ",any(isnan,Error))
println("Checking for Nans in Spectra: ",any(isnan,Spectra))

Noisy_Spectra = Spectra+Error #ADD NOISE
println("Checking for Nans in Noisy Spectra: ",any(isnan,Noisy_Spectra))

#nd=find(isnan())



println("-----------")
println("Spectra Min: ",minimum(Spectra))
println("Spectra Max: ",maximum(Spectra))
println("-----------")
println("Noisy Spectra Min: ",minimum(Noisy_Spectra))
println("Noisy Spectra Max: ",maximum(Noisy_Spectra))
println("-----------")


println("Array Sizes for Output Files:")
println("		Simulated VDM: ", size(vdm_simulated))
println("		Mapping Matrix: ", size(Mats.H))
println("		Simulated Spectr: ", size(Spectra'))


writecsv("Sim_Error.csv",Error')
writecsv("Sim_Spectra.csv",Noisy_Spectra')
