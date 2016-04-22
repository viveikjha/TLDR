include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("DataImportNEW.jl")

include("GenMatrices.jl")


using PyPlot

DATA = Import_Data(2)

Pars = init_Params()
																								#Initial Penalty Parameters
Pars.mu_spec = 5.0#1.0															#Spectral Regularization Weight
Pars.mu_temp = 5.0#0.5															#Temporal Regularization Weight
Pars.mu_smoo = 10.0#0.001															#Smoothing Regularization Weight (TIKHONOV)1
Pars.nits=1000
max_delay=50




#Initialize Matrices
println("-----------")
println("Getting H:")
Mats = Gen_Mats(DATA,Pars)


#Import VDM
vdm = readcsv("simulation/Fim.csv")
Spectra = Mats.H*vdm
writecsv("simulation/Spectra.csv", Spectra)
println("-----------")



dims = size(Spectra)
println("dimensions:", dims)
#Create fake sigmas.
sigma = 1.0
println("Using a sigma of: ",sigma)

n = randn((dims))*sigma #GENERATE NOISE
println("-----------")



Noisy_Spectra = Spectra+n #ADD NOISE


sig_arr = ones(dims)*sigma


writecsv("simulation/SpectraN.csv", Noisy_Spectra')
println("Wrote SpectraN.csv with dimensions: ", size(Noisy_Spectra))

writecsv("simulation/errspectra.csv",sig_arr')
println("Wrote errspectra.csv with dimensions: ", size(sig_arr))

#println("\n Done setting up unit test files! \n")

println("-------------------------------------------")
println("---------  Reconstruction ---------")
println("-------------------------------------------")


Pars = init_Params()
																								#Initial Penalty Parameters
Pars.mu_spec = 5.0#1.0															#Spectral Regularization Weight
Pars.mu_temp = 5.0#0.5															#Temporal Regularization Weight
Pars.mu_smoo = 10.0#0.001															#Smoothing Regularization Weight (TIKHONOV)1
Pars.nits=600
max_delay=50




DATA=Import_DataN("simulation/","new_wavelengths.csv","SpectraN.csv", "errspectra.csv","data/rvm_dates.csv","data/arp151.b.dat")
data_report(DATA)
Mats = Gen_Mats(DATA,Pars)
###############################
#Checking actual chi2
vdm_path = "simulation/Fim.csv"
vdm_act = readcsv(vdm_path)
chi2 = Chi2(Model((vdm_act),Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
println("Chi2: ",chi2)
###############################


#figure()
#imshow(vdm_act,interpolation="None",cmap="Greens",origin="lower")
#show()
tmp,P = TLDR(DATA,Mats,Pars,"True","True","simulation/Fim.csv")
