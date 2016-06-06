include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("DataImportNEW.jl")

include("GenMatrices.jl")
using PyPlot


#Mode = 1 for synthdata Mode = 2 for real data!
#mode = 3 for simulated data.
#mode = 2
#DATA = Import_Data(mode)
println("---------2nd Import----------")
DATA = Import_DataN("","data/rvm_wavelengths_trimmed.csv","data/rvm_fluxes_trimmed.csv", "data/rvm_errfluxes_trimmed.csv","data/rvm_dates.csv","data/arp151.b.dat")
#DATA = Import_DataN("data/","rvm_wavelengths.csv","rvm_fluxes.csv", "rvm_errfluxes.csv","data/rvm_dates.csv","data/arp151.b.dat")





#Initialize ADMM Parameters
Pars = init_Params()
#Initial Penalty Parameters
#flx_scale=1.0/100
#Pars.mu_smoo = 10.0/flx_scale														#Smoothing Regularization Weight (TIKHONOV)1
#Pars.mu_spec = 5.0/flx_scale															#Spectral Regularization Weight
#Pars.mu_temp = 5.0/flx_scale															#Temporal Regularization Weight
#Pars.mu_l1 = 5.0/flx_scale																		#Ell1 Smoothing Regularization WeightPars.nits=800
#Pars.nits=400
#max_delay=50
#Pars.tau=2.0


#Initial Penalty Parameters
maxl=maximum(DATA.L)
Pars.mu_smoo = round(2000.0/maxl)													#Smoothing Regularization Weight (TIKHONOV)1
println("µ initialized as: ", Pars.mu_smoo)
Pars.mu_spec = Pars.mu_smoo/2.0															#Spectral Regularization Weight
Pars.mu_temp = Pars.mu_smoo/2.0															#Temporal Regularization Weight
Pars.mu_l1 = Pars.mu_smoo/2.0																		#Ell1 Smoothing Regularization Weight
Pars.nits=100
max_delay=50

#Initialize Matrices
Mats = Gen_Mats(DATA,Pars)





#Run TLDR
tmp,P = TLDR(DATA,Mats,Pars,"True","True")

#Save Output

#Write_FITS(tmp,P)
#writecsv("vdm.csv",tmp.vdm)

#figure()
#imshow(tmp.vdm,aspect="auto",origin="lower",cmap="Reds",extent=[minimum(DATA.wavelength),maximum(DATA.wavelength),0,max_delay],interpolation="None")
#title("Recovered TDF")
#xlabel("Spectral Channel")
#ylabel("Delay")
#savefig("RecTDF.png",format="png")

#M = Model(tmp.vdm,Mats.H)
#D = DATA.L
#Sigma = DATA.EL
#true_chi2 = Chi2(M,D,Sigma)/(DATA.num_spectra_samples*DATA.num_lines)
#altchi2=alt_chi2(DATA,Mats,tmp.vdm)/(DATA.num_spectra_samples*DATA.num_lines)
#println("------------------------------------------")
#println("---------- FINAL INFORMATION -------------")
#println("------------------------------------------")
#println(Pars.num_tdf_times, " ", DATA.num_lines)
#println("CHI2 on final image: ", true_chi2)
#println("Alt CHI2: ", altchi2)
#println("------------------------------------------")
println("done")
