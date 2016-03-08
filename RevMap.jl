include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("DataImportNEW.jl")

include("GenMatrices.jl")
using PyPlot


#Mode = 1 for synthdata Mode = 2 for real data!
#mode = 3 for simulated data.
mode = 3
DATA = Import_Data(mode)
println("---------2nd Import----------")
DATAN = Import_DataN("simulation/","new_wavelengths.csv","spectra_simulated.csv", "errspectra_simulated.csv","rvm_dates.csv","data/arp151.b.dat")




#Initialize ADMM Parameters
Pars = init_Params()
Pars.mu_spec = 1.00																								#Spectral Regularization Weight
Pars.mu_temp = 0.10																								#Temporal Regularization Weight
Pars.mu_smoo = 0.10   																						#Smoothing Regularization Weight (TIKHONOV)



#Initialize Matrices
Mats = Gen_Mats(DATA,Pars)

#Check for synthetic data mode
if mode ==3
vdm_path = "simulation/simulated_vdm.csv"
vdm_act = readcsv(vdm_path)

M = Model(vdm_act,Mats.H)
D = DATA.L
Sigma = DATA.EL

true_chi2 = Chi2(M,D,Sigma)/(DATA.num_spectra_samples*DATA.num_lines)
altchi2=alt_chi2(DATA,Mats,vdm_act)/(DATA.num_spectra_samples*DATA.num_lines)
println("L: ",size(DATA.L))
println("lines: ",DATA.num_lines)
println("tdf_times: ",Pars.num_tdf_times)
dims = size(DATA.L)
println("> ",sum((DATA.L-M).^2)/(dims[1]*dims[2]))


writecsv("modeln_rec.csv",M)
writecsv("spectran_rec.csv",D)
writecsv("sigma_rec.csv",DATA.EL)
println("------------------------------------------")
println("----------SYNTHETIC DATA MODE-------------")
println("------------------------------------------")
#println(Pars.num_tdf_times, " ", DATA.num_lines)
println("CHI2 on actual image: ", true_chi2)
println("Alt CHI2: ", altchi2)
println("------------------------------------------")
println("Beginning Reconstruction:")
println("")
end


#Run TLDR
#tmp,P = TLDR(DATA,Mats,Pars,"False","True")

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

