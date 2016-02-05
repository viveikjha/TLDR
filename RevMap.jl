include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("GenMatrices.jl")
using PyPlot


#Mode = 1 for synthdata Mode = 2 for real data!
#mode = 3 for simulated data.
mode = 3
DATA = Import_Data(mode)



#Initialize ADMM Parameters
max_delay = 50.0				#Default: 20.0
Pars = init_Params()
Pars.it = 2																													#Current Iteration (Starts at 2)
Pars.nits =200 																											#Number of Iterations Allowed
Pars.tau = 1.2
Pars.sigma=0.75
Pars.G = 10.0
Pars.num_lines = DATA.num_lines
Pars.num_tdf_times = 50																							#Number of points in reconstructed TDF
Pars.tdf_values = zeros(Pars.num_tdf_times)
Pars.tdf_times=collect(1.0:((max_delay-1.0)/(Pars.num_tdf_times-1)):max_delay)
Pars.initial_psi = 0.0																							#Initial value for TDFs
Pars.eps_abs = 0.0	
Pars.eps_rel = 0.01
Pars.alpha = 1.0																										#Muliplier's throttle
Pars.rho0 = 500.0																										#Initial Penalty Parameter
Pars.mu_spec = 1.0																									#Spectral Regularization Weight
Pars.mu_temp = 1.0																									#Temporal Regularization Weight
Pars.mu_smoo = 1.0   																							#Smoothing Regularization Weight (TIKHONOV)



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
println(Pars.num_tdf_times, " ", DATA.num_lines)
println("CHI2 on actual image: ", true_chi2)
println("Alt CHI2: ", altchi2)
println("------------------------------------------")
println("Beginning Reconstruction:")
println("")
end


#Run TLDR
tmp,P = TLDR(DATA,Mats,Pars,"True","True")

#Save Output

Write_FITS(tmp,P)
writecsv("vdm.csv",tmp.vdm)

figure()
imshow(tmp.vdm,aspect="auto",origin="lower",cmap="Reds",extent=[minimum(DATA.wavelength),maximum(DATA.wavelength),0,max_delay],interpolation="None")
title("Recovered TDF")
xlabel("Spectral Channel")
ylabel("Delay")
savefig("RecTDF.png",format="png")
println("done")

