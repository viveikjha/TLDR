include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("DataImportNEW.jl")

include("GenMatrices.jl")

# a=500.0; X=inv(H'*H+a^2*eye(size(H)[2]))*(H'*L); imshow(X.*(X.>0.0),origin="lower",aspect="auto",cmap="Reds")

using PyPlot

DATA = Import_Data(2)

#Create Wavlengths
Ha =6563.0
nlams = 150
lam = collect(1:nlams)+(Ha-nlams/2)
writecsv("UnitTests/UT_wavelengths.csv",lam)
println("Wrote UT_wavelengths.csv with dimensions: ",size(lam))
#Create VDM Vertical Stripe
spread=10.0
ntimes=50

flx_scale=1.0

vdm_vert = zeros(ntimes,nlams)
total_lit=0
for i in 1:ntimes
	for j in 1:nlams
		if i >= (ntimes/2)-(spread/2) && i<(ntimes/2)+spread/2 && j > ((nlams/2)-(spread/2)) && j<=((nlams/2)+(spread/2))
			vdm_vert[i,j]=1.00*flx_scale
			total_lit+=1
		end
	end
end
writecsv("UnitTests/vdm_vert.csv",vdm_vert)
println("Wrote vdm_vert.csv with dimensions: ", size(vdm_vert))
println("For Box, ", total_lit, " pixels lit.")
#Create VDM Horizontal Stripe

#Initialize ADMM Parameters
Pars = init_Params()





#Initialize Matrices
println("-----------")
println("Getting H:")
Mats = Gen_Mats(DATA,Pars)

Spectra_Vert = Mats.H*vdm_vert
writecsv("UnitTests/Spectrac_V.csv", Spectra_Vert)

#Spectra_Horz = Mats.H*vdm_horz
#writecsv("UnitTests/Spectrac_H.csv", Spectra_Horz)

println("-----------")



dims = size(Spectra_Vert)
println("dimensions:", dims)
#Create fake sigmas.
sigma = 1.0*flx_scale
println("Using a sigma of: ",sigma)

n = randn((dims))*sigma #GENERATE NOISE
println("-----------")



Noisy_Spectra_V = Spectra_Vert+n #ADD NOISE
#Noisy_Spectra_H = Spectra_Horz+n

sig_arr = ones(dims)*sigma


writecsv("UnitTests/SpectraN_V.csv", Noisy_Spectra_V')
println("Wrote SpectraN_V.csv with dimensions: ", size(Noisy_Spectra_V'))
#writecsv("UnitTests/SpectraN_H.csv", Noisy_Spectra_H')
#println("Wrote SpectraN_H.csv with dimensions: ", size(Noisy_Spectra_H'))

writecsv("UnitTests/errspectra.csv",sig_arr')
println("Wrote errspectra.csv with dimensions: ", size(sig_arr'))

println("\n Done setting up unit test files! \n")

println("-------------------------------------------")
println("---------  Reconstruction ---------")
println("-------------------------------------------")
DATA=Import_DataN("UnitTests/","UT_wavelengths.csv","SpectraN_V.csv", "errspectra.csv","data/rvm_dates.csv","data/arp151.b.dat")
data_report(DATA)
Mats = Gen_Mats(DATA,Pars)
maxl=maximum(DATA.L)
#Initial Penalty Parameters
Pars.mu_smoo = 2000.0/maxl													#Smoothing Regularization Weight (TIKHONOV)1
Pars.mu_spec = Pars.mu_smoo/2.0															#Spectral Regularization Weight
Pars.mu_temp = Pars.mu_smoo/2.0															#Temporal Regularization Weight
Pars.mu_l1 = Pars.mu_smoo/2.0																		#Ell1 Smoothing Regularization Weight
Pars.nits=100
max_delay=50

###############################
#Checking actual chi2
vdm_path = "UnitTests/vdm_vert.csv"
vdm_act = readcsv(vdm_path)
chi2 = Chi2(Model((vdm_act),Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
println("Chi2: ",chi2)
###############################



tmp,P = TLDR(DATA,Mats,Pars,"True","True","UnitTests/vdm_vert.csv")