push!(LOAD_PATH,"/home/manderson/TLDR/")
#push!(LOAD_PATH,"/home/manderson/Research/TLDR/")

using RMLib
using RMLibMore
using RMTypes
using DataImport
using DataImportNEW
using GenMatrices

using PyPlot

DATA = Import_Data(2)

#Create Wavelengths
Ha =6563.0
nlams = 20
ntimes=50
lam = collect(1:nlams)+(Ha-nlams/2)
writecsv("UnitTests/UT_Wavelengths.csv",lam)
println("Wrote UT_wavelengths.csv with dimensions: ",size(lam))
#Create VDM Vertical Stripe
spread=10.0

lvl=0.1
flx_scale=1.0

vdm_vert = zeros(ntimes,nlams)+0.001
total_lit=0
for i in 1:ntimes
	for j in 1:nlams
		if i >= (ntimes/2)-(spread/2) && i<(ntimes/2)+spread/2 && j > ((nlams/2)-(spread/2)) && j<=((nlams/2)+(spread/2))
			vdm_vert[i,j]=lvl*flx_scale
			total_lit+=1
		end
	end
end
writecsv("UnitTests/UT_vdm.csv",vdm_vert)
println("For Box, ", total_lit, " pixels lit.")
Pars = init_Params()



Mats = Gen_Mats(DATA,Pars)
writecsv("UT_H.csv",Mats.H)
Spectra = Mats.H*vdm_vert
writecsv("UnitTests/Spectrac.csv", Spectra)
V1=var(Spectra)
dims = size(Spectra)
println("dimensions:", dims)


#Create fake sigmas.
DSNR=1.0 #DESIRED SNR
N_AMP=SIG_LVL/DSNR

n = randn((dims))+DSNR #GENERATE NOISE
#GENERATE NOISE BY RANDOMLY SAMPLING REAL DATA ERROR
#nbase=vec(readcsv("data/rvm_errfluxes.csv"))
#n = rand(nbase,length(Spectra)).*rand([-1.0,1.0],length(Spectra))
#n= reshape(n,size(Spectra))
#println("-----------")

Noisy_Spectra = (Spectra+n)' #ADD NOISE
sig_arr = n
writecsv("UnitTests/UT_Spectra.csv", Noisy_Spectra)

Error = sig_arr'
writecsv("UnitTests/UT_Spectra_Error.csv",Error)

println("Max Flux: ",maximum(Noisy_Spectra))
println("Max Flux Error: ",maximum(abs(Error)))
println("Max SNR: ", maximum(Noisy_Spectra./(Error)))
println("------------------------")
println("Median Flux: ",median(Noisy_Spectra))
println("Median Flux Error: ",median(abs(Error)))
println("Median SNR: ", median(Noisy_Spectra./(Error)))
println("------------------------")
println("Mean Flux: ",mean(Noisy_Spectra))
println("Mean Flux Error: ",mean(abs(Error)))
println("Mean SNR: ", mean(Noisy_Spectra./(Error)))
println("------------------------")
println("\n Done setting up unit test files! \n")
