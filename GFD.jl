#push!(LOAD_PATH,"/home/matander/TLDR/")
push!(LOAD_PATH,"/home/manderson/TLDR/")
#push!(LOAD_PATH,"/home/manderson/Research/TLDR/")
#push!(LOAD_PATH,"/Users/manderson/Software/ReverbMap/JuliaVersions/TLDR")
using RMLib
using RMLibMore
using RMTypes
using DataImport
using GenMatrices
using DataImportNEW
using GenMatrices
# a=500.0; X=inv(H'*H+a^2*eye(size(H)[2]))*(H'*L); imshow(X.*(X.>0.0),origin="lower",aspect="auto",cmap="Reds")

using PyPlot

DATA = Import_Data(2)

#Create Wavlengths
Ha =6563.0
nlams = 20
ntimes=50
lam = collect(1:nlams)+(Ha-nlams/2)
writecsv("UnitTests/UT_Wavelengths.csv",lam)

#println("Wrote UT_wavelengths.csv with dimensions: ",size(lam))
#Create VDM Vertical Stripe
spread=10.0

#lvl=0.1
lvl=0.1

#DSNR=100.0
DSNR=50.0 #33.33 is the proper test value

flx_scale=1.0

vdm = zeros(ntimes,nlams)+0.001
total_lit=0
for i in 1:ntimes
	for j in 1:nlams
		if i >= (ntimes/2)-(spread/2) && i<(ntimes/2)+spread/2 && j > ((nlams/2)-(spread/2)) && j<=((nlams/2)+(spread/2))
			vdm[i,j]=lvl*flx_scale
			total_lit+=1
		end
	end
end

#If another pattern is used, load here!
vdm=readcsv("Paper/ring/simulated_ring_vdm.csv")
#vdm=readcsv("Paper/ring/simulated_spiral_vdm.csv")
#vdm=readcsv("checkerboard_20x50/original_checkerboard.csv")
#vdm=readcsv("circle_20x50/original_circle.csv")
#vdm=readcsv("ring_20x50/original_ring.csv")
writecsv("UnitTests/UT_vdm.csv",vdm)
#println("Wrote vdm.csv with dimensions: ", size(vdm))
println("For Box, ", total_lit, " pixels lit.")

#Initialize ADMM Parameters
Pars = init_Params()

#Initialize Matrices
#println("-----------")
#println("Getting H:")
Mats = Gen_Mats(DATA,Pars)
writecsv("UT_H.csv",Mats.H)
Spectra_Vert = Mats.H*vdm
#writecsv("UnitTests/Spectrac_V.csv", Spectra_Vert)


println("-----------")



dims = size(Spectra_Vert)
println("dimensions:", dims)

Spectra_Power = sum(abs.(Spectra_Vert).*abs.(Spectra_Vert))/length(Spectra_Vert)
println("Power level of spectra: ", Spectra_Power)

n =randn((dims))
RN_Power = sum(abs.(n).*abs.(n))/length(n)
println("Power level of raw noise: ", RN_Power)
sf=Spectra_Power/(DSNR*RN_Power)
SN=sqrt(sf)*n
SN_Power= sum(abs.(SN).*abs.(SN))/length(SN)
println("Power level of scaled noise: ", SN_Power)
println("Simulated SNR: ", Spectra_Power/SN_Power)


Noisy_Spectra = (Spectra_Vert+SN)' #ADD NOISE

#sig_arr = ones(dims)*sigma
sig_arr = SN
writecsv("UnitTests/UT_Spectra.csv", Noisy_Spectra)

Error = sig_arr'
writecsv("UnitTests/UT_Spectra_Error.csv",Error)

println("\n Done setting up unit test files! \n")
