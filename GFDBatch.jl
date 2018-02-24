push!(LOAD_PATH,"/home/matander/TLDR/")
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
nlams = 10
ntimes=10
lam = collect(1:nlams)+(Ha-nlams/2)

#Initialize ADMM Parameters
Pars = init_Params()



#DATA.num_spectra_samples=nlams
Pars.num_tdf_times=ntimes

#Initialize Matrices
Mats = Gen_Mats(DATA,Pars)


#DSNR=100.0
DSNR=50.0 #33.33 is the proper test value
flx_scale=1.0


pb="10x10/"
names=["box","checkerboard","circle","diagonal","diagonalinverted","halfbottom","halfleft","halfright","halftop","horizontalstripe","invertedbox","invertedhorizontalstripe","invertedverticalstripe","lowertri","reverseddiagonal","reverseddiagonalinverted","ring","uppertri","verticalstripe"]
for i in range(1,length(names))
	name=names[i]
	fname=string(pb,name,"/",name,".csv")
	vdm=readcsv(fname)


	Spectra_Vert = Mats.H*vdm
	println("-----------")
	println("File: ",name)
	println("-----------")

	dims = size(Spectra_Vert)
	println("dimensions:", dims)

	Spectra_Power = sum(abs.(Spectra_Vert).*abs.(Spectra_Vert))/length(Spectra_Vert)
	#println("Power level of spectra: ", Spectra_Power)

	n =randn((dims))
	RN_Power = sum(abs.(n).*abs.(n))/length(n)
	#println("Power level of raw noise: ", RN_Power)
	sf=Spectra_Power/(DSNR*RN_Power)
	SN=sqrt(sf)*n
	SN_Power= sum(abs.(SN).*abs.(SN))/length(SN)
	#println("Power level of scaled noise: ", SN_Power)
	println("Simulated SNR: ", Spectra_Power/SN_Power)
	Noisy_Spectra = (Spectra_Vert+SN)' #ADD NOISE
	sig_arr = SN
	Error = sig_arr'

	writecsv(string(pb,name,"/","UT_Wavelengths.csv"),lam)
	writecsv(string(pb,name,"/","UT_vdm.csv"),vdm)
	writecsv(string(pb,name,"/","UT_H.csv"),Mats.H)
	writecsv(string(pb,name,"/","UT_Spectra.csv"), Noisy_Spectra)
	writecsv(string(pb,name,"/","UT_Spectra_Error.csv"),Error)

	#println("\n Done setting up unit test files! \n")
end
