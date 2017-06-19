#include("RMLibMore.jl")
#include("RMLib.jl")

push!(LOAD_PATH,"/Users/manderson/Software/ReverbMap/JuliaVersions/TLDR")
using RMTypes
using RMLib
using RMLibMore
using GenMatrices
using DataImportNEW

function GenData(Pars,Mats,lvl,DSNR)
  #Create Wavelengths
  Ha =6563.0
  nlams = 20
  ntimes=50
  lam = collect(1:nlams)+(Ha-nlams/2)
  writecsv("UnitTests/UT_Wavelengths.csv",lam)
  println("Wrote UT_wavelengths.csv with dimensions: ",size(lam))
  #Create VDM Vertical Stripe
  spread=10.0

  #lvl=0.1
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



  #Mats = Gen_Mats(DATA,Pars)
  writecsv("UT_H.csv",Mats.H)
  Spectra = Mats.H*vdm_vert
  writecsv("UnitTests/Spectrac.csv", Spectra)
  V1=var(Spectra)
  dims = size(Spectra)
  println("dimensions:", dims)


  #Create fake sigmas.
  #DESIRED SNR
  N_AMP=lvl/DSNR

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

  #println("Max Flux: ",maximum(Noisy_Spectra))
  #println("Max Flux Error: ",maximum(abs(Error)))
  #println("Max SNR: ", maximum(Noisy_Spectra./(Error)))
  #println("------------------------")
  #println("Median Flux: ",median(Noisy_Spectra))
  #println("Median Flux Error: ",median(abs(Error)))
  #println("Median SNR: ", median(Noisy_Spectra./(Error)))
  #println("------------------------")
  #println("Mean Flux: ",mean(Noisy_Spectra))
  #println("Mean Flux Error: ",mean(abs(Error)))
  #println("Mean SNR: ", mean(Noisy_Spectra./(Error)))
  #println("------------------------")
  #println("\n Done setting up unit test files! \n")
end

lvl=2.0
DSNR=5.0

FILES_ARR=["UnitTests/UT_Wavelengths.csv","UnitTests/UT_Spectra.csv","UnitTests/UT_Spectra_Error.csv","data/rvm_dates.csv","data/arp151.b.dat"] #Data files to load.
wavelengths=FILES_ARR[1];
spectra = FILES_ARR[2];
errspectra = FILES_ARR[3];
dates = FILES_ARR[4];
continuum = FILES_ARR[5];
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum);

Pars=init_Params();
Pars.num_tdf_times=50
Mats=Gen_Mats(DATA,Pars);

GenData(Pars,Mats,lvl,DSNR)
wavelengths=FILES_ARR[1];
spectra = FILES_ARR[2];
errspectra = FILES_ARR[3];
dates = FILES_ARR[4];
continuum = FILES_ARR[5];
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum);


spread=10.0
mus = logspace(0.0,6.0,10)
Res = zeros(length(mus),3)
#Find best mu for Tik initialization under given data.
for i in 1:length(mus)
  println(i)
  m=mus[i]
  tsol=gen_tiksol(Pars,Mats,DATA;scale=1.0,mu_smoo=m,plotting=false,save=false);
  chi2r=Chi2(Mats.H*tsol,DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
  diff= sum(abs(tsol-vdm))/length(vdm) #REDUCED DIFFERENCE ie RESIDUAL
  Res[i,1]=m
  Res[i,2]=chi2r
  Res[i,3]=diff
end
writecsv("Res.csv",Res)
