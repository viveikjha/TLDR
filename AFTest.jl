#include("RMLibMore.jl")
#include("RMLib.jl")
using PyPlot
push!(LOAD_PATH,"/home/manderson/Research/TLDR/")
using RMTypes
using ArrayFire
using GenMatrices
using RMLib
using RMLibMore
using DataImportNEW

function gen_tiksol_AF(Par,Mats,DATA;scale=1.0,mu_smoo=40.0)
  #Par,Mats,DATA=getdata(f)
  #println("µ: ",mu_smoo)
  DATA.L = scale.*DATA.L
  DATA.EL = scale.*DATA.EL
  Pars= init_Params()
  Mats = Gen_Mats(DATA,Pars)

  #end
  Pars.mu_smoo=copy(mu_smoo)

  vdm = zeros(Pars.num_tdf_times,DATA.num_lines)

  W=AFArray(Mats.W)
  HT=AFArray(Mats.HT)
  H=AFArray(Mats.H)
  Gamma=AFArray(Mats.Gammatdf)
  L=AFArray(DATA.L)
  for l = 1:DATA.num_lines
    Ws=reshape(W[l,:,:],size(W[l,:,:])[2],size(W[l,:,:])[3])
    Q=HT*Ws*H+Pars.mu_smoo*Gamma
    B=HT*Ws*L
    G=inv(Q)*B
    vdm[:,l] = Array(G[:,l])
    #vdm =vdm.*(vdm.>0.0) #FILTER OUT NEGATIVE VALUES
    #println(typeof(vdm))
  end
  vdm
end


FILES_ARR=["UnitTests/UT_Wavelengths.csv","UnitTests/UT_Spectra.csv","UnitTests/UT_Spectra_Error.csv","data/rvm_dates.csv","data/arp151.b.dat"] #Data files to load.

wavelengths=FILES_ARR[1];
spectra = FILES_ARR[2];
errspectra = FILES_ARR[3];
dates = FILES_ARR[4];
continuum = FILES_ARR[5];
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum);

nlams = 20
ntimes= 50
spread=10.0
lvl=1.00
mu = 1.0e3
Par=init_Params();
Par.num_tdf_times=ntimes
Mats=Gen_Mats(DATA,Par);
Par.num_tdf_times=ntimes
#tsol=gen_tiksol_AF(Par,Mats,DATA;scale=1.0,mu_smoo=mu);
a=@elapsed gen_tiksol_AF(Par,Mats,DATA;scale=1.0,mu_smoo=mu)
println(a)
#println(size(tsol))
#figure()
#imshow(tsol)
#show()
