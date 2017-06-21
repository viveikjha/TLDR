#include("RMLibMore.jl")
#include("RMLib.jl")

push!(LOAD_PATH,"/Users/manderson/Software/ReverbMap/JuliaVersions/TLDR")
push!(LOAD_PATH,"/home/manderson/TLDR")

using RMTypes
using RMLib
using RMLibMore
using GenMatrices
using DataImportNEW

using PyPlot

FILES_ARR=["UnitTests/UT_Wavelengths.csv","UnitTests/UT_Spectra.csv","UnitTests/UT_Spectra_Error.csv","data/rvm_dates.csv","data/arp151.b.dat"] #Data files to load.
wavelengths=FILES_ARR[1];
spectra = FILES_ARR[2];
errspectra = FILES_ARR[3];
dates = FILES_ARR[4];
continuum = FILES_ARR[5];
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum);
vdm=readcsv("UnitTests/UT_vdm.csv");

Pars=init_Params();
Pars.num_tdf_times=50
Mats=Gen_Mats(DATA,Pars);

mus = logspace(0.0,10.0,20)
Res = zeros(length(mus),3)
println("Dimensions: ", size(DATA.L))
#Find best mu for Tik initialization under given data.
for i in 1:length(mus)
  println(i)
  m=mus[i]
  tsol=gen_tiksol(Pars,Mats,DATA;scale=1.0,mu_smoo=m,plotting=true,save=false);
  chi2r=Chi2(Mats.H*tsol,DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
  diff= sum(abs(tsol-vdm)) #RESIDUAL
  Res[i,1]=m
  Res[i,2]=chi2r
  Res[i,3]=diff
end

println(Res)


lvl=0.1
DSNR=4.0
ax=PyPlot.gca()
ax[:set_ylim](0,200)
semilogx(Res[:,1],Res[:,2])
xlabel("Mu")
ylabel("Chi2")
tstr=string("Lvl= ", string(lvl), " SNR= ", string(DSNR))
title(tstr)
show()
fname=string("TIkSols/Res_L_", string(lvl),"_SNR_",string(DSNR),".csv")
writecsv(fname,Res)
