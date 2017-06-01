#include("RMLibMore.jl")
#include("RMLib.jl")

push!(LOAD_PATH,"/home/manderson/Research/TLDR/")
using RMTypes
using RMLib
using RMLibMore
using DataImportNEW

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
#levels=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05]
levels=[1.0,0.5,0.1,0.05]
#mus = [1.0,10.0,100.0,500.0,1000.0,5000.0,10000.0,5.0e5,1.0e6]
mus = logspace(0.0,8.0,10)
Res = zeros(length(mus),3)
#println("Mu    Chi2r   Diff")
for i in 1:length(mus)
  m=mus[i]
  Par=init_Params();
  Par.num_tdf_times=ntimes
  Mats=Gen_Mats(DATA,Par);
  vdm=gen_UTD(DATA,Par,levels[3],nlams,ntimes) #GENERATES NEW DATA WITH DIFFERENT SIGNAL LEVEL
  Par=init_Params();
  Par.num_tdf_times=ntimes
  Mats=Gen_Mats(DATA,Par);
  tsol=gen_tiksol(Par,Mats,DATA;scale=1.0,mu_smoo=m,plotting=false,save=false);
  chi2r=Chi2(Mats.H*tsol,DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
  diff= sum(abs(tsol-vdm))/length(vdm) #REDUCED DIFFERENCE
  Res[i,1]=m
  Res[i,2]=chi2r
  Res[i,3]=diff
end
writecsv("Res.csv",Res)
