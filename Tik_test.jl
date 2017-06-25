#include("RMLibMore.jl")
#include("RMLib.jl")
using PyPlot
push!(LOAD_PATH,"/home/manderson/TLDR/")
push!(LOAD_PATH,"/home/manderson/Research/TLDR/")
using RMTypes
using GenMatrices
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
levels=[1.0,0.5,0.1,0.05]
mu = 1629750.8
Par=init_Params();
Par.num_tdf_times=ntimes
Mats=Gen_Mats(DATA,Par);
Par.num_tdf_times=ntimes
#tsol=gen_tiksol(Par,Mats,DATA;scale=1.0,mu_smoo=mu,plotting=false,save=false);
a=@elapsed gen_tiksol(Par,Mats,DATA;scale=1.0,mu_smoo=mu,plotting=true,save=false);
println(a)
#figure()
#imshow(tsol)
#show()
