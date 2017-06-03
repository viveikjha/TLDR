#include("RMLibMore.jl")
#include("RMLib.jl")
using PyPlot
using JLD
push!(LOAD_PATH,"/home/manderson/TLDR/")
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
mu = 1000.0
Pars=init_Params();
Pars.num_tdf_times=ntimes
Mats=Gen_Mats(DATA,Pars);


save("TLDR_VARS.jld","Mats",Mats,"Pars",Pars)

function save_vars(fname,Mats,Pars)
  save(fname,"Mats",Mats,"Pars",Pars)
end

function load_vars(fname)
  P,M=load(fname)
  Pars=P[2]
  Mats=M[2]
  Pars,Mats
end
