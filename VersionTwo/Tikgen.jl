push!(LOAD_PATH,"/home/manderson/Research/TLDR/VersionTwo/")
using RMLib
using RMTypes
using DataImportNEW
using RMLibMore
using GenMatrices
using PyPlot

FILES_ARR=["../UnitTests/UT_Wavelengths.csv","../UnitTests/UT_Spectra.csv","../UnitTests/UT_Spectra_Error.csv","../data/rvm_dates.csv","../data/arp151.b.dat"] #Data files to load.
#FILES_ARR=["data/rvm_wavelengths_trimmed.csv","data/rvm_fluxes_trimmed.csv", "data/rvm_errfluxes_trimmed.csv","data/rvm_dates.csv","data/arp151.b.dat"]

wavelengths=FILES_ARR[1]
spectra = FILES_ARR[2]
errspectra = FILES_ARR[3]
dates = FILES_ARR[4]
continuum = FILES_ARR[5]
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)

scale=1.0e0
DATA.L=scale*(DATA.L)
DATA.EL=scale*(DATA.EL)
DATA.continuum_flux=scale*DATA.continuum_flux
DATA.continuum_error_flux=scale*DATA.continuum_error_flux

Pars= init_Params()
Pars.nits=100
Pars.num_tdf_times=50 #This is the default
min=0.0
max=20.0
stepsize=(max-min)/(Pars.num_tdf_times-1)
#collect(1.0:((50-1.0)/(50-1)):50)
Pars.tdf_times=collect(min:stepsize:max)
Mats=Gen_Mats(DATA,Pars)

msmo = 100.0  #was good.
pz=1.0e2

gen_tiksol(Pars,Mats,DATA;mu_l2=msmo,plotting=true,save=false)
