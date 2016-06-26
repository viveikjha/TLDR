include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("DataImportNEW.jl")

include("GenMatrices.jl")
using PyPlot

#FILES_ARR=["UnitTests/UT_Wavelengths.csv","UnitTests/UT_Spectra.csv","UnitTests/UT_Spectra_Error.csv","data/rvm_dates.csv","data/arp151.b.dat"] #Data files to load.
FILES_ARR=["data/rvm_wavelengths_trimmed.csv","data/rvm_fluxes_trimmed.csv", "data/rvm_errfluxes_trimmed.csv","data/rvm_dates.csv","data/arp151.b.dat"]

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
Pars.nits=200
Pars.num_tdf_times=50 #This is the default
Mats=Gen_Mats(DATA,Pars)

msmo = 1000.0
ml1 = 100.0
mspe = 10.0     #GOES WITH V
mtem = 10.0     #GOES WITH T

pz=1.0e8
#pp=1.0e12
pn=1.0e8
pv=1.0e6
pt=1.0e6

K=HOT_LAUNCH(DATA,Mats,Pars;mu_smoo=msmo,mu_spec=mspe,mu_temp=mtem,mu_l1=ml1,scale=1.0,nits=Pars.nits,Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true, rhoN=pn, rhoZ=pz, rhoV=pv,rhoT=pt); #RHOS: rhoZ=pz,rhoN=pn,rhoP=pp, rhoV=pv,rhoT=pt
