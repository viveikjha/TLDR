include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("DataImportNEW.jl")
include("GenMatrices.jl")
using PyPlot

#n_cores = CPU_CORES #NUMBER OF CORES ON MACHINE
#addprocs(n_cores-1) #ADD WORKERS

FILES_ARR=["UnitTests/UT_Wavelengths.csv","UnitTests/UT_Spectra.csv","UnitTests/UT_Spectra_Error.csv","data/rvm_dates.csv","data/arp151.b.dat"] #Data files to load.
#FILES_ARR=["data/rvm_wavelengths_trimmed.csv","data/rvm_fluxes_trimmed.csv", "data/rvm_errfluxes_trimmed.csv","data/rvm_dates.csv","data/arp151.b.dat"]
#FILES_ARR=["data/rvm_wavelengths.csv","data/rvm_fluxes.csv", "data/rvm_errfluxes.csv","data/rvm_dates.csv","data/arp151.b.dat"]

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
Pars.nits=50
Mats=Gen_Mats(DATA,Pars)


msmo = 1.0
mspe =3100.0     #GOES WITH V
ml1 =10.0
mtem =1.0e2     #GOES WITH T

pz=50000.0
pp=16250.0
pn=50000.0
pv=38750.0
pt=16250.0

K=HOT_LAUNCH(DATA,Mats,Pars;mu_smoo=msmo,mu_spec=mspe,mu_temp=mtem,mu_l1=ml1,scale=1.0,nits=50,Tvdm="",Plot_Live=true,Plot_Final=true,RepIt=true,RepF=false,rhoZ=pz,rhoN=pn,rhoP=pp, rhoV=pv,rhoT=pt);
println("Done")
