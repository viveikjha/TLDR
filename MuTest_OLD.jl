#This Script is to run TLDR with various mu values utilizing the hot launch in dev.jl.
#This mode of launch uses a single read on the data and a single setup of the computational
#matrices used in TLDR. This should reduce the setup time between each run of TLDR.



include("dev.jl")

FILES_ARR=["UnitTests/UT_Wavelengths.csv","UnitTests/UT_Spectra.csv","UnitTests/UT_Spectra_Error.csv","data/rvm_dates.csv","data/arp151.b.dat"] #Data files to load.
#FILES_ARR=["data/rvm_wavelengths_trimmed.csv","data/rvm_fluxes_trimmed.csv", "data/rvm_errfluxes_trimmed.csv","data/rvm_dates.csv","data/arp151.b.dat"]

wavelengths=FILES_ARR[1]
spectra = FILES_ARR[2]
errspectra = FILES_ARR[3]
dates = FILES_ARR[4]
continuum = FILES_ARR[5]
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)

scale=1.0
DATA.L=scale*(DATA.L)
DATA.EL=scale*(DATA.EL)
DATA.continuum_flux=scale*DATA.continuum_flux
DATA.continuum_error_flux=scale*DATA.continuum_error_flux

Pars= init_Params()
Pars.nits=50
Mats=Gen_Mats(DATA,Pars)

msmo = 1.0e6
mspe =1.0e1 #GOES WITH V
ml1 =5.0e4
mtem =1.0e1 #GOES WITH T

pz=1.0e10
pp=1.0e5
pn=1.0e9
pv=1.0e5
pt=1.0e5

#μ/ρ=flx_level*(1/20) seems to work quite well for the gradients. Keep the ratio constant as weight changes.

println("For N: μ/ρ= ",ml1/pn)
println("For T: μ/ρ= ",mtem/pt)
println("For V: μ/ρ= ",mspe/pv)
LAUNCH(DATA,Mats,Pars;mu_smoo=msmo,mu_spec=mspe,mu_temp=mtem,mu_l1=ml1,scale=1.0,nits=500,Tvdm="",Plot_Live=false,Plot_Final=true,RepIt=true,RepF=true,rhoZ=pz,rhoN=pn,rhoP=pp, rhoV=pv,rhoT=pt);
