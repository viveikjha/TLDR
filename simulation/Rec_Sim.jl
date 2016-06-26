#RUN THIS SCRIPT WITH julia -p 3 TO USE MULTIPLE PROCESSORS.

include("../RMLib.jl")
include("../RMTypes.jl")
include("../DataImport.jl")
include("../GenMatrices.jl")
include("../DataImportNEW.jl")
include("../dev.jl")
using PyPlot
files=["../data/rvm_wavelengths.csv","../simulation/Sim_Spectra.csv","../simulation/Sim_Error.csv","../data/rvm_dates.csv","../data/arp151.b.dat"]

wavelengths=files[1];
spectra = files[2];
errspectra = files[3];
dates = files[4];
continuum = files[5];
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum);

Pars= init_Params();
Pars.nits=200
Pars.num_tdf_times=50 #This is the default
Mats=Gen_Mats(DATA,Pars);


#println(r)
mu_smoo=1.0e3
mu_spec=0.001
mu_temp=0.001
mu_l1=0.001

pz=1.0e8
#pp=1.0e12
pn=1.0e8
pv=1.0e6
pt=1.0e6

rec=HOT_LAUNCH(DATA,Mats,Pars;mu_smoo=msmo,mu_spec=mspe,mu_temp=mtem,mu_l1=ml1,scale=1.0,nits=Pars.nits,Plot_Live=false,Plot_Final=true,RepIt=false,RepF=false, rhoN=pn, rhoZ=pz, rhoV=pv,rhoT=pt); #RHOS: rhoZ=pz,rhoN=pn,rhoP=pp, rhoV=pv,rhoT=pt


println("Complete.")
