push!(LOAD_PATH,"/home/manderson/Research/TLDR/")
push!(LOAD_PATH,"/home/manderson/TLDR/")

using PyPlot
using JLD
using RMLib
using RMLibMore
using RMTypes
using DataImportNEW
#using DataImport
using GenMatrices
FILES_ARR=["UnitTests/UT_Wavelengths.csv","UnitTests/UT_Spectra.csv","UnitTests/UT_Spectra_Error.csv","data/rvm_dates.csv","data/arp151.b.dat"] #Data files to load.
#FILES_ARR=["data/rvm_wavelengths_trimmed.csv","data/rvm_fluxes_trimmed.csv", "data/rvm_errfluxes_trimmed.csv","data/rvm_dates.csv","data/arp151.b.dat"]

wavelengths=FILES_ARR[1]
spectra = FILES_ARR[2]
errspectra = FILES_ARR[3]
dates = FILES_ARR[4]
continuum = FILES_ARR[5]

#NEW DATASET? IMPORT DATA FILES
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)
save_data("TLDR_data.jld",DATA)
Pars= init_Params()
Mats=Gen_Mats(DATA,Pars)
save_vars("TLDR_vars.jld",Mats,Pars)
#SAME DATA, DIFFERENT RUN? LOAD DATA AND VARIABLES
#DATA = load_data("TLDR_data.jld")
#Pars,Mats=load_vars("TLDR_vars.jld")
print_with_color(:green,"beginning reconstruction\n")
scale=1.0e0

Pars.nits=500

Pars.num_tdf_times=50 #This is the default

#min=0.0
#max=20.0
#stepsize=(max-min)/(Pars.num_tdf_times-1)
#collect(1.0:((50-1.0)/(50-1)):50)
#Pars.tdf_times=collect(min:stepsize:max)
#println(Pars.tdf_times)
rho=1.0e5
Fit=init_fit()
Fit.msmo = 1.0
#Fit.pz=1.0e5
Fit.pz=50000.0

Fit.ml1 = 10.0
#Fit.pn=1.0e5
Fit.pn=50000.0

Fit.mspe = 3100.0    #GOES WITH V
#Fit.pv=1.0e5
Fit.pv=38750.0

Fit.mtem = 100.0   #GOES WITH T
#Fit.pt=1.0e5
Fit.pt = 16250.0

Fit.pp=16250.0
Fit.TI=1629750.8


K=HOT_LAUNCH(DATA,Mats,Pars,Fit;scale=1.0,nits=Pars.nits,Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true,RecD=true); #RHOS: rhoZ=pz,rhoN=pn,rhoP=pp, rhoV=pv,rhoT=pt
vdm=copy(K.vdm)
writecsv("RevMapResult.csv",vdm)
println("wrote result to RevMapResult.csv")
