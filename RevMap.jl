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
println(@elapsed DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum))
save_data("TLDR_data.jld",DATA)
#println( @elapsed DATA = load_data("TLDR_data.jld") )
scale=1.0e0

Pars= init_Params()
println(@elapsed Mats=Gen_Mats(DATA,Pars))

#Save Pars & Mats here.
save_vars("TLDR_vars.jld",Mats,Pars)
writecsv("H.csv",Mats.H)
#println(@elapsed Pars,Mats=load_vars("TLDR_vars.jld"))

Pars.nits=10

Pars.num_tdf_times=50 #This is the default

#min=0.0
#max=20.0
#stepsize=(max-min)/(Pars.num_tdf_times-1)
#collect(1.0:((50-1.0)/(50-1)):50)
#Pars.tdf_times=collect(min:stepsize:max)
#println(Pars.tdf_times)

msmo = 1.0e5
pz=1.0e8


ml1 = 1.0e3
pn=1.0e3


mspe = 1.0e4    #GOES WITH V
pv=1.0e8


mtem = 1.0e3    #GOES WITH T
pt=1.0e4
pp=1.0e2

K=HOT_LAUNCH(DATA,Mats,Pars;mu_smoo=msmo,mu_spec=mspe,mu_temp=mtem,mu_l1=ml1,scale=1.0,nits=Pars.nits,Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true, rhoN=pn, rhoZ=pz, rhoV=pv,rhoT=pt,rhoP=pp); #RHOS: rhoZ=pz,rhoN=pn,rhoP=pp, rhoV=pv,rhoT=pt
vdm=copy(K.vdm)
writecsv("RevMapResult.csv",vdm)
println("wrote result to RevMapResult.csv")
