#push!(LOAD_PATH,"/home/manderson/Research/TLDR/")
#push!(LOAD_PATH,"/home/manderson/TLDR/")
#push!(LOAD_PATH,"/Users/manderson/Software/ReverbMap/JuliaVersions/TLDR")
push!(LOAD_PATH,"/Users/manderson/TLDR")

using PyPlot
using JLD
using RMLib
using RMLibMore
using RMTypes
using DataImportNEW
#using DataImport

using GenMatrices
bpf="box_20x50/" #box 20lams 50tau
prefix=bpf
FILES_ARR=[string(prefix,"Wavelengths.csv"),string(prefix,"Spectra.csv"),string(prefix,"Spectra_Error.csv"),"data/rvm_dates.csv","data/arp151.b.dat"]

#FILES_ARR=["UnitTests/UT_Wavelengths.csv","UnitTests/UT_Spectra.csv","UnitTests/UT_Spectra_Error.csv","data/rvm_dates.csv","data/arp151.b.dat"] #Data files to load.
#FILES_ARR=["data/rvm_wavelengths_trimmed.csv","data/rvm_fluxes_trimmed.csv", "data/rvm_errfluxes_trimmed.csv","data/rvm_dates.csv","data/arp151.b.dat"]

wavelengths=FILES_ARR[1]
spectra = FILES_ARR[2]
errspectra = FILES_ARR[3]
dates = FILES_ARR[4]
continuum = FILES_ARR[5]

#NEW DATASET? IMPORT DATA FILES
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)
save_data(string(prefix,"TLDR_data.jld"),DATA)
Pars= init_Params()
Pars.directory=prefix
Mats=Gen_Mats(DATA,Pars)
save_vars(string(prefix,"TLDR_vars.jld"),Mats,Pars)
#SAME DATA, DIFFERENT RUN? LOAD DATA AND VARIABLES
#DATA = load_data("TLDR_data.jld")
#Pars,Mats=load_vars("TLDR_vars.jld")
print_with_color(:green,"beginning reconstruction\n")
scale=1.0e0


Pars.num_tdf_times=50 #This is the default

#min=0.0
#max=20.0
#stepsize=(max-min)/(Pars.num_tdf_times-1)
#collect(1.0:((50-1.0)/(50-1)):50)
#Pars.tdf_times=collect(min:stepsize:max)
#println(Pars.tdf_times)

Fit=init_fit()
Fit.pp=100.0
Fit.mtem = 10.0   #GOES WITH T
Fit.pt = 100.0
Fit.mspe = 10.0    #GOES WITH V
Fit.pv=100.0
Fit.ml1 = 1000.0
Fit.pn=100.0
Fit.msmo = 1000.0
Fit.pz=100.0

#Fit.TI=1629750.8 #Box
Fit.TI=1000.0 #Box
#Fit.TI=570.0 #Ring

Pars.nits=50000.0
K=HOT_LAUNCH(DATA,Mats,Pars,Fit;scale=1.0,nits=Pars.nits,Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true,RecD=true); #RHOS: rhoZ=pz,rhoN=pn,rhoP=pp, rhoV=pv,rhoT=pt
vdm=copy(K.vdm)
writecsv(string(prefix,"RevMapResult.csv"),vdm)
println("wrote result to RevMapResult.csv")
