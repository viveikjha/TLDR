#push!(LOAD_PATH,"/home/manderson/Research/TLDR/")
push!(LOAD_PATH,"/home/manderson/TLDR/")
#push!(LOAD_PATH,"/Users/manderson/Software/ReverbMap/JuliaVersions/TLDR")
#push!(LOAD_PATH,"/home/matander/TLDR")

using PyPlot
using JLD
using RMLib
using RMLibMore
using RMTypes
using DataImportNEW
#using DataImport

using GenMatrices


prefix="10x10/"
name="box"
bpf=string(prefix,name,"/")
filenames=["UT_Wavelengths.csv","UT_Spectra.csv","UT_Spectra_Error.csv"]
FILES_ARR=[string(bpf,filenames[1]),string(bpf,filenames[2]),string(bpf,filenames[3]),"data/rvm_dates.csv","data/arp151.b.dat"]

wavelengths=FILES_ARR[1]
spectra = FILES_ARR[2]
errspectra = FILES_ARR[3]
dates = FILES_ARR[4]
continuum = FILES_ARR[5]
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)
#save_data(string(prefix,"TLDR_data.jld"),DATA)
Pars= init_Params()
Pars.directory=prefix
Pars.nits=1000
Pars.num_tdf_times=10
Pars.alpha=1.2
Mats=Gen_Mats(DATA,Pars)
#save_vars(string(prefix,"TLDR_vars.jld"),Mats,Pars)
#SAME DATA, DIFFERENT RUN? LOAD DATA AND VARIABLES
#DATA = load_data("TLDR_data.jld")
#Pars,Mats=load_vars("TLDR_vars.jld")
print_with_color(:green,"beginning reconstruction\n")
scale=1.0e0

Pars.threshold=1.0e-4
Fit=init_fit()
Fit.waves=false   #Use wavelets
Fit.mtem = 1.0   #GOES WITH T
Fit.mspe = 1.0    #GOES WITH V
Fit.ml1 = 1.0
Fit.msmo = 1.0

Fit.pz=1.0
Fit.pp=1.0
Fit.pn=1.0
Fit.pv=1.0
Fit.pt = 1.0
#Fit.TI=1629750.8 #Box
Fit.TI=1.0 #Box
#Fit.TI=570.0 #Ring


tic()
K=HOT_LAUNCH(DATA,Mats,Pars,Fit;scale=1.0,nits=Pars.nits,Plot_Live=false,Plot_Final=false,RepIt=true,RepF=false,RecD=false); #RHOS: rhoZ=pz,rhoN=pn,rhoP=pp, rhoV=pv,rhoT=pt
toc()
vdm=copy(K.vdm)
writecsv(string("RevMapResult.csv"),vdm)
")
