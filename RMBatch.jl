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
wavelengths=FILES_ARR[1]
spectra = FILES_ARR[2]
errspectra = FILES_ARR[3]
dates = FILES_ARR[4]
continuum = FILES_ARR[5]

#NEW DATASET? IMPORT DATA FILES
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)
Pars= init_Params()
Pars.nits=100
Pars.num_tdf_times=50 #This is the default
Mats=Gen_Mats(DATA,Pars)

#SAVE DATASET
save_data("TLDR_data.jld",DATA)
save_vars("TLDR_vars.jld",Mats,Pars)


#SAME DATA, DIFFERENT RUN? LOAD DATA AND VARIABLES
DATA = load_data("TLDR_data.jld")
Pars,Mats=load_vars("TLDR_vars.jld")
nps=5
mus=logspace(1.0,10,nps)
ps=logspace(1.0,10,nps)

for i
  for j
    for k
      for

end
Fit=init_fit()
Fit.msmo = 1.0e5
Fit.pz=1.0e8
Fit.ml1 = 1.0e3
Fit.pn=1.0e6
Fit.mspe = 1.0e0    #GOES WITH V
Fit.pv=1.0e0
Fit.mtem = 1.0e0    #GOES WITH T
Fit.pt=1.0e4
Fit.pp=1.0e2
Fit.TI=1629750.8
