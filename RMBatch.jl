push!(LOAD_PATH,"/home/manderson/Research/TLDR/")
push!(LOAD_PATH,"/home/manderson/TLDR/")
push!(LOAD_PATH,"/Users/manderson/Software/ReverbMap/JuliaVersions/TLDR")

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
nps=2
mus=logspace(1.0,10,nps)
ps=logspace(1.0,10,nps)
pp=1.0e1
count=0
for msmo in mus # msmo
  for ml1 in mus #ml1
    for mspe in mus
      for mtem in mus
        for pz in ps
            for pn in ps
              for pv in ps
                for pt in ps
                  count+=1
                  DATA = load_data("TLDR_data.jld")
                  Pars,Mats=load_vars("TLDR_vars.jld")
                  Fit=init_fit()
                  Fit.msmo = msmo
                  Fit.pz=pz
                  Fit.ml1 = ml1
                  Fit.pn=pn
                  Fit.mspe = mspe
                  Fit.pv=pv
                  Fit.mtem =mtem
                  Fit.pt=pt
                  Fit.pp=pp
                  Fit.TI=1629750.8
                  fname=string("Batch/",string(count),".png")
                  K=HOT_LAUNCH(DATA,Mats,Pars,Fit;scale=1.0,nits=Pars.nits,Plot_Live=false,Plot_Final=false,RepIt=false,RepF=false);
                  figure()
                  imshow((K.vdm),aspect="auto",origin="lower",interpolation="None")
                  subplots_adjust(right=0.25)
                  s1=string("Msmo: ",msmo, " RhoS: ",pz )
                  s2=string("Ml1: ",ml1, " RhoN: ",pn )
                  s3=string("Mspe: ",mspe, " RhoV: ",pv )
                  s4=string("Mtem: ",msmo, " RhoT: ",pt )
                  s5=string("RhoP: ", pp, " Chi2: ", Pars.chi2)
                  text(30,0,s1,fontsize=14)
                  text(30,10,s2,fontsize=14)
                  text(30,20,s3,fontsize=14)
                  text(30,30,s4,fontsize=14)
                  text(30,40,s5,fontsize=14)
                  savefig(fname)
                end
              end
            end
          end
        end
      end
    end
  end
end
