push!(LOAD_PATH,"/home/matander/TLDR")

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
Pars.nits=2000
Pars.alpha=1.2

Pars.num_tdf_times=50 #This is the default
Mats=Gen_Mats(DATA,Pars)

#SAVE DATASET
save_data("TLDR_data.jld",DATA)
save_vars("TLDR_vars.jld",Mats,Pars)


#SAME DATA, DIFFERENT RUN? LOAD DATA AND VARIABLES
#DATA = load_data("TLDR_data.jld")
#Pars,Mats=load_vars("TLDR_vars.jld")
nps=5
m1=logspace(2.0,5.0,nps)
#m2=logspace(2.0,2.0,nps)
#m3=logspace(2.0,2.5,nps)


count=0

mspe=1.0
mtem=1.0
ml1=1.0
msmo=1.0
pp=1.0
clim=4.0
pv=100.0
pn=100.0
pz=100.0
pt=100.0

for msmo in m1 #msmo
  for ml1 in m1 #ml1
    for mspe in m1
      for mtem in m1
          count+=1
          DATA = load_data("TLDR_data.jld")
          Pars,Mats=load_vars("TLDR_vars.jld")
          #WILL NEED A NEW IMPORT

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
          Fit.TI=500.0
          K=HOT_LAUNCH(DATA,Mats,Pars,Fit;nits=Pars.nits,Plot_Live=false,Plot_Final=false,RepIt=false,RepF=false);
          if Pars.conflag==1 && Pars.chi2<clim#converged
            col=:green
          elseif Pars.conflag==2 #diverged
            col=:red
          else
            col=:blue #did not converge
          end
          print_with_color(col,string(count, " chi2: ", round(Pars.chi2), " Iterations: ", Pars.it, "\n"))
          if Pars.chi2 < clim
            fname=string("Batch/",string(count),".png")
            figure()
            imshow((K.vdm),aspect="auto",origin="lower",interpolation="None")
            subplots_adjust(right=0.25)
            s1=string("Msmo: ",msmo, " RhoS: ",pz )
            s2=string("Ml1: ",ml1, " RhoN: ",pn )
            s3=string("Mspe: ",mspe, " RhoV: ",pv )
            s4=string("Mtem: ",mtem, " RhoT: ",pt )
            s5=string("RhoP: ", pp, " Chi2: ", Pars.chi2, " Con: ", Pars.conflag)
            s6=string(" Con: ", Pars.conflag)

            text(30,0 ,s1,fontsize=14)
            text(30,10,s2,fontsize=14)
            text(30,20,s3,fontsize=14)
            text(30,30,s4,fontsize=14)
            text(30,40,s5,fontsize=14)
            text(30,50,s6,fontsize=14)

            savefig(fname)
            close("all")
            fname=string("Batch/vdms/",string(count),".csv")
            writecsv(fname,K.vdm)
          end #endif
        end
        end
      end
    end
  #end
