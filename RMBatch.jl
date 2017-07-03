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
Pars.nits=500
Pars.num_tdf_times=50 #This is the default
Mats=Gen_Mats(DATA,Pars)

#SAVE DATASET
save_data("TLDR_data.jld",DATA)
save_vars("TLDR_vars.jld",Mats,Pars)


#SAME DATA, DIFFERENT RUN? LOAD DATA AND VARIABLES
DATA = load_data("TLDR_data.jld")
Pars,Mats=load_vars("TLDR_vars.jld")
nps=4
#m1=logspace(1.5,2.5,nps)
#m2=logspace(0.0,2.0,nps)
#m3=logspace(2.5,3.5,nps)
ps=1.0e5

count=0

ps=linspace(1.0,1.0e4,nps)
mspe=1.0
mtem=1.0
ml1=1.0
msmo=1.0
pv=ps
pt=ps
pz=ps
pn=ps
pp=1.0e5
clim=4.0
#pz=1.0e5
#for msmo in m2 # msmo
#  for ml1 in m2 #ml1
#    for mspe in m3
#      for mtem in m1
for pz in ps # msmo
  for pn in ps #ml1
    for pv in ps
      for pt in ps
        for pp in ps
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

            text(30,0,s1,fontsize=14)
            text(30,10,s2,fontsize=14)
            text(30,20,s3,fontsize=14)
            text(30,30,s4,fontsize=14)
            text(30,40,s5,fontsize=14)
            text(30,50,s6,fontsize=14)

            savefig(fname)
            close("all")
            fname=string("Batch/vdms/",string(count),".csv")
            writecsv(fname,K.vdm)
          end
        end
        end
      end
    end
  end
