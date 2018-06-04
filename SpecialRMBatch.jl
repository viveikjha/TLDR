push!(LOAD_PATH,"/home/manderson/TLDR")

using PyPlot
using JLD
using RMLib
using RMLibMore
using RMTypes
using DataImportNEW
#using DataImport
using GenMatrices

mastercount=0
pb="10x10/"
#names=["box","checkerboard","circle","diagonal","diagonalinverted","halfbottom","halfleft","halfright","halftop","horizontalstripe","invertedbox","invertedhorizontalstripe","invertedverticalstripe","lowertri","reverseddiagonal","reverseddiagonalinverted","ring","uppertri","verticalstripe"]
#names=["reverseddiagonal","reverseddiagonalinverted"]
names=["gradsh"]
nps=6
#m1=logspace(0.0,3.0,nps)
m1=linspace(1.0,150.0,nps)
m2=linspace(18.0,25.0,nps)
m3=linspace(85.0,95.0,nps)
m4=linspace(1.0,10,nps)
#m2=logspace(2.0,2.0,nps)
#m3=logspace(2.0,2.5,nps)
count=0
totalcount=length(names)*nps^4
mspe=88.0
mtem=1.0
ml1=133.0
msmo=143.0
pp=1.0
clim=4.0
pv=1.0
pn=1.0
pz=1.0
pt=1.0

for i in range(1,length(names))
    name=names[i]
    println(name)
	#fname=string(pb,name,"/",name,".csv")
	#vdm=readcsv(fname)
    bpf=string(pb,name,"/")
    filenames=["UT_Wavelengths.csv","UT_Spectra.csv","UT_Spectra_Error.csv"]
    FILES_ARR=[string(bpf,filenames[1]),string(bpf,filenames[2]),string(bpf,filenames[3]),"data/rvm_dates.csv","data/arp151.b.dat"]

    wavelengths=FILES_ARR[1]
    spectra = FILES_ARR[2]
    errspectra = FILES_ARR[3]
    dates = FILES_ARR[4]
    continuum = FILES_ARR[5]

    #NEW DATASET? IMPORT DATA FILES
    DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)
    Pars= init_Params()
    Pars.nits=3000
    Pars.alpha=1.2
    Pars.num_tdf_times=10
    Pars.threshold=1.0e-5

    Mats=Gen_Mats(DATA,Pars)

    #SAVE DATASET
    save_data("TLDR_data.jld",DATA)
    save_vars("TLDR_vars.jld",Mats,Pars)


    #SAME DATA, DIFFERENT RUN? LOAD DATA AND VARIABLES
    #DATA = load_data("TLDR_data.jld")
    #Pars,Mats=load_vars("TLDR_vars.jld")

    for msmo in m1 #msmo
      for ml1 in m1 #ml1
        for mspe in m1
          for mtem in m1
              mastercount+=1
              count+=1
              if count > 684
                  DATA = load_data("TLDR_data.jld")
                  Pars,Mats=load_vars("TLDR_vars.jld")
                  #WILL NEED A NEW IMPORT
                  Fit=init_fit()
                  Fit.waves=true
                  Fit.msmo = msmo
                  Fit.pz=pz
                  Fit.ml1 = ml1
                  Fit.pn=pn
                  Fit.mspe = mspe
                  Fit.pv=pv
                  Fit.mtem =mtem
                  Fit.pt=pt
                  Fit.pp=pp
                  Fit.TI=0.1
                  #Fit.fast=true
                  K,Fit=HOT_LAUNCH(DATA,Mats,Pars,Fit;nits=Pars.nits,Plot_Live=false,Plot_Final=false,RepIt=false,RepF=false);
                  if Pars.conflag==1 && Pars.chi2<clim#converged
                    col=:green
                  elseif Pars.conflag==2 #diverged
                    col=:red
                  else
                    col=:blue #did not converge
                  end
                  print_with_color(col,string("#",mastercount, " of ", totalcount, " chi2: ", round(Pars.chi2,4), " Iterations: ", Pars.it, "\n"))
                  if Pars.chi2 < clim
                    repname=string(bpf,"batch/",string(count),".txt")
                    s1=string("Msmo: ",msmo, " RhoS: ",Fit.pz )
                    s2=string("Ml1: ",ml1, " RhoN: ",Fit.pn )
                    s3=string("Mspe: ",mspe, " RhoV: ",Fit.pv )
                    s4=string("Mtem: ",mtem, " RhoT: ",Fit.pt )
                    s5=string("RhoP: ", pp, " Chi2: ", Pars.chi2, " Con: ", Pars.conflag)
                    s6=string(" Con: ", Pars.conflag)
                    writedlm(repname,[s1,s2,s3,s4,s5,s6],"\n")
                    fname=string(bpf,"batch/",string(count),".csv")
                    writecsv(fname,K.vdm)
                  end #endif
              end
            end
            end
          end
        end
      end
