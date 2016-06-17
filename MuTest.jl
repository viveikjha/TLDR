#This Script is to run TLDR with various mu values utilizing the hot launch in dev.jl.
#This mode of launch uses a single read on the data and a single setup of the computational
#matrices used in TLDR. This should reduce the setup time between each run of TLDR.

using PyPlot
using LaTeXStrings
using PyCall
include("RMLib.jl")

H=readcsv("H.csv")
FILES_ARR=["UnitTests/UT_Wavelengths.csv","UnitTests/UT_Spectra.csv","UnitTests/UT_Spectra_Error.csv","data/rvm_dates.csv","data/arp151.b.dat"] #Data files to load.
#FILES_ARR=["data/rvm_wavelengths_trimmed.csv","data/rvm_fluxes_trimmed.csv", "data/rvm_errfluxes_trimmed.csv","data/rvm_dates.csv","data/arp151.b.dat"]

wavelengths=FILES_ARR[1]
spectra = FILES_ARR[2]
errspectra = FILES_ARR[3]
dates = FILES_ARR[4]
continuum = FILES_ARR[5]
DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)

scale=1.0
DATA.L=scale*(DATA.L)
DATA.EL=scale*(DATA.EL)
DATA.continuum_flux=scale*DATA.continuum_flux
DATA.continuum_error_flux=scale*DATA.continuum_error_flux

Pars= init_Params()
Mats=Gen_Mats(DATA,Pars)
vdm_act=readcsv("UnitTests/UT_vdm.csv")
num=5




msmo = 1000000.0
ml1 = 1000.0
mspe = 100.0     #GOES WITH V
mtem = 100.0     #GOES WITH T


pz=1.0e8
#pp=1.0e12
pn=1.0e8
pv=1.0e6
pt=1.0e6


MUX = msmo
MUN = logspace(0,4,num)
MUT = logspace(0,3,num)
MUV = logspace(0,3,num)
PZ =  logspace(7,9,num)
PP =  1.0e5
PN =  logspace(7,9,num)
PV =  logspace(4,6,num)
PT =  logspace(4,6,num)

count=0
for muX in MUX
  for muN in MUN
    for muT in MUT
      for muV in MUV
        for pz in PZ
          for pn in PN
            for pv in PV
              for pt in PT
                for pp in PP
                  Pars= init_Params() #REINITIALIZE PARAMETERS
                  save=false
                  K=HOT_LAUNCH(DATA,Mats,Pars;mu_smoo=muX,mu_spec=muV,mu_temp=muT,mu_l1=muN,nits=100,Tvdm="",Plot_Live=true,Plot_Final=false,RepIt=true,RepF=false,rhoZ=pz,rhoN=pn,rhoP=pp, rhoV=pv,rhoT=pt);
                  if Pars.l1N_state==true && Pars.l1T_state==true && Pars.l1V_state==true && Pars.chi2 < 500 #ONLY SAVE IF ALL THREE HAVE NOT FAILED
                    save=true
                    count+=1
                    chi2 = Chi2(Model(K.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
                    resn = ell2norm(H*K.vdm-DATA.L)
                    line = string(count,",",Pars.conflag,",",muX,",",muV,",",muT,",",muN,",",pz,",",pp,",",pn,",",pv,",",pt,",",chi2,",",resn,"\n")
                    title =  string("MuTests/",count,".csv")
                    writecsv(title,K.vdm)
                  end
                  if save==true
                    f=open("MuTests/test.txt","a")
                    write(f,line)
                    close(f)
                    println(string(Pars.conflag,",",muX,",",muV,",",muT,",",muN,",",pz,",",pp,",",pn,",",pv,",",pt,",",muT/pt, ",",muV/pv))
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end
