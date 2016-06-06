#This Script is to run TLDR with various mu values utilizing the hot launch in dev.jl.
#This mode of launch uses a single read on the data and a single setup of the computational
#matrices used in TLDR. This should reduce the setup time between each run of TLDR.


using PyPlot
include("RMLib.jl")

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
Pars.nits=50
Mats=Gen_Mats(DATA,Pars)

msmo = 1.0e6
mspe =1.0e-5     #GOES WITH V
ml1 =1.0e-5
mtem =1.0e-5     #GOES WITH T

pz=1.0e10
pp=1.0e5
pn=1.0e-5
pv=1.0e-5
pt=1.0e-5

#μ/ρ=flx_level*(1/20) seems to work quite well for the gradients. Keep the ratio constant as weight changes.

#println("For N: μ/ρ= ",ml1/pn)
#println("For T: μ/ρ= ",mtem/pt)
#println("For V: μ/ρ= ",mspe/pv)
vdm_act=readcsv("UnitTests/UT_vdm.csv")

mu_l2 = [1.0e3,5.0e3,1.0e4,2.0e4,4.0e4,5.0e4,6.0e4,8.0e4,1.0e5,2.0e5,4.0e5,]

chi2=zeros(size(mu_l2))
reg = zeros(size(mu_l2))
res = zeros(size(mu_l2))
con = zeros(size(mu_l2))
for i in 1:size(mu_l2)[1]
  Pars.conflag=false
  Pars.mu_smoo = mu_l2[i]

  K=HOT_LAUNCH(DATA,Mats,Pars;mu_smoo=mu_l2[i],mu_spec=mspe,mu_temp=mtem,mu_l1=ml1,scale=1.0,nits=2000,Tvdm="",Plot_Live=false,Plot_Final=false,RepIt=false,RepF=true,rhoZ=pz,rhoN=pn,rhoP=pp, rhoV=pv,rhoT=pt);
  chi2[i] = Chi2(Model(K.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
  reg[i]=regX(K,Pars)
  res[i] = sum(abs(K.vdm-vdm_act))
  if Pars.conflag == true
    con[i]=true
  end
end


outarr=zeros(size(mu_l2)[1],5)
outarr[:,1] = mu_l2
outarr[:,2] = chi2
outarr[:,3] = reg
outarr[:,4] = res
outarr[:,5] = con
writecsv("MuTests/Results.csv",outarr)

figure()
loglog(mu_l2,chi2,"k")
for i in 1:size(mu_l2)[1]
  if con[i] == true
    loglog(mu_l2[i],chi2[i],"g.")
  else
    loglog(mu_l2[i],chi2[i],"r.")
  end
end
xlabel("mu_l2")
ylabel("chi2")
show()
figure()
plot(mu_l2,reg)
xlabel("mu_l2")
ylabel("reg")
show()
figure()
loglog(res,chi2,"k")
for i in 1:size(mu_l2)[1]
  if con[i] == true
    loglog(res[i],chi2[i],"g.")
  else
    loglog(res[i],chi2[i],"r.")
  end
end
ylabel("chi2")
xlabel("res")
show()
