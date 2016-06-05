include("RMLib.jl")
include("RMLibMore.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("DataImportNEW.jl")
include("GenMatrices.jl")

#using PyPlot

function getdata(FILES_ARR)
  #FILES_ARR = ["data/rvm_wavelengths_trimmed.csv","data/rvm_fluxs_trimmed.csv","data/rvm_errfluxes_trimmed.csv","data/rvm_dates.csv","data/arp151.b.dat"]
  wavelengths=FILES_ARR[1]
  spectra = FILES_ARR[2]
  errspectra = FILES_ARR[3]
  dates = FILES_ARR[4]
  continuum = FILES_ARR[5]
  DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)

  #data_report(DATA)
  Pars = init_Params()
  Mats=Gen_Mats(DATA,Pars)

  Pars,Mats,DATA
end

#function gen_tiksol(DATA,Pars,Mats;scale=1.0,mu_smoo=40.0)
function gen_tiksol(f;scale=1.0,mu_smoo=40.0,plotting=true,save=false)
  Par,Mats,DATA=getdata(f)
  #println("µ: ",mu_smoo)
  DATA.L = scale.*DATA.L
  DATA.EL = scale.*DATA.EL
  Pars= init_Params()
  Mats = Gen_Mats(DATA,Pars)

  #end
  Pars.mu_smoo=mu_smoo

  vdm = zeros(Pars.num_tdf_times,DATA.num_lines)
  for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
    W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
    Q = Mats.HT * W_slice * Mats.H +  (Pars.mu_smoo)*Mats.Gammatdf #INCLUCES L1 NORM ON X
    B = Mats.HT* W_slice * DATA.L
    vdm[:,l] = Q\B[:,l]
  end
  vdm =vdm.*(vdm.>0.0) #FILTER OUT NEGATIVE VALUES
  if save==true
    writecsv("devsol.csv",vdm)
  end
  if plotting==true
    figure()
    imshow(vdm,origin="lower",cmap="Greens",interpolation="None",aspect="auto")
    xlabel("Spectral Channel")
    ylabel("Delay Time")
    titlestring=string("Tikhonov Image for µ=",Pars.mu_smoo, " and scale=",scale)
    title(titlestring)
    show()
  end
  vdm
end


#=--------------------------------------------------=#
#================= TLDR HOT LAUNCHER ====================#
#=--------------------------------------------------=#
#1. FILES REQUIRED - FILES_ARR = [WAVELENGTHS,SPECTRA,ERRSPECTRA,DATES,CONTINUUM]
#2. IS THIS A TEST? TRUE VDM FILE? -> OPTIONAL
#3. mus: mu_smoo required. others optional?

#Tvdm file: Optional file imput of the real tdf. Used for testing puroposes.
#mu_spec, mu_l1, and mu_temp flags are options for providing individual regularization weights
#otherwise, the values are based on those for mu_smoo.
#Plot_Live option shows the active reconstruction every so many iterations.
#Plot_Final option shows the final plot that displays reconstruction images X,Z,V,T.
function LAUNCH(DATA,Mats,Pars;mu_smoo=40.0,mu_spec=false,mu_temp=false,mu_l1=false,scale=1.0,nits=50,Tvdm="",Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true,rhoZ=8000.0,rhoN=800.0,rhoP=800.0, rhoV=800.0,rhoT=800.0)
  #IMPORT DATA FROM FILES_ARR
	#wavelengths=FILES_ARR[1]
	#spectra = FILES_ARR[2]
	#errspectra = FILES_ARR[3]
	#dates = FILES_ARR[4]
	#continuum = FILES_ARR[5]
	#DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)

	DATA.L=scale*(DATA.L)
	DATA.EL=scale*(DATA.EL)
	DATA.continuum_flux=scale*DATA.continuum_flux
	DATA.continuum_error_flux=scale*DATA.continuum_error_flux

	#data_report(DATA)
	#Pars = init_Params()
	#Pars.num_tdf_times=50

	Mats=Gen_Mats(DATA,Pars)

	Pars.nits=nits

	#SET RECONSTRUCTION PARAMETERS
	scale=1.0
  if mu_temp != false
    Pars.mu_temp = mu_temp
  else
    Pars.mu_temp = 0.25*mu_smoo/scale
  end
  if mu_spec != false
    Pars.mu_spec = mu_spec
  else
    Pars.mu_spec = 0.25*mu_smoo/scale
  end
  if mu_l1 != false
    Pars.mu_l1 = mu_l1
  else
    Pars.mu_l1 = 0.25*mu_smoo/scale
  end
  Pars.mu_smoo=mu_smoo/scale^2
  tmp,P = TLDR(1.0,DATA,Mats,Pars;Plot_A=Plot_Live,Plot_F=Plot_Final,vdmact=Tvdm,RepIt=RepIt,RepF=RepF,rhoZ=rhoZ,rhoN=rhoN,rhoP=rhoP,rhoT=rhoT,rhoV=rhoV)
  tmp.vdm
end
