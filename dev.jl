include("RMLib.jl")
include("RMLibMore.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("DataImportNEW.jl")
include("GenMatrices.jl")

using PyPlot

function getdata(FILES_ARR)
  #FILES_ARR = ["data/rvm_wavelengths_trimmed.csv","data/rvm_fluxs_trimmed.csv","data/rvm_errfluxes_trimmed.csv","data/rvm_dates.csv","data/arp151.b.dat"]
  wavelengths=FILES_ARR[1]
  spectra = FILES_ARR[2]
  errspectra = FILES_ARR[3]
  dates = FILES_ARR[4]
  continuum = FILES_ARR[5]
  DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)

  #data_report(DATA)
  #Pars = init_Params()
  #Mats=Gen_Mats(DATA,Pars)

  #Pars,Mats,DATA
  DATA
end

#function gen_tiksol(DATA,Pars,Mats;scale=1.0,mu_smoo=40.0)
function gen_tiksol(f;scale=1.0,mu_smoo=40.0,plotting=true,save=false)
  DATA=getdata(f)
  println("Scale: ",scale)
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
