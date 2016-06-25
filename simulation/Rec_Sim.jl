#RUN THIS SCRIPT WITH julia -p 3 TO USE MULTIPLE PROCESSORS.

include("../RMLib.jl")
include("../RMTypes.jl")
include("../DataImport.jl")
include("../GenMatrices.jl")
include("../DataImportNEW.jl")
include("../dev.jl")
using PyPlot




vdmfile="Spiral/simulated_vdm.csv"
truevdm=readcsv(vdmfile)
H = readcsv("H.csv")
files=["../data/rvm_wavelengths.csv","../simulation/Sim_Spectra.csv","../simulation/Sim_Error.csv","../data/rvm_dates.csv","../data/arp151.b.dat"]
Error = (readcsv("Sim_Error.csv"))'
Flux = (readcsv("Sim_Spectra.csv"))'
off=0.000001



#tic()
#for i=1:length(µ)
#  vdm=gen_tiksol(files;scale=1.0,mu_smoo=µ[i],plotting=false,save=false)
#  M=H*vdm
#  residual = sqrt(sum(M-Flux).^2)
#  flux = (sqrt(sum(vdm.^2)))
#  chi2 = sum(((M-Flux)./(Error)).^2)/length(Flux)
#  r[i,:]=[µ[i],residual,flux,chi2]

#  println("µ " , i, " of ", length(µ)," complete.")
#end
#toc()

#println(r)
mu_smoo=1.0e3
mu_spec=0.001
mu_temp=0.001
mu_l1=0.001
rec = COLD_LAUNCH(files;mu_smoo=mu_smoo,mu_spec=mu_spec,mu_temp=mu_temp,mu_l1=mu_l1,nits=50,Plot_Live=false,Plot_Final=false,RepIt=false,RepF=false)
println("Complete.")
