include("../RMLib.jl")
include("../RMTypes.jl")
include("../DataImport.jl")
include("../GenMatrices.jl")
include("../DataImportNEW.jl")
include("../dev.jl")
using PyPlot

function res(A,B)
  residual=sqrt(sum((A-B).^2))#./(sqrt(sum(B.^2)))
end

vdmfile="Sim_VDM.csv"
truevdm=readcsv(vdmfile)

files=["../data/rvm_wavelengths.csv","Sim_Spectra.csv","Sim_Error.csv","../data/rvm_dates.csv","../data/arp151.b.dat"]
#files=["../data/rvm_wavelengths.csv","../data/rvm_fluxes.csv","../data/rvm_errfluxes.csv","../data/rvm_dates.csv","../data/arp151.b.dat"]

off=0.000001

µ=[1.0,10.0,1.0e2,500.0,1.0e3,5000.0,1.0e4,5.0e4,1.0e5,5.0e5,1.0e6,5.0e6,1.0e7,5.0e7,1.0e8,5.0e8,1.0e9,5.0e9,1.0e10,5.0e10,1.0e11,5.0e11,1.0e12,5.0e12,1.0e13]
#µ=[1.0,10.0,1.0e2]
flx=sqrt(sum(truevdm.^2))
r=zeros(length(µ))
println("Beginning minimizations with µs.")
for i=1:length(µ)
  Rec,p=LAUNCH(files;mu_smoo=µ[i],mu_spec=off,mu_temp=off,mu_l1=off,scale=1.0,nits=50,Tvdm="",Plot_Live=false,Plot_Final=false,RepIt=false,RepF=false)
  #Rec=gen_tiksol(files;scale=1.0,mu_smoo=µ[i],plotting=false,save=false)
  #r[i]=res(Rec,truevdm)
  r[i]=res(Rec.vdm,truevdm)
  println("µ " , i, " of ", length(µ)," complete.")
end

#figure()
#loglog(µ,r,"k.")
#show()

writecsv("TikMus.csv",µ)
writecsv("TikRes.csv",r)
writecsv("TikFlx.csv"l,flx)
println("Complete.")
