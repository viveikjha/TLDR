#RUN THIS SCRIPT WITH julia -p 3 TO USE MULTIPLE PROCESSORS.

include("../RMLib.jl")
include("../RMTypes.jl")
include("../DataImport.jl")
include("../GenMatrices.jl")
include("../DataImportNEW.jl")
@everywhere include("../dev.jl")
using PyPlot


function pmap(f, lst)
  np = nprocs()   #GETS THE NUMBER OF PROCESSES AVAILABLE
  n = length(lst)
  results = cell(n,4)
  i = 1
  #function to produce the next work item from the queue
  #in this case it's just an index
  nextidx() = (idx=i; i+=1; idx)
  @sync begin
    for p=1:np
      if p != myid() || np == 1
        @async begin
          while true
            idx = nextidx()
            if idx > n
              break
            end #IF
            results[idx,:] = remotecall_fetch(p,f,lst[idx]) #CALLS FUNCTION f ON PROCESSOR p WITH PARAMETER lst[idx]
            println("µ " , idx, " of ", length(lst)," complete.")
          end #WHILE
        end #@ASYNC
      end #IF
    end #FOR
  end #@SYNC
  results
end #FUNCTION


@everywhere vdmfile="Sim_VDM.csv"
@everywhere truevdm=readcsv(vdmfile)
@everywhere H = readcsv("H.csv")
@everywhere files=["../data/rvm_wavelengths.csv","../simulation/Sim_Spectra.csv","../simulation/Sim_Error.csv","../data/rvm_dates.csv","../data/arp151.b.dat"]
@everywhere Error = (readcsv("Sim_Error.csv"))'
@everywhere Flux = (readcsv("Sim_Spectra.csv"))'
@everywhere off=0.000001

µ=[1.0,10.0,1.0e2,500.0,1.0e3,5000.0,1.0e4,5.0e4,1.0e5,5.0e5,1.0e6,5.0e6,1.0e7,5.0e7,1.0e8,5.0e8,1.0e9,5.0e9,1.0e10,5.0e10,1.0e11,5.0e11,1.0e12,5.0e12,1.0e13]
#µ=[1.0,10.0,1.0e2]
flx=sqrt(sum(truevdm.^2))
#r=cell(length(µ),4)
#println("Beginning minimizations with µs.")


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

nextidx() = (idx=i; i+=1; idx)

@everywhere function f(mu)
  vdm = gen_tiksol(files;mu_smoo=mu,plotting=false,save=false,nits=5)
  M=H*vdm
  residual = sqrt(sum(M-Flux).^2)
  flux = (sqrt(sum(vdm.^2)))
  chi2 = sum((((M)-Flux)./(Error)).^2)/length(Flux)
  [mu,residual,flux,chi2]
end
tic()
out = pmap(f,µ)
toc()
println(out)
writecsv("RecSol.csv",out)
# In output file: [Mu, residual, flux, chi2]
#writecsv("TikRes.csv",r)
println("Complete.")
