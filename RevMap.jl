using Optim
#using PyPlot
include("RMLib.jl")
#path="data/"
path ="synth/"
clear()
#=IMPORTING FILES=#
#SPECTRA
wavelength_filename = "rvm_wavelengths.csv"
wavelength_path = string(path,wavelength_filename)
wavelength = readcsv(wavelength_path)                    #List of measured wavelengths

#spectra_filename= "rvm_fluxes.csv"
spectra_filename= "rvm_flxS.csv"
spectrapath = string(path,spectra_filename)
L = readcsv(spectrapath)                           #SPECTRAL FLUXES (L)
num_lines = size(L,1)                                #NUMBER OF SPECTRAL LINES

#spectra_error_filename = "rvm_errfluxes.csv"
spectra_error_filename = "rvm_flx_errS.csv"
spectra_error_path = string(path,spectra_error_filename)
EL = readcsv(spectra_error_path)                 #SPECTRAL FLUX ERRORS

spectra_dates_filename = "rvm_dates.csv"
spectra_dates_path = string(path,spectra_dates_filename)
spectra_dates = readcsv(spectra_dates_path)                   #SPECTRAL SAMPLING DATES
num_spectra_samples = length(spectra_dates)      				  #NUMBER OF DATA POINTS IN SPECTRA

#CONTINUUM
#continuum_array_filename = "arp151.b.dat"
continuum_array_filename = "continuumS.csv"
continuum_array_path =string(path,continuum_array_filename)
#continuum_array = readdlm(continuum_array_path)
continuum_array = readcsv(continuum_array_path)
println(size(continuum_array))
#CONTINUUM_ARRAY CONTAINTS THE CONTINUUM DATES, FLUX, AND FLUX ERRORS.
#IN THAT ORDER.

continuum_dates = continuum_array[:,1]
continuum_scale = 1.0
continuum_flux = continuum_array[:,2]*continuum_scale
continuum_error_flux = continuum_array[:,3]*continuum_scale
num_continuum_dates = length(continuum_dates)
#=DONE IMPORTING FILES=#

#=SETTING UP TIME DELAY FUNCTION=#
num_tdf_times = 50
tdf_values = zeros(num_tdf_times)
tdf_times = linspace(1,20,num_tdf_times)
#=DONE SETTING UP TIME DELAY FUNCTION=#

#= COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS =#
interpolation_points = zeros(num_spectra_samples,num_tdf_times)
ICF = zeros(num_spectra_samples,num_tdf_times)
ICFE= zeros(num_spectra_samples,num_tdf_times)
for date in 1:num_spectra_samples
	for delay in 1:num_tdf_times
		interpolation_points[date,delay]=spectra_dates[date]-tdf_times[delay]
  end
  P = interpolation_points[date,:]
  #ICF[date,:] = interp(P[:],continuum_dates,continuum_flux)
  ICF[date,:] =interp(interpolation_points[date,:],continuum_dates,continuum_flux)
  #println(interp(interpolation_points[date,:],continuum_dates,continuum_flux))
  ICFE[date,:] = interp(interpolation_points[date,:],continuum_dates,continuum_error_flux)
end


#=
plot(continuum_dates,continuum_flux)
plot(interpolation_points[1,:],ICF[1,:],"r*")
show()
=#

#for i in 1:10
#  println(interpolation_points[i,1],"     ",ICF[i,1])
#end

#=    PRECOMPUTING TIKHONOV MATRICES     =#

#Build Mapping Matrix
num_spectra_dates=size(spectra_dates)[1]
H = zeros((num_spectra_dates, num_tdf_times))
for lns in 1:num_spectra_dates
  for tdf in 1:num_tdf_times
    if lns >= tdf
        H[lns,tdf] = ICF[lns,tdf]
    end
  end
end

#Build Covariance Matrix
W = zeros((size(wavelength)[1],size(L)[2],size(L)[2]))
for lam in 1:num_lines
  T = eye(num_spectra_dates)
  for i in 1:num_spectra_dates
    for j in 1:num_spectra_dates
      #println("!")
      if i == j
	println(EL[lam,i])
        T[i,j] =1.0/ EL[lam,i].^2
      end
    end
  end
  W[lam,:,:] = T
end

HT = H'
Gamma = eye(num_tdf_times)
GammaT = Gamma'

#= INITIALIZING ADMM PARAMETERS =#
nits = 30 #Number of ADMM Iterations (An upper limit if convergence is not reached)
final_it = nits
initial_psi = 0.1  #Initial value for PSI


	#= REGULARIZATION WEIGHTS =#
mu = 10.0

  #= CONVERGANCE PARAMETERS =#
G = 10.0 # only for first iteration 1.5 otherwise
tau = 1.2
Con_Arr = zeros(num_lines)

X =   zeros((nits+1, size(L,1),num_tdf_times))
X_s = zeros((nits+1, size(L,1),num_tdf_times))
Z =   zeros((nits+1, size(L,1),num_tdf_times))
Z_s = zeros((nits+1, size(L,1),num_tdf_times))
rho = zeros((nits+1, size(L,1)))
U =   zeros((nits+1, size(L,1),num_tdf_times))
tau_prim = zeros((nits,num_lines))
tau_dual = zeros((nits,num_lines))
value = 1.0
initial_rho = 2.0
fill!(X,initial_psi)
fill!(X_s,initial_psi)
fill!(Z,initial_psi)
fill!(Z_s,initial_psi)
fill!(rho,initial_rho)
#=  Calculate Initial Multipliers  =#
for p in 1:num_lines
#for p in 1:2
  U[1,p,:]=HT * squeeze(W[p,:,:],1) * ( H * vec(Z[1,p,:]) - vec(L[p,:]))
end

#= Calculate Initial Rho =#
#Holding off on this due to minimization;
#NEEDS GRADIENT TO WORK!
#f(V) = Chi2(  Model((Z[1,1,:]-U[1,1,:]./V),ICF) ,vec(L[1,:]),vec(EL[1,:]))
#println(optimize(f,[1.0], method= :l_bfgs))

converged = 0
it = 2
nits = 20

while it <= nits && converged==0        #ADMM ITERATION LOOP
  println("Iteration: ", it)
  for l in 1:num_lines        #SPECTAL CHANNEL LOOP
  #for l in 1:1
    if Con_Arr[l] != 1
    #Step 1:
      X_s[it,l,:] = Z[it-1,l,:]-U[it-1,l,:]/rho[it,l]
      X[it,l,:] = Prox(X_s[it,l,:],mu,rho[it,l])
      #converged = 1
    #Step 1 done!
    #Step 2:
      Z_s[it,l,:] = X[it,l,:]-U[it-1,l,:]/rho[it,l]
      slice = reshape(W[l,:,:],size(W[l,:,:])[2],size(W[l,:,:])[3])
      #slice = squeeze(W[l,:,:],1)	#RESHAPE IS SUPPOSED TO BE FASTER THAN SQUEEZE!

      A =inv(((HT*slice)*H)+(rho[it,l]/2.0)*(GammaT*Gamma))
      B = (HT*slice)*vec(L[l,:])+(rho[it,l])*GammaT*Gamma*vec(Z_s[it,l,:])
      Z[it,l,:] = A * B

    #Step 2 done!
    #Step 3:
      eq =U[it-1,l,:]+(rho[it,l]/2.0)*(X[it,l,:]-Z[it,l,:])
      #U[it,l,:] = U[it-1,l,:]+(rho[it,l]/2.0)*(X[it,l,:]-Z[it,l,:])
      #= OR =#
      U[it,l,:] = HT * slice * ( H * vec(Z[it,l,:]) - vec(L[l,:]))
      
    #Step 3 done!
    # Evaluation:
      eps_abs = 0.00     #Must be >= 0
      eps_rel = 0.001  #Must be between 0 and 1
      S = rho[it,l]*(Z[it,l,:]-Z[it-1,l,:])
      R = X[it,l,:]-Z[it,l,:]

      tau_prim[it,l] = sqrt(num_tdf_times)*eps_abs + eps_rel*(maximum([ell2norm(X[it,l,:]),ell2norm(Z[it,l,:])]))
      tau_dual[it,l] = sqrt(num_tdf_times)*eps_abs + eps_rel*ell2norm(U[it,l,:])

      #eta = ell2norm(R)*tau_dual[it,l] / (ell2norm(S)*tau_prim[it,l])		#OPTION 1
      eta = ell2norm(R)*tau_dual[it-1,l] / (ell2norm(S)*tau_prim[it-1,l])	#OPTION 2


      if l == 250
	e1 = H*vec(X[it,l,:]).-vec(L[l,:])
	ext =  e1'*squeeze(W[l,:,:],1)*e1
        #println("Line = ", l, "Chi2 = ",Chi2(Model((X[it,l,:]),ICF),L[l,:],EL[l,:]), " Mat = ",ext)
	println("Line = ", l, "Chi2 = ",Chi2(Model((X[it,l,:]),H),L[l,:],EL[l,:]), " Mat = ",ext)        
	println("Rho: ",rho[it,l]," U: ",mean(U[it-1,l,:]))
	println("S: ", ell2norm(S), "    R: ", ell2norm(R))
	println("tau_prim: ", tau_prim[it,l], "    tau_dual: ", tau_dual[it,l])
	println("--------------------------------------")
     end
	
      #Check for convergance
      if ell2norm(R) < tau_prim[it,l] && ell2norm(S) < tau_dual[it,l]
        println("Line Converged")
        Con_Arr[l]=1
      end
      if sum(Con_Arr) == num_lines
        converged = 1
        final_it = it
        println("All lines Converged!!!") 
      end


    #Evaluation Done
    #Update Rho
      if it == 2
        G = 1.5
      end
      if eta > tau
        rho[it+1,l] = G * rho[it,l]
      elseif eta < 1.0/tau
        rho[it+1,l] = rho[it,l]/G
      else
        rho[it+1,l] = rho[it,l]
    #Rho done
      end
    else
      X[it,l,:]=X[it-1,l,:]
      X_s[it,l,:]=X_s[it-1,l,:]
      Z[it,l,:]=Z[it-1,l,:]
      Z_s[it,l,:]=Z_s[it-1,l,:]
    end
  end
  it = it+1
end
Save_Array = reshape(X[final_it,:,:],size(X)[2],size(X)[3])
output_filename = "tdf.csv"
writecsv(output_filename,Save_Array)


println("done")
