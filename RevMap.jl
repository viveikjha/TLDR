include("RMLib.jl")
clear()


#mode = 1 	#SYNTHETIC DATA MODE
mode = 2	#REAL DATA MODE
if mode == 1
	path ="synth/"
	wavelength_filename = "wavelengthS.csv"
	spectra_filename= "rvm_flxS.csv"
	spectra_error_filename = "rvm_flx_errS.csv"
	spectra_dates_filename = "rvm_dateS.csv"
	continuum_array_filename = "continuumS.csv"
	continuum_array_path =string(path,continuum_array_filename)
	continuum_array = readcsv(continuum_array_path)
elseif mode == 2
	path="data/"
	wavelength_filename = "rvm_wavelengths.csv"
	spectra_filename= "rvm_fluxes.csv"
	spectra_error_filename = "rvm_errfluxes.csv"
	spectra_dates_filename = "rvm_dates.csv"
	continuum_array_filename = "arp151.b.dat"
	continuum_array_path =string(path,continuum_array_filename)
	continuum_array = readdlm(continuum_array_path)
end
#=IMPORTING FILES=#
#SPECTRA
wavelength_path = string(path,wavelength_filename)
wavelength = readcsv(wavelength_path)                    #List of measured wavelengths
spectrapath = string(path,spectra_filename)
L = readcsv(spectrapath)                           #SPECTRAL FLUXES (L)
num_lines = size(L,1)                                #NUMBER OF SPECTRAL LINES
spectra_error_path = string(path,spectra_error_filename)
EL = readcsv(spectra_error_path)                 #SPECTRAL FLUX ERRORS
spectra_dates_path = string(path,spectra_dates_filename)
spectra_dates = readcsv(spectra_dates_path)                   #SPECTRAL SAMPLING DATES
num_spectra_samples = length(spectra_dates)      				  #NUMBER OF DATA POINTS IN SPECTRA
#CONTINUUM
println("Continuum: ",size(continuum_array))
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
#= COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS =#
interpolation_points = zeros(num_spectra_samples,num_tdf_times)
H = zeros(num_spectra_samples,num_tdf_times)
HE= zeros(num_spectra_samples,num_tdf_times)
for date in 1:num_spectra_samples
	for delay in 1:num_tdf_times
		interpolation_points[date,delay]=spectra_dates[date]-tdf_times[delay]
  end
  P = interpolation_points[date,:]
  H[date,:] =interp(interpolation_points[date,:],continuum_dates,continuum_flux)
  HE[date,:] = interp(interpolation_points[date,:],continuum_dates,continuum_error_flux)
end

#= 		FINITE DIFFERENCES MATRICES			=#
Ds = zeros(num_tdf_times,num_tdf_times)
for i in 1:num_tdf_times
	for j in 1:num_tdf_times
		if (i+1) == j
			Ds[i,j] = 1
		end
		if (i-1) == j
			Ds[i,j] = -1
		end
	end
end
DsT = Ds'
Dv = zeros(num_lines,num_lines)
for i in 1:num_lines
	for j in 1:num_lines
		if (i+1) == j
			Dv[i,j] = 1
		end
		if (i-1) == j
			Dv[i,j] = -1
		end
	end
end
DvT = Dv'


#=    PRECOMPUTING TIKHONOV MATRICES     =#
num_spectra_dates=size(spectra_dates)[1]
W = zeros((size(wavelength)[1],size(L)[2],size(L)[2]))
for lam in 1:num_lines
  T = eye(num_spectra_dates)
  for i in 1:num_spectra_dates
    for j in 1:num_spectra_dates
      if i == j
        T[i,j] =1.0/ EL[lam,i].^2
      end
    end
  end
  W[lam,:,:] = T
end
HT = H'
Gammatdf = eye(num_tdf_times)
GammatdfT = Gammatdf'
Gammaspe = eye(num_lines)
GammaspeT = Gammaspe'



function TLDR(mu_smoo,mu_spec,mu_temp,L, num_lines, EL, spectra_dates, num_spectra_samples,continuum_dates,continuum_scale,continuum_flux,continuum_error_flux,num_continuum_dates)
	L = L'
	EL = EL'
	#= INITIALIZING ADMM PARAMETERS =#
	nits = 20 #Number of ADMM Iterations (An upper limit if convergence is not reached)
	final_it = nits
	initial_psi = 0.0  #Initial value for PSI
		#= REGULARIZATION WEIGHT =#
	#mu = 100.0 #Begins to work around 48 fails around 115000.0
	  #= CONVERGANCE PARAMETERS =#
	G = 10.0 # only for first iteration 1.5 otherwise
	tau = 1.1
	sigma = 0.9
	#ADMM ARRAYS:
	Con_Arr = zeros(num_lines)
	#IMAGE ARRAYS
	X = zeros((num_tdf_times,num_lines))
	fill!(X,initial_psi)
	X_s = zeros((num_tdf_times,num_lines))
	fill!(X_s,initial_psi)
 
	Z = zeros((num_tdf_times,num_lines))
	fill!(Z,initial_psi)
	Z_previous = zeros((num_tdf_times,num_lines)) 	#-> UNECESSARY????????? <-
	Z_s = zeros((num_tdf_times,num_lines))
	fill!(Z_s,initial_psi)

	P = zeros((num_tdf_times,num_lines))
	fill!(P,initial_psi)
	P_s = zeros((num_tdf_times,num_lines))
	fill!(P_s,initial_psi)

	T = zeros((num_tdf_times,num_lines))
	fill!(T,initial_psi)
	T_s = zeros((num_tdf_times,num_lines))
	fill!(T_s,initial_psi)

	V = zeros((num_tdf_times,num_lines))
	fill!(V,initial_psi)
	V_s = zeros((num_tdf_times,num_lines))
	fill!(V_s,initial_psi)


	#Aug. Lagrangian Parameters
	rho_Z = 0.0
	rho_Z_min = 0.0
	rho_Z_max = Inf

	rho_P = 0.0
	rho_P_min =0.0
	rho_P_max = Inf

	rho_T = 0.0
	rho_T_min = 0.0
	rho_T_max = Inf

	rho_V = 0.0
	rho_V_min = 0.0
	rho_V_max = Inf

#LAGRANGE MULTIPLIERS
	U_Z = zeros((num_tdf_times,num_lines))
	U_P = zeros((num_tdf_times,num_lines))
	U_T = zeros((num_tdf_times,num_lines))
	U_V = zeros((num_tdf_times,num_lines))

#CONVERGANCE ARRAYS
	phi = zeros((size(L,1)))
	phi_previous = zeros((size(L,1)))
	tau_prim = zeros((num_lines))
	tau_prim_previous = zeros((num_lines))
	tau_dual = zeros((num_lines))
	tau_dual_previous = zeros((num_lines))

#ACTUAL TDF ARRAY FOR PLOTTING
	X_act = zeros((50))
	fill!(X_act,initial_psi)
	
	#=  Calculate Initial Multipliers & RHO =#
	println("H: ", size(H))
	println("HT: ", size(HT))
	println("W: ",size(squeeze(W[1,:,:],1)))
	println("Z: ", size(Z))
	println("L: ", size((L)))
	#println("U: ", size(U))


	for p in 1:num_lines
	  U_Z[:,p]=HT * squeeze(W[p,:,:],1) * ( H * vec(Z[:,p]) - vec(L[:,p]))
	end
	U_T = U_Z
	U_V = U_Z
	U_P = U_Z
#	rho_Z = U_Z'*HT*H*U_Z*inv(U_Z'*U_Z)
	rho_Z = 1000.0
	rho_T = rho_Z
	rho_V = rho_Z
	rho_P = rho_Z
#NO IDEA WHAT TO DO FOR INITIALIZATOIN OF NEW MULTIPLIER'S AND WEIGHTS..........
	Qinv = zeros(num_tdf_times,num_tdf_times)
	B = zeros(num_tdf_times,num_lines)
	converged = 0
	it = 2
	while it <= nits && converged==0        #ADMM ITERATION LOOP
	 println("Iteration: ", it, "  ", sum(Con_Arr), "/", num_lines, " Converged.")

	#Step 1: MINIMIZATION W.R.T. X
	  for l in 1:num_lines        #SPECTAL CHANNEL LOOP
			if it == 3
				println(l)
			end		
		#MINIMIZATION ON X_v INDIVIDUALLY WITHIN THIS SPECTRAL LOOP. (THIS PART CAN BE DONE IN PARALLEL!!!!)			
			W_slice = reshape(W[l,:,:],size(W[l,:,:])[2],size(W[l,:,:])[3])
			#println("##########")
			#println(size(rho_T))
			i = HT * W_slice*H
			j = rho_T*DsT*Ds
			k = (mu_smoo+rho_Z+rho_T)*Gammatdf
			#Qinv = inv(HT *W_slice* H + rho_T * DsT*Ds + (mu_smoo+rho_Z+rho_T)*Gammatdf) #GAMMA MAY NOT BE CORRECT DIMS
			#Qinv = inv(i+j+k)
			Q = i+j+k			
			#println("&&&&&&&&&&")
			B = HT* W_slice * L + DsT*(U_T+rho_T*T)+U_P+rho_P*P+U_Z+rho_Z*Z
			#println("@@@@@@@@@@")
			X[:,l] = Q\B[:,l] 
		 
		end
	#Step 1 done!
	println("STEP 1 COMPLETE")
	#Step 2: UPDATE REGULARIZATION TERMS
		#A: POSITIVITY: P
		P_s = X-U_P/rho_P
		P = pos_prox_op(P_s)
		#B: THRESHOLD: T
		T_s = Ds*X-U_T/rho_T
		T = ell1_prox_op(T_s,mu_temp,rho_T)
		#C: FREQUENCY: V
		V_s = Z*Dv - U_V/rho_V
		V = ell1_prox_op(V_s,mu_spec,rho_V)
	#Step 2: done!
	println("STEP 2 COMPLETE")
	#Step 3: MINIMIZATION W.R.T Z
		#TIKHONOV SOLUTION
		Z = (1.0/(rho_V+rho_Z))*(U_V+rho_V*V)*DvT-U_Z+rho_Z*X 
	#Step 3: done!
		println("STEP 3 COMPLETE")
	#Step 4: UPDATE LAGRANGE MULTIPLIERS
		#A: U_P
		U_P = U_P+(rho_P/2.0)*(X-P)
		#B: U_T
		U_T = U_T+(rho_T/2.0)*(X-T)
		#C: U_V
		U_V = U_V+(rho_V/2.0)*(X-V)
		#D: U_S
		U_Z = U_Z+(rho_Z/2.0)*(X-Z)
	#Step 4: done!

#Step 5: CONVERGANCE TESTING!!!!!!!!!!!
   # Evaluation:
    eps_abs = 0.0    #Must be >= 0
    eps_rel = 0.1  #Must be between 0 and 1
    S = rho_Z*(Z-Z_previous)
    R = X-Z
		tau_prim_previous = tau_prim
    tau_prim = sqrt(num_tdf_times)*eps_abs + eps_rel*(maximum([ell2norm(X),ell2norm(Z)]))
		tau_dual_previous = tau_dual
    tau_dual = sqrt(num_tdf_times)*eps_abs + eps_rel*ell2norm(U_Z)
    eta = ell2norm(R)*tau_dual / (ell2norm(S)*tau_prim)		#OPTION 1\
	  #eta = ell2norm(R)*tau_dual_previous[l] / (ell2norm(S)*tau_prim_previous[l])	#OPTION 2
		phi_previous = phi
		phi = max(ell2norm(R/tau_prim),ell2norm(S)/tau_dual)
		if phi <= 1.0
			Con_Arr = 1
		end
		if it == 2 
			G = 1.5
		end
		#println("tau: ", (tau[1]))
		#println("eta: ", (eta[1]))
		#println("Phi: ", (phi[1]))
		#println("sigma: ", sigma[1])
		
		if (1.0/tau[1] <= eta[1]) && (eta[1] <= tau[1]) || (phi[1] < sigma[1]*phi_previous[1])
			#NOTHING HAPPENS IN HERE	
		elseif eta < (1.0/tau)
			rho_max = rho_Z
			if rho_Z_min > 0.0
				rho_Z = sqrt(rho_Z_min*rho_Z_max)
			else rho_Z = rho_Z_max/G 
			end
		elseif eta > tau
			rho_Z_min = rho_Z
			if rho_Z_max < Inf
				rho_Z = sqrt(rho_Z_min*rho_Z_max)
			else rho_Z = rho_Z_min*G
			end
		end
		rho_T = rho_Z
		rho_V = rho_Z
		rho_P = rho_Z
    if sum(Con_Arr) == num_lines
	    converged = 1
      #println("All lines Converged!!!") 
    end 
		it = it+1
  end
  
	
	Mod = Model(X,H) 
	#outarr = [it-1, Chi2(Mod,L[1,:],EL[1,:])/num_tdf_times]
	if converged != 1
		println("Mu = ", mu," did not converge")
	end
	#outarr #NEED A NEW OUTPUT METHOD....FILES I SUSPECT
	println("ADMM Call Completed")
end


num_mus=1
#mu = linspace(40,120000.0,num_mus)   #A GOOD RANGE FOR SIGMA = 0.75 On SYNTH DATA : TDF == 0.04
mu = 10.0

mu_spec = 10.0		#Spectral Regularization Weight
mu_temp = 10.0		#Temporal Regularization Weight
mu_smoo = 10.0		#Smoothing Regularization Weight

chi2 = zeros((num_mus))
its = zeros((num_mus))


tmp = TLDR(mu_smoo,mu_spec,mu_temp,L, num_lines, EL, spectra_dates, num_spectra_samples,continuum_dates,continuum_scale,continuum_flux,continuum_error_flux,num_continuum_dates)
	#chi2[i] = tmp[2]
	#its[i] = tmp[1]
	#println(i/num_mus)

#outarr = [mu chi2 its]

#writecsv("mu_test_mu.csv",outarr)

println("Done")














