include("RMLib.jl")
clear()


mode = 1 	#SYNTHETIC DATA MODE
#mode = 2	#REAL DATA MODE
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
Dv = zeros(num_spectra_samples,num_spectra_samples)
for i in 1:num_spectra_samples
	for j in 1:num_spectra_samples
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
Gamma = eye(num_tdf_times)
GammaT = Gamma'
DT = D'


function TLDR(mu_smoo,mu_spec,mu_temp,L, num_lines, EL, spectra_dates, num_spectra_samples,continuum_dates,continuum_scale,continuum_flux,continuum_error_flux,num_continuum_dates)

	
	#= INITIALIZING ADMM PARAMETERS =#
	nits = 500 #Number of ADMM Iterations (An upper limit if convergence is not reached)
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
	X =   zeros((size(L,1),num_tdf_times))
	fill!(X,initial_psi)
	X_s = zeros((size(L,1),num_tdf_times))
	fill!(X_s,initial_psi)

	Z =   zeros((size(L,1),num_tdf_times))
	fill!(Z,initial_psi)
	Z_previous = zeros((size(L,1),num_tdf_times)) 	#-> UNECESSARY????????? <-
	Z_s = zeros((size(L,1),num_tdf_times))
	fill!(Z_s,initial_psi)

	P = zeros((size(L,1),num_tdf_times))
	fill!(P,initial_psi)
	P_s = zeros((size(L,1),num_tdf_times))
	fill!(P_s,initial_psi)

	T = zeros((size(L,1),num_tdf_times))
	fill!(T,initial_psi)
	T_s = zeros((size(L,1),num_tdf_times))
	fill!(T_s,initial_psi)

	V = zeros((size(L,1),num_tdf_times))
	fill!(V,initial_psi)
	V_s = zeros((size(L,1),num_tdf_times))
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
	U_Z =   zeros((size(L,1),num_tdf_times))
	U_P =   zeros((size(L,1),num_tdf_times))
	U_T =   zeros((size(L,1),num_tdf_times))
	U_V =   zeros((size(L,1),num_tdf_times))

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
	for p in 1:num_lines
	  U[p,:]=HT * squeeze(W[p,:,:],1) * ( H * vec(Z[p,:]) - vec(L[p,:]))
		A =(vec(U[p,:])'*HT*squeeze(W[p,:,:],1)*H*vec(U[p,:]))
		B =(vec(U[p,:])'*vec(U[p,:]))
		rho[p] = A[1]/B[1]
	end
	
#NO IDEA WHAT TO DO FOR INITIALIZATOIN OF NEW MULTIPLIER'S AND WEIGHTS..........


	converged = 0
	it = 2
	while it <= nits && converged==0        #ADMM ITERATION LOOP
	  #println("Iteration: ", it, "  ", sum(Con_Arr), "/", num_lines, " Converged.")

	#Step 1: MINIMIZATION W.R.T. X
		Q = HT * H + rho_T * DsT*Ds + (mu_smoo+rho_S+rho_T)*Gamma
		B = HT * L + DsT*(U_T+rho_T*T)+U_P+rho_P*P+U_S+rho_S*S
	  for l in 1:num_lines        #SPECTAL CHANNEL LOOP
			b=
		#MINIMIZATION ON X_v INDIVIDUALLY WITHIN THIS SPECTRAL LOOP. (THIS PART CAN BE DONE IN PARALLEL!!!!)			
					
		end
	#Step 1 done!
	#Step 2: UPDATE REGULARIZATION TERMS
		#A: POSITIVITY: P
		P_s = X-U_p/rho_P
		P = pos_prox_op(P_s)
		#B: THRESHOLD: T
		T_s = D*X-U_T/rho_T
		T = ell1_prox_op(T_s,mu_temp,rho_T)
		#C: FREQUENCY: V
		V_s = Z*D - U_V/rho_V
		V = ell1_prox_op(V_s,mu_spec,rho_V)
	#Step 2: done!

	#Step 3: MINIMIZATION W.R.T Z
		#TIKHONOV SOLUTION
		Rinv = inv(rho_V*Dv*DvT+rho_S*Gamma)
		C =	(U_V+rho_V*V)*DvT-U_S+rho_S*X 
		Z = Rinv*C
	#Step 3: done!

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









 #Step 2: MINIMIZATION W.R.T Z
#      Z_s[l,:] = X[l,:]-U[l,:]/rho[l]
#      slice = reshape(W[l,:,:],size(W[l,:,:])[2],size(W[l,:,:])[3])
	      #slice = squeeze(W[l,:,:],1)	#RESHAPE IS SUPPOSED TO BE FASTER THAN SQUEEZE!
#	      A =inv(HT*slice*H+rho[l]/2.0*GammaT*Gamma)
	      #B = HT*slice*vec(L[l,:])+rho[it,l]*GammaT*Gamma*vec(Z_s[it,l,:])
#	      B = HT*slice*vec(L[l,:])
#		    Z[l,:] = A * B
    #Step 2 done!


    #Step 3: UPDATE LAGRANGE MULTIPLIERS
      #eq =U[l,:]+(rho[l]/2.0)*(X[l,:]-Z[l,:])
      U[l,:] = U[l,:]+(rho[l]/2.0)*(X[l,:]-Z[l,:])
     #= OR =#   
			#U[l,:] = HT * slice * ( H * vec(Z[l,:]) - vec(L[l,:]))		#UPDATE MULTIPLIERS
    
			#Step 3 done!

		#CONVERGANCE TESTING!!!!!!!!!!!
    # Evaluation:
      eps_abs = 0.0    #Must be >= 0
      eps_rel = 0.1  #Must be between 0 and 1
      S = rho[l]*(Z[l,:]-Z_previous[l,:])
      R = X[l,:]-Z[l,:]
			tau_prim_previous = tau_prim[l]
      tau_prim[l] = sqrt(num_tdf_times)*eps_abs + eps_rel*(maximum([ell2norm(X[l,:]),ell2norm(Z[l,:])]))
			tau_dual_previous[l] = tau_dual[l]
      tau_dual[l] = sqrt(num_tdf_times)*eps_abs + eps_rel*ell2norm(U[l,:])
      eta = ell2norm(R)*tau_dual[l] / (ell2norm(S)*tau_prim[l])		#OPTION 1\
		  #eta = ell2norm(R)*tau_dual_previous[l] / (ell2norm(S)*tau_prim_previous[l])	#OPTION 2
			phi_previous[l] = phi[l]
			phi[l] = max(ell2norm(R/tau_prim[l]),ell2norm(S)/tau_dual[l])
			if phi[l] <= 1.0
				Con_Arr[l] = 1
			end
			if it == 2 
				G = 1.5
			end
			if (1.0/tau <= eta) && (eta <= tau) || (phi[l] < sigma*phi_previous[l])
				#NOTHING HAPPENS IN HERE	
			elseif eta < (1.0/tau)
				rho_max[l] = rho[l]
				if rho_min[l] > 0.0
					rho[l] = sqrt(rho_min[l]*rho_max[l])
				else rho[l] = rho_max[l]/G 
				end
			elseif eta > tau
				rho_min[l] = rho[l]
				if rho_max[l] < Inf
					rho[l] = sqrt(rho_min[l]*rho_max[l])
				else rho[l] = rho_min[l]*G
				end
			end
      if sum(Con_Arr) == num_lines
        converged = 1
        #println("All lines Converged!!!") 
      end 
  end
  it = it+1
	end
	Mod = Model(X[1,:],H) 
	outarr = [it-1, Chi2(Mod,L[1,:],EL[1,:])/num_tdf_times]
	if converged != 1
		println("Mu = ", mu," did not converge")
	end
	outarr
end


num_mus=500
mu = linspace(40,120000.0,num_mus)   #A GOOD RANGE FOR SIGMA = 0.75 On SYNTH DATA : TDF == 0.04


mu_spec = 1.0		#Spectral Regularization Weight
mu_temp = 1.0		#Temporal Regularization Weight
mu_smoo = 1.0		#Smoothing Regularization Weight

chi2 = zeros((num_mus))
its = zeros((num_mus))

for i in 1:num_mus
	tmp = TLDR(mu[i],L, num_lines, EL, spectra_dates, num_spectra_samples,continuum_dates,continuum_scale,continuum_flux,continuum_error_flux,num_continuum_dates)
	chi2[i] = tmp[2]
	its[i] = tmp[1]
	println(i/num_mus)
end

outarr = [mu chi2 its]

writecsv("mu_test_mu.csv",outarr)

println("Done")














