include("RMLib.jl")
include("RMTypes.jl")
using PyPlot
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


#figure()
#imshow(L)
#show()


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
			Ds[i,j] = 1.
		end
		if (i-1) == j
			Ds[i,j] = -1.
		end
	end
end
DsT = Ds'
Dv = zeros(num_lines,num_lines)
for i in 1:num_lines
	for j in 1:num_lines
		if (i+1) == j
			Dv[i,j] = 1.
		end
		if (i-1) == j
			Dv[i,j] = -1.
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
	nits = 50 #Number of ADMM Iterations (An upper limit if convergence is not reached)
	final_it = nits
	initial_psi = 0.0  #Initial value for PSI
		#= REGULARIZATION WEIGHT =#
	#mu = 100.0 #Begins to work around 48 fails around 115000.0
	  #= CONVERGANCE PARAMETERS =#
	G = 10.0 # only for first iteration 1.5 otherwise
	tau = 1.2
	sigma = 0.9
	#ADMM ARRAYS:
	Con_Arr = zeros(num_lines)
	#IMAGE ARRAYS
	
	Pars = init_Params()
	Pars.it = 2
	Pars.nits =100
	Pars.tau = 1.2
	Pars.sigma=0.9
	Pars.G = 10.0
	Pars.num_lines = num_lines
	Pars.num_tdf_times = num_tdf_times
	Pars.initial_psi = 0.5
	Pars.mu_smoo = mu_smoo
	Pars.mu_spec = mu_spec
	Pars.mu_temp = mu_temp
	Pars.eps_abs = 0.1	
	Pars.eps_rel = 0.1
	Pars.alpha = 8.0

	rho0=1000.0
	X = init_IMAGE(rho0)
	X.vdm = zeros((num_tdf_times,num_lines))
	fill!(X.vdm,initial_psi)
	X.vdm_squiggle = zeros((num_tdf_times,num_lines))
	fill!(X.vdm_squiggle,initial_psi)
 
	Z = init_IMAGE(rho0)
	Z.vdm = zeros((num_tdf_times,num_lines))
	fill!(Z.vdm,initial_psi)
	Z.vdm_squiggle = zeros((num_tdf_times,num_lines))
	fill!(Z.vdm_squiggle,initial_psi)
	Z.U=zeros((num_tdf_times,num_lines))
	
	T = init_IMAGE(rho0)
	T.vdm = zeros((num_tdf_times,num_lines))
	fill!(T.vdm,initial_psi)
	T.vdm_squiggle = zeros((num_tdf_times,num_lines))
	fill!(T.vdm_squiggle,initial_psi)
	T.U=zeros((num_tdf_times,num_lines))
	
	V = init_IMAGE(rho0)
	V.vdm = zeros((num_tdf_times,num_lines))
	fill!(V.vdm,initial_psi)
	V.vdm_squiggle = zeros((num_tdf_times,num_lines))
	fill!(V.vdm_squiggle,initial_psi)
	V.U=zeros((num_tdf_times,num_lines))
	


	#=  Calculate Initial Multipliers & RHO =#
	
	for p in 1:num_lines
	  Z.U[:,p]=HT * squeeze(W[p,:,:],1) * ( H * vec(Z.vdm[:,p]) - vec(L[:,p]))
	end
	V.U = Z.U
	T.U = Z.U
	
#NO IDEA WHAT TO DO FOR INITIALIZATOIN OF NEW MULTIPLIER'S AND WEIGHTS..........
	Qinv = zeros(num_tdf_times,num_tdf_times)
	B = zeros(num_tdf_times,num_lines)
	converged = 0
	while Pars.it <= Pars.nits && converged==0        #ADMM ITERATION LOOP
		println("Iteration: ", Pars.it-1, "  ")
		X.vdm_previous = X.vdm
	#Step 1: MINIMIZATION W.R.T. X
	  for l in 1:num_lines        #SPECTAL CHANNEL LOOP
		#MINIMIZATION ON X_v INDIVIDUALLY WITHIN THIS SPECTRAL LOOP. (THIS PART CAN BE DONE IN PARALLEL!!!!)			
			W_slice = reshape(W[l,:,:],size(W[l,:,:])[2],size(W[l,:,:])[3])
			Q = HT * W_slice*H + T.rho*DsT*Ds + (mu_smoo+Z.rho+T.rho)*Gammatdf 
			B = HT* W_slice * L + DsT*(T.U+T.rho*T.vdm)+Z.U+Z.rho*Z.vdm
			X.vdm[:,l] = Q\B[:,l] 
		end

	#Step 2: UPDATE REGULARIZATION TERMS
		T.vdm_squiggle = Ds*X.vdm-T.U/T.rho
		T.vdm_previous = T.vdm
		T.vdm = ell2_prox_op(T.vdm_squiggle,mu_temp,T.rho)

		V.vdm_squiggle = Z.vdm*Dv - V.U/V.rho
		V.vdm_previous = V.vdm
		V.vdm = ell2_prox_op(V.vdm_squiggle,mu_spec,V.rho)

		Rinv = inv(V.rho*Dv*DvT+Z.rho*Gammaspe)
		C = (V.U+V.rho*V.vdm)*DvT-Z.U+Z.rho*X.vdm
		Z.vdm_previous = Z.vdm
		Z.vdm = 	C*Rinv

		V.vdm_squiggle = Z.vdm*Dv - V.U/V.rho
		V.vdm_previous = V.vdm
		V.vdm = ell2_prox_op(V.vdm_squiggle,mu_spec,V.rho)
	
	#Step 4: UPDATE LAGRANGE MULTIPLIERS	
		T = LG_update(X,T,Pars.alpha)
		V = LG_update(X,V,Pars.alpha)
		Z = LG_update(X,Z,Pars.alpha)
		println(" U_T: ", mean(T.U), " U_V: " ,mean(V.U), " U_Z: ", mean(Z.U))

  	if Pars.it == 2 
			Pars.G = 1.5
		end
	#Step 5: Update Penalty Parameters
		Z = Pen_update(X,Z,Pars)
		T = Pen_update(X,T,Pars)
		V = Pen_update(X,V,Pars)
		Pars.it = Pars.it+1
	Mod = Model(X.vdm,H) 
	println("Chi2: ", Chi2(Mod,L,EL)/(num_tdf_times*num_lines))
	#imshow(X)
	#draw()
  end
  #ioff()
	
	
	#outarr = [it-1, Chi2(Mod,L[1,:],EL[1,:])/num_tdf_times]
	if converged != 1
		println("Mu = ", mu," did not converge")
	end
	#outarr #NEED A NEW OUTPUT METHOD....FILES I SUSPECT
	println("ADMM Call Completed")
	figure(figsize=(20,10))
	subplot(221)
	imshow(Z.vdm)
	title("Z")
	colorbar()
	subplot(222)
	imshow(X.vdm)
	title("X")
	colorbar()
	subplot(223)
	imshow(T.vdm)
	title("T")
	colorbar()
	subplot(224)
	imshow(V.vdm)
	title("V")
	colorbar()

	show()
end


num_mus=1
#mu = linspace(40,120000.0,num_mus)   #A GOOD RANGE FOR SIGMA = 0.75 On SYNTH DATA : TDF == 0.04
mu = 10.0

mu_spec = 10.0		#Spectral Regularization Weight
mu_temp = 100.0		#Temporal Regularization Weight
mu_smoo = 0.000001    #Smoothing Regularization Weight (TIKHONOV)

chi2 = zeros((num_mus))
its = zeros((num_mus))


tmp = TLDR(mu_smoo,mu_spec,mu_temp,L, num_lines, EL, spectra_dates, num_spectra_samples,continuum_dates,continuum_scale,continuum_flux,continuum_error_flux,num_continuum_dates)
	#chi2[i] = tmp[2]
	#its[i] = tmp[1]
	#println(i/num_mus)

#outarr = [mu chi2 its]

#writecsv("mu_test_mu.csv",outarr)

println("Done")














