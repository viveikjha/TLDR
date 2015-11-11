include("RMLib.jl")
using PyPlot
clear()

function ell2norm(X)
	sqrt(sum(vec(X).^2))   #Array Returned!
	#sqrt(sum(X).^2)
end

function Model(X,H)
  #L = size(ICF)[1]
  MF = zeros(L)
  MF = H*vec(X)
	#MF = H*X
  MF #Array Returned!
end

function Chi2(M,D,Sigma)
	  #sum(   (vec(M)-vec(D)).^2  ./vec(Sigma) .^2)   #Value Returned!
#		sum(   (M-D).^2 ./(Sigma).^2)
	sum(((vec(M)-vec(D))/vec(Sigma))^2)
end

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


true_vdm = readcsv("synth/TDF.csv")


#SPECTRA
wavelength_path = string(path,wavelength_filename)
wavelength = readcsv(wavelength_path)     	              #List of measured wavelengths
spectrapath = string(path,spectra_filename)
L = readcsv(spectrapath)                          			 	#SPECTRAL FLUXES (L)
num_lines = size(L,1)                             		  	#NUMBER OF SPECTRAL LINES
spectra_error_path = string(path,spectra_error_filename)
EL = readcsv(spectra_error_path)                 					#SPECTRAL FLUX ERRORS
spectra_dates_path = string(path,spectra_dates_filename)
spectra_dates = readcsv(spectra_dates_path)               #SPECTRAL SAMPLING DATES
num_spectra_samples = length(spectra_dates)      				  #NUMBER OF DATA POINTS IN SPECTRA
#CONTINUUM
#println("Continuum: ",size(continuum_array))
#CONTINUUM_ARRAY CONTAINTS THE CONTINUUM DATES, FLUX, AND FLUX ERRORS.
#IN THAT ORDER.
continuum_dates = continuum_array[:,1]
continuum_scale = 1.0
continuum_flux = continuum_array[:,2]*continuum_scale
continuum_error_flux = continuum_array[:,3]*continuum_scale
num_continuum_dates = length(continuum_dates)





function TLDR(mu,L, num_lines, EL, spectra_dates, num_spectra_samples,continuum_dates,continuum_scale,continuum_flux,continuum_error_flux,num_continuum_dates)
	#=DONE IMPORTING FILES=#
	#=SETTING UP TIME DELAY FUNCTION=#
	num_tdf_times = 50
	tdf_values = zeros(num_tdf_times)
	tdf_times = collect(1.0:((20.0-1.0)/(num_tdf_times-1)):20.0)

	#= COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS =#
	interpolation_points = zeros(num_spectra_samples,num_tdf_times)
	H = zeros(num_spectra_samples,num_tdf_times)
	HE= zeros(num_spectra_samples,num_tdf_times)
	for date in collect(1:num_spectra_samples)
		for delay in collect(1:num_tdf_times)
			interpolation_points[date,delay]=spectra_dates[date]-tdf_times[delay]
	  end
	  P = interpolation_points[date,:]
	  H[date,:] =interp(interpolation_points[date,:],continuum_dates,continuum_flux)
	  HE[date,:] = interp(interpolation_points[date,:],continuum_dates,continuum_error_flux)
	end







	#=    PRECOMPUTING TIKHONOV MATRICES     =#
	#println("L: ", size(L))
	#println("H: ", size(H))

	#Build Mapping Matrix
	num_spectra_dates=size(spectra_dates)[1]
	W = zeros((size(wavelength)[1],size(L)[2],size(L)[2]))

	
	for lam in collect(1:num_lines)
	  T = eye(num_spectra_dates)
	  for i in collect(1:num_spectra_dates)
	    for j in collect(1:num_spectra_dates)
	      if i == j
					#println("!!", i, " " , j)
	        T[i,j] =1.0/ EL[lam,i].^2
	      end
	    end
	  end
	  W[lam,:,:] = T
	end
	HT = H'
	Gamma = eye(num_tdf_times)
	GammaT = Gamma'
	println("size of W:" ,size(W))
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
	X =   zeros((size(L,1),num_tdf_times))

	
	X_s = zeros((size(L,1),num_tdf_times))
	Z =   zeros((size(L,1),num_tdf_times))
	
	Z_previous = zeros((size(L,1),num_tdf_times))
	Z_s = zeros((size(L,1),num_tdf_times))
	rho = zeros((size(L,1)))
	rho_min = zeros((size(L,1)))
	rho_max = zeros((size(L,1)))
	fill!(rho_min,0.0)
	fill!(rho_max,Inf)
	phi = zeros((size(L,1)))
	phi_previous = zeros((size(L,1)))
	U =   zeros((size(L,1),num_tdf_times))
	tau_prim = zeros((num_lines))
	tau_prim_previous = zeros((num_lines))
	tau_dual = zeros((num_lines))
	tau_dual_previous = zeros((num_lines))
	fill!(X,initial_psi)
	fill!(X_s,initial_psi)
	fill!(Z,initial_psi)
	fill!(Z_s,initial_psi)
	
	X_act = zeros((50))
	fill!(X_act,initial_psi)
	
	#=  Calculate Initial Multipliers & RHO =#
	for p in collect(1:num_lines)
	  U[p,:]=HT * squeeze(W[p,:,:],1) * ( H * vec(Z[p,:]) - vec(L[p,:]))
		A =(vec(U[p,:])'*HT*squeeze(W[p,:,:],1)*H*vec(U[p,:]))
		B =(vec(U[p,:])'*vec(U[p,:]))
		rho[p] = A[1]/B[1]
	end
	println("#############################################################################################")
	println("INITIAL RHO: ", rho)
	println("#############################################################################################")
	println("INITIAL U- ", " min: ",minimum(U)," max: ", maximum(U), " mean: ", mean(U) )
	println("#############################################################################################")
	#= Calculate Initial Rho =#
	converged = 0
	it = 2	#DON'T CHANGE THIS
	while it <= nits && converged==0        #ADMM ITERATION LOOP
	  #println("Iteration: ", it, "  ", sum(Con_Arr), "/", num_lines, " Converged.")
	 for l in collect(1:num_lines)        #SPECTAL CHANNEL LOOP
		#for l in collect(1:2)
	    if Con_Arr[l] != 1
	    #Step 1:
	      X_s[l,:] = Z[l,:]-U[l,:]/rho[l]
	      X[l,:] = ell1_prox_op(X_s[l,:],mu,rho[l])
	    #Step 1 done!
	    #Step 2:
	      Z_s[l,:] = X[l,:]-U[l,:]/rho[l]
	      slice = reshape(W[l,:,:],size(W[l,:,:])[2],size(W[l,:,:])[3])
	      A =inv(HT*slice*H+rho[l]/2.0*GammaT*Gamma)
	      B = HT*slice*vec(L[l,:])
		    Z[l,:] = A * B
	    #Step 2 done!
	    #Step 3:
	      U[l,:] = U[l,:]+(rho[l]/2.0)*(X[l,:]-Z[l,:])
	    #Step 3 done!
	    # Evaluation:
	

				Mod = Model((X[l,:]),H) 
				
			
	      eps_abs = 0.0    #Must be >= 0
	      eps_rel = 0.1  #Must be between 0 and 1
				
				S = rho[l]*(Z[l,:]-Z_previous[l,:])
	      R = X[l,:]-Z[l,:]
				
				S_previous =  S	      
				R_previous =  R
				tau_prim_previous = tau_prim[l]
	      tau_prim[l] = sqrt(num_tdf_times)*eps_abs + eps_rel*(maximum([ell2norm(X[l,:]),ell2norm(Z[l,:])]))
				tau_dual_previous[l] = tau_dual[l]
	      tau_dual[l] = sqrt(num_tdf_times)*eps_abs + eps_rel*ell2norm(U[l,:])
	      eta = ell2norm(R)*tau_dual[l] / (ell2norm(S)*tau_prim[l])		#OPTION 1\
			  #eta = ell2norm(R)*tau_dual_previous[l] / (ell2norm(S)*tau_prim_previous[l])	#OPTION 2

				#if (ell2norm(R)/tau_dual[l]) < (ell2norm(R_previous)/tau_dual_previous[l])
					Acol = :red
				#else
					#Acol = :green
				#end

				#if (ell2norm(S)/tau_prim[l]) < (ell2norm(S_previous)/tau_prim_previous[l])
					Bcol = :red
				#else
					Bcol = :green
				#end
				#Diagnostic Print Statements:
				println("Chi2: ",Chi2(Mod,L[1,:],EL[1,:])/num_tdf_times)
				println("Tau_prim: ", tau_prim[l],"   tau_dual: ", tau_dual[l])
				println("       R: ", ell2norm(R),             "        S: ", ell2norm(S))
				print("R/Tau_pr: ")
				print_with_color(Acol,string(ell2norm(R)/tau_prim[l]))
				print( " S/Tau_du: ")
				print_with_color(Bcol, string(ell2norm(S)/tau_dual[l]))
				print("\n")
				println("Phi: ", phi[l], "  Rho: ", rho[l])
				println("--------------------------------------------------------------------")
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
	  end
	  it = it+1
	end


	X
end 


mu = 2000.0
tmp = TLDR(mu,L, num_lines, EL, spectra_dates, num_spectra_samples,continuum_dates,continuum_scale,continuum_flux,continuum_error_flux,num_continuum_dates)

println(tmp[1,:])





figure(1)
imshow(tmp',aspect="auto")
xlabel("lines")
ylabel("delay")
colorbar()
show()












