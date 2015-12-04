using FITSIO
#=--------------------------------------------------=#
#====================== CLEAR =======================#
#=--------------------------------------------------=#
function clear()
  for i in collect(1:30)
    println("\n")
  end
end


function print_rgb(r, g, b, t)
           println("\e[1m\e[38;2;$r;$g;$b;249m",t)
end
#=--------------------------------------------------=#
#==================INTERPOLATION=====================#
#=--------------------------------------------------=#
#interp takes x and y data and desired x points
#returns linearally interpolated y desired points
#extrapolates off either end of x data using the
#slope of the first and last sets of data points given.
function interp(x_DESIRED,x_DATA,y_DATA)
  y_DESIRED = zeros(length(x_DESIRED))
	index_DATA = collect(1:length(x_DATA))
  index_DESIRED = collect(1:length(x_DESIRED))
  minx = minimum(x_DATA)
  maxx = maximum(x_DATA)
  for i in index_DESIRED
    yset=0
    y_NEW =0
    for y in index_DATA
      if yset ==0 && x_DESIRED[i] > x_DATA[y] && x_DESIRED[i] < x_DATA[y+1]  y < (length(index_DATA)-1) && i < (length(index_DATA)-1)
        rise = y_DATA[y+1]-y_DATA[y]
        run = x_DATA[y+1]-x_DATA[y]
        slope = rise / run
        y_change = slope*(x_DESIRED[i]-x_DATA[y])
        y_NEW = y_DATA[y]+y_change
        yset = 1
      elseif yset ==0 && x_DESIRED[i] < minx #EXTRAPOLATION DOWNWARD
        #println("EXTRAP LOW!!")
        rise =y_DATA[2]-y_DATA[1]
        run = x_DATA[2]-x_DATA[1]
        first_slope= rise / run
        y_change = first_slope*(x_DATA[1]-x_DESIRED[i])
        y_NEW = y_DATA[1]-y_change
        yset = 1
      elseif yset ==0 && x_DESIRED[i] > maxx #EXTRAPOLATION UPWARD
        #println("EXTRAP HIGH!!")
        rise=y_DATA[length(y_DATA)]-y_DATA[length(y_DATA)-1]
        run = x_DATA[length(x_DATA)]-x_DATA[length(x_DATA)-1]
        last_slope = rise / run
        y_change = last_slope*(x_DESIRED[i]-x_DATA[length(x_DATA)])
        y_NEW = y_DATA[length(y_DATA)]+ y_change
        yset = 1
      end
      y_DESIRED[i] = y_NEW
    end
  end
  y_DESIRED #This is the returned value. In Julia, return statements are not required.
end
#=--------------------------------------------------=#
#============== Chi Squared Gradient ================#
#=--------------------------------------------------=#
#Z is the current TDF
#D is the spectral data
#Sigma is the error on D
#Rho is the regularization weight
#ICF is the interpolated continuum flux
function chigrad(Z,D,Sigma,rho,ICF)
  M = Model(Z,ICF)
  chi2grad = zeros(size(Z))
  for i in collect(1:size(Z)[1])
    #chi2grad[i] = 2.0*sum((M-D)/Sigma^2. .* ICF[:,i])
    chi2grad[i] = sum(2.0*sum(D - sum(vec(Z).*vec(ICF[i,:])))-sum(ICF,2))
		chi2grad[i] = sum(2.0*sum(D - sum(Z.* H))-sum(ICF))
  end
  chi2grad #Array Returned!
end


#=--------------------------------------------------=#
#====================== Model =======================#
#=--------------------------------------------------=#
#X is the current TDF being modeled.
#ICF is the interpolated continuum flux.
#function Model(X,ICF)
function Model(X,H)
  #L = size(ICF)[1]
  #MF = zeros(L)
  #MF = H*vec(X)
	MF = H*X
  MF #Array Returned!
end

#=--------------------------------------------------=#
#====================== Chi^2 =======================#
#=--------------------------------------------------=#
# M should be the model from the Model function
# D should be the spectral data
# Sigma should be the error on the spectral data
function Chi2(M,D,Sigma)
	  #sum(   (vec(M)-vec(D)).^2  ./vec(Sigma) .^2)   #Value Returned!
		sum(   (M-D).^2 ./(Sigma).^2)
		#sum(((vec(M)-vec(D))/vec(Sigma))^2)
end

#=--------------------------------------------------=#
#============= ELL1 Proximal Operator ===============#
#=--------------------------------------------------=#
#Proximal operator for the L1 norm with positivity.
#X is the current model TDF
#mu is the regularization weight
#rho is the current hyperparameter
function ell1_prox_op(X,mu,rho)
  for i in collect(1:length(X))
    if X[i] > mu/rho
      X[i] = X[i] - mu/rho
    else
      X[i] = 0
    end
  end	
	X #array returned
end

#=--------------------------------------------------=#
#============= ELL2 Proximal Operator ===============#
#=--------------------------------------------------=#
#Proximal operator for the L2 norm.
#X is the current model TDF
#mu is the regularization weight
function ell2_prox_op(X,mu,rho)
	X_returned = 1.0/(1.0+2.0*mu)*X
end

#=--------------------------------------------------=#
#========= Positivity Proximal Operator  ============#
#=--------------------------------------------------=#
function pos_prox_op(X)
	for i in collect(1:size(X)[1])
		for j in collect(1:size(X)[2])
			if X[i,j] < 0.0
				X[i,j] = 0.0
			end
		end	
	end
	X
end

#=--------------------------------------------------=#
#=================== Ell 2 Norm =====================#
#=--------------------------------------------------=#
function ell2norm(X)
	sqrt(sum(X).^2)
end

#=--------------------------------------------------=#
#=============== Update Multipliers =================#
#=--------------------------------------------------=#
#X and Y are IMAGE structs
#alpha is a throttling term on multiplier update.
function LG_update(U,Z,X,rho,alpha)
	U = U + (rho/alpha)*(Z - X)
	Z
end

#=--------------------------------------------------=#
#================= Update Penalty ===================#
#=--------------------------------------------------=#
function Pen_update(X,Y,P)
	S = Y.rho*(Y.vdm-Y.vdm_previous)
	R = X.vdm-Y.vdm
	tau_prim_previous = Y.tau_prim
  Y.tau_prim = sqrt(P.num_tdf_times)*P.eps_abs + P.eps_rel*(maximum([ell2norm(X.vdm),ell2norm(Y.vdm)]))
	tau_dual_previous = Y.tau_dual
  Y.tau_dual = sqrt(P.num_tdf_times)*P.eps_abs + P.eps_rel*ell2norm(Y.U)
	eta = ell2norm(R)*Y.tau_dual / (ell2norm(S)*Y.tau_prim)		
	#eta = ell2norm(R)*tau_dual_previous[l] / (ell2norm(S)*tau_prim_previous[l])	#OPTION 2
	phi_previous = Y.phi
	Y.phi = max(ell2norm(R/Y.tau_prim),ell2norm(S)/Y.tau_dual)
	if (1.0/P.tau[1] <= eta[1]) && (eta[1] <= P.tau[1]) || (Y.phi[1] < P.sigma[1]*phi_previous[1])
		#NOTHING HAPPENS IN HERE	
	elseif eta < (1.0/P.tau)
		Y.rho_max = Y.rho
		if Y.rho_min > 0.0
			Y.rho = sqrt(Y.rho_min*Y.rho_max)
		else Y.rho = Y.rho_max/P.G 
		end
	elseif eta > P.tau
		Y.rho_min = Y.rho
		if Y.rho_max < Inf
			Y.rho = sqrt(Y.rho_min*Y.rho_max)
		else Y.rho = Y.rho_min*P.G
		end
	end
	Y
end

#=--------------------------------------------------=#
#============= Generate ADMM Variable ===============#
#=--------------------------------------------------=#
#Intializes and fills the regularization terms for ADMM
function Gen_Var(rhoi, num_tdf_times,num_lines,psi)
	x = init_IMAGE(rhoi)
	x.vdm = zeros((num_tdf_times,num_lines))
	fill!(x.vdm,psi)
	x.vdm_squiggle = zeros((num_tdf_times,num_lines))
	fill!(x.vdm_squiggle,psi)
	x.U=zeros((num_tdf_times,num_lines))
	x
end

#=--------------------------------------------------=#
#================ Write FITS FILE ===================#
#=--------------------------------------------------=#
#Writes VDM and relevent info to FITS file & Header
function Write_FITS(X,P)
	Date = ["Date",string(Dates.today()),""]
	it = ["it", P.it, "Iterations to convergance"]
	nit = ["max_its", P.nits, "Max number of iterations allowed."]
	mu_smoo = ["mu_smo", P.mu_smoo,"Weight of smoothing regularizer."]
	mu_spec = ["mu_spe", P.mu_spec, "Weight of spectral regularizer."]
	mu_temp = ["mu_tem", P.mu_temp, "Weight of the temporal regularizer."]
	eps_abs = ["eps_abs",P.eps_abs,""]
	eps_rel = ["eps_rel",P.eps_rel,""]
	sigma = ["sigma", P.sigma,""]
	G = ["G",P.G, ""]
	alpha = ["alpha",P.alpha,"Mulitplier Throttle"]
	rho0 = ["in_rho", P.rho0,"Initial Penalty"]
	tau = ["tau", P.tau,""]
	chi2 = ["Chi2", P.chi2,"Final Chi Squared"]
	keys = [Date[1],it[1],nit[1],mu_smoo[1],mu_spec[1],mu_temp[1],eps_abs[1],eps_abs[1],sigma[1],G[1],alpha[1],rho0[1],tau[1],chi2[1]]
	vals = [Date[2],it[2]-1,nit[2],mu_smoo[2],mu_spec[2],mu_temp[2],eps_abs[2],eps_abs[2],sigma[2],G[2],alpha[2],rho0[2],tau[2],chi2[2]]
	coms = [Date[3],it[3],nit[3],mu_smoo[3],mu_spec[3],mu_temp[3],eps_abs[3],eps_abs[3],sigma[3],G[3],alpha[3],rho0[3],tau[3],chi2[3]]
	head = FITSHeader(keys,vals,coms)
	name=string(now(),"vdm.fits")
	f = FITS(name,"w")
	write(f,X.vdm,header=head)
	println("VDM saved as ",name)
end

#=--------------------------------------------------=#
#================ ADMM Algorithm ====================#
#=--------------------------------------------------=#
function TLDR(DATA,Mats,Pars,Plot_F="True",Plot_A="False")
	if Plot_A == "True"
		ion()
		figure()

		title("Active Reconstruction")
		xlabel("Wavelength")
		ylabel("Delay")
	end
	Pars.it = 2
	chibest=Inf
	chibit=1
	diffbest=10.0
	difbit=1
	#= INITIALIZING ADMM PARAMETERS =#
	initial_psi = 0.0  #Initial value for PSI
	#ADMM ARRAYS:
	Con_Arr = zeros(DATA.num_lines)
	#IMAGE ARRAYS
	init_vdm = Tik_Init(1.0,Mats,DATA)
	X = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	X.vdm = init_vdm
	Z = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	Z.vdm = init_vdm
	T = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	T.vdm = Mats.Ds*init_vdm
	V = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	V.vdm = init_vdm*Mats.Dv
	P = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	P.vdm = pos_prox_op(init_vdm)
	#=  Calculate Initial Multipliers & RHO =#
	for p in 1:DATA.num_lines
	  Z.U[:,p]=Mats.HT * squeeze(Mats.W[p,:,:],1) * ( Mats.H * vec(Z.vdm[:,p]) - vec(DATA.L[:,p]))
#		Z.U[:,p]=Mats.HT * ( Mats.H * vec(Z.vdm[:,p]) - vec(DATA.L[:,p]))
	end
	V.U = Z.U
	T.U = Z.U
	P.U = Z.U
	Qinv = zeros(Pars.num_tdf_times,Pars.num_tdf_times)
	B = zeros(Pars.num_tdf_times,DATA.num_lines)
	converged = 0
	while Pars.it <= Pars.nits && converged==0        #ADMM ITERATION LOOP
		X.vdm_previous = X.vdm
	#Step 1: MINIMIZATION W.R.T. X
	  for l in 1:DATA.num_lines        #SPECTAL CHANNEL LOOP
			W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
			Q = Mats.HT * W_slice * Mats.H + T.rho*Mats.DsT*Mats.Ds + (Pars.mu_smoo+Z.rho+T.rho)*Mats.Gammatdf 
#			Q = Mats.HT*Mats.H + T.rho*Mats.DsT*Mats.Ds + (Pars.mu_smoo+Z.rho+T.rho)*Mats.Gammatdf 

			B = Mats.HT* W_slice * DATA.L + Mats.DsT*(T.U+T.rho*T.vdm)+P.U+P.rho*P.vdm+Z.U+Z.rho*Z.vdm
#			B = Mats.HT* W_slice * DATA.L + Mats.DsT*(T.U+T.rho*T.vdm)+Z.U+Z.rho*Z.vdm
#			B = Mats.HT*DATA.L + Mats.DsT*(T.U+T.rho*T.vdm)+Z.U+Z.rho*Z.vdm

			X.vdm[:,l] = Q\B[:,l] 
		end
		X.vdm = pos_prox_op(X.vdm)
	#Step 2: UPDATE REGULARIZATION TERMS
		P.vdm_squiggle = X.vdm - P.U/P.rho	
		P.vdm_previous = P.vdm
		P.vdm = pos_prox_op(P.vdm_squiggle)

		T.vdm_squiggle = Mats.Ds*X.vdm-T.U/T.rho
		T.vdm_previous = T.vdm
		T.vdm = ell2_prox_op(T.vdm_squiggle,Pars.mu_temp,T.rho)

		V.vdm_squiggle = Z.vdm*Mats.Dv - V.U/V.rho
		V.vdm_previous = V.vdm
		V.vdm = ell2_prox_op(V.vdm_squiggle,Pars.mu_spec,V.rho)
	
		Rinv = inv(V.rho*Mats.Dv*Mats.DvT+Z.rho*Mats.Gammaspe)
		C = (V.U+V.rho*V.vdm)*Mats.DvT-Z.U+(Z.rho*X.vdm)
		Z.vdm_previous = Z.vdm
		Z.vdm = 	C*Rinv
	#Step 4: UPDATE LAGRANGE MULTIPLIERS	
		T.U = LG_update(T.U,T.vdm,Mats.Ds*X.vdm,T.rho,Pars.alpha)
		V.U = LG_update(V.U,V.vdm,Z.vdm*Mats.Dv,V.rho,Pars.alpha)
		P.U = LG_update(P.U,P.vdm,X.vdm,P.rho,Pars.alpha)
		Z.U = LG_update(Z.U,Z.vdm,X.vdm,Z.rho,Pars.alpha)
		
  	if Pars.it == 2 
			Pars.G = 1.5
		end
	#Step 5: Update Penalty Parameters
		Z = Pen_update(X,Z,Pars)
		T = Pen_update(X,T,Pars)
		V = Pen_update(X,V,Pars)
		P = Pen_update(X,P,Pars)
	#Step 6: Check for Convergence
		Pars.it = Pars.it+1
		Mod = Model(X.vdm,Mats.H) 
		chiprev=Pars.chi2
		Pars.chi2 = Chi2(Mod,DATA.L,DATA.EL)/(Pars.num_tdf_times*DATA.num_lines)
		#sqdiff = sum((X.vdm-vdm_act).^2)
		#println("Iteration: ",Pars.it-1, " Chi2: ", chi2, " SQDiff: ", sqdiff)
		println("Iteration: ",Pars.it-1, " Chi2: ", Pars.chi2)
		#if sqdiff < diffbest
			#diffbest = sqdiff
			#difbit = Pars.it
		#end
		if abs(Pars.chi2 - chiprev) < (0.0001)
			converged = 1
			#println("CONVERGED")
		end
		if Plot_A == "True"
			imshow(log(X.vdm),extent=[minimum(DATA.wavelength),maximum(DATA.wavelength),0.0,60.0],aspect="auto",origin="lower",cmap="YlOrRd")
			draw()
		end		
  end
	if Plot_A == "True"
		ioff()
	end
	#Final Plot
	if Plot_F == "True"
		figure(figsize=(20,10))
		title("Final Reconstruction")
		set_cmap("YlOrRd")
		subplot(221)
		#figure("Z")	
		imshow((Z.vdm),aspect="auto",origin="lower")
		title("Z")
		colorbar()
		subplot(222)
		#figure("X")
		imshow((X.vdm),aspect="auto",origin="lower")
		title("X")
		colorbar()
		subplot(223)
		#figure("T")
		imshow((T.vdm),aspect="auto",origin="lower")
		title("T")
		colorbar()
		subplot(224)
		#figure("V")
		imshow((V.vdm),aspect="auto",origin="lower")
		title("V")
		colorbar()
		show()
	end
	X,Pars
end


#=--------------------------------------------------=#
#================== Initial VDM =====================#
#=--------------------------------------------------=#
function Tik_Init(beta,M,D)
	gamma = eye(size(M.HT)[1])
	W_slice = reshape(Mats.W[1,:,:],size(Mats.W[1,:,:])[2],size(Mats.W[1,:,:])[3])
	x = inv(M.HT*W_slice*M.H+beta.*gamma)*M.HT*D.L
	x
end


