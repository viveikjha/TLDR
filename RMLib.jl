using FITSIO
#=--------------------------------------------------=#
#===================== ALTCHI2 ======================#
#=--------------------------------------------------=#
function alt_chi2(DATA,MATS,x)
	ls = MATS.H*x
	total = 0.0
	for l in 1:DATA.num_lines        #SPECTAL CHANNEL LOOP
			W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
			res = (ls-DATA.L)[:,l]
			bit = res'*W_slice*res
			total += sum(bit)
	end
	total
end




#=--------------------------------------------------=#
#====================== CLEAR =======================#
#=--------------------------------------------------=#
function clear()
  for i in collect(1:30)
    println("\n")
  end
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
		sum(   ((M-D)./(Sigma)).^2)
end

#=--------------------------------------------------=#
#============= ELL1 Proximal Operator ===============#
#=--------------------------------------------------=#
#Proximal operator for the L1 norm with positivity.
#X is the current model TDF
#mu is the regularization weight
#rho is the current hyperparameter
function ell1_prox_op(X,XS,mu,rho)
  for i in collect(1:length(X))
    if abs(XS[i]) > mu/rho
			if XS[i] > mu/rho
				X[i] = XS[i] - mu/rho
			else
				X[i] = XS[i]+mu/rho
			end
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
	for i in 1:length(X)
		if X[i] < 0.0
			X[i] = 0.0
		end
	end
	X
end

#=--------------------------------------------------=#
#=================== Ell 1 Norm =====================#
#=--------------------------------------------------=#
function ell1norm(X)
	sum(abs(X))
end
#=--------------------------------------------------=#
#=================== Ell 2 Norm =====================#
#=--------------------------------------------------=#
function ell2norm(X)
	sqrt(sum(X.^2))
end
#=--------------------------------------------------=#
#=================== Ell 2 Norm Squared=====================#
#=--------------------------------------------------=#
function ell2normsquared(X)
	sum(X.^2)
end

#=--------------------------------------------------=#
#=============== Minimization wrt X =================#
#=--------------------------------------------------=#
function min_wrt_x(X,T,P,Z,Pars,DATA,Mats)
	#println("!")
	s = size(X.vdm)
	#println("size: ", s)
	vdm = SharedArray(Float64,s[1],s[2])
	#println(size(vdm))
	#tic()
	for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
			W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
			Q = Mats.HT * W_slice * Mats.H + T.rho*Mats.DsT*Mats.Ds + (Pars.mu_smoo+Z.rho+T.rho)*Mats.Gammatdf #ORIGINAL PAPER VERSION

			B = Mats.HT* W_slice * DATA.L + Mats.DsT*(T.U+T.rho*T.vdm)+P.U+P.rho*P.vdm+Z.U+Z.rho*Z.vdm
			vdm[:,l] = Q\B[:,l]
	end
	#toc()
	X.vdm=sdata(vdm) #sdata() pulls the underlying shared array

	X.vdm


#	QC = T.rho*Mats.DsT*Mats.Ds + (Pars.mu_smoo+Z.rho+T.rho)*Mats.Gammatdf

#  @parallel for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
#			W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
#			Q = Mats.HT * W_slice * Mats.H + QC
#			L_slice =    DATA.L[:,l]
#			DsT_slice = 	reshape(Mats.DsT[l,:],size(Mats.DsT[l,:])[2])
#			T_slice = 	T.vdm[:,l]
#			UT_slice = 	vec(T.U[:,l])
#			P_slice = 	vec(P.vdm[:,l])
#			UP_slice = 	vec(P.U[:,l])
#			Z_slice = 	vec(Z.vdm[:,l])
#			UZ_slice = 	vec(Z.U[:,l])

#			Bl = Mats.HT*W_slice*L_slice + DsT_slice.*(UT_slice+T.rho*T_slice)+UP_slice+P.rho*P_slice+UZ_slice+Z.rho*Z_slice
#			vdm[:,l] = Q\Bl
#	end
#	X.vdm=sdata(vdm) #sdata() pulls the underlying shared array

#	X.vdm
end

function min_wrt_z(X,V,Z,Pars,DATA,Mats)
  	Rinv = inv(V.rho*Mats.Dv*Mats.DvT+Z.rho*Mats.Gammaspe)
    C = (V.U+V.rho*V.vdm)*Mats.DvT-Z.U+(Z.rho*X.vdm)	#ORIGINAL PAPER VERSION
    Z.vdm = 	C*Rinv
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
#function Pen_update(X,Y,P)
#	scale=4.5
#	#println(size(Y.vdm),"-",size(Y.vdm_previous))
#	S = Y.rho*(Y.vdm-Y.vdm_previous)
#	R = X.vdm-Y.vdm
#	tau_prim_previous = Y.tau_prim
#  Y.tau_prim = sqrt(P.num_tdf_times)*P.eps_abs + P.eps_rel*(maximum([ell2norm(X.vdm),ell2norm(Y.vdm)]))
#	tau_dual_previous = Y.tau_dual
#  Y.tau_dual = sqrt(P.num_tdf_times)*P.eps_abs + P.eps_rel*ell2norm(Y.U)
#	eta = ell2norm(R)*Y.tau_dual / (ell2norm(S)*Y.tau_prim)
	#eta = ell2norm(R)*tau_dual_previous[l] / (ell2norm(S)*tau_prim_previous[l])	#OPTION 2
#	phi_previous = Y.phi
#	Y.phi = max(ell2norm(R/Y.tau_prim),ell2norm(S)/Y.tau_dual)
#	if ell2norm(R) > scale*ell2norm(S)
#		Y.rho=Y.rho*P.tau
#	elseif ell2norm(R)*scale < ell2norm(S)
#		Y.rho=Y.rho/P.tau
#	else
		#nothin.
#	end
#	Y
#end


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
#	eta = ell2norm(R)*tau_dual_previous[l] / (ell2norm(S)*tau_prim_previous[l])	#OPTION 2
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
function TLDR(DATA,Mats,Pars,Plot_F="True",Plot_A="False",vdmact="None")
	Pars.tau=2.0
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

	#When loading a known TDF
	if vdmact != "None"
		vdm_path = vdmact
		vdm_act = readcsv(vdm_path)
		init_vdm =vdm_act
	else
			init_vdm=ones(Pars.num_tdf_times,DATA.num_lines)*initial_psi
	end

	#println("!!!!!",sum(init_vdm))
	#init_vdm = Tik_Init(1.0,Mats,DATA)
	#init_vdm =vdm_act

	init_vdm=randn(size(init_vdm)) #Start from Random

	X = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	X.vdm = init_vdm
	Z = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	Z.vdm = init_vdm
	T = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	T.vdm = Mats.Ds*init_vdm
	V = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	V.vdm = init_vdm*Mats.Dv
	P = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	#println("Initial J: ", J(X,P,T,V,DATA,Mats,Pars))
	#	P.vdm = pos_prox_op(init_vdm)
	#Z.rho=1.0
	#T.rho=5.0
	#V.rho=2.0
	#P.rho=1.0

	init_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
	Pars.chi2=init_chi2
	println("initialized chi2: ",init_chi2,"\tInitial J: ", J(X,P,T,V,DATA,Mats,Pars),"\t L2x: ",regX(X,Pars),"\tL1T: ",regT(T,Pars),"\tL1V: ",regV(V,Pars))


	#=  Calculate Initial Multipliers & RHO =#
	for p in 1:DATA.num_lines
	  Z.U[:,p]=Mats.HT * squeeze(Mats.W[p,:,:],1) * ( Mats.H * vec(Z.vdm[:,p]) - vec(DATA.L[:,p]))
	  #Z.U[:,p]=Mats.HT * squeeze(Mats.W[p,:,:],1) * ( vec(DATA.L[:,p])-Mats.H * vec(Z.vdm[:,p]))
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
		X.vdm = min_wrt_x(X,T,P,Z,Pars,DATA,Mats)

	#Step 2: UPDATE REGULARIZATION TERMS
		P.vdm_squiggle = X.vdm - P.U/P.rho
		P.vdm_previous = P.vdm
		P.vdm = pos_prox_op(P.vdm_squiggle)

		T.vdm_squiggle = Mats.Ds*X.vdm-T.U/T.rho
		T.vdm_previous = T.vdm
		T.vdm = ell1_prox_op(T.vdm,T.vdm_squiggle,Pars.mu_temp,T.rho)

		V.vdm_squiggle = Z.vdm*Mats.Dv - V.U/V.rho
		V.vdm_previous = V.vdm
		V.vdm = ell1_prox_op(V.vdm,V.vdm_squiggle,Pars.mu_spec,V.rho)

		Z.vdm_previous = Z.vdm
		Z.vdm = min_wrt_z(X,V,Z,Pars,DATA,Mats)



	#Step 4: UPDATE LAGRANGE MULTIPLIERS
		T.U = LG_update(T.U,T.vdm,Mats.Ds*X.vdm,T.rho,Pars.alpha)
		V.U = LG_update(V.U,V.vdm,Z.vdm*Mats.Dv,V.rho,Pars.alpha)
		P.U = LG_update(P.U,P.vdm,X.vdm,P.rho,Pars.alpha)
		Z.U = LG_update(Z.U,Z.vdm,X.vdm,Z.rho,Pars.alpha)

  	if Pars.it == 2
			Pars.G = 1.5
		end
	#Step 5: Update Penalty Parameters
#		Z = Pen_update(X,Z,Pars)
		Z.rho=200.0
#		T = Pen_update(X,T,Pars)
		T.rho=2000.0 #5.0
#		V = Pen_update(X,V,Pars)
		V.rho=2000.0
#		P = Pen_update(X,P,Pars)
		P.rho=2000.0

	#Step 6: Check for Convergence
		Pars.it = Pars.it+1

		true_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
		chiprev = Pars.chi2
		Pars.chi2= true_chi2



		#Reporting
		chiz = Chi2(Model((Z.vdm),Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)

		Jstring=@sprintf "\tJ: %0.3f" J(X,P,T,V,DATA,Mats,Pars)
		L2xstring=@sprintf "\tL2x: %0.0f" regX(X,Pars)
		L1Tstring=@sprintf "\tL1T: %0.1f" regT(T,Pars)
		L1Vstring=@sprintf "\tL1V: %0.1f" regV(V,Pars)
		chi2xstring=@sprintf "\tChi2x: %0.3f" Pars.chi2
		chi2zstring=@sprintf "\tChi2z: %0.3f" chiz
		sstring=@sprintf "\ts: %0.3f" ell2norm(X.rho*(X.vdm-X.vdm_previous))


		rastring=@sprintf "\trP: %0.1f" ell2norm(P.vdm-X.vdm)
		rbstring=@sprintf "\trZ: %0.1f" ell2norm(Z.vdm-X.vdm)
		rcstring=@sprintf "\trV: %0.1f" ell2norm(V.vdm-Z.vdm*Mats.Dv)
		rdstring=@sprintf "\trT: %0.1f" ell2norm(T.vdm-Mats.Ds*X.vdm)

		rhozstring=@sprintf "\tZro: %0.1f" Z.rho
		rhopstring=@sprintf "\tPro: %0.1f" P.rho
		rhotstring=@sprintf "\tTro: %0.1f" T.rho
		rhovstring=@sprintf "\tVro: %0.1f" V.rho






		println("It: ",Pars.it-2,Jstring,L2xstring,L1Tstring,L1Vstring,chi2xstring,sstring,rastring, rbstring,rcstring,rdstring,rhozstring,rhopstring,rhotstring,rhovstring)
		#println("It: ",Pars.it-2,Jstring,L2xstring,L1Tstring,chi2xstring,sstring,rastring, rbstring,rdstring)
		#if sqdiff < diffbest
			#diffbest = sqdiff
			#difbit = Pars.it
		#end
		if abs(Pars.chi2 - chiprev) < (0.001)
			#converged = 1
			#println("CONVERGED")
		end
		#Plotting
		if Plot_A == "True" && (Pars.it%20 == 3)
			imshow((X.vdm),extent=[minimum(DATA.wavelength),maximum(DATA.wavelength),0.0,50.0],aspect="auto",origin="lower",cmap="Reds",interpolation="None")
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
		imshow((Z.vdm),aspect="auto",origin="lower",interpolation="None")
		title("Z")
		colorbar()
		subplot(222)
		imshow((X.vdm),aspect="auto",origin="lower",interpolation="None")
		title("X")
		colorbar()
		subplot(223)
		imshow((T.vdm),aspect="auto",origin="lower",interpolation="None")
		title("T")
		colorbar()
		subplot(224)
		imshow((V.vdm),aspect="auto",origin="lower",interpolation="None")
		title("V")
		colorbar()
		show()
		#writecsv("UnitTests/Ds.csv",Mats.Ds)
		#writecsv("UnitTests/T.csv",T.vdm)
		#writecsv("UnitTests/X.csv",T.vdm)
		#savefig("Spiral_rec.png")
	end
	X,Pars
end


#=--------------------------------------------------=#
#================== Initial VDM =====================#
#=--------------------------------------------------=#
#function Tik_Init(beta,M,D)
#	gamma = eye(size(M.HT)[1])
#	W_slice = reshape(Mats.W[1,:,:],size(Mats.W[1,:,:])[2],size(Mats.W[1,:,:])[3])
#	x = inv(M.HT*W_slice*M.H+beta.*gamma)*M.HT*D.L
#	x
#end
#=--------------------------------------------------=#
#================== DATA REPORT =====================#
#=--------------------------------------------------=#
function data_report(d)
	println("-----------------------------")
	println("--------- DATA INFO ---------")
	println("-----------------------------")
	println("	Wavelengths: ", size(d.wavelength))
	println("	L: ", size(d.L))
	println("	EL: ", size(d.EL))
	println("	numlines: ", d.num_lines)
	println("	spectra_samples: ", d.num_spectra_samples)
	#println("	spectra dates: ", size(d.spectra_dates))
	println("	continuum dates: ", size(d.continuum_dates))
	#println("	continuum flux: ", size(d.continuum_flux))
	#println("	continuum erflux: ", size(d.continuum_error_flux))
	#println("	num continuum dates: ", d.num_continuum_dates)
	println("-----------------------------")
	println("--------- END  INFO ---------")
	println("-----------------------------")
end

#=--------------------------------------------------=#
#================== THE FUNCTION ====================#
#=--------------------------------------------------=#
function J(X,P,T,V,DATA,Mats,Pars)
	0.5*(Pars.chi2*DATA.num_spectra_samples*DATA.num_lines)+Pars.mu_smoo*0.5*ell2norm(X.vdm).^2+Pars.mu_temp*ell1norm(T.vdm)+Pars.mu_spec*ell1norm(V.vdm)
#	0.5*(Pars.chi2*DATA.num_spectra_samples*DATA.num_lines)+Pars.mu_smoo*0.5*ell2norm(X.vdm).^2+Pars.mu_temp*ell1norm(T.vdm)
end

function regX(X,Pars)
	Pars.mu_smoo*0.5*ell2normsquared(X.vdm)
end

function regT(T,Pars)
	Pars.mu_temp*ell1norm(T.vdm)
end

function regV(V,Pars)
	Pars.mu_spec*ell1norm(V.vdm)
end
