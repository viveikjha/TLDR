using FITSIO
include("RMLibMore.jl")

#=--------------------------------------------------=#
#================ ADMM Algorithm ====================#
#=--------------------------------------------------=#
function TLDR(flx_scale,DATA,Mats,Pars,Plot_F="True",Plot_A="False",vdmact="None")
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

	#INITIALIZATION FROM TIKHONOV SOLUTION
println("a: ",flx_scale^2*Pars.mu_smoo)

vdm = zeros(Pars.num_tdf_times,DATA.num_lines)
for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
		W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
		Q = Mats.HT * W_slice * Mats.H +  (Pars.mu_smoo)*Mats.Gammatdf #INCLUCES L1 NORM ON X
		B = Mats.HT* W_slice * DATA.L
		vdm[:,l] = Q\B[:,l]
end
Ini=vdm #sdata() pulls the underlying shared array




	#Ini = inv(Mats.H'*Mats.H+(flx_scale^2*Pars.mu_smoo)^2*eye(size(Mats.H)[2]))*(Mats.H'*DATA.L) #INITIALIZATION FROM TIKHONOV SOLUTION
	init_vdm =Ini.*(Ini.>0.0) #FILTER OUT NEGATIVE VALUES
	if Plot_A == "True"
		imshow(init_vdm,aspect="auto",origin="lower",cmap="Reds")
		#colorbar()
		#show()
	end

	#init_vdm=randn(size(init_vdm)) #Start from Random
	#init_vdm=0.0*randn(size(init_vdm)) #Start from Random

	X = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	Z = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	T = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	V = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	P = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	N = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)

	#When loading a known TDF
	if vdmact != "None"
		vdm_path = vdmact
		vdm_act = readcsv(vdm_path)
		X.vdm = vdm_act
		Z.vdm = vdm_act
		T.vdm = Mats.Ds*vdm_act
		V.vdm = vdm_act*Mats.Dv
		P.vdm = pos_prox_op(vdm_act)
		N.vdm = vdm_act
		act_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
		Pars.chi2=act_chi2
		Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=true,s=false,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,Msg=" -True_VDM-")
	end

	X.vdm = init_vdm
	Z.vdm = init_vdm
	T.vdm = Mats.Ds*init_vdm
	V.vdm = init_vdm*Mats.Dv
	P.vdm = pos_prox_op(init_vdm)
	N.vdm = init_vdm

	#Initiailize Penalty Parameters.
	#flx_scale= 5.0
	#flx_scale=1.0
	println("Z: ",Z.rho,"P: ",P.rho,"T: ",T.rho,"V: ",V.rho,"N: ",N.rho)
	Z.rho=4000000.0/flx_scale^2#20.0*Pars.mu_smoo#2000.0/flx_scale#2000.0
	P.rho=200.00/flx_scale^2
	T.rho=200.00/flx_scale^2
	V.rho=200.00/flx_scale^2
	N.rho=200.00/flx_scale^2
#	T.rho=2000#20.0*Pars.mu_temp#2000.0/flx_scale#2000.0 #5.0
#	V.rho=2000#20.0*Pars.mu_spec#2000.0/flx_scale#2000.0
#	P.rho=2000#2000.0/(flx_scale)#2000.0
#	N.rho=2000#20.0*Pars.mu_l1#2000.0/flx_scale

#Diagnostics on VDM initialized with Tihonov Solution.
	init_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
	Pars.chi2=init_chi2
	Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=true,s=false,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,Msg=" -Tik_Init-")

	#=  Calculate Initial Multipliers & RHO =#
	for p in 1:DATA.num_lines
	  Z.U[:,p]=Mats.HT * squeeze(Mats.W[p,:,:],1) * ( Mats.H * vec(Z.vdm[:,p]) - vec(DATA.L[:,p]))
	end
	#println("INITIAL MULTIPLIERS: ",mean(Z.U))
	Z.U=zeros(size(Z.U))
	V.U = Z.U
	T.U = Z.U
	P.U = Z.U
	N.U = Z.U
	Qinv = zeros(Pars.num_tdf_times,Pars.num_tdf_times)
	B = zeros(Pars.num_tdf_times,DATA.num_lines)
	converged = 0
	while Pars.it <= Pars.nits && converged==0        #ADMM ITERATION LOOP
		X.vdm_previous = X.vdm
	#Step 1: MINIMIZATION W.R.T. X
		X.vdm = min_wrt_x(X,T,P,N,Z,Pars,DATA,Mats)

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

		N.vdm_squiggle = X.vdm - N.U/N.rho
		N.vdm_previous = N.vdm
		N.vdm = ell1_prox_op(N.vdm,N.vdm_squiggle,Pars.mu_l1,N.rho)

		Z.vdm_previous = Z.vdm
		Z.vdm = min_wrt_z(X,V,Z,Pars,DATA,Mats)

		#Step 4: UPDATE LAGRANGE MULTIPLIERS
		T.U = LG_update(T.U,T.vdm,Mats.Ds*X.vdm,T.rho,Pars.alpha)
		V.U = LG_update(V.U,V.vdm,Z.vdm*Mats.Dv,V.rho,Pars.alpha)
		N.U = LG_update(N.U,N.vdm,X.vdm,N.rho,Pars.alpha)
		P.U = LG_update(P.U,P.vdm,X.vdm,P.rho,Pars.alpha)
		Z.U = LG_update(Z.U,Z.vdm,X.vdm,Z.rho,Pars.alpha)

  	if Pars.it == 2
			Pars.G = 1.5
		end
	#Step 5: Update Penalty Parameters
#		Z = Pen_update(X,Z,Pars)
#		T = Pen_update(X,T,Pars)
#		V = Pen_update(X,V,Pars)
#		P = Pen_update(X,P,Pars)
Z.rho=800.0/flx_scale^2#20.0*Pars.mu_smoo#2000.0/flx_scale#2000.0
T.rho=8000.0/flx_scale^2#20.0*Pars.mu_temp#2000.0/flx_scale#2000.0 #5.0
V.rho=8000.0/flx_scale^2#20.0*Pars.mu_spec#2000.0/flx_scale#2000.0
P.rho=8000.0/flx_scale^2#2000.0/(flx_scale)#2000.0
N.rho=8000.0/flx_scale^2#20.0*Pars.mu_l1#2000.0/flx_scale

	#Step 6: Check for Convergence
		Pars.it = Pars.it+1

		true_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
		chiprev = Pars.chi2
		Pars.chi2= true_chi2

		#Reporting
		Report(X,Z,P,T,V,N,DATA,Mats,Pars,Jf=true,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,s=true,Pres=true,Zres=true,Tres=true,Vres=true,Nres=true)

		#println("It: ",Pars.it-2,Jstring,L2xstring,L1Tstring,L1Vstring,chi2xstring,sstring,rastring, rbstring,rcstring,rdstring,restring,rhozstring,rhopstring,rhotstring,rhovstring,rhonstring)
		#println("It: ",Pars.it-2,Jstring,L2xstring,L1Tstring,chi2xstring,sstring,rastring, rbstring,rdstring)
		#if sqdiff < diffbest
			#diffbest = sqdiff
			#difbit = Pars.it
		#end
		if abs(Pars.chi2 - chiprev) < (0.000001)
			#converged = 1
			#println("CONVERGED")
		end
		#Plotting
		if Plot_A == "True" && (Pars.it%10 == 0)
			imshow((X.vdm),extent=[minimum(DATA.wavelength),maximum(DATA.wavelength),0.0,50.0],aspect="auto",origin="lower",cmap="Reds",interpolation="None")
			draw()
		end
  end
	if Plot_A == "True"
		ioff()
	end
	#Final Plot
	plotfin(Plot_F,X,Z,T,V)

	X,Pars
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
#=============== Minimization wrt X =================#
#=--------------------------------------------------=#
function min_wrt_x(X,T,P,N,Z,Pars,DATA,Mats)
	s = size(X.vdm)
	vdm = SharedArray(Float64,s[1],s[2])
	for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
			W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
			Q = Mats.HT * W_slice * Mats.H + T.rho*Mats.DsT*Mats.Ds + (Pars.mu_smoo+Z.rho+T.rho+N.rho)*Mats.Gammatdf #INCLUCES L1 NORM ON X
			B = Mats.HT* W_slice * DATA.L + Mats.DsT*(T.U+T.rho*T.vdm)+P.U+P.rho*P.vdm+Z.U+Z.rho*Z.vdm+N.U+N.rho*N.vdm #INCLUDES L1 NORM ON X
			vdm[:,l] = Q\B[:,l]
	end
	X.vdm=sdata(vdm) #sdata() pulls the underlying shared array
	X.vdm
end

#=--------------------------------------------------=#
#=============== Minimization wrt Z =================#
#=--------------------------------------------------=#
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
