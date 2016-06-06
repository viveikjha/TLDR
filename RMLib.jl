using FITSIO
using PyPlot
include("RMLibMore.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("DataImportNEW.jl")
include("GenMatrices.jl")
#=--------------------------------------------------=#
#================= TLDR HOT LAUNCHER ====================#
#=--------------------------------------------------=#
#1. FILES REQUIRED - FILES_ARR = [WAVELENGTHS,SPECTRA,ERRSPECTRA,DATES,CONTINUUM]
#2. IS THIS A TEST? TRUE VDM FILE? -> OPTIONAL
#3. mus: mu_smoo required. others optional?

#Tvdm file: Optional file imput of the real tdf. Used for testing puroposes.
#mu_spec, mu_l1, and mu_temp flags are options for providing individual regularization weights
#otherwise, the values are based on those for mu_smoo.
#Plot_Live option shows the active reconstruction every so many iterations.
#Plot_Final option shows the final plot that displays reconstruction images X,Z,V,T.
function HOT_LAUNCH(Data,Mats,Pars;mu_smoo=40.0,mu_spec=false,mu_temp=false,mu_l1=false,scale=1.0,nits=50,Tvdm="",Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true,rhoZ=8000.0,rhoN=800.0,rhoP=800.0, rhoV=800.0,rhoT=800.0)

	Data.L=scale*(Data.L)
	Data.EL=scale*(Data.EL)
	Data.continuum_flux=scale*Data.continuum_flux
	Data.continuum_error_flux=scale*Data.continuum_error_flux
	Pars.nits=nits

	#SET RECONSTRUCTION PARAMETERS
	scale=1.0
  if mu_temp != false
    Pars.mu_temp = mu_temp
  else
    Pars.mu_temp = 0.25*mu_smoo/scale
  end
  if mu_spec != false
    Pars.mu_spec = mu_spec
  else
    Pars.mu_spec = 0.25*mu_smoo/scale
  end
  if mu_l1 != false
    Pars.mu_l1 = mu_l1
  else
    Pars.mu_l1 = 0.25*mu_smoo/scale
  end
  Pars.mu_smoo=mu_smoo/scale^2
  tmp,P = TLDR(1.0,DATA,Mats,Pars;Plot_A=Plot_Live,Plot_F=Plot_Final,vdmact=Tvdm,RepIt=RepIt,RepF=RepF,rhoZ=rhoZ,rhoN=rhoN,rhoP=rhoP,rhoT=rhoT,rhoV=rhoV)
  tmp;
end
#=--------------------------------------------------=#
#================= TLDR COLD LAUNCHER ====================#
#=--------------------------------------------------=#
#1. FILES REQUIRED - FILES_ARR = [WAVELENGTHS,SPECTRA,ERRSPECTRA,DATES,CONTINUUM]
#2. IS THIS A TEST? TRUE VDM FILE? -> OPTIONAL
#3. mus: mu_smoo required. others optional?

#Tvdm file: Optional file imput of the real tdf. Used for testing puroposes.
#mu_spec, mu_l1, and mu_temp flags are options for providing individual regularization weights
#otherwise, the values are based on those for mu_smoo.
#Plot_Live option shows the active reconstruction every so many iterations.
#Plot_Final option shows the final plot that displays reconstruction images X,Z,V,T.
function COLD_LAUNCH(FILES_ARR;mu_smoo=40.0,mu_spec=false,mu_temp=false,mu_l1=false,scale=1.0,nits=50,Tvdm="",Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true,rhoZ=8000.0,rhoN=800.0,rhoP=800.0, rhoV=800.0,rhoT=800.0)
  #IMPORT DATA FROM FILES_ARR
	wavelengths=FILES_ARR[1]
	spectra = FILES_ARR[2]
	errspectra = FILES_ARR[3]
	dates = FILES_ARR[4]
	continuum = FILES_ARR[5]
	 DATA = Import_DataN("",wavelengths,spectra,errspectra,dates,continuum)

	DATA.L=scale*(DATA.L)
	DATA.EL=scale*(DATA.EL)
	DATA.continuum_flux=scale*DATA.continuum_flux
	DATA.continuum_error_flux=scale*DATA.continuum_error_flux

	#data_report(DATA)
	 Pars = init_Params()
	Pars.num_tdf_times=50

	 Mats=Gen_Mats(DATA,Pars)

	Pars.nits=nits

	#SET RECONSTRUCTION PARAMETERS
	scale=1.0
  if mu_temp != false
    Pars.mu_temp = mu_temp
  else
    Pars.mu_temp = 0.25*mu_smoo/scale
  end
  if mu_spec != false
    Pars.mu_spec = mu_spec
  else
    Pars.mu_spec = 0.25*mu_smoo/scale
  end
  if mu_l1 != false
    Pars.mu_l1 = mu_l1
  else
    Pars.mu_l1 = 0.25*mu_smoo/scale
  end
  Pars.mu_smoo=mu_smoo/scale^2
  tmp,P = TLDR(1.0,DATA,Mats,Pars;Plot_A=Plot_Live,Plot_F=Plot_Final,vdmact=Tvdm,RepIt=RepIt,RepF=RepF,rhoZ=rhoZ,rhoN=rhoN,rhoP=rhoP,rhoT=rhoT,rhoV=rhoV)
end


#=--------------------------------------------------=#
#================ ADMM Algorithm ====================#
#=--------------------------------------------------=#
function TLDR(flx_scale,DATA,Mats,Pars;Plot_F=true,Plot_A=false,vdmact="",RepIt=true,RepF=true,rhoZ=8000.0,rhoN=800.0,rhoP=800.0,rhoT=800.0,rhoV=800.0,savefits=false)
	Pars.tau=2.0
	threshold = 1.0e-4 #CONVERGANCE THRESHOLD
	CX=false
	CN=false
	CT=false
	CV=false
	if Plot_A == true
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
	#println("a: ",flx_scale^2*Pars.mu_smoo)

	vdm = zeros(Pars.num_tdf_times,DATA.num_lines)
	for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
			W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
			Q = Mats.HT * W_slice * Mats.H +  (Pars.mu_smoo)*Mats.Gammatdf #INCLUCES L1 NORM ON X
			B = Mats.HT* W_slice * DATA.L
			vdm[:,l] = Q\B[:,l]
		end
	Ini=vdm.*(vdm .>= 0.0) #sdata() pulls the underlying shared array
	#Ini = inv(Mats.H'*Mats.H+(flx_scale^2*Pars.mu_smoo)^2*eye(size(Mats.H)[2]))*(Mats.H'*DATA.L) #INITIALIZATION FROM TIKHONOV SOLUTION
	init_vdm =Ini #FILTER OUT NEGATIVE VALUES
	if Plot_A == true
		imshow(init_vdm,aspect="auto",origin="lower",interpolation="None",cmap="Blues")
		#colorbar()
		#show()
	end
	writecsv("tiksol.csv",init_vdm)
	#init_vdm=randn(size(init_vdm)) #Start from Random
	#init_vdm=0.0*randn(size(init_vdm)) #Start from Random

	 X = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 Z = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 T = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 V = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 P = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 N = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)



	#When loading a known TDF
	if vdmact != ""
		vdm_path = vdmact
		vdm_act = readcsv(vdm_path)
		X.vdm = vdm_act
		Z.vdm = vdm_act
		T.vdm = Mats.Ds*vdm_act
		V.vdm = vdm_act*Mats.Dv
		P.vdm = pos_prox_op(vdm_act)
		N.vdm = vdm_act
		act_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
		#println("actual chi2: ",act_chi2)
		Pars.chi2=act_chi2
		if RepIt==true
			Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=true,s=false,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,Msg=" -True_VDM-")
		end
	end

	X.vdm = init_vdm
	Z.vdm = init_vdm
	T.vdm = Mats.Ds*init_vdm
	V.vdm = init_vdm*Mats.Dv
	P.vdm = pos_prox_op(init_vdm)
	N.vdm = init_vdm

	#Initiailize Penalty Parameters.
	#flx_scale= 5.0
	#flx_scale=1.0/10000.0
	#println("Z: ",Z.rho,"P: ",P.rho,"T: ",T.rho,"V: ",V.rho,"N: ",N.rho)
	Z.rho= rhoZ #4000000.0/flx_scale^2#20.0*Pars.mu_smoo#2000.0/flx_scale#2000.0
	P.rho=rhoP#200.00/flx_scale^2
	T.rho=rhoT#200.00/flx_scale^2
	V.rho=rhoV#200.00/flx_scale^2
	#N.rho=4000000.0/flx_scale^2
	N.rho=rhoN#2000.00/flx_scale^2

	siglvl=abs(median(DATA.L))
	#println("Signal Level: ",siglvl, " - ", Pars.mu_smoo/Z.rho)
	#println("ρZ:", Z.rho," μsmoo/ρZ: ",Pars.mu_smoo/Z.rho, " ρN:", N.rho, " μ l1: ",Pars.mu_l1, " μl1/ρN: ", Pars.mu_l1/N.rho)
#	T.rho=2000#20.0*Pars.mu_temp#2000.0/flx_scale#2000.0 #5.0
#	V.rho=2000#20.0*Pars.mu_spec#2000.0/flx_scale#2000.0
#	P.rho=2000#2000.0/(flx_scale)#2000.0
#	N.rho=2000#20.0*Pars.mu_l1#2000.0/flx_scale

#Diagnostics on VDM initialized with Tihonov Solution.
	init_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
	Pars.chi2=init_chi2
	if RepIt==true
		Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=true,s=false,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,Msg=" -Tik_Init-")
	end
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
	converged = false
	while Pars.it <= Pars.nits && Pars.conflag==false        #ADMM ITERATION LOOP
		X.vdm_previous = copy(X.vdm)
	#Step 1: MINIMIZATION W.R.T. X
		X.vdm = min_wrt_x(X,T,P,N,Z,Pars,DATA,Mats)
		#X.vdm = X.vdm.*(X.vdm.>0.0)
	#Step 2: UPDATE REGULARIZATION TERMS
		P.vdm_squiggle = X.vdm - P.U/P.rho
		P.vdm_previous = copy(P.vdm)
		P.vdm = pos_prox_op(P.vdm_squiggle)

		T.vdm_squiggle = Mats.Ds*X.vdm-T.U/T.rho
		T.vdm_previous = copy(T.vdm)

		T.vdm = ell1_prox_op(T.vdm,T.vdm_squiggle,Pars.mu_temp,T.rho)
		V.vdm_squiggle = Z.vdm*Mats.Dv - V.U/V.rho
		V.vdm_previous = copy(V.vdm)
		V.vdm = ell1_prox_op(V.vdm,V.vdm_squiggle,Pars.mu_spec,V.rho)

		N.vdm_squiggle = X.vdm - N.U/N.rho
		N.vdm_previous = copy(N.vdm)
		N.vdm = ell1_prox_op(N.vdm,N.vdm_squiggle,Pars.mu_l1,N.rho)

		Z.vdm_previous = copy(Z.vdm)
		Z.vdm = min_wrt_z(X,V,Z,Pars,DATA,Mats)
		#Z.vdm = Z.vdm.*(Z.vdm.>0.0)
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
#Z.rho=#800.0/flx_scale^2#20.0*Pars.mu_smoo#2000.0/flx_scale#2000.0
#T.rho=#8000.0/flx_scale^2#20.0*Pars.mu_temp#2000.0/flx_scale#2000.0 #5.0
#V.rho=#8000.0/flx_scale^2#20.0*Pars.mu_spec#2000.0/flx_scale#2000.0
#P.rho=#8000.0/flx_scale^2#2000.0/(flx_scale)#2000.0
#N.rho=#8000.0/flx_scale^2#20.0*Pars.mu_l1#2000.0/flx_scale

	#Step 6: Check for Convergence
		if Pars.it != 2
			if abs(regX(X,Pars)-PregX) < threshold
				CX=true
			end
			if abs(regN(N,Pars)-PregN) < threshold
				CN=true
			end
			if abs(regT(T,Pars)-PregT) < threshold
				CT=true
			end
			if abs(regV(V,Pars)-PregV) < threshold
				CV=true
			end
		end
		if CX == true && CN == true && CT == true && CV == true
			Pars.conflag = true
			print_with_color(:blue,"TLDR CONVERGED \n")
		end
		PregX=regX(X,Pars)
		PregN=regN(N,Pars)
		PregT=regT(T,Pars)
		PregV=regV(V,Pars)
		Pars.it = Pars.it+1

		true_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
		chiprev = Pars.chi2
		Pars.chi2= true_chi2

		#Reporting
		if RepIt==true
			Report(X,Z,P,T,V,N,DATA,Mats,Pars,Jf=true,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,s=true,Pres=true,Zres=true,Tres=true,Vres=true,Nres=true)
		end
		#println("It: ",Pars.it-2,Jstring,L2xstring,L1Tstring,L1Vstring,chi2xstring,sstring,rastring, rbstring,rcstring,rdstring,restring,rhozstring,rhopstring,rhotstring,rhovstring,rhonstring)
		#println("It: ",Pars.it-2,Jstring,L2xstring,L1Tstring,chi2xstring,sstring,rastring, rbstring,rdstring)
		#if sqdiff < diffbest
			#diffbest = sqdiff
			#difbit = Pars.it
		#end

		#Plotting
		if Plot_A == true && (Pars.it%10 == 0)
			clf()
			imshow((X.vdm),extent=[minimum(DATA.wavelength),maximum(DATA.wavelength),0.0,50.0],cmap="Reds",aspect="auto",origin="lower",interpolation="None")
			colorbar()
			draw()
		end
  end
	if RepF==true
		Report(X,Z,P,T,V,N,DATA,Mats,Pars,Jf=true,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,s=true,Pres=true,Zres=true,Tres=true,Vres=true,Nres=true,Msg=" -Final-")
	end

	if Plot_A == true
		ioff()
	end
	#Final Plot
	plotfin(Plot_F,X,Z,T,V)
	if Pars.it==Pars.nits-1 && Pars.conflag==false
		print_with_color(:red,"TLDR FAILED TO CONVERGE \n")
	end
	if savefits==true
		Write_FITS(X,P);
	end
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
