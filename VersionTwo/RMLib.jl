using FITSIO
using PyPlot
include("RMLibMore.jl")
include("RMTypes.jl")
#include("DataImport.jl")
include("DataImportNEW.jl")
include("GenMatrices.jl")
#=--------------------------------------------------=#
#================= TLDR HOT LAUNCHER ====================#
#=--------------------------------------------------=#
#1. FILES REQUIRED - FILES_ARR = [WAVELENGTHS,SPECTRA,ERRSPECTRA,DATES,CONTINUUM]
#2. IS THIS A TEST? TRUE VDM FILE? -> OPTIONAL
#3. mus: mu_l2 required. others optional?

#Tvdm file: Optional file imput of the real tdf. Used for testing puroposes.
#mu_spec, mu_l1, and mu_temp flags are options for providing individual regularization weights
#otherwise, the values are based on those for mu_l2.
#Plot_Live option shows the active reconstruction every so many iterations.
#Plot_Final option shows the final plot that displays reconstruction images X,Z,V,T.
function HOT_LAUNCH(Data,Mats,Pars;mu_l2=40.0,mu_spec=false,mu_temp=false,mu_l1=false,mu_p=false,nits=50,Tvdm="",Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true,rhoZ=8000.0,rhoN=800.0,rhoP=800.0, rhoV=800.0,rhoT=800.0)

	Data.L=scale*(Data.L)
	Data.EL=scale*(Data.EL)
	Data.continuum_flux=scale*Data.continuum_flux
	Data.continuum_error_flux=scale*Data.continuum_error_flux
	Pars.nits=nits

	#SET RECONSTRUCTION PARAMETERS

  if mu_temp != false
    Pars.mu_temp = mu_temp
		Tt=true
	else
		Tt=false
  end
  if mu_spec != false
    Pars.mu_spec = mu_spec
		Vt=true
	else
		Vt=false
  end
  if mu_l1 != false
    Pars.mu_l1 = mu_l1
		Nt=true
	else
		Nt=false
  end
	if mu_p == true
		Pt=true
	else
		Pt=false
	end
  tmp,P = TLDR(1.0,DATA,Mats,Pars;Plot_A=Plot_Live,Plot_F=Plot_Final,vdmact=Tvdm,RepIt=RepIt,RepF=RepF,rhoZ=rhoZ,rhoN=rhoN,rhoP=rhoP,rhoT=rhoT,rhoV=rhoV,T_toggle=Tt,V_toggle=Vt,N_toggle=Nt,P_toggle=Pt)
  tmp;
end
#=--------------------------------------------------=#
#================= TLDR COLD LAUNCHER ====================#
#=--------------------------------------------------=#
#1. FILES REQUIRED - FILES_ARR = [WAVELENGTHS,SPECTRA,ERRSPECTRA,DATES,CONTINUUM]
#2. IS THIS A TEST? TRUE VDM FILE? -> OPTIONAL
#3. mus: mu_l2 required. others optional?

#Tvdm file: Optional file imput of the real tdf. Used for testing puroposes.
#mu_spec, mu_l1, and mu_temp flags are options for providing individual regularization weights
#otherwise, the values are based on those for mu_l2.
#Plot_Live option shows the active reconstruction every so many iterations.
#Plot_Final option shows the final plot that displays reconstruction images X,Z,V,T.
function COLD_LAUNCH(FILES_ARR;mu_l2=40.0,mu_spec=false,mu_temp=false,mu_l1=false,mu_p=false,scale=1.0,nits=50,Tvdm="",Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true,rhoZ=8000.0,rhoN=800.0,rhoP=800.0, rhoV=800.0,rhoT=800.0)
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
	#NOTE: mu_pos WOULD BE MEANINGLESS

  if mu_temp != false
    Pars.mu_temp = mu_temp
  	Tt=true
  end
  if mu_spec != false
    Pars.mu_spec = mu_spec
		Vt = true
  end
  if mu_l1 != false
    Pars.mu_l1 = mu_l1
		Nt = true
  end
	if mu_p == true
		Pt = true
  end
  tmp,P = TLDR(1.0,DATA,Mats,Pars;Plot_A=Plot_Live,Plot_F=Plot_Final,vdmact=Tvdm,RepIt=RepIt,RepF=RepF,rhoZ=rhoZ,rhoN=rhoN,rhoP=rhoP,rhoT=rhoT,rhoV=rhoV,T_toggle=Tt,V_toggle=Vt,N_toggle=Nt,P_toggle=Pt)
end


#=--------------------------------------------------=#
#================ ADMM Algorithm ====================#
#=--------------------------------------------------=#
function TLDR(flx_scale,DATA,Mats,Pars;Plot_F=true,Plot_A=false,vdmact="",RepIt=true,RepF=true,rhoZ=8000.0,rhoN=8000.0,rhoP=8000.0,rhoT=8000.0,rhoV=8000.0,savefits=false,P_toggle=false,T_toggle=false,V_toggle=false,N_toggle=false)
	Pars.tau=2.0
	threshold = 1.0e-4 #CONVERGANCE THRESHOLD


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
		vdm = zeros(Pars.num_tdf_times,DATA.num_lines)
	for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
			Wslice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
			Q = Mats.HT * Wslice * Mats.H +  (Pars.mu_l2)*Mats.Gammatdf #INCLUCES L2 NORM ON X
			B = Mats.HT* Wslice * DATA.L
			vdm[:,l] = Q\B[:,l]
		end
	Ini=vdm.*(vdm .>= 0.0) #sdata() pulls the underlying shared array
	#Ini = inv(Mats.H'*Mats.H+(flx_scale^2*Pars.mu_l2)^2*eye(size(Mats.H)[2]))*(Mats.H'*DATA.L) #INITIALIZATION FROM TIKHONOV SOLUTION
	init_vdm =Ini #FILTER OUT NEGATIVE VALUES
	if Plot_A == true
		imshow(init_vdm,aspect="auto",origin="lower",interpolation="None",cmap="Blues")
		#colorbar()
		#show()
	end
	writecsv("tiksol.csv",init_vdm)
	#init_vdm=randn(size(init_vdm)) #Start from Random
	#init_vdm=0.0*randn(size(init_vdm)) #Start from Random
	#GENERATE REG TERMS
	X = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)	#REQUIRED
	X.vdm = copy(init_vdm)

	Z = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi) 	#REQUIRED
	Z.vdm = copy(init_vdm)
	Z.rho= rhoZ
	Wslice=0
	min_x_eqn_string_nt="Mats.HT * Wslice * DATA.L + Z.U + Z.rho .* Z.vdm"
	min_x_eqn_string_it="Mats.HT * Wslice * Mats.H + Pars.mu_l2 * Mats.Gammatdf"

	min_z_eqn_string_nt="-Z.U+Z.rho*X.vdm"
	min_z_eqn_string_it="Z.rho*Mats.Gammaspe"

	if P_toggle==true #CHECK IF USING POSITIVITY
		P = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)	#GENERATE VDM ARRAYS
		P.vdm = pos_prox_op(init_vdm)																					#POPULATE VDM
		P.rho=copy(rhoP)																											#POPULATE REG. HYPERPARAMETER
		min_x_eqn_string_nt=string(min_x_eqn_string_nt,"+","P.U+P.rho.*P.vdm")#ADD TERMS TO EQUATION NORMAL TERM
		min_x_eqn_string_it=string(min_x_eqn_string_it,"+","P.rho.*Mats.Gammatdf")#ADD TERMS TO EQUATION INVERSE TERM
	else
		P=false
	end
	if T_toggle==true #CHECK IF USING TEMPORAL GTV
		T = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)	#GENERATE VDM ARRAYS
		T.vdm = Mats.Ds*init_vdm																							#POPULATE VDM
		T.rho=copy(rhoT)																											#POPULATE REG. HYPERPARAMETER
		min_x_eqn_string_nt=string(min_x_eqn_string_nt,"+","Mats.DsT*(T.U+T.rho.*T.vdm)")#ADD TERMS TO EQUATION NORMAL TERM
		min_x_eqn_string_it=string(min_x_eqn_string_it,"+","T.rho.*Mats.DsT*Mats.Ds")#ADD TERMS TO EQUATION INVERSE TERM
	else
		T=false
	end
	if V_toggle==true	#CHECK IF USING FREQUENCY GTV
		V = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)	#GENERATE VDM ARRAYS
		V.vdm = init_vdm*Mats.Dv																							#POPULATE VDM
		V.rho=copy(rhoV)																											#POPULATE REG. HYPERPARAMETER
		min_z_eqn_string_nt=string(min_z_eqn_string_nt,"+","(V.U+V.rho.*V.vdm)*Mats.DvT")#ADD TERMS TO EQUATION NORMAL TERM
		min_z_eqn_string_it=string(min_z_eqn_string_it,"+","V.rho*Mats.Dv*Mats.DvT")#ADD TERMS TO EQUATION INVERSE TERM
	else
		V=false
	end
	if N_toggle==true	#CHECK IF USING L1-NORM
		N = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)	#GENERATE VDM ARRAYS
		N.vdm = copy(init_vdm)																								#POPULATE VDM
		N.rho=copy(rhoN)																											#POPULATE REG. HYPERPARAMETER
		min_x_eqn_string_nt=string(min_x_eqn_string_nt,"+","N.U+N.rho.*N.vdm")#ADD TERMS TO EQUATION NORMAL TERM
		min_x_eqn_string_it=string(min_x_eqn_string_it,"+","N.rho*Mats.Mats.Gammatdf")#ADD TERMS TO EQUATION INVERSE TERM
	else
		N=false
	end

	#GENERATE MINIMIZATION FUNCTIONS FROM EQUATION STRINGS     what could go wrong?
	@generated function minx_nt(X,T,P,N,Z,Pars,DATA,Mats,Wslice)
			return parse(min_x_eqn_string_nt)
	end

	@generated function minx_it(X,T,P,N,Z,Pars,DATA,Mats,Wslice)
			return parse(min_x_eqn_string_it)
	end

	@generated function minz_nt(X,T,P,N,Z,Pars,DATA,Mats)
			return parse(min_z_eqn_string_nt)
	end

	@generated function minz_it(X,T,P,N,Z,Pars,DATA,Mats)
			return parse(min_z_eqn_string_it)
	end
	

	siglvl=abs(median(DATA.L))

#Diagnostics on VDM initialized with Tihonov Solution.
	init_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
	Pars.chi2=copy(init_chi2)
	if RepIt==true
		Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=true,s=false,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,Msg=" -Tik_Init-")
	end
	#=  Calculate Initial Multipliers & RHO =#
	for p in 1:DATA.num_lines
	  Z.U[:,p]=Mats.HT * squeeze(Mats.W[p,:,:],1) * ( Mats.H * vec(Z.vdm[:,p]) - vec(DATA.L[:,p]))  #MAY NOT BE CORRECT FOR ALL VDMs.
	end
	#println("INITIAL MULTIPLIERS: ",mean(Z.U))
	Z.U=zeros(size(Z.U))
	if V_toggle==true; V.U=copy(Z.U); end
	if T_toggle==true; T.U=copy(Z.U); end
	if P_toggle==true; P.U=copy(Z.U); end
	if N_toggle==true; N.U=copy(Z.U); end

	#Qinv = zeros(Pars.num_tdf_times,Pars.num_tdf_times)  #I DON'T THINIK THIS IS USED OUTSIDE OF FUNCTIONS
	#B = zeros(Pars.num_tdf_times,DATA.num_lines)					#I DON'T THINIK THIS IS USED OUTSIDE OF FUNCTIONS
	converged = false
	CX=false
	if N_toggle==true; CN=false; end
	if P_toggle==true; CP=false; end
	if T_toggle==true; CT=false; end
	if V_toggle==true; CV=false; end

	while Pars.it <= Pars.nits && Pars.conflag==false        #ADMM ITERATION LOOP
		println(Pars.it)
		X.vdm_previous = copy(X.vdm)

	#Step 1: MINIMIZATION W.R.T. X

		s = size(X.vdm) #THIS LINE IS A WASTE OF TIME. ITS ALSO A RHYME.
		vdm = zeros(s[1],s[2])#Array(Float64,s[1],s[2])
		for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
				Wslice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
				Q = inv(minx_it(X,T,P,N,Z,Pars,DATA,Mats,Wslice)) #INVERSE TERM
				#Q = Mats.HT * Wslice * Mats.H + T.rho*Mats.DsT*Mats.Ds + (Pars.mu_l2+Z.rho+T.rho+N.rho)*Mats.Gammatdf #INVERSE TERM
				B = minx_nt(X,T,P,N,Z,Pars,DATA,Mats,Wslice) #NORMAL TERM
				#B = Mats.HT* Wslice * DATA.L + Mats.DsT*(T.U+T.rho.*T.vdm)+P.U+P.rho.*P.vdm+Z.U+Z.rho.*Z.vdm+N.U+N.rho*N.vdm #NORMAL TERM

				G=Q*B

				#vdm[:,l] = Q\B[:,l]
				vdm[:,l] = G[:,l]

		end
		X.vdm=copy(vdm)

		#X.vdm = X.vdm.*(X.vdm.>0.0)
	#Step 2: UPDATE REGULARIZATION TERMS ONLY IF THEY EXIST
		if P_toggle == true
			P.vdm_squiggle = X.vdm - P.U./P.rho
			P.vdm_previous = copy(P.vdm)
			P.vdm = pos_prox_op(P.vdm_squiggle)
		end

		if T_toggle == true
			T.vdm_squiggle = Mats.Ds*X.vdm-T.U./T.rho
			T.vdm_previous = copy(T.vdm)
			T.vdm = ell1_prox_op(T.vdm,T.vdm_squiggle,Pars.mu_temp,T.rho)
			#T.vdm=ell2_prox_op(T.vdm_squiggle,Pars.mu_temp,T.rho)
		end

		if V_toggle==true
			V.vdm_squiggle = Z.vdm*Mats.Dv - V.U./V.rho
			V.vdm_previous = copy(V.vdm)
			V.vdm = ell1_prox_op(V.vdm,V.vdm_squiggle,Pars.mu_spec,V.rho)
			#V.vdm=ell2_prox_op(V.vdm_squiggle,Pars.mu_spec,V.rho)
		end

		if N_toggle==true
			N.vdm_squiggle = X.vdm - N.U./N.rho
			N.vdm_previous = copy(N.vdm)
			N.vdm = ell1_prox_op(N.vdm,N.vdm_squiggle,Pars.mu_l1,N.rho)
		end

		Z.vdm_previous = copy(Z.vdm)
		Z.vdm = minz_nt(X,T,P,N,Z,Pars,DATA,Mats)*inv(minz_it(X,T,P,N,Z,Pars,DATA,Mats))

		#Z.vdm = Z.vdm.*(Z.vdm.>0.0)
		#Step 4: UPDATE LAGRANGE MULTIPLIERS ONLY IF THEY EXIST
		if T_toggle==true; T.U = LG_update(T.U,T.vdm,Mats.Ds*X.vdm,T.rho,Pars.alpha); end
		if V_toggle==true; V.U = LG_update(V.U,V.vdm,Z.vdm*Mats.Dv,V.rho,Pars.alpha); end
		if N_toggle==true; N.U = LG_update(N.U,N.vdm,X.vdm,N.rho,Pars.alpha); end
		if P_toggle==true; P.U = LG_update(P.U,P.vdm,X.vdm,P.rho,Pars.alpha); end
		Z.U = LG_update(Z.U,Z.vdm,X.vdm,Z.rho,Pars.alpha)

  	if Pars.it == 2
			Pars.G = 1.5
		end
	#Step 5: Update Penalty Parameters
		#Z = Pen_update(X,Z,Pars)
		#T = Pen_update_E(X,T,Pars; Ds=Mats.Ds)
		#V = Pen_update_E(Z,V,Pars; Dv=Mats.Dv)
		#P = Pen_update(X,P,Pars)


	#Step 6: Check for Convergence
		if Pars.it != 2
			if abs(regX(X,Pars)-PregX) < threshold; CX=true; end
			if N_toggle == true && abs(regN(N,Pars)-PregN) < threshold; CN=true; end
			if T_toggle == true && abs(regT(T,Pars)-PregT) < threshold; CT=true; end
			if V_toggle == true && abs(regV(V,Pars)-PregV) < threshold;	CV=true; end
			if P_toggle == true && abs(regP(P,Pars)-PregP) < threshold; CP=true; end
			if N_toggle == true && sum(N.vdm) == 0		#CHECK IF L1(N) HAS FAILED
				Pars.l1N_state=false
				Pars.conflag = true
			end
			if T_toggle == true && sum(T.vdm) == 0 		#CHECK IF L1(T) HAS FAILED
				Pars.l1T_state=false
				Pars.conflag = true
			end
			if V_toggle == true && sum(V.vdm) == 0		#CHECK IF L1(V) HAS FAILED
				Pars.l1V_state=false
				Pars.conflag = true
			end
		end
		if CX == true && CN == true && CT == true && CV == true
			Pars.conflag = true
			print_with_color(:blue,"TLDR CONVERGED \n")
		end
		#STILL NEED CHECKS.

		PregX=regX(X,Pars)
		if N_toggle==true; PregN=regN(N,Pars); end
		if T_toggle==true; PregT=regT(T,Pars); end
		if V_toggle==true; PregV=regV(V,Pars); end
		if P_toggle==true; PregP=regP(P,Pars); end

		Pars.it = Pars.it+1

		true_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
		chiprev = copy(Pars.chi2)
		Pars.chi2= copy(true_chi2)

		#Reporting
		if RepIt==true
			Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=true,L2x=true,L1T=true,L1V=true,L1N=true,s=false,Chi2x=true,Ppen=false,Zpen=true,Tpen=true,Vpen=true,Npen=true)
		end

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
function min_wrt_x(X,T,P,N,Z,Pars,DATA,Mats,eqns)
	s = size(X.vdm)
	vdm = Array(Float64,s[1],s[2])
	for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
			Wslice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
			println(eqns[2])
			Q = inv(eval(eqns[2])) #INVERSE TERM
			#Q = Mats.HT * Wslice * Mats.H + T.rho*Mats.DsT*Mats.Ds + (Pars.mu_l2+Z.rho+T.rho+N.rho)*Mats.Gammatdf #INVERSE TERM
			B = eval(eqns[1]) #NORMAL TERM
			#B = Mats.HT* Wslice * DATA.L + Mats.DsT*(T.U+T.rho.*T.vdm)+P.U+P.rho.*P.vdm+Z.U+Z.rho.*Z.vdm+N.U+N.rho*N.vdm #NORMAL TERM

			G=Q*B

			#vdm[:,l] = Q\B[:,l]
			vdm[:,l] = G[:,l]

	end
	X.vdm=copy(vdm) #sdata() pulls the underlying shared array
	X.vdm
end


#=--------------------------------------------------=#
#=============== Update Multipliers =================#
#=--------------------------------------------------=#
#X and Y are IMAGE structs
#alpha is a throttling term on multiplier update.
function LG_update(U,Z,X,rho,alpha)
	U = U + (rho/alpha).*(Z - X)
end

#=--------------------------------------------------=#
#================= Update Penalty ===================#
#=--------------------------------------------------=#
function Pen_update(X,Y,P)

	S = Y.rho*(Y.vdm-Y.vdm_previous)
	R = X.vdm-Y.vdm
	tau_prim_previous = copy(Y.tau_prim)
  Y.tau_prim = sqrt(P.num_tdf_times)*P.eps_abs + P.eps_rel*(maximum([ell2norm(X.vdm),ell2norm(Y.vdm)]))
	tau_dual_previous = copy(Y.tau_dual)
  Y.tau_dual = sqrt(P.num_tdf_times)*P.eps_abs + P.eps_rel*ell2norm(Y.U)
	eta = ell2norm(R)*Y.tau_dual / (ell2norm(S)*Y.tau_prim)
#	eta = ell2norm(R)*tau_dual_previous[l] / (ell2norm(S)*tau_prim_previous[l])	#OPTION 2
	phi_previous = copy(Y.phi)
	Y.phi = max(ell2norm(R/Y.tau_prim),ell2norm(S)/Y.tau_dual)
	if (1.0/P.tau[1] <= eta[1]) && (eta[1] <= P.tau[1]) || (Y.phi[1] < P.sigma[1]*phi_previous[1])
		#NOTHING HAPPENS IN HERE
	elseif eta < (1.0/P.tau)
		Y.rho_max = copy(Y.rho)
		if Y.rho_min > 0.0
			Y.rho = sqrt(Y.rho_min*Y.rho_max)
		else Y.rho = Y.rho_max/P.G
		end
	elseif eta > P.tau
		Y.rho_min = copy(Y.rho)
		if Y.rho_max < Inf
			Y.rho = sqrt(Y.rho_min*Y.rho_max)
		else Y.rho = Y.rho_min*P.G
		end
	end
	Y
end

#=--------------------------------------------------=#
#================= Update Penalty ===================#
#=--------------------------------------------------=#
function Pen_update_E(X,Y,P;Dv=false,Ds=false)
	if Dv==false
		Dv = eye(size(X.vdm)[2])
	end
	if Ds == false
		Ds = eye(size(X.vdm)[1])
	end
	S = Y.rho*(Y.vdm-Y.vdm_previous)
	R = Ds*X.vdm*Dv-Y.vdm
	tau_prim_previous = copy(Y.tau_prim)
  Y.tau_prim = sqrt(P.num_tdf_times)*P.eps_abs + P.eps_rel*(maximum([ell2norm(X.vdm),ell2norm(Y.vdm)]))
	tau_dual_previous = copy(Y.tau_dual)
  Y.tau_dual = sqrt(P.num_tdf_times)*P.eps_abs + P.eps_rel*ell2norm(Y.U)
	eta = ell2norm(R)*Y.tau_dual / (ell2norm(S)*Y.tau_prim)
#	eta = ell2norm(R)*tau_dual_previous[l] / (ell2norm(S)*tau_prim_previous[l])	#OPTION 2
	phi_previous = copy(Y.phi)
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
		Y.rho_min = copy(Y.rho)
		if Y.rho_max < Inf
			Y.rho = sqrt(Y.rho_min*Y.rho_max)
		else Y.rho = Y.rho_min*P.G
		end
	end
	Y
end
