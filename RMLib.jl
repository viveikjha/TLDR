
module RMLib
using JLD
using FITSIO
using PyPlot
using RMLibMore
using RMTypes
using DataImportNEW
using GenMatrices

#include("RMLibMore.jl")
#include("RMTypes.jl")
#include("DataImport.jl")
#include("DataImportNEW.jl")
#include("GenMatrices.jl")
export HOT_LAUNCH,COLD_LAUNCH
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
#RecD option stores convergence data.
function HOT_LAUNCH(Data,Mats,Pars,Fit;scale=1.0,nits=50,Tvdm="",Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true,RecD=false)

	#Data.L=scale*(Data.L)
	#Data.EL=scale*(Data.EL)
	#Data.continuum_flux=scale*Data.continuum_flux
	#Data.continuum_error_flux=scale*Data.continuum_error_flux
	Pars.nits=nits

	#SET RECONSTRUCTION PARAMETERS
	scale=1.0
  Fit.msmo=Fit.msmo

  tmp,P = TLDR(1.0,Data,Mats,Pars,Fit;Plot_A=Plot_Live,Plot_F=Plot_Final,vdmact=Tvdm,RepIt=RepIt,RepF=RepF,RecD=RecD)
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
function COLD_LAUNCH(FILES_ARR,Fit;scale=1.0,nits=50,Tvdm="",Plot_Live=true,Plot_Final=true,RepIt=true,RepF=true)
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
    tmp,P = TLDR(1.0,DATA,Mats,Pars,Fit;Plot_A=Plot_Live,Plot_F=Plot_Final,vdmact=Tvdm)
end


#=--------------------------------------------------=#
#================ ADMM Algorithm ====================#
#=--------------------------------------------------=#
function TLDR(flx_scale,DATA,Mats,Pars,Fit;Plot_F=true,Plot_A=false,vdmact="",RepIt=true,RepF=true,savefits=false,RecD=false)
	if RecD==true
		Record=Init_Record_Data(Pars.nits)
	end
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


	#sol=gen_tiksol(Pars,Mats,DATA;scale=1.0,mu_smoo=Fit.TI,plotting=false,save=false)
	#Ini=sol.*(sol .>= 0.0) #sdata() pulls the underlying shared array
	#Ini = inv(Mats.H'*Mats.H+(flx_scale^2*Fit.msmo)^2*eye(size(Mats.H)[2]))*(Mats.H'*DATA.L) #INITIALIZATION FROM TIKHONOV SOLUTION
	#init_vdm =Ini #FILTER OUT NEGATIVE VALUES
	init_vdm=zeros(Pars.num_tdf_times,DATA.num_lines)
	Ini=0 #Memory Release
	if Plot_A == true
		imshow(init_vdm,aspect="auto",origin="lower",interpolation="None")
		#colorbar()
		#show()
	end
	writecsv(string(Pars.directory,"tiksol_init.csv"),init_vdm)
	#init_vdm=randn(size(init_vdm)) #Start from Random
	#init_vdm=0.0*randn(size(init_vdm)) #Start from Random
	tmax=5.0
	 X = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 X.tau_max=tmax
	 Z = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 Z.tau_max=tmax
	 T = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 T.tau_max=tmax
	 V = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 V.tau_max=tmax
	 P = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 P.tau_max=tmax
	 N = Gen_Var(Pars.rho0,Pars.num_tdf_times,DATA.num_lines,initial_psi)
	 N.tau_max=tmax


	#When loading a known TDF
	if vdmact != ""
		vdm_path = vdmact
		vdm_act = readcsv(vdm_path)
		X.vdm = copy(vdm_act)
		Z.vdm = copy(vdm_act)
		T.vdm = Mats.Ds*vdm_act
		V.vdm = vdm_act*Mats.Dv
		P.vdm = pos_prox_op(vdm_act)
		N.vdm = copy(vdm_act)
		act_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
		#println("actual chi2: ",act_chi2)
		Pars.chi2=copy(act_chi2)
		if RepIt==true
			Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=true,s=false,Chi2x=true,Ppen=false,Zpen=true,Tpen=true,Vpen=true,Npen=true,Msg=" -True_VDM-")
		end
	end

	 X.vdm = copy(init_vdm)
	#zinis=gen_tiksol(Pars,Mats,DATA;scale=1.0,mu_smoo=1000.0,plotting=false,save=false)
	Z.vdm = copy(init_vdm)
	T.vdm = Mats.Ds*init_vdm
	V.vdm = init_vdm*Mats.Dv
	P.vdm = pos_prox_op(init_vdm)
	N.vdm = copy(init_vdm)
	init_vdm=0 #Memory Release
	#Initiailize Penalty Parameters.
	Z.rho= Fit.pz
	P.rho=copy(Fit.pp)
	T.rho=copy(Fit.pt)
	V.rho=copy(Fit.pv)
	N.rho=copy(Fit.pn)

	siglvl=abs(median(DATA.L))

#Diagnostics on VDM initialized with Tihonov Solution.
	#println("-------")
	init_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
	Pars.chi2=copy(init_chi2)
	init_chi2=0 #Memory Release
	if RepIt==true
			Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=true,s=false,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,Msg=" -Tik_Init-")
	end
	#=  Calculate Initial Multipliers & RHO =#
	Z.U=Z.U+1
	#slice_size=size(Mats.W)
	#for p in 1:DATA.num_lines
	#	W_slice = reshape(Mats.W[p,:,:],slice_size[2],slice_size[3])
	  #Z.U[:,p]=Mats.HT * squeeze(Mats.W[p,:,:],1) * ( Mats.H * vec(Z.vdm[:,p]) - vec(DATA.L[:,p]))
	#	Z.U[:,p]=Mats.HT * W_slice * ( Mats.H * vec(Z.vdm[:,p]) - vec(DATA.L[:,p]))
	#end

	#println("INITIAL MULTIPLIERS: ",mean(Z.U))
	#Z.U=ones(size(Z.U))

	V.U = copy(Z.U)
	T.U = copy(Z.U)
	P.U = copy(Z.U)
	N.U = copy(Z.U)

	Qinv = zeros(Pars.num_tdf_times,Pars.num_tdf_times)
	B = zeros(Pars.num_tdf_times,DATA.num_lines)
	converged = false
	#Pars.conflag=true #should never begin main loop.
	gc() #GARBAGE CLEANUP
	fcn_vals=zeros(Pars.nits+1)
	while Pars.it <= Pars.nits && Pars.conflag==0        #ADMM ITERATION LOOP

		X.vdm_previous = copy(X.vdm)
	#Step 1: MINIMIZATION W.R.T. X
		X.vdm = min_wrt_x(X,T,P,N,Z,Pars,DATA,Mats,Fit)
		#X.vdm = X.vdm.*(X.vdm.>0.0)
	#Step 2: UPDATE REGULARIZATION TERMS
		P.vdm_squiggle = X.vdm - P.U./P.rho
		P.vdm_previous = copy(P.vdm)
		P.vdm = pos_prox_op(P.vdm_squiggle)

		T.vdm_squiggle = Mats.Ds*X.vdm-T.U./T.rho
		T.vdm_previous = copy(T.vdm)
		T.vdm = ell1_prox_op(T.vdm,T.vdm_squiggle,Fit.mtem,T.rho)
		#T.vdm=ell2_prox_op(T.vdm_squiggle,Fit.mtem,T.rho)

		V.vdm_squiggle = Z.vdm*Mats.Dv - V.U./V.rho
		V.vdm_previous = copy(V.vdm)
		V.vdm = ell1_prox_op(V.vdm,V.vdm_squiggle,Fit.mspe,V.rho)
		#V.vdm=ell2_prox_op(V.vdm_squiggle,Fit.mspe,V.rho)

		N.vdm_squiggle = X.vdm - N.U./N.rho
		N.vdm_previous = copy(N.vdm)
		N.vdm = ell1_prox_op(N.vdm,N.vdm_squiggle,Fit.ml1,N.rho)

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
		if Pars.it < 5
			update_penN(X,N,Pars,Fit)
			update_penP(X,P,Pars,Fit)
			update_penT(X,T,Pars,Mats,Fit)
			update_penV(Z,V,Pars,Mats,Fit)
			update_penZ(X,Z,Pars,Fit)
		end #endif
		#println("Z: ",Z.rho," ",Z.tau," T: ",T.rho," ",T.tau," V: ",V.rho," ",V.tau," P: ",P.rho," ",N.tau)

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
			if sum(N.vdm) == 0		#CHECK IF L1(N) HAS FAILED
				Pars.l1N_state=false
				#Pars.conflag = true
			end
			if sum(T.vdm) == 0 		#CHECK IF L1(T) HAS FAILED
				Pars.l1T_state=false
				#Pars.conflag = true
			end
			if sum(V.vdm) == 0		#CHECK IF L1(V) HAS FAILED
				Pars.l1V_state=false
				#Pars.conflag = true
			end
		end
		if CX == true && CN == true && CT == true && CV == true
			Pars.conflag = 1
			#print_with_color(:blue,"TLDR CONVERGED \n")
		end
		PregX=regX(X,Pars)
		PregN=regN(N,Pars)
		PregT=regT(T,Pars)
		PregV=regV(V,Pars)
		Pars.it = Pars.it+1
		#println(DATA.num_spectra_samples, " ", DATA.num_lines, size(DATA.L))
		true_chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
		chiprev = copy(Pars.chi2)
		Pars.chi2= copy(true_chi2)
		fcn_vals[Pars.it]=J(X,P,T,V,N,DATA,Mats,Pars)
		if J(X,P,T,V,N,DATA,Mats,Pars) > 1.0e20
			Pars.conflag = 2
			#print_with_color(:red,"!!! FUNCTION DIVERGING ABORTING !!!\n")
			RepF=false
			savefits=false
			plot_A=false
		end
		#Reporting
		if RepIt==true
			#Report(X,Z,P,T,V,N,DATA,Mats,Pars,Jf=true,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,s=true,Pres=true,Zres=true,Tres=true,Vres=true,Nres=true)
			Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=true,L2x=true,L1T=true,L1V=true,L1N=true,s=false,Chi2x=true,Ppen=false,Zpen=false,Tpen=false,Vpen=false,Npen=false)
		end
		#Recording
		if RecD==true
			Record.ConFlag[Pars.it-1]=Pars.conflag
			Record.Chi2[Pars.it-1]=Pars.chi2
			Record.J[Pars.it-1]=J(X,P,T,V,N,DATA,Mats,Pars)
			Record.Z_res[Pars.it-1]=ell2norm(abs(X.vdm-Z.vdm))
			Record.T_res[Pars.it-1]=ell2norm(abs(Mats.Ds*X.vdm-T.vdm))
			Record.N_res[Pars.it-1]=ell2norm(abs(X.vdm-N.vdm))
			Record.V_res[Pars.it-1]=ell2norm(abs(X.vdm*Mats.Dv-V.vdm))
			Record.P_res[Pars.it-1]=ell2norm(abs(X.vdm-P.vdm))
		end
		#Plotting
		if Plot_A == true && (Pars.it%10 == 0)
			clf()
			imshow((X.vdm),extent=[minimum(DATA.wavelength),maximum(DATA.wavelength),0.0,50.0],aspect="auto",origin="lower",interpolation="None")
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
	#if Pars.it==Pars.nits-1 && Pars.conflag==0 || Pars.conflag==2
	#	print_with_color(:red,"TLDR FAILED TO CONVERGE \n")
	#end
	if savefits==true
		Write_FITS(X,P);
	end
	if RecD ==true
		#save_data("TLDR_Con_data.jld",Record)
		writecsv(string(Pars.directory,"Reconstruction/Conflag.csv"),Record.ConFlag)
		writecsv(string(Pars.directory,"Reconstruction/Chi2.csv"),Record.Chi2)
		writecsv(string(Pars.directory,"Reconstruction/J.csv"),Record.J)
		writecsv(string(Pars.directory,"Reconstruction/Z_res.csv"),Record.Z_res)
		writecsv(string(Pars.directory,"Reconstruction/T_res.csv"),Record.T_res)
		writecsv(string(Pars.directory,"Reconstruction/V_res.csv"),Record.V_res)
		writecsv(string(Pars.directory,"Reconstruction/N_res.csv"),Record.N_res)
		writecsv(string(Pars.directory,"Reconstruction/P_res.csv"),Record.P_res)
	end
	#writecsv("J_vals.csv", fcn_vals)
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
function min_wrt_x(X,T,P,N,Z,Pars,DATA,Mats,Fit)
	s = size(X.vdm)
	vdm = Array(Float64,s[1],s[2])
	slice_size=size(Mats.W)
	for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
			#W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
			W_slice = reshape(Mats.W[l,:,:],slice_size[2],slice_size[3])
			Q = Mats.HT * W_slice * Mats.H + T.rho*Mats.DsT*Mats.Ds + (Fit.msmo+Z.rho+T.rho+N.rho)*Mats.Gammatdf #INCLUCES L1 NORM ON X
			B = Mats.HT* W_slice * DATA.L + Mats.DsT*(T.U+T.rho.*T.vdm)+P.U+P.rho.*P.vdm+Z.U+Z.rho.*Z.vdm+N.U+N.rho*N.vdm #INCLUDES L1 NORM ON X

			G=inv(Q)*B
			#vdm[:,l] = Q\B[:,l]
			vdm[:,l] = G[:,l]

	end
	X.vdm=copy(vdm) #sdata() pulls the underlying shared array
	X.vdm
end

#=--------------------------------------------------=#
#=============== Minimization wrt Z =================#
#=--------------------------------------------------=#
function min_wrt_z(X,V,Z,Pars,DATA,Mats)
  	Rinv = inv(V.rho.*Mats.Dv*Mats.DvT+Z.rho.*Mats.Gammaspe)
    C = (V.U+V.rho.*V.vdm)*Mats.DvT-Z.U+(Z.rho.*X.vdm)	#ORIGINAL PAPER VERSION
		Z.vdm = 	C*Rinv
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
#================= Update PenN ======================#
#=--------------------------------------------------=#
function update_penN(X,N,Pars,Fit)
	rc=ell2norm(X.vdm-N.vdm)
	sc=ell2norm(N.rho*(N.vdm-N.vdm_previous))
	if rc > Fit.ml1*sc
		N.rho=N.rho*Pars.tau
	elseif  sc > Fit.ml1*rc
		N.rho=N.rho/Pars.tau
	else
	end #endif
end

#=--------------------------------------------------=#
#================= Update PenP ======================#
#=--------------------------------------------------=#
function update_penP(X,P,Pars,Fit)
rc= ell2norm(X.vdm-P.vdm)
sc=ell2norm(P.rho*(P.vdm-P.vdm_previous))
	if rc > sc
		P.rho=P.rho*Pars.tau
	elseif  sc > rc
		P.rho=P.rho/Pars.tau
	else
	end #endif
end

#=--------------------------------------------------=#
#================= Update PenT ======================#
#=--------------------------------------------------=#
function update_penT(X,T,Pars,Mats,Fit)
rc=ell2norm(T.vdm-Mats.Ds*X.vdm)
sc=ell2norm(T.rho*(T.vdm-T.vdm_previous))
if rc > sc
	T.rho=T.rho*Pars.tau
elseif  sc > rc
	T.rho=T.rho/Pars.tau
else
end #endif

end

#=--------------------------------------------------=#
#================= Update PenV ======================#
#=--------------------------------------------------=#
function update_penV(Z,V,Pars,Mats,Fit)
rc=ell2norm(V.vdm-(Z.vdm*Mats.Dv))
sc=ell2norm(V.rho*(V.vdm-V.vdm_previous))
if rc > sc
	V.rho=V.rho*Pars.tau
elseif  sc > rc
	V.rho=V.rho/Pars.tau
else
end #endif

end
#=--------------------------------------------------=#
#================= Update PenZ ======================#
#=--------------------------------------------------=#
function update_penZ(X,Z,Pars,Fit)
rc=ell2norm(X.vdm-Z.vdm)
sc=ell2norm(Z.rho*(Z.vdm-Z.vdm_previous))
if rc > sc
	Z.rho=Z.rho*Pars.tau
elseif  sc > rc
	Z.rho=Z.rho/Pars.tau
else
end #endif

end
end #endmodule
