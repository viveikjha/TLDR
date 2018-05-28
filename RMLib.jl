
module RMLib
using JLD
using FITSIO
using PyPlot
using RMLibMore
using RMTypes
using DataImportNEW
using GenMatrices
using Wavelets
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


	Pars.nits=nits

	#SET RECONSTRUCTION PARAMETERS

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


	#data_report(DATA)
	Pars = init_Params()
	Pars.num_tdf_times=10 #SHOULD NOT BE HARDcoded

	 Mats=Gen_Mats(DATA,Pars)

	Pars.nits=nits

	#SET RECONSTRUCTION PARAMETERS

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

	dims=(Pars.num_tdf_times,DATA.num_lines)
	Pars.wavelet_bases=[WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db7,WT.db8,WT.haar];
	Pars.DvM=Array{Float64}(dims[1],dims[2],length(Pars.wavelet_bases))
	Pars.DsM=Array{Float64}(dims[1],dims[2],length(Pars.wavelet_bases))

	#IWx=Array{Floa64}(nx*ny,length(wavelet_bases)) #nx and ny will be num_tdf_times and num???
	Pars.IDvM=Array{Float64}(dims[1],dims[2],length(Pars.wavelet_bases))
	Pars.IDsM=Array{Float64}(dims[1],dims[2],length(Pars.wavelet_bases))


	X= Gen_Var(Pars.rho0,dims,initial_psi)
	Z = Gen_Var(Pars.rho0,dims,initial_psi)
	P = Gen_Var(Pars.rho0,dims,initial_psi)
	N = Gen_Var(Pars.rho0,dims,initial_psi)

	#THESE WILL NEED TO CHANGE FOR WAVELETS. 3D OR JUST EXTENDED...
	wdims=(Pars.num_tdf_times,DATA.num_lines,length(Pars.wavelet_bases)) #WAVELETS
	T = Gen_Var(Pars.rho0,wdims,initial_psi)
	V = Gen_Var(Pars.rho0,wdims,initial_psi)

	#if isnan(Fit.TI) == false && isnan(Fit.TI2) == false
		#if Fit.TI != Fit.TI2
			#sol1=gen_tiksol(Pars,Mats,DATA;scale=1.0,mu_smoo=Fit.TI,plotting=false,save=false)
			#writecsv(string(Pars.directory,"tiksol_init.csv"),sol1)
			#X.vdm=copy(sol1)
			#sol1=0
			#sol2=gen_tiksol(Pars,Mats,DATA;scale=1.0,mu_smoo=Fit.TI2,plotting=false,save=false)
			#writecsv(string(Pars.directory,"tiksol_init2.csv"),sol2)
			#Z.vdm=copy(sol2)
			#sol2=0
		#else
			#sol1=gen_tiksol(Pars,Mats,DATA;scale=1.0,mu_smoo=Fit.TI,plotting=false,save=false)
			#writecsv(string(Pars.directory,"tiksol_init.csv"),sol1)
			#X.vdm=copy(sol1)
			#Z.vdm=copy(sol1)
			#sol1=0
		#end
	#elseif isnan(Fit.TI) == false && isnan(Fit.TI2) == true
		#sol1=gen_tiksol(Pars,Mats,DATA;scale=1.0,mu_smoo=Fit.TI,plotting=false,save=false)
		#writecsv(string(Pars.directory,"tiksol_init.csv"),sol1)
		#X.vdm=copy(sol1)
		#Z.vdm=copy(sol1)
		#sol1=0
	#else
		#println("You have made an error with Tikhonov initialization. If only one Tikhonov initialization is to be used, Fit.TI should be used, not Fit.TI2")
		#X.vdm=zeros(Pars.num_tdf_times,DATA.num_lines)
		#Z.vdm=zeros(Pars.num_tdf_times,DATA.num_lines)
	#end
	#if Plot_A == true
		#imshow(init_vdm,aspect="auto",origin="lower",interpolation="None",cmap="Reds")
		#colorbar()
		#show()

	#end
	#writecsv(string(Pars.directory,"tiksol_init.csv"),init_vdm)
	#init_vdm=randn(size(init_vdm)) #Start from Random
	#init_vdm=0.0*randn(size(init_vdm)) #Start from Random
	tmax=5.0




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

	#X.vdm = copy(init_vdm)
	#zinis=gen_tiksol(Pars,Mats,DATA;scale=1.0,mu_smoo=1000.0,plotting=false,save=false)
	#Z.vdm = copy(init_vdm)
	#T.vdm = Mats.Ds*copy(X.vdm)
	#T.vdm = Ds(X.vdm) #Wavelets
	#V.vdm = copy(Z.vdm)*Mats.Dv
	#V.vdm = Dv(Z.vdm) #WAVELETS
	#P.vdm = pos_prox_op(copy(X.vdm))
	#N.vdm = copy(X.vdm)
	init_vdm=0 #Memory Release

	#Initiailize Penalty Parameters.
	Z.rho= Fit.pz
	P.rho=copy(Fit.pp)
	N.rho=copy(Fit.pn)
	T.rho=copy(Fit.pt)
	V.rho=copy(Fit.pv)

	#println(typeof(DATA.L))
	siglvl=abs(median(DATA.L))

#Diagnostics on VDM initialized with Tikhonov Solution.
	#println("-------")
	chi2 = Chi2(Model(X.vdm,Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
	Pars.chi2=copy(chi2)

	if RepIt==true
			Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=true,s=false,L2x=true,L1T=true,L1V=true,L1N=true,Chi2x=true,Msg=" -Tik_Init-")
	end
	#=  Calculate Initial Multipliers & RHO =#
	IU=1.0
	fill!(Z.U,IU)
	fill!(V.U,IU)
	fill!(T.U,IU)
	fill!(P.U,IU)
	fill!(N.U,IU)

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

		#T.vdm_squiggle = Mats.Ds*X.vdm-T.U./T.rho
		T.vdm_squiggle = Ds(X.vdm,Pars)-T.U./T.rho	#WAVELETS
		T.vdm_previous = copy(T.vdm)
		T.vdm = ell1_prox_op(T.vdm,T.vdm_squiggle,Fit.mtem,T.rho)
		#T.vdm=ell2_prox_op(T.vdm_squiggle,Fit.mtem,T.rho)

		#V.vdm_squiggle = Z.vdm*Mats.Dv - V.U./V.rho
		V.vdm_squiggle = Dv(Z.vdm,Pars) - V.U./V.rho	#WAVELETS
		V.vdm_previous = copy(V.vdm)
		V.vdm = ell1_prox_op(V.vdm,V.vdm_squiggle,Fit.mspe,V.rho)
		#V.vdm=ell2_prox_op(V.vdm_squiggle,Fit.mspe,V.rho)

		N.vdm_squiggle = X.vdm - N.U./N.rho
		N.vdm_previous = copy(N.vdm)
		N.vdm = ell1_prox_op(N.vdm,N.vdm_squiggle,Fit.ml1,N.rho)

		Z.vdm_previous = copy(Z.vdm)
		Z.vdm = min_wrt_z(X,V,Z,Pars,DATA,Mats)

		#Step 4: UPDATE LAGRANGE MULTIPLIERS
		#CHANGE
		#T.U = LG_update(T.U,T.vdm,Mats.Ds*X.vdm,T.rho,Pars.alpha)
		T.U = LG_update(T.U,T.vdm,Ds(X.vdm,Pars),T.rho,Pars.alpha)	#WAVELETS

		#V.U = LG_update(V.U,V.vdm,Z.vdm*Mats.Dv,V.rho,Pars.alpha)
		V.U = LG_update(V.U,V.vdm,Dv(Z.vdm,Pars),V.rho,Pars.alpha) #WAVELETS


		N.U = LG_update(N.U,N.vdm,X.vdm,N.rho,Pars.alpha)
		P.U = LG_update(P.U,P.vdm,X.vdm,P.rho,Pars.alpha)
		Z.U = LG_update(Z.U,Z.vdm,X.vdm,Z.rho,Pars.alpha)

		#Step 5: UPDATE INTERMEDIATE MULTIPLIERS
		T.UI_previous=copy(T.UI)
		#T.UI = IM_update(T.UI,T.vdm,Mats.Ds*X.vdm_previous,T.rho)
		T.UI = IM_update(T.UI,T.vdm,Ds(X.vdm_previous,Pars),T.rho) #Wavelets
		V.UI_previous=copy(V.UI)
		#V.UI = IM_update(V.UI,V.vdm,Z.vdm_previous*Mats.Dv,V.rho)
		V.UI = IM_update(V.UI,V.vdm,Dv(Z.vdm_previous,Pars),V.rho)	#WAVELETS


		N.UI_previous=copy(N.UI)
		N.UI = IM_update(N.UI,N.vdm,X.vdm_previous,N.rho)
		P.UI_previous=copy(P.UI)
		P.UI = IM_update(P.UI,P.vdm,X.vdm_previous,P.rho)
		Z.UI_previous=copy(Z.UI)
		Z.UI = IM_update(Z.UI,Z.vdm,X.vdm_previous,Z.rho)

		#Step 6: UPDATE RHO
		if (Pars.it%10 == 0 || Pars.it==2) #UPUDATE EVERY ## ITERATIONS
			#adapt(T,Mats.Ds*X.vdm,Mats.Ds*X.vdm_previous,Z.UI_previous,Z.UI,Pars.it)
			adapt(T,Ds(X.vdm,Pars),Ds(X.vdm_previous,Pars),Ds(Z.UI_previous,Pars),Ds(Z.UI,Pars),Pars.it) #WAVELETS #TROUBLE HERE!

			#adapt(V,Z.vdm*Mats.Dv,Z.vdm_previous*Mats.Dv,Z.UI_previous,Z.UI,Pars.it)
			adapt(V,Dv(Z.vdm,Pars),Dv(Z.vdm_previous,Pars),Dv(Z.UI_previous,Pars),Dv(Z.UI,Pars),Pars.it) #WAVELETS #TROUBLE HERE!
#
			adapt(N,X.vdm,X.vdm_previous,Z.UI_previous,Z.UI,Pars.it)
			adapt(P,X.vdm,X.vdm_previous,Z.UI_previous,Z.UI,Pars.it)
			adapt(Z,X.vdm,X.vdm_previous,Z.UI_previous,Z.UI,Pars.it)
		end

  	if Pars.it == 2

			Pars.G = 1.5
		end

		#Step 6: Check for Convergence
		if Pars.it != 2
			if abs.(regX(X,Pars)-PregX) < Pars.threshold
				CX=true
			end
			if abs.(regN(N,Pars)-PregN) < Pars.threshold
				CN=true
			end
			if abs.(regT(T,Pars)-PregT) < Pars.threshold
				CT=true
			end
			if abs.(regV(V,Pars)-PregV) < Pars.threshold
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
			print_with_color(:blue,"TLDR CONVERGED \n")
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
			print_with_color(:red,"!!! FUNCTION DIVERGING ABORTING !!!\n")
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

			#Record.Z_res[Pars.it-1]=ell2norm(abs.(X.vdm-Z.vdm))
			Record.Z_res[Pars.it-1]=vecnorm(abs.(X.vdm-Z.vdm)[:],2)
			#Record.T_res[Pars.it-1]=ell2norm(abs.(Mats.Ds*X.vdm-T.vdm))
			Record.T_res[Pars.it-1]=vecnorm(abs.(Mats.Ds*X.vdm-T.vdm)[:],2)
			#Record.N_res[Pars.it-1]=ell2norm(abs.(X.vdm-N.vdm))
			Record.N_res[Pars.it-1]=vecnorm(abs.(X.vdm-N.vdm)[:],2)
			#Record.V_res[Pars.it-1]=ell2norm(abs.(X.vdm*Mats.Dv-V.vdm))
			Record.V_res[Pars.it-1]=vecnorm(abs.(X.vdm*Mats.Dv-V.vdm)[:],2)
			#Record.P_res[Pars.it-1]=ell2norm(abs.(X.vdm-P.vdm))
			Record.P_res[Pars.it-1]=vecnorm(abs.(X.vdm-P.vdm)[:],2)
		end
		#Plotting
		if Plot_A == true && (Pars.it%50 == 0)
			clf()
			imshow((X.vdm),extent=[minimum(DATA.wavelength),maximum(DATA.wavelength),0.0,50.0],aspect="auto",origin="lower",interpolation="None",cmap="Reds")
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
		save_data("TLDR_Con_data.jld",Record)
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
 function Gen_Var(rhoi, dims,psi)
	#dims=num_tdf_times,num_lines
	x = init_IMAGE(rhoi,dims)
	x.vdm = zeros(dims)
	fill!(x.vdm,psi)
	x.vdm_squiggle = zeros(dims)
	fill!(x.vdm_squiggle,psi)
	x.U=zeros(dims)
	x.UI=zeros(dims)
	x.UI_previous=zeros(dims)
	x
end


#=--------------------------------------------------=#
#=============== Minimization wrt X =================#
#=--------------------------------------------------=#
function min_wrt_x(X,T,P,N,Z,Pars,DATA,Mats,Fit)
	s = size(X.vdm)
	#vdm = Array(Float64,s[1],s[2])
	vdm = Array{Float64}(s[1],s[2]) #UPDATE FOR JV0.6.2
	for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
			#Q = Mats.HT * Mats.W[:,:,l] * Mats.H + T.rho*Mats.DsT*Mats.Ds + (Fit.msmo+Z.rho+T.rho+N.rho)*Mats.Gammatdf #INCLUCES L1 NORM ON X
			#B = Mats.HT* Mats.W[:,:,l] * DATA.L + Mats.DsT*(T.U+T.rho.*T.vdm)+P.U+P.rho.*P.vdm+Z.U+Z.rho.*Z.vdm+N.U+N.rho*N.vdm #INCLUDES L1 NORM ON X
			DsTDs=length(Pars.wavelet_bases)*eye(s[1],s[2]) #not sure where the length( should end......
			Q = Mats.HT * Mats.W[:,:,l] * Mats.H + T.rho*DsTDs+(Fit.msmo+Z.rho+T.rho+N.rho)*Mats.Gammatdf #INCLUCES L1 NORM ON X
			B = Mats.HT* Mats.W[:,:,l] * DATA.L + DsT(T.U+T.rho.*T.vdm,Pars)+P.U+P.rho.*P.vdm+Z.U+Z.rho.*Z.vdm+N.U+N.rho*N.vdm #INCLUDES L1 NORM ON X

			G=inv(Q)*B
			vdm[:,l] = G[:,l]
	end
	X.vdm=copy(vdm) #sdata() pulls the underlying shared array
	X.vdm
end

#=--------------------------------------------------=#
#=============== Minimization wrt Z =================#
#=--------------------------------------------------=#
function min_wrt_z(X,V,Z,Pars,DATA,Mats)
  	#Rinv = inv(V.rho.*Mats.Dv*Mats.DvT+Z.rho.*Mats.Gammaspe)
    #C = (V.U+V.rho.*V.vdm)*Mats.DvT-Z.U+(Z.rho.*X.vdm)	#ORIGINAL PAPER VERSION
	dims=size(X.vdm)
	DvTDv=length(Pars.wavelet_bases)*eye(dims[1],dims[2]) #not sure where the length( should end......
	Rinv = inv(V.rho.*DvTDv+Z.rho.*Mats.Gammaspe)
    C = DvT(V.U+V.rho.*V.vdm,Pars)-Z.U+(Z.rho.*X.vdm)	#ORIGINAL PAPER VERSION
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
#========= Update Intermediate Multipliers ==========#
#=--------------------------------------------------=#
#FOR ADAPTIVE GCADMM
#X and Y are IMAGE structs
function IM_update(UI,Z_prev,X,rho)
	#UI = UI + (rho).*(X-Z_prev)
	UI = UI + rho.*(Z_prev-X) #According to the paper.
end

#=--------------------------------------------------=#
#====================== ADAPT =======================#
#=--------------------------------------------------=#
#Pass in the data structs
#Calculates the stepsize automatically according to Xu et al 2017
#Must operate on each minimization pair individually
#J will be the regularization Variable
#K will be comparision
#Changes rho in data struct J
function adapt(J,Kvdm,Kvdm_previous,KUI_previous,KUI,k)

	tau_previous=J.rho
	DJUI=J.UI.-J.UI_previous
	DJ=J.vdm.-J.vdm_previous

	DKUI=KUI-KUI_previous
	DK=Kvdm-Kvdm_previous

	#SPECTRAL STEP SIZES
	Alpha_SD=vecdot(DJUI,DJUI)/vecdot(DJ,DJUI)
	Alpha_MG=vecdot(DJ,DJUI)/vecdot(DJ,DJ)

	if 2.0*Alpha_MG > Alpha_SD
		Alpha=Alpha_MG
	else
		Alpha = Alpha_SD-(Alpha_MG/2.0)
	end #END IF ELSE

	Beta_SD=dot(DKUI,DKUI)/dot(DK,DJUI)
	Beta_MG=dot(DK,DKUI)/dot(DK,DK)
	#Beta_SD=vecdot(DKUI,DKUI)/vecdot(DK,DJUI)
	#Beta_MG=vecdot(DK,DKUI)/vecdot(DK,DK)

	if 2.0*Beta_MG > Beta_SD
		Beta=Beta_MG
	else
		Beta = Beta_SD-(Beta_MG/2.0)
	end #END IF ELSE
	#CORRELATIONS
	Alpha_corr=dot(DJ,DJUI)/(BLAS.asum(DJ)*BLAS.asum(DJUI))
	Beta_corr=dot(DK,DKUI)/(BLAS.asum(DK)*BLAS.asum(DKUI))
	#Alpha_corr=vecdot(DJ,DJUI)/(BLAS.asum(DJ)*BLAS.asum(DJUI))
	#Beta_corr=vecdot(DK,DKUI)/(BLAS.asum(DK)*BLAS.asum(DKUI))
	#println("Alpha_Corr: ", Alpha_corr, " Beta_Corr: ", Beta_corr)
	eps_corr=0.05 #Correlation 0.2 recommended
	CCG=1.0e10	 #10e10 recommended
	if Alpha_corr > eps_corr && Beta_corr > eps_corr
		#println("SSS")
		tau_hat=sqrt(Alpha*Beta)
	elseif Alpha_corr > eps_corr && Beta_corr <= eps_corr
		#println("XXX")
		tau_hat=Alpha
	elseif Alpha_corr <= eps_corr && Beta_corr > eps_corr
		#println("ZZZ")
		tau_hat=Beta
	else
		#println("===")
			tau_hat=tau_previous
	end #END IF ELSE BLOCK
	tau=max(min(tau_hat,(1.0+(CCG/k^2)*tau_previous)),tau_previous/(1.0+(CCG/k^2)))
	#if J.rho != tau
	#	println("Rho changed from ", J.rho, " to ", tau)
	#end
	J.rho=tau
end

#=--------------------------------------------------=#
#===============  Wavelets  =================#
#=--------------------------------------------------=#
	#wavelet_bases=[WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db7,WT.db8,WT.haar]; #Set once at beginning
	#Wx=Array{Float64}(nx*ny,length(wavelet_bases)) #nx and ny will be num_tdf_times and num???
	#DvM=Array{Float64}(nx,ny,length(wavelet_bases))
	#DsM=Array{Float64}(nx,ny,length(wavelet_bases))

	#IWx=Array{Floa64}(nx*ny,length(wavelet_bases)) #nx and ny will be num_tdf_times and num???
	#IDvM=Array{Float64}(nx,ny,length(wavelet_bases))
	#IDsM=Array{Float64}(nx,ny,length(wavelet_bases))


	#X[:,#] gives columns X[#,:] gives rows
	function Dv(mat,Params)
		for i = 1:length(Params.wavelet_bases)
			for j = 1:size(mat)[1]
				Params.DvM[j,:,i]=vec(dwt(vec(mat[j,:]),wavelet(Params.wavelet_bases[i])))
			end
		end
		Params.DvM
	end

	function Ds(mat,Params)
		for i = 1:length(Params.wavelet_bases)
			for j = 1:size(mat)[2]
				Params.DsM[:,j,i]=vec(dwt(vec(mat[:,j]),wavelet(Params.wavelet_bases[i])))
			end
		end
		Params.DsM
	end


#	function Wt(mat)
#		for i:length(wavelet_bases)
#			IWu[:,i]=vec(idwt(reshape(mat[:,i],(nx,ny)), wavelet(wavelet_bases[i])))
#		end
#		sum(IWu,2);
#	end

	function DvT(mat,Params)
		for i = 1:length(Params.wavelet_bases)
			for j=1:size(mat)[1]
				Params.IDvM[j,:,i]=vec(idwt(vec(mat[j,:,i]),wavelet(Params.wavelet_bases[i])))
			end
		end
		reshape(sum(Params.IDvM,3)/length(Params.wavelet_bases),(size(mat)[1],size(mat)[2]));
	end

	function DsT(mat,Params)
		for i = 1:length(Params.wavelet_bases)
			for j=1:size(mat)[2]
				Params.IDsM[:,j,i]=vec(idwt(vec(mat[:,j]),wavelet(Params.wavelet_bases[i])))
			end
		end
		reshape(sum(Params.IDsM,3)/length(Params.wavelet_bases),(size(mat)[1],size(mat)[2]));
	end



	#DvTDv=length(wavelet_bases*eye(nx,ny))
	#DsTDs=length(wavelet_bases*eye(nx,ny))
end #endmodule
