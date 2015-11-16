include("RMLib.jl")
include("RMTypes.jl")
include("DataImport.jl")
include("GenMatrices.jl")
using PyPlot
#clear()
ion()
#Mode = 1 for synthdata Mode = 2 for real data!
DATA = Import_Data(2)

vdm_path = "synth/TDF.csv"
vdm_act = readcsv(vdm_path)

Pars = init_Params()
Pars.it = 2
Pars.nits =50
Pars.tau = 1.2
Pars.sigma=0.75
Pars.G = 10.0
Pars.num_lines = DATA.num_lines
Pars.num_tdf_times = 50
Pars.initial_psi = 0.5
Pars.eps_abs = 0.0	
Pars.eps_rel = 0.1
Pars.alpha = 1.0

rho0=500.0

#=SETTING UP TIME DELAY FUNCTION=#
Pars.num_tdf_times = 50
Pars.tdf_values = zeros(Pars.num_tdf_times)
num_tdf_times = Pars.num_tdf_times
Pars.tdf_times=collect(1.0  :((20.0-1.0)/(num_tdf_times-1)):20.0)
println(size(Pars.tdf_times))
Mats = Gen_Mats(DATA,Pars)

println("L Pre: ", size(DATA.L))
#DATA.L = DATA.L'
#DATA.EL = DATA.EL'
println("L: ", size(DATA.L))


#println("vdm: ", size(vdm_act))
#println("H: ", size(Mats.H))
#Mod = Model(vdm_act,Mats.H) 
#println("Mod: ", size(Mod))
#chi2_act = Chi2(Mod,DATA.L,DATA.EL)/(DATA.num_lines*Pars.num_tdf_times)
#println("#####################################")
#println("Chi2 on actual vdm: ", chi2_act)
#println("#####################################")

figure()

function TLDR(mu_smoo,mu_spec,mu_temp,DATA,Mats)
	
	#= INITIALIZING ADMM PARAMETERS =#
		initial_psi = 0.5  #Initial value for PSI
		#ADMM ARRAYS:
	Con_Arr = zeros(DATA.num_lines)
	#IMAGE ARRAYS
	X = Gen_Var(rho0,num_tdf_times,DATA.num_lines,initial_psi)
	Z = Gen_Var(rho0,num_tdf_times,DATA.num_lines,initial_psi)
	T = Gen_Var(rho0,num_tdf_times,DATA.num_lines,initial_psi)
	V = Gen_Var(rho0,num_tdf_times,DATA.num_lines,initial_psi)
	#=  Calculate Initial Multipliers & RHO =#
	for p in 1:DATA.num_lines
	  Z.U[:,p]=Mats.HT * squeeze(Mats.W[p,:,:],1) * ( Mats.H * vec(Z.vdm[:,p]) - vec(DATA.L[:,p]))
	end
	V.U = Z.U
	T.U = Z.U
	Qinv = zeros(num_tdf_times,num_tdf_times)
	B = zeros(num_tdf_times,DATA.num_lines)
	converged = 0
	while Pars.it <= Pars.nits && converged==0        #ADMM ITERATION LOOP
		X.vdm_previous = X.vdm
	#Step 1: MINIMIZATION W.R.T. X
	  for l in 1:DATA.num_lines        #SPECTAL CHANNEL LOOP
	#MINIMIZATION ON X_v INDIVIDUALLY WITHIN THIS SPECTRAL LOOP. (THIS PART CAN BE DONE IN PARALLEL!!!!)			
			W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
			Q = Mats.HT * W_slice * Mats.H + T.rho*Mats.DsT*Mats.Ds + (mu_smoo+Z.rho+T.rho)*Mats.Gammatdf 
			B = Mats.HT* W_slice * DATA.L + Mats.DsT*(T.U+T.rho*T.vdm)+Z.U+Z.rho*Z.vdm
			X.vdm[:,l] = Q\B[:,l] 
		end

	#Step 2: UPDATE REGULARIZATION TERMS
		T.vdm_squiggle = Mats.Ds*X.vdm-T.U/T.rho
		T.vdm_previous = T.vdm
		T.vdm = ell2_prox_op(T.vdm_squiggle,mu_temp,T.rho)

		V.vdm_squiggle = Z.vdm*Mats.Dv - V.U/V.rho
		V.vdm_previous = V.vdm
		V.vdm = ell2_prox_op(V.vdm_squiggle,mu_spec,V.rho)
	
		Rinv = inv(V.rho*Mats.Dv*Mats.DvT+Z.rho*Mats.Gammaspe)
		C = (V.U+V.rho*V.vdm)*Mats.DvT-Z.U+(Z.rho*X.vdm)
		Z.vdm_previous = Z.vdm
		Z.vdm = 	C*Rinv
	#Step 4: UPDATE LAGRANGE MULTIPLIERS	
		T.U = LG_update(T.U,T.vdm,Mats.Ds*X.vdm,T.rho,Pars.alpha)
		V.U = LG_update(V.U,V.vdm,Z.vdm*Mats.Dv,V.rho,Pars.alpha)
		Z.U = LG_update(Z.U,Z.vdm,X.vdm,Z.rho,Pars.alpha)
		
  	if Pars.it == 2 
			Pars.G = 1.5
		end
	#Step 5: Update Penalty Parameters
		Z = Pen_update(X,Z,Pars)
		T = Pen_update(X,T,Pars)
		V = Pen_update(X,V,Pars)
	#Step 6: Check for Convergence
		Pars.it = Pars.it+1
		Mod = Model(X.vdm,Mats.H) 
		println("Iteration: ",Pars.it-1, " Chi2: ", Chi2(Mod,DATA.L,DATA.EL)/(num_tdf_times*DATA.num_lines))
		imshow(X.vdm,aspect="auto",cmap="YlOrRd",origin="lower")
		draw()
  end
	iof()
	#outarr #NEED A NEW OUTPUT METHOD....FILES I SUSPECT
	println("ADMM Call Completed")
	
	#Final Plot
	figure(figsize=(20,10))
	set_cmap("YlOrRd")
	subplot(221)
	#figure("Z")	
	imshow((Z.vdm),aspect="auto")
	title("Z")
	colorbar()
	subplot(222)
	#figure("X")
	imshow((X.vdm),aspect="auto")
	title("X")
	colorbar()
	subplot(223)
	#figure("T")
	imshow((T.vdm),aspect="auto")
	title("T")
	colorbar()
	subplot(224)
#	figure("V")
	imshow((V.vdm),aspect="auto")
	title("V")
	colorbar()

	show()
	X
end


mu_spec = 1000.0		#Spectral Regularization Weight
mu_temp = 1.0		#Temporal Regularization Weight
mu_smoo = 1.0    #Smoothing Regularization Weight (TIKHONOV)


tmp = TLDR(mu_smoo,mu_spec,mu_temp,DATA,Mats)


#println("#####################################")
#println("Chi2 on actual vdm: ", chi2_act)
#println("#####################################")


show()

println("Done")


                            











