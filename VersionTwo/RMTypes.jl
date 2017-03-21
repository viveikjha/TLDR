type ADMM_IMAGE
	vdm::Array{Float64}										#1
	vdm_previous::Array{Float64}					#2
	vdm_squiggle::Array{Float64}					#3
	U::Array{Float64}											#4
	rho::Real															#5
	rho_max::Real													#6
	rho_min::Real													#7
	phi::Real															#8
	tau_prim::Real												#9
	tau_dual::Real												#10
end


type Params
	it::Int 				#Iteration Counter											#1
	conflag::Bool		#Convergence Flag												#2
	nits::Int				#Total Iterations Allowed								#3
	tau::Real																								#4
	num_lines::Int																					#5
	num_tdf_times::Int																			#6
	initial_psi::Real																				#7
	mu_l2::Real		#Tik Smoothing Weight										#8
	mu_l1::Real			#Ell1 Reg weight												#9
	mu_spec::Real		#Spectral Reg Weight										#10
	mu_temp::Real		#Temporal Reg Weight										#11
	eps_abs::Real																						#12
	eps_rel::Real																						#13
	sigma::Real																							#14
	G::Real																									#15
	alpha::Real																							#16
	rho0::Real																							#17
	chi2::Real																							#18
	tdf_times::Array{Float64}																#19
	tdf_values::Array{Float64}															#20
	l1N_state::Bool																					#21
	l1T_state::Bool																					#22
	l1V_state::Bool																					#23
end

type Data
	wavelength::Array{Float64}
	L::Array{Float64}
	EL::Array{Float64}
	num_lines::Int
	num_spectra_samples::Int
	spectra_dates::Array{Float64}
	continuum_dates::Array{Float64}
	continuum_flux::Array{Float64}
	continuum_error_flux::Array{Float64}
	num_continuum_dates::Int
end

type Matrices_S
	H::Array{Float64}
	HE::Array{Float64}
	HT::Array{Float64}
	Ds::Array{Float64}
	DsT::Array{Float64}
	Dv::Array{Float64}
	DvT::Array{Float64}
	W::Array{Float64}
	Gammatdf::Array{Float64}
	GammatdfT::Array{Float64}
	Gammaspe::Array{Float64}
	GammaspeT::Array{Float64}
end

function init_IMAGE(rho)
	ADMM_IMAGE([],[],[],[],rho,Inf,0.0,0.0,0.0,0.0)
end

function init_Params()
  Params(2,false,50,1.2,0.0,50,0.1,1.0,1.0,1.0,1.0,0.0,0.01,0.75,10.0,1.0,50.0,Inf,collect(1.0:((50-1.0)/(50-1)):50),zeros(50),true,true,true)
end

function init_Data()
	Data([],[],[],0,0,[],[],[],[],0)
end

function init_Mats()
	Matrices_S([],[],[],[],[],[],[],[],[],[],[],[])
end
