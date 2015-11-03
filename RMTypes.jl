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
	nits::Int				#Total Iterations Allowed								#2
	tau::Real																								#3
	num_lines::Int																					#4
	num_tdf_times::Int																			#5
	initial_psi::Real																				#6
	mu_smoo::Real		#Tik Smoothing Weight										#7
	mu_spec::Real		#Spectral Reg Weight										#8
	mu_temp::Real		#Temporal Reg Weight										#9
	eps_abs::Real																						#10
	eps_rel::Real																						#11
	sigma::Real																							#12
	G::Real																									#13
	alpha::Real																							#14
end

function init_IMAGE(rho)
	ADMM_IMAGE([],[],[],[],rho,Inf,0.0,0.0,0.0,0.0)
end

function init_Params()
	Params(0,0,0.0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0.0)
end
