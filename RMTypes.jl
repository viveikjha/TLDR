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
	rho0::Real																							#15
	chi2::Real																							#16
	tdf_times::Array{Float64}																#17
	tdf_values::Array{Float64}															#18
	

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
	Params(0,0,0.0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0.0,0.0,Inf,[],[])
end

function init_Data()
	Data([],[],[],0,0,[],[],[],[],0)
end

function init_Mats()
	Matrices_S([],[],[],[],[],[],[],[],[],[],[],[])
end

