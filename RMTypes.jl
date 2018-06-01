module RMTypes
using Wavelets
export ADMM_IMAGE,Data,Params,Fit_Params,DATA,Matrices_S,init_IMAGE,init_Params,init_Data,init_Mats,init_fit,Init_Record_Data
#export ADMM_IMAGE,Fit_Params,Params,Matrices_S,Data,Rec_Data
mutable struct ADMM_IMAGE
	vdm::Array{Float64}										#1
	vdm_previous::Array{Float64}					#2
	vdm_squiggle::Array{Float64}					#3
	U::Array{Float64}											#4
	UI::Array{Float64}										#5
	UI_previous::Array{Float64}						#6
	rho::Real															#7
	rho_max::Real													#8
	rho_min::Real													#9
	phi::Real															#10
	tau_prim::Real												#11
	tau_dual::Real												#12
	tau_max::Real													#13
	tau::Real															#14
end

mutable struct Fit_Params
	msmo::Real
	pz::Real
	ml1::Real
	pn::Real
	mspe::Real
	pv::Real
	mtem::Real
	pt::Real
	pp::Real
	TI::Real
	TI2::Real
	waves::Bool
	fast::Bool
end

mutable struct Params
	it::Int 				#Iteration Counter										#1
	conflag::Int		#Convergence Flag											#2
	nits::Int				#Total Iterations Allowed								#3
	tau::Real																		#4
	num_lines::Int																	#5
	num_tdf_times::Int																#6
	initial_psi::Real																#7
	mu_smoo::Real		#Tik Smoothing Weight										#8
	mu_l1::Real			#Ell1 Reg weight											#9
	mu_spec::Real		#Spectral Reg Weight										#10
	mu_temp::Real		#Temporal Reg Weight										#11
	eps_abs::Real																	#12
	eps_rel::Real																	#13
	sigma::Real																		#14
	G::Real																			#15
	alpha::Real																		#16
	rho0::Real																		#17
	chi2::Real																		#18
	tdf_times::Array{Float64}														#19
	tdf_values::Array{Float64}														#20
	l1N_state::Bool																	#21
	l1T_state::Bool																	#22
	l1V_state::Bool																	#23
	directory::String
	threshold::Real
	wavelet_bases::AbstractArray{}
	DvM::Array{Float64}
	DsM::Array{Float64}
	IDvM::Array{Float64}
	IDsM::Array{Float64}
end

mutable struct Data
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

mutable struct Matrices_S
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
	HTWH::Array{Float64}
	HTWL::Array{Float64}
	Qinv::Array{Float64}
	Rinv::Array{Float64}
end

mutable struct Rec_Data
	ConFlag::Array{Float64}
	J::Array{Float64}
	Chi2::Array{Float64}
	Z_res::Array{Float64}
	N_res::Array{Float64}
	T_res::Array{Float64}
	V_res::Array{Float64}
	P_res::Array{Float64}
end

function Init_Record_Data(nits)
	Rec_Data(ones(nits)+5,zeros(nits),zeros(nits),zeros(nits),zeros(nits),zeros(nits),zeros(nits),zeros(nits))
end

function init_IMAGE(rho,dims)
	#           1, 2, 3, 4, 5, 6,  7,  8,  9, 10, 11, 12,   13, 14
	ADMM_IMAGE(Array{Float64}(dims),Array{Float64}(dims),Array{Float64}(dims),Array{Float64}(dims),Array{Float64}(dims),Array{Float64}(dims),rho,Inf,0.0,0.0,0.0,0.0,1.0e3,2.0)
end

function init_Params()
  Params(2,0,50,2.0,0.0,50,0.1,1.0,1.0,1.0,1.0,0.0,0.01,0.75,10.0,1.0,50.0,Inf,collect(1.0:((50-1.0)/(50-1)):50),zeros(50),true,true,true,"",1.0e-7,[],[],[],[],[])
end

function init_Data()
	Data([],[],[],0,0,[],[],[],[],0)
end

function init_Mats()
	Matrices_S([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[])
end

function init_fit()
	Fit_Params(1.0e3,1.0e3,1.0e3,1.0e3,1.0e3,1.0e3,1.0e3,1.0e3,1.0e3,NaN,NaN,true,false)
end

end #end module
