#using ArrayFire
module RMTypesAF
using ArrayFire
export ADMM_IMAGE,Params,DATA,Matrices_S,init_IMAGE,init_Params,init_Data,init_Mats
	type ADMM_IMAGE
		vdm::AFArray{Float64,2}										#1
		vdm_previous::AFArray{Float64,2}					#2
		vdm_squiggle::AFArray{Float64,2}					#3
		U::AFArray{Float64,2}											#4
		rho::Real																	#5
		rho_max::Real															#6
		rho_min::Real															#7
		phi::Real																	#8
		tau_prim::Real														#9
		tau_dual::Real														#10
	end


	type Params
		it::Int 				#Iteration Counter											#1
		conflag::Bool		#Convergence Flag												#2
		nits::Int				#Total Iterations Allowed								#3
		tau::Real																								#4
		num_lines::Int																					#5
		num_tdf_times::Int																			#6
		initial_psi::Real																				#7
		mu_smoo::Real		#Tik Smoothing Weight										#8
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
		tdf_times::AFArray{Float64,1}														#19
		tdf_values::AFArray{Float64,1}													#20
		l1N_state::Bool																					#21
		l1T_state::Bool																					#22
		l1V_state::Bool																					#23
		AF::Bool
	end


	type Data
		wavelength::ArrayFire.AFArray{Float64,1}
		L::AFArray{Float64,2}
		EL::AFArray{Float64,2}
		num_lines::Int
		num_spectra_samples::Int
		spectra_dates::AFArray{Float64,1}
		continuum_dates::AFArray{Float64,1}
		continuum_flux::AFArray{Float64,1}
		continuum_error_flux::AFArray{Float64,1}
		num_continuum_dates::Int
	end

	type Matrices_S
		H::AFArray{Float64,2}
		HE::AFArray{Float64,2}
		HT::AFArray{Float64,2}
		Ds::AFArray{Float64,2}
		DsT::AFArray{Float64,2}
		Dv::AFArray{Float64,2}
		DvT::AFArray{Float64,2}
		W::AFArray{Float64,3}
		Gammatdf::AFArray{Float64,2}
		GammatdfT::AFArray{Float64,2}
		Gammaspe::AFArray{Float64,2}
		GammaspeT::AFArray{Float64,2}
	end

	function init_IMAGE(rho,dims)
		ADMM_IMAGE(zeros(dims),zeros(dims),zeros(dims),zeros(dims),rho,Inf,0.0,0.0,0.0,0.0)
	end

	function init_Params()
	  Params(2,false,50,1.2,0.0,50,0.1,1.0,1.0,1.0,1.0,0.0,0.01,0.75,10.0,1.0,50.0,Inf,collect(1.0:((50-1.0)/(50-1)):50),zeros(50),true,true,true,false)
	end

	function init_Data(ws,ls,els,sds,cds,cfs,cefs)
		Data(zeros(ws),zeros(ls),zeros(els),0,0,zeros(sds),zeros(cds),zeros(cfs),zeros(cefs),0)
	end

	function init_Mats(W,H,HT,HE,Gt,GtT,Gs,GsT,Dv,DvT,Ds,DsT)
		Matrices_S(zeros(H),zeros(HE),zeros(HT),zeros(Ds),zeros(DsT),zeros(Dv),zeros(DvT),zeros(W),zeros(Gt),zeros(GtT),zeros(Gs),zeros(GsT))
	end
end #end module
