include("RMLib.jl")
include("RMTypes.jl")

function Gen_Mats(DATA,Params)
	Mat=init_Mats()
	#= COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS =#
	interpolation_points = zeros(DATA.num_spectra_samples,Params.num_tdf_times)
	println("interpolation points:",size(interpolation_points))
	H = zeros(DATA.num_spectra_samples,Params.num_tdf_times)
	HE= zeros(DATA.num_spectra_samples,Params.num_tdf_times)
	for date in collect(1:DATA.num_spectra_samples)
		for delay in collect(1:Params.num_tdf_times)
			interpolation_points[date,delay]=DATA.spectra_dates[date]-Params.tdf_times[delay]
	  end
	  P = interpolation_points[date,:]
	  H[date,:] =interp(interpolation_points[date,:],DATA.continuum_dates,DATA.continuum_flux)
	  HE[date,:] = interp(interpolation_points[date,:],DATA.continuum_dates,DATA.continuum_error_flux)
	end

	Mat.H = H
	Mat.HE = HE
	#= 		FINITE DIFFERENCES MATRICES			=#
	Ds = zeros(Params.num_tdf_times,Params.num_tdf_times)
	for i in collect(1:Params.num_tdf_times)
		for j in collect(1:Params.num_tdf_times)
			#CHANGED TO 2ND ORDER CENTRAL DIFFERENCES
			if i == j
				Ds[i,j]==-2.0
			end

			if (i+1) == j
				Ds[i,j] = 1.0
			end
			if (i-1) == j
				Ds[i,j] = 1.0
			end
		end
	end
	Mat.Ds = Ds
	Mat.DsT = Ds'
	Dv = zeros(DATA.num_lines,DATA.num_lines)
	for i in collect(1:DATA.num_lines)
		for j in collect(1:DATA.num_lines)
			if i == j
				Dv[i,j]==-2.0
			end
			if (i+1) == j
				Dv[i,j] = 1.0
			end
			if (i-1) == j
				Dv[i,j] = 1.0
			end
		end
	end
	Mat.Dv = Dv
	Mat.DvT = Dv'

	#=    PRECOMPUTING TIKHONOV MATRICES     =#
	#println("DATA.L: ", size(DATA.L))
	#println("DATA.EL: ", size(DATA.EL))
	#println(collect(1:DATA.num_lines))
	num_spectra_dates=size(DATA.spectra_dates)[1]
	#println("# spec dates: ", num_spectra_dates)
	#println("num_lines: ", DATA.num_lines)
	W = zeros((DATA.num_lines,size(DATA.L)[1],size(DATA.L)[1]))
	#println(size(W))
	for lam in collect(1:DATA.num_lines)
	  T = eye(num_spectra_dates)
	  for i in collect(1:num_spectra_dates)
	    for j in collect(1:num_spectra_dates)
				#println(lam, "  ",i, "  ",j)
	      if i == j
	        T[i,j] =1.0/ DATA.EL[i,lam].^2
	      end
	    end
	  end
	  W[lam,:,:] = T
	end

	

	Mat.W= W
	Mat.HT = H'
	Mat.Gammatdf = eye(Params.num_tdf_times)
	Mat.GammatdfT = Mat.Gammatdf'
	Mat.Gammaspe = eye(DATA.num_lines)
	Mat.GammaspeT = Mat.Gammaspe'
	Mat
end




