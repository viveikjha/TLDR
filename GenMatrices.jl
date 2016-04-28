#include("RMLib.jl")
include("RMTypes.jl")

function Gen_Mats(DATA,Params)
	Mat=init_Mats()
	#= COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS =#
	interpolation_points = zeros(DATA.num_spectra_samples,Params.num_tdf_times)
#	println("interpolation points:",size(interpolation_points))
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
#	Mat.H = H./0.5
	Mat.H = H
	#println("		H: ",size(Mat.H))
	writecsv("H.csv",Mat.H)
	Mat.HE = HE

	#= 		FINITE DIFFERENCES MATRICES			=#
	#Ds = zeros(Params.num_tdf_times,Params.num_tdf_times)
	#for i in collect(1:Params.num_tdf_times)
	#	for j in collect(1:Params.num_tdf_times)
	#		#CHANGED TO 2ND ORDER CENTRAL DIFFERENCES
	#		if i == j
	#			Ds[i,j]= -2.0
	#		end
	#		if (i+1) == j
	#			Ds[i,j] = 1.0
	#		end
	#		if (i-1) == j
	#			Ds[i,j] = 1.0
	#		end
	#	end
	#end

	Ds = zeros(Params.num_tdf_times,Params.num_tdf_times)
	for i in 1:Params.num_tdf_times
		if i == 1
			Ds[i,i]=-1
			Ds[i,i+1]=1
		elseif i == Params.num_tdf_times
			Ds[i,i] = -1
			Ds[i,i-1]= 1
		else
			Ds[i,i+1]=-0.5
			Ds[i,i-1]= 0.5
		end
	end
	s=size(Ds)
	Mat.Ds=eye(s[1],s[2])
	Mat.Ds = Ds
	Mat.DsT = Ds'

	#Dv = zeros(DATA.num_lines,DATA.num_lines)
	#for i in collect(1:DATA.num_lines)
	#	for j in collect(1:DATA.num_lines)
	#		if i == j
	#			Dv[i,j]= -2.0
	#		end
	#		if (i+1) == j
	#			Dv[i,j] = 1.0
	#		end
	#		if (i-1) == j
	#			Dv[i,j] = 1.0
	#		end
	#	end
	#end

	Dv = zeros(DATA.num_lines,DATA.num_lines)
	for i in 1:DATA.num_lines
		if i == 1
			Dv[i,i]=-1
			Dv[i,i+1]=1
		elseif i == DATA.num_lines
			Dv[i,i] = -1
			Dv[i,i-1]= 1
		else
			Dv[i,i+1]=-0.5
			Dv[i,i-1]= 0.5
		end
	end
	Dv=Dv'

	Mat.Dv = Dv
	Mat.DvT = Dv'



	#=    PRECOMPUTING TIKHONOV MATRICES     =#
	num_spectra_dates=size(DATA.spectra_dates)[1]
	W = zeros((DATA.num_lines,size(DATA.L)[1],size(DATA.L)[1]))
	#println("AT W Generation: ", size(W))
	for lam in collect(1:DATA.num_lines)
	  T = eye(num_spectra_dates)
	  for i in collect(1:num_spectra_dates)
       T[i,i] =1.0/ DATA.EL[i,lam].^2
	  end
	  W[lam,:,:] = T
	end
	Mat.W= W
	Mat.HT = H'
	Mat.Gammatdf = eye(Params.num_tdf_times)
	Mat.GammatdfT = Mat.Gammatdf'
	Mat.Gammaspe = eye(DATA.num_lines)
	Mat.GammaspeT = Mat.Gammaspe'
	Mat #Return Data Structure
end
