#include("RMLib.jl")
module GenMatrices
using RMTypes
using RMLibMore
using PyPlot
using JLD
export Gen_Mats
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
	loc=string(Params.directory,"H.csv")
	writecsv(loc,Mat.H)
	Mat.HE = HE


	o=ones(Params.num_tdf_times)
	Ds=spdiagm((-o[1:Params.num_tdf_times-1],o), (-1.0,0.0),Params.num_tdf_times,Params.num_tdf_times)
	s=size(Ds)
	Mat.Ds=eye(s[1],s[2])
	Mat.Ds = Ds
	writecsv(string(Params.directory,"Ds.csv"),Ds)
	Mat.DsT = Ds'

	o=ones(DATA.num_lines)
	Dv=spdiagm((-o[1:DATA.num_lines-1],o), (-1.0,0.0),DATA.num_lines,DATA.num_lines)


	Dv=Dv'#THE SWITCHEROO
	Mat.Dv = Dv
	writecsv(string(Params.directory,"Dv.csv"),Dv)
	Mat.DvT = Dv'


	num_spectra_dates=size(DATA.spectra_dates)[1]

	#println(DATA.num_lines)
	#println(num_spectra_dates)
	Mat.HT = H'

	W = zeros((size(DATA.L)[1],size(DATA.L)[1],DATA.num_lines)) 	#New Version
	HTWH=zeros((size(Mat.HT)[1],size(H)[2],DATA.num_lines))
	HTWL=zeros((size(Mat.HT)[1],size(DATA.L)[2],DATA.num_lines))
	#println("AT W Generation: ", size(W))
		for lam in collect(1:DATA.num_lines)
		  T = eye(num_spectra_dates)
		  for i in collect(1:num_spectra_dates)
	       T[i,i] =1.0/ DATA.EL[i,lam].^2
		  end
		  #W[lam,:,:] = T #Old Version
			W[:,:,lam] = T # New  Version
			tmp1=Mat.HT*T*Mat.H
			tmp2=Mat.HT*T*DATA.L
			HTWH[:,:,lam]=tmp1
			HTWL[:,:,lam]=tmp2
		end


	Mat.W= W
	Mat.Gammatdf = eye(Params.num_tdf_times)
	Mat.GammatdfT = Mat.Gammatdf'
	Mat.Gammaspe = eye(DATA.num_lines)
	Mat.GammaspeT = Mat.Gammaspe'
	Mat.HTWH=HTWH
	writecsv("HTWL.csv",HTWL)
	Mat.HTWL=HTWL
	gc()
	Mat #Return Data Structure
end #endfunction
end #endmodule
