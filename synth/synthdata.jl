#using Optim
#using PyPlot
#using Gadfly
include("synthlib.jl")
path="../data/"
#path ="synth/"
clear()



#=IMPORTING FILES=#
#SPECTRA
wavelength_filename = "rvm_wavelengths.csv"
wavelength_path = string(path,wavelength_filename)
wavelength = readcsv(wavelength_path)                    #List of measured wavelengths

spectra_filename= "rvm_fluxes.csv"
spectrapath = string(path,spectra_filename)
L = readcsv(spectrapath)                           #SPECTRAL FLUXES (L)
num_lines = size(L,1)                                #NUMBER OF SPECTRAL LINES

spectra_error_filename = "rvm_errfluxes.csv"
spectra_error_path = string(path,spectra_error_filename)
EL = readcsv(spectra_error_path)                 #SPECTRAL FLUX ERRORS

spectra_dates_filename = "rvm_dates.csv"
spectra_dates_path = string(path,spectra_dates_filename)
spectra_dates = readcsv(spectra_dates_path)                   #SPECTRAL SAMPLING DATES
num_spectra_samples = length(spectra_dates)      				  #NUMBER OF DATA POINTS IN SPECTRA

#CONTINUUM
continuum_array_filename = "arp151.b.dat"
continuum_array_path =string(path,continuum_array_filename)
continuum_array = readdlm(continuum_array_path)
#continuum_array = readcsv(continuum_array_path)
println(size(continuum_array))
#CONTINUUM_ARRAY CONTAINTS THE CONTINUUM DATES, FLUX, AND FLUX ERRORS.
#IN THAT ORDER.

continuum_dates = continuum_array[:,1]
continuum_scale = 1.0
continuum_flux = continuum_array[:,2]*continuum_scale
continuum_error_flux = continuum_array[:,3]*continuum_scale
num_continuum_dates = length(continuum_dates)


########################
#=DONE IMPORTING FILES=#
########################


#=SETTING UP TIME DELAY FUNCTION=#

num_tdf_times = 50
tdf_values = zeros(num_tdf_times)
tdf_times = linspace(1,20,num_tdf_times)

X =   zeros((num_tdf_times))
initial_x = 0.04																#INITIAL VALUE FOR THE TDF
fill!(X,initial_x)															#FILL TDF WITH INITIAL VALUE
#=DONE SETTING UP TIME DELAY FUNCTION=#



#= COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS =#

interpolation_points = zeros(num_spectra_samples,num_tdf_times)
H = zeros(num_spectra_samples,num_tdf_times)
HE= zeros(num_spectra_samples,num_tdf_times)
for date in 1:num_spectra_samples
	for delay in 1:num_tdf_times
		interpolation_points[date,delay]=spectra_dates[date]-tdf_times[delay]
  end
  P = interpolation_points[date,:]
  #ICF[date,:] = interp(P[:],continuum_dates,continuum_flux)
  H[date,:] =interp(interpolation_points[date,:],continuum_dates,continuum_flux)
  #println(interp(interpolation_points[date,:],continuum_dates,continuum_flux))
  HE[date,:] = interp(interpolation_points[date,:],continuum_dates,continuum_error_flux)
end



println("Size of L: ", size(L))

synd = Model(X,H)		

#= No have the synthetic spectral line (L) and the TDF (x).
Now need to build the output arrays to match the arp151 data. =#


#OUTPUT ARRAYS
#= SETTING UP OUTPUT ARRAYS =#

println("NUM LINES: ",num_lines)
println("NUM SPECTRA: ", num_spectra_samples)
flx_arr = zeros(num_lines,num_spectra_samples)
err_arr = zeros(num_lines,num_spectra_samples)
TDF = zeros(num_tdf_times)

println("Data Array: ", size(flx_arr))



#FILL OUPUT ARRAYS!
percent_err = 0.0															#FOR ADDING ERROR
for j in 1:num_lines
	for h in 1:num_spectra_samples
		flx_arr[j,:] = synd
 		err_arr[j,h] = 1.0											#NO ERROR HERE YET.
	end
end

#WRITE OUPUT FILES



#writecsv("filename.csv",array)
writecsv("rvm_flxS.csv", flx_arr)
writecsv("rvm_flx_errS.csv",err_arr)
writecsv("continuumS.csv",continuum_array)
writecsv("TDF.csv",X)
println(size(flx_arr),size(err_arr),size(continuum_array),size(TDF))
println(L[1,:])

#figure(1)
#plot(synd,"r--")
#show()


println("Done")
