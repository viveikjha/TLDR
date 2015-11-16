#using Optim
using PyPlot
#using Gadfly
include("synthlib.jl")
path="../data/"
#path ="synth/"
#clear()




#=IMPORTING FILES=#
#SPECTRA
wavelength_filename = "rvm_wavelengths.csv"
wavelength_path = string(path,wavelength_filename)
wavelength = readcsv(wavelength_path)                    #List of measured wavelengths

spectra_filename= "rvm_fluxes.csv"
spectrapath = string(path,spectra_filename)
L = readcsv(spectrapath)                           #SPECTRAL FLUXES (L)
num_lines = size(L,1)                                #NUMBER OF SPECTRAL LINES
println("L: ", size(L))
spectra_error_filename = "rvm_errfluxes.csv"
spectra_error_path = string(path,spectra_error_filename)
EL = readcsv(spectra_error_path)                 #SPECTRAL FLUX ERRORS

spectra_dates_filename = "rvm_dates.csv"
spectra_dates_path = string(path,spectra_dates_filename)
spectra_dates = readcsv(spectra_dates_path)                   #SPECTRAL SAMPLING DATES
println(size(spectra_dates))
num_spectra_samples = length(spectra_dates)      				  #NUMBER OF DATA POINTS IN SPECTRA

#CONTINUUM
continuum_array_filename = "arp151.b.dat"
continuum_array_path =string(path,continuum_array_filename)
continuum_array = readdlm(continuum_array_path)
#continuum_array = readcsv(continuum_array_path)
#println(size(continuum_array))
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

tdf_times = collect(1.0:((20.0-1.0)/(num_tdf_times-1)):20.0)

tdf_values = abs(sin(pi*(tdf_times+5.0)/20.0))
#lines
#tdf_values = abs(sin(pi*(lines)/num_lines))

#= Creating Artificial Data Arrays =#
num_lines = 10
num_spectra_samples = 100


X =   zeros((num_tdf_times,num_lines))
#= CONSTANT =#
#initial_x = 1.0																#INITIAL VALUE FOR THE TDF
#fill!(X,initial_x)															#FILL TDF WITH INITIAL VALUE


#for i in collect(1:num_lines)
#	X[:,i] = tdf_values
#end

#= 2D SINE =#
for i in collect(1:num_tdf_times)
	for j in collect(1:num_lines)
		X[i,j] = 1.0*sin( (sin(j/num_lines*pi)*i)/num_tdf_times*pi)
	end
end


#= STRIPE 1 =#
#for i in collect(1:num_tdf_times)
#	for j in collect(1:num_lines)
#		if (j>=3 && j<=7)
#			X[i,j] = 1.0
#		else X[i,j]=0.00001
#		end
#	end
#end

#= STRIPE 2 =#
#for i in collect(1:num_tdf_times)
#	for j in collect(1:num_lines)
#		if (i > 20 && i<30)
#			X[i,j] = 1.0
#		else X[i,j]=0.0000
#		end
#	end
#end



#=DONE SETTING UP TIME DELAY FUNCTION=#


#For Wavlength:
low = 450.0
high = 1000.0
Lrange = high - low
wavelengths = Lrange * sort(rand(num_lines))+low


#For Spectra Dates/Samples:
beginning = minimum(spectra_dates)
ending = maximum(spectra_dates)
range  = ending -beginning
spectra_dates = range * sort(rand(num_spectra_samples))+beginning




#= COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS =#

interpolation_points = zeros(num_spectra_samples,num_tdf_times)
println("interpolation points:",size(interpolation_points))
H = zeros(num_spectra_samples,num_tdf_times)
HE= zeros(num_spectra_samples,num_tdf_times)

for date in collect(1:num_spectra_samples)
	for delay in collect(1:num_tdf_times)
		interpolation_points[date,delay]=spectra_dates[date]-tdf_times[delay]
  end
  P = interpolation_points[date,:]
  H[date,:] =interp(interpolation_points[date,:],continuum_dates,continuum_flux)
  HE[date,:] = interp(interpolation_points[date,:],continuum_dates,continuum_error_flux)
end


#H2 = zeros(num_spectra_samples,num_tdf_times)
#for lns in 1:num_spectra_samples
#  for tdf in 1:num_tdf_times
#    if lns >= tdf
#        H2[lns,tdf] = H[lns,tdf]
#    end
#  end
#end





figure(1)
title("H Matrix")
imshow(H,aspect="auto")




synthetic_data = Model(X,H)		
figure(2)
title("synth data")
imshow(synthetic_data,aspect="auto")

figure(3)
title("VDM")
imshow(X,aspect="auto")




#= Now have the synthetic spectral line (L) and the TDF (x).
Now need to build the output arrays to match the arp151 data. =#


#OUTPUT ARRAYS
#= SETTING UP OUTPUT ARRAYS =#

println("NUM LINES: ",num_lines)
println("NUM SPECTRA: ", num_spectra_samples)
flx_arr = zeros(num_spectra_samples,num_lines)
err_arr = zeros(num_spectra_samples,num_lines)


println("synth_data: ",size(synthetic_data))

#FILL OUPUT ARRAYS!
sigma = 0.05													#FOR ADDING ERROR
for j in 1:num_lines
	for h in 1:num_spectra_samples
		#flx_arr[h,j] = synthetic_data[h,j]+sigma*randn(1)[1]		
		#println(synthetic_data[h])
		#err_arr[h,j] = sigma
		flx_arr[h,j] = synthetic_data[h,j]										#NO ERROR
 		err_arr[h,j] = 1.0																	#NO ERROR 
	end
end





#WRITE OUPUT FILES

#writecsv("filename.csv",array)
writecsv("wavelengthS.csv",wavelengths)
writecsv("rvm_dateS.csv",spectra_dates)
#println("spectra_dates: ", size(spectra_dates))
println("flx_arr: ", size(flx_arr))
#println("err_arr: ", size(err_arr))
#println("continuumS: ", size(continuum_array))
#println("TDF: ", size(X))

println("H: ", size(H))
Modl = Model(X,H)
println("Model: ", size(Modl))
chi2 = Chi2(Modl,flx_arr,err_arr)
println("Chi2 on the VDM: ",chi2 )


writecsv("rvm_flxS.csv", flx_arr)
writecsv("rvm_flx_errS.csv",err_arr)
writecsv("continuumS.csv",continuum_array)
writecsv("TDF.csv",X)
writecsv("Hsy.csv",H)
#println(H[1:20])
#println(X[1:20])

#println(size(flx_arr),size(err_arr),size(continuum_array),size(TDF))
#println(L[1,:])

#figure(1)
#plot(synthetic_data,"b") #PLOT MODEL

#plot(flx_arr[1,:]',"r*")		 #PLOT snyth noise
#plot(flx_arr[2,:]',"y*")		 #PLOT snyth noise
#plot(flx_arr[3,:]',"g*")		 #PLOT snyth noise
#plot(flx_arr[4,:]',"k*")		 #PLOT snyth noise
#plot(flx_arr[5,:]',"m*")		 #PLOT snyth noise
#plot(flx_arr[6,:]',"b*")		 #PLOT snyth noise



#figure(2)
#imshow(X,aspect="auto",cmap="YlOrRd")
#xlabel("Lines")
#ylabel("Delay")
#colorbar()
figure(4)
title("L")
imshow(flx_arr,aspect="auto")



show()


println("Done")
