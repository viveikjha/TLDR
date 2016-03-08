include("RMTypes.jl")

function Import_DataN(path="data/",waves="rvm_wavelengths.csv",spectra="rvm_fluxes.csv",spectraerr="rvm_errfluxes.csv",spectradat="rvm_dates.csv",continuum="arp151.b.dat")
	Data_Arrs = init_Data()

#Import_DataN("simulation/","simulated_vdm.csv","spectra_simulated.csv", "errspectra_simulated.csv","rvm_dates.csv","arp151.b.dat")
	#=IMPORTING FILES=#
	path = path
	wavelength_filename = waves
	spectra_filename= spectra
	spectra_error_filename = spectraerr
	spectra_dates_filename = spectradat
	continuum_array_filename = continuum

	#CONTINUUM
	continuum_array_path = continuum_array_filename#string(path,continuum_array_filename)
	continuum_array = readdlm(continuum_array_path)
	println("Read: ", continuum_array_path, " size: ", size(continuum_array))


	#SPECTRA
	wavelength_path = string(path,wavelength_filename)
	Data_Arrs.wavelength = readcsv(wavelength_path)           			#List of measured wavelengths
	println("Read: ", wavelength_path, " size: ", size(Data_Arrs.wavelength))
	spectrapath = string(path,spectra_filename)

	Data_Arrs.L = (readcsv(spectrapath))'                      			#SPECTRAL FLUXES (L)
	println("Read: ",spectrapath, " size: ",size(Data_Arrs.L))
	Data_Arrs.num_lines = size(Data_Arrs.L,2)                  			#NUMBER OF SPECTRAL LINES
#	Data_Arrs.num_lines = size(Data_Arrs.L,1)                  			#NUMBER OF SPECTRAL LINES
	
	spectra_error_path = string(path,spectra_error_filename)
	Data_Arrs.EL = (readcsv(spectra_error_path))'            			#SPECTRAL FLUX ERRORS
	println("Read: ",spectra_error_path, " size: ",size(Data_Arrs.EL))


	spectra_dates_path = spectra_dates_filename#string(path,spectra_dates_filename)
	Data_Arrs.spectra_dates = readcsv(spectra_dates_path)           #SPECTRAL SAMPLING DATES
	println("Read: ",spectra_dates_path, " size: ",size(Data_Arrs.spectra_dates))


	Data_Arrs.num_spectra_samples = length(Data_Arrs.spectra_dates) #NUMBER OF DATA POINTS IN SPECTRA
	#CONTINUUM
	#CONTINUUM_ARRAY CONTAINTS THE CONTINUUM DATES, FLUX, AND FLUX ERRORS.
	#IN THAT ORDER.
	Data_Arrs.continuum_dates = continuum_array[:,1]
	continuum_scale = 1.0
	Data_Arrs.continuum_flux = continuum_array[:,2]*continuum_scale
	Data_Arrs.continuum_error_flux = continuum_array[:,3]*continuum_scale
	Data_Arrs.num_continuum_dates = length(Data_Arrs.continuum_dates)
	
	println("Import Summary:")
	println( "		Wavelengths: ", length(Data_Arrs.wavelength))
	println( "		Spectra: ", size(Data_Arrs.L))
	println( "		Continuum: ", size(Data_Arrs.continuum_flux))
	Data_Arrs
end

