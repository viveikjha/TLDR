include("RMTypes.jl")
function Import_Data(mode)
		#mode = 1 	#SYNTHETIC DATA MODE
		#mode = 2	#REAL DATA MODE
		#mode =3 #SIMULATED DATA
	if mode == 1
		path ="synth/"
		wavelength_filename = "wavelengthS.csv"
		spectra_filename= "rvm_flxS.csv"
		spectra_error_filename = "rvm_flx_errS.csv"
		spectra_dates_filename = "rvm_dateS.csv"
		continuum_array_filename = "continuumS.csv"
		continuum_array_path =string(path,continuum_array_filename)
		continuum_array = readcsv(continuum_array_path)
	elseif mode == 2
		path="data/"
		wavelength_filename = "rvm_wavelengths.csv"
		spectra_filename= "rvm_fluxes.csv"
		spectra_error_filename = "rvm_errfluxes.csv"
		spectra_dates_filename = "rvm_dates.csv"
		continuum_array_filename = "arp151.b.dat"
		continuum_array_path =string(path,continuum_array_filename)
		continuum_array = readdlm(continuum_array_path)
	elseif mode == 3 
		path="simulation/"
#		wavelength_filename = "rvm_wavelengths.csv"
		wavelength_filename = "simulated_vdm.csv"
		spectra_filename= "spectra_simulated.csv"
#		spectra_error_filename = "rvm_errfluxes.csv"
		spectra_error_filename = "errspectra_simulated.csv"
		spectra_dates_filename = "rvm_dates.csv"
		continuum_array_filename = "arp151.b.dat"
		continuum_array_path =string(path,continuum_array_filename)
		continuum_array = readdlm(continuum_array_path)
	end



	Data_Arrs = init_Data()

	#=IMPORTING FILES=#
	#SPECTRA
	wavelength_path = string(path,wavelength_filename)
	Data_Arrs.wavelength = readcsv(wavelength_path)           			#List of measured wavelengths
	spectrapath = string(path,spectra_filename)
	Data_Arrs.L = (readcsv(spectrapath))'                      			#SPECTRAL FLUXES (L)
	println("AT IMPORT: ",size(Data_Arrs.L))
	Data_Arrs.num_lines = size(Data_Arrs.L,2)                  			#NUMBER OF SPECTRAL LINES
#	Data_Arrs.num_lines = size(Data_Arrs.L,1)                  			#NUMBER OF SPECTRAL LINES
	
	spectra_error_path = string(path,spectra_error_filename)
	Data_Arrs.EL = (readcsv(spectra_error_path))'            			#SPECTRAL FLUX ERRORS
	spectra_dates_path = string(path,spectra_dates_filename)
	Data_Arrs.spectra_dates = readcsv(spectra_dates_path)           #SPECTRAL SAMPLING DATES
	Data_Arrs.num_spectra_samples = length(Data_Arrs.spectra_dates) #NUMBER OF DATA POINTS IN SPECTRA
	#CONTINUUM
	#CONTINUUM_ARRAY CONTAINTS THE CONTINUUM DATES, FLUX, AND FLUX ERRORS.
	#IN THAT ORDER.
	Data_Arrs.continuum_dates = continuum_array[:,1]
	continuum_scale = 1.0
	Data_Arrs.continuum_flux = continuum_array[:,2]*continuum_scale
	Data_Arrs.continuum_error_flux = continuum_array[:,3]*continuum_scale
	Data_Arrs.num_continuum_dates = length(Data_Arrs.continuum_dates)
	
println("Imported Data Arrays:")
println( "		Wavelengths: ", length(Data_Arrs.wavelength))
println( "		Spectra: ", size(Data_Arrs.L))
println( "		Continuum: ", size(Data_Arrs.continuum_flux))
	Data_Arrs
end

