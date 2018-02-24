module DataImportNEW

using RMTypes
export Import_DataN

#include("RMTypes.jl")
function Import_DataN(path="data/",waves="rvm_wavelengths.csv",spectra="rvm_fluxes.csv",spectraerr="rvm_errfluxes.csv",spectradate="rvm_dates.csv",continuum="arp151.b.dat";cont_scale=1.0,spec_scale=1.0,Reports=false)


#Import_DataN("simulation/","simulated_vdm.csv","spectra_simulated.csv", "errspectra_simulated.csv","rvm_dates.csv","arp151.b.dat")
	#=IMPORTING FILES=#
	path = path
	wavelength_filename = waves
	spectra_filename= spectra
	spectra_error_filename = spectraerr
	spectra_dates_filename = spectradate
	continuum_array_filename = continuum

	#CONTINUUM
	continuum_array_path = string(path,continuum_array_filename)
	continuum_array = readdlm(continuum_array_path)
	#println("Read: ", continuum_array_path, " size: ", size(continuum_array))


	#SPECTRA
	spec_scale=1.0
	wavelength_path = string(path,wavelength_filename)
	wavelength = readcsv(wavelength_path)         #List of measured wavelengths
	#println("Read: ", wavelength_path, " size: ", size(wavelength))
	spectrapath = string(path,spectra_filename)

	L = spec_scale.*(readcsv(spectrapath))'                      			#SPECTRAL FLUXES (L)
#	println("Read: ",spectrapath, " size: ",size(L))
	num_lines = size(L,2)                  			#NUMBER OF SPECTRAL LINES
#	num_lines = size(L,1)                  			#NUMBER OF SPECTRAL LINES

	spectra_error_path = string(path,spectra_error_filename)
	EL = spec_scale.*(readcsv(spectra_error_path))'            			#SPECTRAL FLUX ERRORS
	#println("Read: ",spectra_error_path, " size: ",size(EL))


	spectra_dates_path = string(path,spectra_dates_filename)
	spectra_dates = readcsv(spectra_dates_path)           #SPECTRAL SAMPLING DATES
	#println("Read: ",spectra_dates_path, " size: ",size(spectra_dates))


	num_spectra_samples = length(spectra_dates) #NUMBER OF DATA POINTS IN SPECTRA
	#CONTINUUM
	#CONTINUUM_ARRAY CONTAINTS THE CONTINUUM DATES, FLUX, AND FLUX ERRORS.
	#IN THAT ORDER.
	continuum_dates = continuum_array[:,1]
	continuum_scale = 1.0
	continuum_flux = cont_scale.*continuum_array[:,2]*cont_scale
	continuum_error_flux = cont_scale.*continuum_array[:,3]
	num_continuum_dates = length(continuum_dates)
	Data_Arrs = Data(wavelength,L,EL,num_lines,num_spectra_samples,spectra_dates,continuum_dates,continuum_flux,continuum_error_flux,num_continuum_dates)

	if Reports == true
		println("Import Summary:")
		println( "		Wavelengths: ", length(wavelength))
		println( "		Spectra: ", size(L))
		println( "		Continuum: ", size(continuum_flux))
	end
	Data_Arrs
end
end #endmodule
