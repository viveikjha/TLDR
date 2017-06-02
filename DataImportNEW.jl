module DataImportNEW
if "ArrayFire" in keys(Pkg.installed())
  AF=true
  using ArrayFire
	using RMTypesAF
else
	using RMTypes
end

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
  w_size=size(wavelength)
	#println("Read: ", wavelength_path, " size: ", size(Data_Arrs.wavelength))
	spectrapath = string(path,spectra_filename)

	L = spec_scale.*(readcsv(spectrapath))'                      			#SPECTRAL FLUXES (L)
  L_size=size(L)
	#println("Read: ",spectrapath, " size: ",size(Data_Arrs.L))
	num_lines = size(L,2)                  			#NUMBER OF SPECTRAL LINES
#	Data_Arrs.num_lines = size(Data_Arrs.L,1)                  			#NUMBER OF SPECTRAL LINES

	spectra_error_path = string(path,spectra_error_filename)
	EL = spec_scale.*(readcsv(spectra_error_path))'            			#SPECTRAL FLUX ERRORS
  EL_size=size(EL)
	#println("Read: ",spectra_error_path, " size: ",size(Data_Arrs.EL))


	spectra_dates_path = string(path,spectra_dates_filename)
	spectra_dates = readcsv(spectra_dates_path)           #SPECTRAL SAMPLING DATES
  sd_size=size(spectra_dates)
	#println("Read: ",spectra_dates_path, " size: ",size(Data_Arrs.spectra_dates))


	num_spectra_samples = length(spectra_dates) #NUMBER OF DATA POINTS IN SPECTRA
	#CONTINUUM
	#CONTINUUM_ARRAY CONTAINTS THE CONTINUUM DATES, FLUX, AND FLUX ERRORS.
	#IN THAT ORDER.
	continuum_dates = continuum_array[:,1]
  cd_size=size(continuum_dates)
	continuum_scale = 1.0
	continuum_flux = cont_scale.*continuum_array[:,2]*cont_scale
  cf_size=size(continuum_flux)
	continuum_error_flux = cont_scale.*continuum_array[:,3]
  cef_size=size(continuum_error_flux)
	num_continuum_dates = length(continuum_dates)
	if Reports == true
		println("Import Summary:")
		println( "		Wavelengths: ", length(wavelength))
		println( "		Spectra: ", size(L))
		println( "		Continuum: ", size(continuum_flux))
	end

	#ws,ls,els,sds,cds,cfs,cefs)

  Data_Arrs = init_Data(w_size[1],L_size,EL_size,sd_size[1],cd_size,cf_size,cef_size)
	if AF == true
		println("Hello")
		#println(typeof(Data_Arrs))
		#println(typeof(wavelength))

		Data_Arrs.L =AFArray(L)
		Data_Arrs.num_lines=num_lines
		Data_Arrs.EL=AFArray(EL)
    Data_Arrs.wavelength=AFArray(vec(wavelength))
		Data_Arrs.spectra_dates=AFArray(vec(spectra_dates))
    Data_Arrs.num_spectra_samples=num_spectra_samples
		Data_Arrs.continuum_flux=AFArray(continuum_flux)
		Data_Arrs.continuum_error_flux=AFArray(continuum_error_flux)
		Data_Arrs.continuum_dates=AFArray(continuum_dates)
		Data_Arrs.num_continuum_dates = num_continuum_dates
	else
		Data_Arrs.wavelength=wavelenth
		Data_Arrs.L =L
		Data_Arrs.num_lines=num_lines
		Data_Arrs.EL=EL
		Data_Arrs.spectra_dates=spectra_dates
    Data_Arrs.num_spectra_samples=num_spectra_samples
		Data_Arrs.continuum_flux=continuum_flux
		Data_Arrs.continuum_error_flux=continuum_error_flux
		Data_Arrs.continuum_dates=continuum_dates
		Data_Arrs.num_continuum_dates = num_continuum_dates
end
		Data_Arrs #RETURNED

end
end #endmodule
