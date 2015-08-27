;arp151.b.dat is the continuum lightcurve
;HJD,FLUX,FLUXERR
;arp151.formatted.trs is the line data. 1st cube is data 2nd is error

;IMPORT CONTINUUM LIGHTCURVE
LIGHTCURVE = 'arp151.b.dat'
READCOL,LIGHTCURVE,HJD,FLUX,FLUXERR

;NORMALIZE THE LIGHTCURVE
NORM_LC = FLUX/MAX(FLUX)

;PLOT THE NORMALIZED LIGHT CURVE
;PLOT,HJD,NORM_LC

;IMPORT JD
DATAFILE = 'rvm_dates.csv'
HJD =  READ_CSV(DATAFILE)

;IMPORT FlUX
DATAFILE = 'rvm_fluxes.csv'
FLUX = READ_CSV(DATAFILE)

;IMPORT FLUX ERROR
DATAFILE = 'rvm_errfluxes.csv'
ERRFLUX = READ_CSV(DATAFILE)

;IMPORT THE WAVELENGTHS
DATAFILE = 'rvm_wavelength.csv'
;IMPORT THE HJD DATA
DATAFILE = 'DataHJD.csv'
HJD = READ_CSV(DATAFILE)



end
