
#dates = readcsv("data/rvm_dates.csv")
flux = readcsv("data/rvm_fluxes.csv")
errflux =readcsv("data/rvm_errfluxes.csv")
waves = readcsv("data/rvm_wavelengths.csv")

ind = find(6300.<waves.<6800)
#println(ind)
println("F: ",size(flux))
#datesT = dates[:,ind]
fluxT = flux[ind,:]
errfluxT = errflux[ind,:]
wavesT = waves[ind]

writecsv("data/rvm_fluxes_trimmed.csv",fluxT)
writecsv("data/rvm_errfluxes_trimmed.csv",errfluxT)
writecsv("data/rvm_wavelengths_trimmed.csv",wavesT)



#println(size(fluxT))
