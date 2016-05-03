
#dates = readcsv("data/rvm_dates.csv")
flux = readcsv("data/rvm_fluxes.csv")
errflux =readcsv("data/rvm_errfluxes.csv")
waves = readcsv("data/rvm_wavelengths.csv")


#ind = find(6300.<waves.<6800) #Halpha
#ind = find(4650.<waves.<5150) #Hbeta
#ind = find(4100.<waves.<4600)#Hgamma
ind = find(4600.<waves.<5100)
#ind = find(5200.<waves.<5600) #ZERO REGION

println(minimum(waves))
#println(ind)
println("F: ",size(flux))
#datesT = dates[:,ind]
fluxT = flux[ind,:]
errfluxT = errflux[ind,:]
wavesT = waves[ind]

writecsv("data/rvm_fluxes_trimmed.csv",fluxT)
writecsv("data/rvm_errfluxes_trimmed.csv",errfluxT)
writecsv("data/rvm_wavelengths_trimmed.csv",wavesT)
println(wavesT)


#println(size(fluxT))
