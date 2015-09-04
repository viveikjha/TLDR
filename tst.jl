include("RMLib.jl")

path ="synth/"
clear()
#=IMPORTING FILES=#
#SPECTRA
wavelength_filename = "rvm_wavelengths.csv"
wavelength_path = string(path,wavelength_filename)
wavelength = readcsv(wavelength_path)                    #List of measured wavelengths

#spectra_filename= "rvm_fluxes.csv"
spectra_filename= "rvm_flxS.csv"
spectrapath = string(path,spectra_filename)
L = readcsv(spectrapath)                           #SPECTRAL FLUXES (L)
num_lines = size(L,1)                                #NUMBER OF SPECTRAL LINES

#spectra_error_filename = "rvm_errfluxes.csv"
spectra_error_filename = "rvm_flx_errS.csv"
spectra_error_path = string(path,spectra_error_filename)
EL = readcsv(spectra_error_path)                 #SPECTRAL FLUX ERRORS

spectra_dates_filename = "rvm_dates.csv"
spectra_dates_path = string(path,spectra_dates_filename)
spectra_dates = readcsv(spectra_dates_path)                   #SPECTRAL SAMPLING DATES
num_spectra_samples = length(spectra_dates)      				  #NUMBER OF DATA POINTS IN SPECTRA

#CONTINUUM
#continuum_array_filename = "arp151.b.dat"
continuum_array_filename = "continuumS.csv"
continuum_array_path =string(path,continuum_array_filename)
#continuum_array = readdlm(continuum_array_path)
continuum_array = readcsv(continuum_array_path)
println(size(continuum_array))
#CONTINUUM_ARRAY CONTAINTS THE CONTINUUM DATES, FLUX, AND FLUX ERRORS.
#IN THAT ORDER.

continuum_dates = continuum_array[:,1]
continuum_scale = 1.0
continuum_flux = continuum_array[:,2]*continuum_scale
continuum_error_flux = continuum_array[:,3]*continuum_scale
num_continuum_dates = length(continuum_dates)
#=DONE IMPORTING FILES=#

#=SETTING UP TIME DELAY FUNCTION=#
num_tdf_times = 50
tdf_values = zeros(num_tdf_times)
tdf_times = linspace(1,20,num_tdf_times)
#=DONE SETTING UP TIME DELAY FUNCTION=#

#= COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS =#
interpolation_points = zeros(num_spectra_samples,num_tdf_times)
#println(size(interpolation_points))
ICF = zeros(num_spectra_samples,num_tdf_times)
ICFE= zeros(num_spectra_samples,num_tdf_times)
for date in 1:num_spectra_samples
	for delay in 1:num_tdf_times
		interpolation_points[date,delay]=spectra_dates[date]-tdf_times[delay]
  end
  P = interpolation_points[date,:]
  #ICF[date,:] = interp(P[:],continuum_dates,continuum_flux)
  ICF[date,:] =interp(interpolation_points[date,:],continuum_dates,continuum_flux)
  #println(interp(interpolation_points[date,:],continuum_dates,continuum_flux))
  ICFE[date,:] = interp(interpolation_points[date,:],continuum_dates,continuum_error_flux)
end
#println("------------------------------------------")
#println(ICF)
#println("------------------------------------------")


#=
plot(continuum_dates,continuum_flux)
plot(interpolation_points[1,:],ICF[1,:],"r*")
show()
=#

#for i in 1:10
#  println(interpolation_points[i,1],"     ",ICF[i,1])
#end

#=    PRECOMPUTING TIKHONOV MATRICES     =#

#Build Mapping Matrix
num_spectra_dates=size(spectra_dates)[1]
H = zeros((num_spectra_dates, num_tdf_times))
for lns in 1:num_spectra_dates
  for tdf in 1:num_tdf_times
    if lns >= tdf
        H[lns,tdf] = ICF[lns,tdf]
    end
  end
end
#println("H!: ", H)
#Build Covariance Matrix
W = zeros((size(wavelength)[1],size(L)[2],size(L)[2]))
for lam in 1:num_lines
  T = eye(num_spectra_dates)
  for i in 1:num_spectra_dates
    for j in 1:num_spectra_dates
      #println("!")
      if i == j
        T[i,j] =1.0/ EL[lam,i].^2
      end
    end
  end
  W[lam,:,:] = T
end

HT = H'
Gamma = eye(num_tdf_times)
GammaT = Gamma'

#= INITIALIZING ADMM PARAMETERS =#
nits = 20 #Number of ADMM Iterations (An upper limit if convergence is not reached)
final_it = nits
initial_psi = 0.04  #Initial value for PSI

	#= REGULARIZATION WEIGHTS =#
mu = 5.0

  #= CONVERGANCE PARAMETERS =#
G = 10.0 # only for first iteration 1.5 otherwise
tau = 1.2
Con_Arr = zeros(num_lines)

X =   zeros((nits+1, size(L,1),num_tdf_times))
X_s = zeros((nits+1, size(L,1),num_tdf_times))
Z =   zeros((nits+1, size(L,1),num_tdf_times))
Z_s = zeros((nits+1, size(L,1),num_tdf_times))
rho = zeros((nits+1, size(L,1)))
U =   zeros((nits+1, size(L,1),num_tdf_times))
tau_prim = zeros((nits,num_lines))
tau_dual = zeros((nits,num_lines))
value = 1.0
initial_rho = 2.0
fill!(X,initial_psi)
fill!(X_s,initial_psi)
fill!(Z,initial_psi)
fill!(Z_s,initial_psi)
fill!(rho,initial_rho)


X_act = zeros((50))
fill!(X_act,initial_psi)

#=  Calculate Initial Multipliers  =#
for p in 1:num_lines
#for p in 1:2
  U[1,p,:]=HT * squeeze(W[p,:,:],1) * ( H * vec(Z[1,p,:]) - vec(L[p,:]))
end

X_act = zeros((50))
fill!(X_act,0.04)




#Need to look at Model:println("X_act: ",X_act)
Mod_act = Model(X_act,H)
mod_sct = Model(X_act,ICF)



println("modact(with H): ", Mod_act[1])
println("modact(with ICF): ", Model(X_act,ICF)[1])
println("initial chi2: ", Chi2(mod_sct,L[1,:],EL[1,:])/num_tdf_times)
println("--------------------------------------")
println("Model: ", mod_sct)
println("--------------------------------------")
println("Data: ", L[1,:])
println("--------------------------------------")
diff = mod_sct - vec(L[1,:])
println("Summed Difference: ",sum(diff) )
println("--------------------------------------")
println("Error: ", EL[1,:])
println("--------------------------------------")
println("Chi2: ", sum((diff./vec(EL[1,:])).^2)/num_tdf_times)



















