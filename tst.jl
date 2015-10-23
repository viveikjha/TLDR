#Building Difference Matrices
num_spectra_samples = 10
num_tdf_times = 10
D = zeros(num_spectra_samples,num_tdf_times)

for i in 1:num_spectra_samples
	for j in 1:num_tdf_times
		if (i+1) == j
			D[i,j] = 1
		end
		if (i-1) == j
			D[i,j] = -1
		end
	end
end
println(D)

function ell2_prox_op(X,mu)
	X_returned = 1.0/(1.0*2.0*mu)*X
end


