=--------------------------------------------------=#
#=============== Minimization wrt X =================#
#=--------------------------------------------------=#
function min_wrt_x(X,T,P,N,Z,Pars,DATA,Mats)
	s = size(X.vdm)
	vdm = SharedArray(Float64,s[1],s[2])
	for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
			W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
			Q = Mats.HT * W_slice * Mats.H + T.rho*Mats.DsT*Mats.Ds + (Pars.mu_smoo+Z.rho+T.rho+N.rho)*Mats.Gammatdf #INCLUCES L1 NORM ON X
			B = Mats.HT* W_slice * DATA.L + Mats.DsT*(T.U+T.rho*T.vdm)+P.U+P.rho*P.vdm+Z.U+Z.rho*Z.vdm+N.U+N.rho*N.vdm #INCLUDES L1 NORM ON X
			vdm[:,l] = Q\B[:,l]
	end
	X.vdm=sdata(vdm) #sdata() pulls the underlying shared array
	X.vdm
end


=---------------------------------------------------=#
#========== Parallel Minimization wrt X =============#
#=--------------------------------------------------=#
function P_min_wrt_x(X,T,P,N,Z,Pars,DATA,Mats)
	l=collect(1:Pars.num_lines) #WAVELENGTH INDICES
	s = size(X.vdm)
	@everywhere vdm = Array(Float64,s[1],s[2])
	tic()
	vdm = pmap(f,l)
	toc()
	vdm
end


@everywhere function f(l) #THE GUTS OF THE MINIMIZATION FUNCTION
	W_slice = reshape(Mats.W[l,:,:],size(Mats.W[l,:,:])[2],size(Mats.W[l,:,:])[3])
	Q = Mats.HT * W_slice * Mats.H + T.rho*Mats.DsT*Mats.Ds + (Pars.mu_smoo+Z.rho+T.rho+N.rho)*Mats.Gammatdf #INCLUCES L1 NORM ON X
	B = Mats.HT* W_slice * DATA.L + Mats.DsT*(T.U+T.rho*T.vdm)+P.U+P.rho*P.vdm+Z.U+Z.rho*Z.vdm+N.U+N.rho*N.vdm #INCLUDES L1 NORM ON X
	vdm[:,l] = Q\B[:,l]
end

=---------------------------------------------------=#
#================= Parallel Maping ==================#
#=--------------------------------------------------=#

function pmap(f, lst)
    np = nprocs()  # determine the number of processes available
    n = length(lst)
    results = Vector{Any}(n)
    i = 1
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if idx > n
                            break
                        end
                        results[idx] = remotecall_fetch(f, p, lst[idx])
                    end
                end
            end
        end
    end
    results
end
