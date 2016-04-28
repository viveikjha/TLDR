#THIS FILE CONTAINS THE BASE LEVEL FUNCTIONS FOR TLDR THAT ARE NOT PART OF THE CORE ADMM ALGORITHM
#Inputs include all the reconstruction variables. Options are for printed lines.


#=--------------------------------------------------=#
#=========== THE Minimization FUNCTION ==============#
#=--------------------------------------------------=#
#Gives the value for the minimization problem given parameters.
function J(X,P,T,V,N,DATA,Mats,Pars)
	(0.5*Pars.chi2*DATA.num_spectra_samples*DATA.num_lines)+regX(X,Pars)+regT(T,Pars)+regV(V,Pars)+regN(N,Pars)
end
#=--------------------------------------------------=#
#=========== Regularization term on X ===============#
#=--------------------------------------------------=#
function regX(X,Pars)
	Pars.mu_smoo*0.5*ell2normsquared(X.vdm)
end
#=--------------------------------------------------=#
#========== Regularization term on T,V,N ============#
#=--------------------------------------------------=#
#These functions give the value of the given regularization term.
#The functions are identical, which is silly.
function regT(T,Pars)
	Pars.mu_temp*ell1norm(T.vdm)
end

function regV(V,Pars)
	Pars.mu_spec*ell1norm(V.vdm)
end

function regN(N,Pars)
	Pars.mu_l1*ell1norm(N.vdm)
end



#=--------------------------------------------------=#
#===================== REPORT =======================#
#=--------------------------------------------------=#
#This function prints a numerical report of the desired parameters listed
#for the provided images X,Z,P,T,V,N, and data in DATA, Mats, Pars.
function Report(X,Z,P,T,V,N,DATA,Mats,Pars;Jf=false,L2x=false,L1T=false,L1V=false,L1N=false,Chi2x=false,Chi2z=false,s=false,Pres=false,Zres=false,Tres=false,Vres=false,Nres=false,Ppen=false,Zpen=false,Tpen=false,Vpen=false,Npen=false,Msg="")
#Values
  if Jf==true
    Jstring=@sprintf "\tJ: %0.3f" J(X,P,T,V,N,DATA,Mats,Pars)
  else Jstring=""
  end
  if L2x==true
    L2xstring=@sprintf "\tL2x: %0.0f" regX(X,Pars)
  else L2xstring =""
  end
  if L1T==true
    L1Tstring=@sprintf "\tL1T: %0.1f" regT(T,Pars)
  else L1Tstring=""
  end
  if L1V==true
    L1Vstring=@sprintf "\tL1V: %0.1f" regV(V,Pars)
  else L1Vstring=""
  end
  if L1N==true
    L1Nstring=@sprintf "\tL1N: %0.1f" regN(N,Pars)
  else L1Nstring=""
  end
  if Chi2x==true
    chi2xstring=@sprintf "\tChi2x: %0.3f" Pars.chi2
  else chi2xtring=""
  end
  if Chi2z==true
    chiz = Chi2(Model((Z.vdm),Mats.H),DATA.L,DATA.EL)/(DATA.num_spectra_samples*DATA.num_lines)
    chi2zstring=@sprintf "\tChi2z: %0.3f" chiz
  else chi2zstring=""
  end

	rx=ell2norm(X.vdm)
  if s==true sstring=@sprintf "\ts: %0.3e" ell2norm(X.rho*(X.vdm-X.vdm_previous))/rx
  else sstring=""
  end
#Residuals

  if Pres==true
    PRstring=@sprintf "\trP: %0.1e" ell2norm(P.vdm-X.vdm)/rx							#P residual
  else PRstring=""
  end
  if Zres==true
    ZRstring=@sprintf "\trZ: %0.1e" ell2norm(Z.vdm-X.vdm)/rx							#Z residual
  else ZRstring=""
  end
  if Vres==true
    VRstring=@sprintf "\trV: %0.1e" ell2norm(V.vdm-Z.vdm*Mats.Dv)/rx			#V residual
  else VRstring=""
  end
  if Tres==true
    TRstring=@sprintf "\trT: %0.1e" ell2norm(T.vdm-Mats.Ds*X.vdm)/rx			#T residual
  else TRstring=""
  end
  if Nres==true
    NRstring=@sprintf "\trN: %0.1e" ell2norm(N.vdm-X.vdm)/rx							#N residual
  else NRstring=""
  end
  if Zpen==true
    ZPstring=@sprintf "\tZro: %0.1f" Z.rho                            #Z Penalty
  else ZPstring=""
  end
  if Ppen==true
    PPstring=@sprintf "\tPro: %0.1f" P.rho                            #P Penalty
  else PPstring=""
  end
  if Tpen==true
    TPstring=@sprintf "\tTro: %0.1f" T.rho                            #T Penalty
  else TPstring=""
  end
  if Vpen==true
    VPstring=@sprintf "\tVro: %0.1f" V.rho                            #V Penalty
  else VPstring=""
  end
  if Npen==true
    NPstring=@sprintf "\tNro: %0.1f" N.rho                            #N Penalty
  else NPstring=""
  end
  println(" It: ",Pars.it-2,Jstring,L2xstring,L1Tstring,L1Vstring,L1Nstring,chi2xstring,sstring,PRstring, ZRstring,TRstring,VRstring,NRstring,PPstring,ZPstring,TPstring,VPstring,NPstring,Msg)
end

#=--------------------------------------------------=#
#================== DATA REPORT =====================#
#=--------------------------------------------------=#
#This frunction prints information related to the data files read in.
function data_report(d)
	println("-----------------------------")
	println("--------- DATA INFO ---------")
	println("-----------------------------")
	println("	Wavelengths: ", size(d.wavelength))
	println("	L: ", size(d.L))
	println("	EL: ", size(d.EL))
	println("	numlines: ", d.num_lines)
	println("	spectra_samples: ", d.num_spectra_samples)
	println("	continuum dates: ", size(d.continuum_dates))
	println("-----------------------------")
	println("--------- END  INFO ---------")
	println("-----------------------------")
end

#=--------------------------------------------------=#
#================ Write FITS FILE ===================#
#=--------------------------------------------------=#
#Writes VDM and relevent info to FITS file & Header
function Write_FITS(X,P)
	Date = ["Date",string(Dates.today()),""]
	it = ["it", P.it, "Iterations to convergance"]
	nit = ["max_its", P.nits, "Max number of iterations allowed."]
	mu_smoo = ["mu_smo", P.mu_smoo,"Weight of smoothing regularizer."]
	mu_spec = ["mu_spe", P.mu_spec, "Weight of spectral regularizer."]
	mu_temp = ["mu_tem", P.mu_temp, "Weight of the temporal regularizer."]
	eps_abs = ["eps_abs",P.eps_abs,""]
	eps_rel = ["eps_rel",P.eps_rel,""]
	sigma = ["sigma", P.sigma,""]
	G = ["G",P.G, ""]
	alpha = ["alpha",P.alpha,"Mulitplier Throttle"]
	rho0 = ["in_rho", P.rho0,"Initial Penalty"]
	tau = ["tau", P.tau,""]
	chi2 = ["Chi2", P.chi2,"Final Chi Squared"]
	keys = [Date[1],it[1],nit[1],mu_smoo[1],mu_spec[1],mu_temp[1],eps_abs[1],eps_abs[1],sigma[1],G[1],alpha[1],rho0[1],tau[1],chi2[1]]
	vals = [Date[2],it[2]-1,nit[2],mu_smoo[2],mu_spec[2],mu_temp[2],eps_abs[2],eps_abs[2],sigma[2],G[2],alpha[2],rho0[2],tau[2],chi2[2]]
	coms = [Date[3],it[3],nit[3],mu_smoo[3],mu_spec[3],mu_temp[3],eps_abs[3],eps_abs[3],sigma[3],G[3],alpha[3],rho0[3],tau[3],chi2[3]]
	head = FITSHeader(keys,vals,coms)
	name=string(now(),"vdm.fits")
	f = FITS(name,"w")
	write(f,X.vdm,header=head)
	println("VDM saved as ",name)
end

#=--------------------------------------------------=#
#==================INTERPOLATION=====================#
#=--------------------------------------------------=#
#interp takes x and y data and desired x points
#returns linearally interpolated y desired points
#extrapolates off either end of x data using the
#slope of the first and last sets of data points given.
function interp(x_DESIRED,x_DATA,y_DATA)
  y_DESIRED = zeros(length(x_DESIRED))
	index_DATA = collect(1:length(x_DATA))
  index_DESIRED = collect(1:length(x_DESIRED))
  minx = minimum(x_DATA)
  maxx = maximum(x_DATA)
  for i in index_DESIRED
    yset=0
    y_NEW =0
    for y in index_DATA
      if yset ==0 && x_DESIRED[i] > x_DATA[y] && x_DESIRED[i] < x_DATA[y+1]  y < (length(index_DATA)-1) && i < (length(index_DATA)-1)
        rise = y_DATA[y+1]-y_DATA[y]
        run = x_DATA[y+1]-x_DATA[y]
        slope = rise / run
        y_change = slope*(x_DESIRED[i]-x_DATA[y])
        y_NEW = y_DATA[y]+y_change
        yset = 1
      elseif yset ==0 && x_DESIRED[i] < minx #EXTRAPOLATION DOWNWARD
        #println("EXTRAP LOW!!")
        rise =y_DATA[2]-y_DATA[1]
        run = x_DATA[2]-x_DATA[1]
        first_slope= rise / run
        y_change = first_slope*(x_DATA[1]-x_DESIRED[i])
        y_NEW = y_DATA[1]-y_change
        yset = 1
      elseif yset ==0 && x_DESIRED[i] > maxx #EXTRAPOLATION UPWARD
        #println("EXTRAP HIGH!!")
        rise=y_DATA[length(y_DATA)]-y_DATA[length(y_DATA)-1]
        run = x_DATA[length(x_DATA)]-x_DATA[length(x_DATA)-1]
        last_slope = rise / run
        y_change = last_slope*(x_DESIRED[i]-x_DATA[length(x_DATA)])
        y_NEW = y_DATA[length(y_DATA)]+ y_change
        yset = 1
      end
      y_DESIRED[i] = y_NEW
    end
  end
  y_DESIRED #This is the returned value. In Julia, return statements are not required.
end

#=--------------------------------------------------=#
#====================== Model =======================#
#=--------------------------------------------------=#
#X is the current TDF being modeled.
#ICF is the interpolated continuum flux.
#function Model(X,ICF)
function Model(X,H)
	MF = H*X
  MF #Array Returned!
end

#=--------------------------------------------------=#
#====================== Chi^2 =======================#
#=--------------------------------------------------=#
# M should be the model from the Model function
# D should be the spectral data
# Sigma should be the error on the spectral data
function Chi2(M,D,Sigma)
		sum(   ((M-D)./(Sigma)).^2)
end

#=--------------------------------------------------=#
#============= ELL1 Proximal Operator ===============#
#=--------------------------------------------------=#
#Proximal operator for the L1 norm with positivity.
#X is the current model TDF
#mu is the regularization weight
#rho is the current hyperparameter
function ell1_prox_op(X,XS,mu,rho)
#  for i in collect(1:length(X))
#    if abs(XS[i]) > mu/rho
#			if XS[i] >= mu/rho
#				X[i] = XS[i] - mu/rho
#			else
#				X[i] = XS[i]+mu/rho
#			end
#    else
		#	println("!!!!!!!")
#      X[i] = 0
#    end
#  end
#	X #array returned
	for i in collect(1:length(X))
		Lam = mu/rho
  	if XS[i] >= Lam
			X[i] = XS[i] - Lam
		elseif XS[i] <= -Lam
			X[i] = XS[i] + Lam
		else
			X[i] = 0.0
		end
	end
	X
end

#=--------------------------------------------------=#
#============= ELL2 Proximal Operator ===============#
#=--------------------------------------------------=#
#Proximal operator for the L2 norm.
#X is the current model TDF
#mu is the regularization weight
function ell2_prox_op(X,mu,rho)
	X_returned = 1.0/(1.0+2.0*mu)*X
end

#=--------------------------------------------------=#
#========= Positivity Proximal Operator  ============#
#=--------------------------------------------------=#
function pos_prox_op(X)
	for i in 1:length(X)
		if X[i] < 0.0
			X[i] = 0.0
		end
	end
	X
end

#=--------------------------------------------------=#
#=================== Ell 1 Norm =====================#
#=--------------------------------------------------=#
function ell1norm(X)
	sum(abs(X))
end
#=--------------------------------------------------=#
#=================== Ell 2 Norm =====================#
#=--------------------------------------------------=#
function ell2norm(X)
	sqrt(sum(X.^2))
end
#=--------------------------------------------------=#
#============== Ell 2 Norm Squared ==================#
#=--------------------------------------------------=#
function ell2normsquared(X)
	sum(X.^2)
end

#=--------------------------------------------------=#
#=================== Final Plot =====================#
#=--------------------------------------------------=#
function plotfin(Plot_F,X,Z,T,V)
  if Plot_F == true
  		figure(figsize=(20,10))
  		title("Final Reconstruction")
  		set_cmap("YlOrRd")
  		subplot(221)
  		imshow((Z.vdm),aspect="auto",origin="lower",interpolation="None")
  		title("Z")
  		colorbar()
  		subplot(222)
  		imshow((X.vdm),aspect="auto",origin="lower",interpolation="None")
  		title("X")
  		colorbar()
  		subplot(223)
  		imshow((T.vdm),aspect="auto",origin="lower",interpolation="None")
  		title("T")
  		colorbar()
  		subplot(224)
  		imshow((V.vdm),aspect="auto",origin="lower",interpolation="None")
  		title("V")
  		colorbar()
  		show()
  		#figure()
  		#imshow((N.vdm),aspect="auto",origin="lower",interpolation="None")
  		#colorbar()
  		#show()
  		#writecsv("UnitTests/Ds.csv",Mats.Ds)
  		#writecsv("UnitTests/T.csv",T.vdm)
  		#writecsv("UnitTests/X.csv",T.vdm)
  		#savefig("Spiral_rec.png")
  	end
  end
