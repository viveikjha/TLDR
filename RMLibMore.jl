module RMLibMore
using RMTypes
using GenMatrices
using PyPlot
using JLD
export J,regX,regT,regV,regN,regP,Report,data_report,Write_FITS,interp,gen_UTD,gen_tiksol,plotfin,Model,Chi2,ell1_prox_op,ell2_prox_op,pos_prox_op,ell1norm,ell2norm,ell2normsquared,save_vars,load_vars,save_data,load_data,load_record
#THIS FILE CONTAINS THE BASE LEVEL FUNCTIONS FOR TLDR THAT ARE NOT PART OF THE CORE ADMM ALGORITHM
#Inputs include all the reconstruction variables. Options are for printed lines.


function save_vars(fname,Mats,Pars)
  save(fname,"Mats",Mats,"Pars",Pars)
end

function load_vars(fname)
  P,M=load(fname)
  Pars=P[2]
  Mats=M[2]
  Pars,Mats
end

function save_data(fname,Data)
	save(fname,"Data",Data)
end

function load_data(fname)
	D=load(fname)
	Data=D["Data"]
end
function load_record(fname)
	D=load(fname)
	Data=D["Record"]
end
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
	#writecsv("Arrays/Xv.csv",X.vdm)
	#A=Pars.mu_smoo*0.5*ell2normsquared(X.vdm)
	G=copy(X.vdm)
	A=Pars.mu_smoo*0.5*	sum(G.^2)
	A
end
#=--------------------------------------------------=#
#========== Regularization term on T,V,N ============#
#=--------------------------------------------------=#
#These functions give the value of the given regularization term.
#The functions are identical, which is silly.
function regT(T,Pars)
	#Pars.mu_temp*ell1norm(T.vdm)
  Pars.mu_temp*vecnorm(T.vdm,1)
end

function regV(V,Pars)
	#Pars.mu_spec*ell1norm(V.vdm)
  Pars.mu_spec*vecnorm(V.vdm,1)
end

function regN(N,Pars)
	#Pars.mu_l1*ell1norm(N.vdm)
  Pars.mu_l1*vecnorm(N.vdm,1)
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
    L2xstring=@sprintf "\tL2x: %0.2f" regX(X,Pars)
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

	#rx=ell2norm(X.vdm)
  rx=vecnorm(X.vdm,2)
  if s==true
			#println("dims of x:", size(X.vdm))
			#println("dims of xp:", size(X.vdm_previous))
			#sstring=@sprintf "\ts: %0.3e" ell2norm(X.rho*(X.vdm-X.vdm_previous))/rx
      sstring=@sprintf "\ts: %0.3e" vecnorm(X.rho*(X.vdm-X.vdm_previous),2)/rx
  else sstring=""
  end
#Residuals

  if Pres==true
    #PRstring=@sprintf "\trP: %0.1e" ell2norm(P.vdm-X.vdm)/rx							#P residual
    PRstring=@sprintf "\trP: %0.1e" vecnorm(P.vdm-X.vdm,2)/rx							#P residual
  else PRstring=""
  end
  if Zres==true
    #ZRstring=@sprintf "\trZ: %0.1e" ell2norm(Z.vdm-X.vdm)/rx							#Z residual
    ZRstring=@sprintf "\trZ: %0.1e" vecnorm(Z.vdm-X.vdm,2)/rx							#Z residual
  else ZRstring=""
  end
  if Vres==true
    #VRstring=@sprintf "\trV: %0.1e" ell2norm(V.vdm-Z.vdm*Mats.Dv)/rx			#V residual
    VRstring=@sprintf "\trV: %0.1e" vecnorm(V.vdm-Z.vdm*Mats.Dv,2)/rx			#V residual
  else VRstring=""
  end
  if Tres==true
    #TRstring=@sprintf "\trT: %0.1e" ell2norm(T.vdm-Mats.Ds*X.vdm)/rx			#T residual
    TRstring=@sprintf "\trT: %0.1e" vecnorm(T.vdm-Mats.Ds*X.vdm,2)/rx			#T residual
  else TRstring=""
  end
  if Nres==true
    #NRstring=@sprintf "\trN: %0.1e" ell2norm(N.vdm-X.vdm)/rx							#N residual
    NRstring=@sprintf "\trN: %0.1e" vecnorm(N.vdm-X.vdm,2)/rx							#N residual
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
	it = ["it", P.it, "Iterations to convergence"]
	conflag = ["ConFlag", P.conflag, "Convergence Flag"]
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
	keys = [Date[1],conflag[1],it[1],nit[1],mu_smoo[1],mu_spec[1],mu_temp[1],eps_abs[1],eps_abs[1],sigma[1],G[1],alpha[1],rho0[1],tau[1],chi2[1]]
	vals = [Date[2],conflag[2],it[2]-1,nit[2],mu_smoo[2],mu_spec[2],mu_temp[2],eps_abs[2],eps_abs[2],sigma[2],G[2],alpha[2],rho0[2],tau[2],chi2[2]]
	coms = [Date[3],conflag[3],it[3],nit[3],mu_smoo[3],mu_spec[3],mu_temp[3],eps_abs[3],eps_abs[3],sigma[3],G[3],alpha[3],rho0[3],tau[3],chi2[3]]
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
      y_DESIRED[i] = copy(y_NEW)
    end
  end
  y_DESIRED #This is the returned value. In Julia, return statements are not required.
end

#=--------------------------------------------------=#
#====================== Model =======================#
#=--------------------------------------------------=#
#X is the current TDF being modeled.
#ICF is the interpolated continuum flux.
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

	#X=(abs(XS)-Lam).*(abs(XS).>=Lam).*((1.0.*(XS.>0.0)).-(1.0.*(XS.<0.0))) #THIS IS SLOWER THAN FOR LOOP.
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
	X=X.*(X.>0.0)
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
			savefig("Rec.png")
  		show()
  		#figure()
  		#imshow((N.vdm),aspect="auto",origin="lower",interpolation="None")
  		#colorbar()
  		#show()
  		#writecsv("UnitTests/Ds.csv",Mats.Ds)
  		#writecsv("UnitTests/T.csv",T.vdm)
  		#writecsv("UnitTests/X.csv",T.vdm)

  	end
  end

	#function gen_tiksol(DATA,Pars,Mats;scale=1.0,mu_smoo=40.0)


	function gen_tiksol(Par,Mats,DATA;scale=1.0,mu_smoo=40.0,plotting=false,save=false)
	  #Par,Mats,DATA=getdata(f)
	  #println("µ: ",mu_smoo)
	  Pars= init_Params()
	  Mats = Gen_Mats(DATA,Pars)

	  #end
	  Pars.mu_smoo=copy(mu_smoo)
	  vdm = zeros(Pars.num_tdf_times,DATA.num_lines)
    #slice_size=size(Mats.W)
	  for l=1:DATA.num_lines        #SPECTAL CHANNEL LOOP
      Q = Mats.HT * Mats.W[:,:,l] * Mats.H +  (Pars.mu_smoo)*Mats.Gammatdf #INCLUCES L1 NORM ON X
      B = Mats.HT* Mats.W[:,:,l] * DATA.L
			G=inv(Q)*B
			vdm[:,l] = G[:,l]
	  end
	  vdm =vdm.*(vdm.>0.0) #FILTER OUT NEGATIVE VALUES
	  if save==true
	    writecsv("devsol.csv",vdm)
	  end
	  if plotting==true
	    figure()
	    imshow(vdm,origin="lower",cmap="Greens",interpolation="None",aspect="auto")
	    xlabel("Spectral Channel")
	    ylabel("Delay Time")
	    titlestring=string("Tikhonov Image for µ=",Pars.mu_smoo, " and scale=",scale)
	    title(titlestring)
	    show()
	  end
	  vdm
	end

	function gen_UTD(DATA,Pars,lvl,nlams,ntimes;CW=6563.0,spread=10)
		lams = collect(1:nlams)+(CW-nlams/2.0)
		vdm = zeros(ntimes,nlams)
		for i in 1:ntimes
			for j in 1:nlams
				if i>= (ntimes/2)-(spread/2) && i < (ntimes/2)+(spread/2) && j> (nlams/2)-(spread/2) && j <= (nlams/2)+(spread/2)
					vdm[i,j]=lvl
				end
			end
		end
		Par = init_Params()
		Mats = Gen_Mats(DATA,Pars)
		Spectra=Mats.H*vdm
		nbase=vec(DATA.EL)
		n=reshape(nbase,size(Spectra))
		Noisy_Spectra = (Spectra+n)
		Error = n
		DATA.L=Noisy_Spectra
		DATA.EL=Error
		DATA.num_lines=size(DATA.L,2)
		DATA.num_spectra_samples=length(DATA.spectra_dates)
		DATA.num_continuum_dates=length(DATA.continuum_dates)
		Mats=Gen_Mats(DATA,Pars)
		vdm
	end

end #END MODULE
