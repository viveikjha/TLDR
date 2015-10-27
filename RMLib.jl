#=--------------------------------------------------=#
#====================== CLEAR =======================#
#=--------------------------------------------------=#
function clear()
  for i in collect(1:30)
    println("\n")
  end
end



function print_rgb(r, g, b, t)
           println("\e[1m\e[38;2;$r;$g;$b;249m",t)
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
#  index_DATA = linspace(1,length(x_DATA)-1,length(x_DATA)-1)
	index_DATA = collect(1:length(x_DATA))
  index_DESIRED = collect(1:length(x_DESIRED))
  #println([1:length(x_DESIRED)])
  #println(index_DESIRED)
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
#============== Chi Squared Gradient ================#
#=--------------------------------------------------=#
#Z is the current TDF
#D is the spectral data
#Sigma is the error on D
#Rho is the regularization weight
#ICF is the interpolated continuum flux
function chigrad(Z,D,Sigma,rho,ICF)
  M = Model(Z,ICF)
  chi2grad = zeros(size(Z))
  for i in collect(1:size(Z)[1])
    #chi2grad[i] = 2.0*sum((M-D)/Sigma^2. .* ICF[:,i])
    chi2grad[i] = sum(2.0*sum(D - sum(vec(Z).*vec(ICF[i,:])))-sum(ICF,2))
		chi2grad[i] = sum(2.0*sum(D - sum(Z.* H))-sum(ICF))
  end
  chi2grad #Array Returned!
end




#=--------------------------------------------------=#
#====================== Model =======================#
#=--------------------------------------------------=#
#X is the current TDF being modeled.
#ICF is the interpolated continuum flux.
#function Model(X,ICF)
function Model(X,H)
  #L = size(ICF)[1]
  #MF = zeros(L)
  #MF = H*vec(X)
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
	  #sum(   (vec(M)-vec(D)).^2  ./vec(Sigma) .^2)   #Value Returned!
		sum(   (M-D).^2 ./(Sigma).^2)
		#sum(((vec(M)-vec(D))/vec(Sigma))^2)
end

#=--------------------------------------------------=#
#============= ELL1 Proximal Operator ===============#
#=--------------------------------------------------=#
#Proximal operator for the L1 norm with positivity.
#X is the current model TDF
#mu is the regularization weight
#rho is the current hyperparameter
function ell1_prox_op(X,mu,rho)
  #if mu/rho == 0.0
  #  println("Error in Proximal Operator: mu/rho = 0")
  #elseif mu/rho == inf || mu/rho == inf
  #  println("Error in Proximal Operator: mu/rho = infinity")
  #end
  for i in collect(1:length(X))
    if X[i] > mu/rho
      X[i] = X[i] - mu/rho
    else
      X[i] = 0
    end
  end
  X #Array Returned!
end

#=--------------------------------------------------=#
#============= ELL2 Proximal Operator ===============#
#=--------------------------------------------------=#
#Proximal operator for the L2 norm.
#X is the current model TDF
#mu is the regularization weight
function ell2_prox_op(X,mu)
	X_returned = 1.0/(1.0*2.0*mu).*X
end

#=--------------------------------------------------=#
#========= Positivity Proximal Operator  ============#
#=--------------------------------------------------=#
function pos_prox_op(X)
	for i in collect(1:size(X)[1])
		for j in collect(1:size(X)[2])
			if X[i,j] < 0.0
				X[i,j] = 0.0
			end
		end	
	end
	X
end

#=--------------------------------------------------=#
#=================== Ell 2 Norm =====================#
#=--------------------------------------------------=#
function ell2norm(X)
	#sqrt(sum(vec(X).^2))   #Array Returned!
	sqrt(sum(X).^2)
end
