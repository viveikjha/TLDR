using PyPlot


x=zeros(10,10)
nps=10
lvl=1.0
grad=linspace(0.0,lvl,5)
gradrev=linspace(lvl,0.0,5)
for i in 1:10
    if i <=nps/2
        x[:,i]=grad[i]
    elseif i>nps/2
        x[:,i]=gradrev[i-5]
    end
end
#println(x)

writecsv("gradsv.csv",x)
