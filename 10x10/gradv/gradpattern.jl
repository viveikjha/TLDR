using PyPlot


x=zeros(10,10)
x2=zeros(10,10)
nps=10
lvl=1.0
grad=linspace(0.0,lvl,nps)
for i in 1:10
    x[:,i]=grad
    x2[i,:]=grad
end
#println(x)

writecsv("gradv.csv",x)
writecsv("gradh.csv",x2)
