using PyPlot

delays = 50
waves = 100
synthdata=zeros((delays,waves))

for i in 1:delays
	for j in 1:waves
		synthdata[i,j] = sin( (sin(j/100*pi)*i)/50*pi)
	end
end
#println(synthdata) 

PyPlot.figure(1)
PyPlot.imshow(synthdata)
PyPlot.show()



