using PyPlot


vdm = readcsv("simulated_vdm.csv")


figure()
imshow(log(vdm),cmap="Reds",aspect="auto",origin="lower")
show()
