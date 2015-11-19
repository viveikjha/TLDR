using PyPlot

cits = readcsv("mu_test/cits.csv")
chis = readcsv("mu_test/chis.csv")
dits = readcsv("mu_test/dits.csv")
diffs = readcsv("mu_test/diffs.csv")
mu = readcsv("mu_test/mu_spec.csv")






figure()
plot(cits,chis,label="Chi")
plot(dits,diffs, label="Dif")
xlim(15,25)
xlabel("iterations")
ylabel("value")
legend()
#show()

#figure()
#plot(cits,mu,"r.")
#plot(dits,mu,"b.")
#title("iterations to best")
#xlabel("Iterations")
#ylabel("Mu")
#xlim(10,25)
#show()

