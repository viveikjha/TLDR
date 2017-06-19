#using PyPlot



vdm = readcsv("vdm.csv")
lam = readcsv("data/rvm_wavelengths.csv")


nlines =size(lam)

lam_target = 4861.3
vmin = -4000.0
c=299792.458
lammin = (vmin/c*lam_target)+lam_target

vmax = 4000.0
lammax = (vmax/c*lam_target)+lam_target
target_index = indmin(abs(lam-lam_target))

println(lam>lammin)