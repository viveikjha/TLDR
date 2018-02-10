#actual=readcsv("home/matander/TLDR/UnitTests/UT_vdm.csv");
actual=readcsv("../../UnitTests/UT_vdm.csv");
#actual=readcsv("../../Paper/ring/simulated_ring_vdm.csv");

mindif=100.0
BPSNR=0.0
minf=""
minp=""
for f in filter(x -> ismatch(r"\.csv", x),readdir())
  solution = readcsv(f)
  absdif=(sum(abs(actual-solution))) #abs() or abs.() will depend on Julia version.
  MSE=(1.0/length(f))*sum((actual-solution).^2)
  MAX=255.0
  PSNR=10.0*log10(MAX^2/MSE)
  if absdif < mindif
    mindif=absdif
    minf=f
  end
  if PSNR > BPSNR
    BPSNR=PSNR
    minp=f
  end

end
println("Minimum Difference: ", mindif, " File: ", minf)
println("Best PSNR: ", BPSNR, " dB File: ", minp)
