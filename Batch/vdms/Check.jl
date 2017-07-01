#actual=readcsv("home/manderson/TLDR/UnitTests/UT_vdm.csv");
actual=readcsv("../../UnitTests/UT_vdm.csv");
mindif=100.0
minf=""
for f in filter(x -> ismatch(r"\.csv", x),readdir())
  solution = readcsv(f)
  absdif=(sum(abs(actual-solution)))
  if absdif < mindif
    mindif=absdif
    minf=f
  end
end
println("Minimum Difference: ", mindif, " File: ", minf)
