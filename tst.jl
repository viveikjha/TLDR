using FITSIO




x = readcsv("vdm.csv")



Keys = ["Date"]
vals = [string(Dates.today())]
coms = ["date produced"]

head = FITSHeader(Keys,vals,coms)


#println(get_comment(head,"Comment"))

f = FITS("vdm.fits","w")
write(f,x, header=head)
