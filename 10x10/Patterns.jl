#using PyPlot

nx=10
ny=10


checkerboard=zeros(nx,ny); n=2; #n sets the size of the checks.
circle=zeros(nx,ny); center=[nx/2,ny/2]; radius=(min(nx,ny)/4);
ring=zeros(nx,ny); center=[nx/2,ny/2]; radius=(min(nx,ny)/4); width=2
lowertri=zeros(nx,ny)
uppertri=ones(nx,ny)
diag=zeros(nx,ny)
diaginv=ones(nx,ny)
for i in collect(1:nx)
  for j in collect(1:ny)
    if i%(2*n)==0 && j%(2*n) ==0 && i!=nx && j!=ny
      for p in collect(0:n-1)
        for q in collect(0:n-1)
          checkerboard[i-p,j-q]=1
        end
      end
    end
    if sqrt((i-center[1])^2+(j-center[2])^2)<=radius
      circle[i,j]=1
    end

    if sqrt((i-center[1])^2+(j-center[2])^2)<=(radius+width) && sqrt((i-center[1])^2+(j-center[2])^2)>=(radius-width)
      ring[i,j]=1
    end

    if i < j
        lowertri[i,j]=1
        uppertri[i,j]=0
    end
    if i == j
        diag[i,j]=1
        diaginv[i,j]=0
    end
  end
end
#figure()
#imshow(ring,interpolation="None",cmap="Greys",origin="Lower")
#colorbar()
#show()
diagrev=copy(diag[end:-1:1,end:-1:1])
diaginvrev=copy(diaginv[end:-1:1,end:-1:1])

box=zeros(nx,ny);box[Int(nx/2-2):Int(nx/2+2),Int(ny/2-2):Int(ny/2+2)]=1.0
ibox=ones(nx,ny);ibox[Int(nx/2-2):Int(nx/2+2),Int(ny/2-2):Int(ny/2+2)]=0.0

vl=zeros(nx,ny);vl[:,1:Int(nx/2)]=1
vr=zeros(nx,ny);vr[:,Int(nx/2):nx]=1
vm=ones(nx,ny);vm[:,Int(nx/2)]=0
vmi=zeros(nx,ny);vmi[:,Int(nx/2)]=1

hm=ones(nx,ny);hm[Int(nx/2),:]=0
hmi=zeros(nx,ny);hmi[Int(nx/2),:]=1
ht=zeros(nx,ny);ht[1:Int(nx/2),:]=1
hb=zeros(nx,ny);hb[Int(nx/2):nx,:]=1

writecsv("checkerboard.csv",checkerboard)
writecsv("circle.csv",circle)
writecsv("ring.csv",ring)
writecsv("lowertri.csv",lowertri)
writecsv("uppertri.csv",uppertri)
writecsv("diagonal.csv",diag)
writecsv("diagonalinverted.csv",diaginv)
writecsv("reverseddiagonal.csv",diagrev)
writecsv("reverseddiagonalinverted.csv",diaginvrev)
writecsv("box.csv",box)
writecsv("invertedbox.csv",ibox)
writecsv("halfleft.csv",vl)
writecsv("halfright.csv",vr)
writecsv("verticalstripe.csv",vm)
writecsv("invertedverticalstripe.csv",vmi)
writecsv("horizontalstripe.csv",hm)
writecsv("invertedhorizontalstripe.csv",hmi)
writecsv("halftop.csv",ht)
writecsv("halfbottom.csv",hb)


#writecsv("checkerboard_20x50/original_checkerboard.csv",checkerboard)
#writecsv("circle_20x50/original_circle.csv",circle)
#writecsv("ring_20x50/original_ring.csv",ring)
