using PyPlot

nx=50
ny=20


checkerboard=zeros(nx,ny); n=2; #n sets the size of the checks.
circle=zeros(nx,ny); center=[nx/2,ny/2]; radius=(min(nx,ny)/4);
ring=zeros(nx,ny); center=[nx/2,ny/2]; radius=(min(nx,ny)/4); width=3
for i in collect(1:nx)
  for j in collect(1:ny)
    if i%(2*n)==0 && j%(2*n) ==0 && i!=nx && j!=ny
      for p in collect(0:n-1)
        for q in collect(0:n-1)
          checkerboard[i-p,j-q]=1.0
        end
      end
    end
    if sqrt((i-center[1])^2+(j-center[2])^2)<=radius
      circle[i,j]=1.0
    end

    if sqrt((i-center[1])^2+(j-center[2])^2)<=(radius+width) && sqrt((i-center[1])^2+(j-center[2])^2)>=(radius-width)
      ring[i,j]=1.0
    end

  end
end
figure()
imshow(ring,interpolation="None",cmap="Greys",origin="Lower")
colorbar()
show()
