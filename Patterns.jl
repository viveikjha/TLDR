using PyPlot

nx=50
ny=20


checkerboard=zeros(nx,ny)
circle=zeros(nx,ny); center=[nx/2,ny/2]; radius=(min(nx,ny)/4);
for i in collect(1:nx)
  for j in collect(1:ny)
    if i%2==0 && j %2 ==0 && i != nx && j!=ny
      checkerboard[i,j]=1.0
    end
    if sqrt((i-center[1])^2+(j-center[2])^2)<=radius
      circle[i,j]=1.0
    end
  end
end
figure()
imshow(circle,interpolation="None",cmap="Greys",origin="Lower")
show()
