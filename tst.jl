


A = [1,2,3,4,5]

W = eye(5)

chiA = A'*W*A
println(chiA)

function chi2(M,D,Sigma)
 sum( ( (vec(M)-vec(D))  ./vec(Sigma) ).^2)   #Value Returned!
end

chiB = chi2([0,0,0,0,0],A,[1.0,1.0,1.0,1.0,1.0])
println(chiB)
