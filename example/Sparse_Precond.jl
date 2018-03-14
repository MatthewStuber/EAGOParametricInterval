using EAGOParametricInterval

nx = 5
A = sprand(Float64,nx,nx,0.9)
c = rand(Float64,nx)
P = sprand(Float64,nx,nx,0.9)
lufct = lufact(P)
At = copy(A)
ct = copy(c)
Pt = copy(P)
preA = inv(lufct)*A
prec = inv(lufct)*c
SSto = SparseInSto()
SSto.Xh = zeros(Float64,nx)
SSto.Yh = zeros(Float64,nx)
SSto.Zh = zeros(Float64,nx)
SSto.Xj = spzeros(Float64,nx,nx)
SSto.Yj = spzeros(Float64,nx,nx)
SSto.Zj = spzeros(Float64,nx,nx)
SSto.nx = nx
Sparse_Precondition!(c,A,P,SSto)
