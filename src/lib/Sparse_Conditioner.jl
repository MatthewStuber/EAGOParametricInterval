"""
     SparseInSto

Storage type for in-place sparse calculations
"""
type SparseInSto
     Xh
     Xj
     Yh
     Yj
     Zh
     Zj
     nx
end
SparseInSto() = SparseInSto([],[],[],[],[],[],[])

"""
     Preconditioner(h,X,P;jac="User")

Directly applies inverse preconditioning matrix.
"""
function Preconditioner(h,X,P;jac="User")
  J = h(X,P)
  #println("J:   ",J)
  if (length(X)>1)
    #println("mid.(J)", mid.(J))
    Y = inv(mid.(J))
  else
    Y = [1.0/(mid(J[1]))]
  end
  return Y
end

"""
     Sparse_Forward_Elimination!(x::Vector{Tq},L::SparseMatrixCSC{Tv,Ti},
                                 b::Vector{Tq},nx::Int64) where {Tv,Tq,Ti}

Solves for `Lx=b` via forward elimination. A must be a lower triangular sparse
matrix of CSC format.
"""
function Sparse_Forward_Elimination!(x::Vector{Tq},L::SparseMatrixCSC{Tv,Ti},
                                     b::Vector{Tq},nx::Int64) where {Tv,Tq,Ti}
     # converts to CSR (expensive if dense ... much cheaper calc if sparse)
     A::SparseMatrixCSC{Tv,Ti} = L'
     # standard row-oriented forward elimination
     x[1,:] = b[1,:]/A.nzval[1]
     for i=2:(nx)
          for k=(A.colptr[i]):(A.colptr[i+1]-2)
               b[i,:] = b[i,:] - A.nzval[k]*x[A.rowval[k],:]
          end
          x[i,:] = b[i,:]/A.nzval[A.colptr[i+1]-1]
     end
end

"""
     Sparse_Forward_Elimination!(x::SparseMatrixCSC{Tq,Ti},U::SparseMatrixCSC{Tv,Ti},
                                  b::SparseMatrixCSC{Tq,Ti},nx::Int64)

Solves for `Lx=b` via forward elimination. A must be a lower triangular sparse
matrix of CSC format.
"""
# TO DO: Improve row access for b
function Sparse_Forward_Elimination!(x::SparseMatrixCSC{Tq,Ti},L::SparseMatrixCSC{Tv,Ti},
                                     b::SparseMatrixCSC{Tq,Ti},nx::Int64) where {Tv,Tq,Ti}
     # converts to CSR (expensive if dense ... much cheaper calc if sparse)
     A::SparseMatrixCSC{Tv,Ti} = L'
     # standard row-oriented forward elimination
     x[1,:] = b[1,:]/A.nzval[1]
     for i=2:(nx)
          for k=(A.colptr[i]):(A.colptr[i+1]-2)
               b[i,:] = b[i,:] - A.nzval[k]*x[A.rowval[k],:]
          end
          x[i,:] = b[i,:]/A.nzval[A.colptr[i+1]-1]
     end
end

"""
     Sparse_Back_Substitution!(x::Vector{Tq},U::SparseMatrixCSC{Tv,Ti},
                              b::Vector{Tq},nx::Int64)

Solves for `Ux=b` via backsubstitution. A must be a upper triangular sparse
matrix of CSC format.
"""
function Sparse_Back_Substitution!(x::Vector{Tq},U::SparseMatrixCSC{Tv,Ti},
                                   b::Vector{Tq},nx::Int64) where {Tv,Tq,Ti}
     # converts to CSR (expensive if dense ... much cheaper calc if sparse)
     A::SparseMatrixCSC{Tv,Ti} = U'
     # standard row-oriented back substituion
     x[end,:] = b[end,:]/A.nzval[end]
     for i=(nx-1):-1:1
          for k=(A.colptr[i+1]-1):-1:(A.colptr[i]+1)
               b[i,:] = b[i,:] - A.nzval[k]*x[A.rowval[k],:]
          end
          x[i,:] = b[i,:]/A.nzval[A.colptr[i]]
     end
end

"""
     Sparse_Back_Substitution!(x::SparseMatrixCSC{Tq,Ti},U::SparseMatrixCSC{Tv,Ti},
                               b::SparseMatrixCSC{Tq,Ti},nx::Int64)

Solves for `Ux=b` via backsubstitution. A must be a upper triangular sparse
matrix of CSC format.
"""
# TO DO: Improve row access for b
function Sparse_Back_Substitution!(x::SparseMatrixCSC{Tq,Ti},U::SparseMatrixCSC{Tv,Ti},
                                   b::SparseMatrixCSC{Tq,Ti},nx::Int64) where {Tv,Tq,Ti}
     # converts to CSR (expensive if dense ... much cheaper calc if sparse)
     A::SparseMatrixCSC{Tv,Ti} = U'
     # standard row-oriented back substituion
     x[end,:] = b[end,:]/A.nzval[end]
     for i=(nx-1):-1:1
          for k=(A.colptr[i+1]-1):-1:(A.colptr[i]+1)
               b[i,:] = b[i,:] - A.nzval[k]*x[A.rowval[k],:]
          end
          x[i,:] = b[i,:]/A.nzval[A.colptr[i]]
     end
end

"""
     Sparse_Precondition!(H::Vector{Ta},J::SparseMatrixCSC{Ta,Ti},
                          P::SparseMatrixCSC{Tp,Ti},st::SparseInSto)

Preconditions the H & J to inv(P)H and inv(P)J using a sparse LU factorization
method with full pivoting. J and P must be of size nx-by-nx and H must be of
size nx. st is the inplace storage type.
"""
function Sparse_Precondition!(H::Vector{Ta},J::SparseMatrixCSC{Ta,Ti},
                              P::SparseMatrixCSC{Tp,Ti},st::SparseInSto) where {Ta,Tp,Ti}

     # generate LU-PDQ factorization
     lu = lufact(P)

     # solves Lz = PDH for z
     Sparse_Forward_Elimination!(st.Zh,lu[:L],(lu[:Rs].*H)[lu[:p]],st.nx)
     Sparse_Forward_Elimination!(st.Zj,lu[:L],(lu[:Rs].*J)[lu[:p],:],st.nx)

     # solves Uy = z for y
     Sparse_Back_Substitution!(st.Yh,lu[:U],st.Zh,st.nx)
     Sparse_Back_Substitution!(st.Yj,lu[:U],st.Zj,st.nx)

     # solves x = Qy
     st.Xh = st.Yh[lu[:q]]
     st.Xj = st.Yj[lu[:q],:]

     # stores the preconditioned matrices back in place
     H[:] = st.Xh[:,1]
     J[:] = st.Xj[:,1:(st.nx)]
end

"""
     Dense_Precondition!(H::Vector{Ta},J::Array{Ta,2},P::Array{Tp,2})

Preconditions the H & J to inv(P)H and inv(P)J using a sparse LU factorization
method with full pivoting. J and P must be of size nx-by-nx and H must be of
size nx. st is the inplace storage type.
"""

function LowMatMult!(a::Array{U,2},b::Array{T,2},m::Int64,n::Int64) where {T,U<:AbstractFloat}
    @inbounds for i=1:m
        @inbounds for j = 1:n
            if (b[j,i] != zero(T))
                temp::T = b[j,i]
                @inbounds for k=(j+1):(n)
                    b[k,i] = b[k,i]-A[k,j]*temp
                end
            end
        end
    end
end

function UppMatMult!(a::Array{U,2},b::Array{T,2},m::Int64,n::Int64) where {T,U<:AbstractFloat}
     @inbounds for i = 1:m
          @inbounds for j = n:-1:1
               if (b[j,i] != zero(T))
                    temp::T = b[j,i]
                    b[j,i] = temp/A[j,j]
                    @inbounds for k=1:(j-1)
                    b[k,i] = b[k,i]-A[k,j]*temp
                    end
               end
          end
     end
end

function Dense_Precondition!(H::Vector{S},J::Array{S,2},Y::Array{U,2},nx::Int64) where {S,U<:AbstractFloat}
         lu = lufact(Y)
         HJ::Array{S,2} = [H, J][lu[:p],:]
         LowMatMult!(lu[:L],HJ,nx+1,nx)
         UppMatMult!(lu[:U],HJ,nx+1,nx)
         H[:] = HJ[:,1]
         J[:] = HJ[:,2:(nx+1)]
end
