"""
    PI_NewtonGS

Uses the interval arithmetic presented in ValidatedNumerics to bound unique
implicit functions x:P->X defined by a system of equations h(z,p)=0 via
Gauss-Siedel Newton. Inputs:
* `X0::Vector{Interval{T}}`: Bounds for dependent variables
* `P::Vector{Interval{T}}`: Bounds for independent variables
* `hj::Function`: Jacobian of h(x,p) with respect to x
* `h::Function`: Equations defining potential implicit function
* `opt::Vector{Any}`: `opt[1]::Int64` is the number of iterations, `opt[2]::Float64`
                       is the equality tolerance for intervals, `opt[3]::Float64`
                       is the amount added to M[i,i] when processing extended
                       interval division.
* `Eflag::Bool`: True if no implicit function exists in box
* `Iflag::Bool`: True if unique implicit function exists in box
* `eDflag::Bool`: True if extended division occured in iteration
Returns a tuple `(X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh)` where
* `X::Vector{Interval{T}}`: Output interval box
* `Xtemp::Vector{Interval{T}}`: Extra box formed from extended division
* `Eflag::Bool`: True if no implicit function exists in box
* `Iflag::Bool`: True if unique implicit function exists in box
* `inclusionLow::Vector{Bool}`: Each component is true if corresponding variable
                                contracted from lower bound
* inclusionHigh::Vector{Bool}`: Each component is true if corresponding variable
                                contracted from upper bound
"""
function PI_NewtonGS(X0::Vector{Interval{T}},P::Vector{Interval{T}},
                     hj::Function,h::Function,
                     opt::Vector{Any},Eflag::Bool,
                     Iflag::Bool,eDflag::Bool) where {T<:AbstractFloat}
  # unpacks option file
  kmax::Int64 = opt[1]
  etol::Float64 = opt[2]
  rtol::Float64 = opt[3]

  # Initializes variables
  nx::Int64 = length(X0)
  S1::Interval{T} = Interval(0.0)
  S2::Interval{T} = Interval(0.0)
  inclusion::Vector{Bool} = [false for i=1:nx]
  inclusionHigh::Vector{Bool} = [false for i=1:nx]
  inclusionLow::Vector{Bool} = [false for i=1:nx]
  eD::Int64 = 0
  exclusion::Bool = false
  Iflag::Bool = false
  Eflag::Bool = false
  eDflag::Bool = false

  X::Vector{Interval{T}} = copy(X0)
  Xi::Vector{Interval{T}} = copy(X0)
  N::Vector{Interval{T}} = copy(X)
  x_mid::Vector{T} = mid.(X)
  k::Int64 = 1
  H::Union{Array{Interval{T},2},Vector{Interval{T}}} = h(x_mid,P)
  J::Union{Array{Interval{T},2},Vector{Interval{T}}} = hj(X,P)
  Y::Union{Array{T,2},Vector{T}} = Preconditioner(hj,X,P,jac="User")
  if (nx == 1)
    B::Vector{Interval{T}} = Y.*H
    M::Union{Array{Interval{T},2},Vector{Interval{T}}} = eye(nx)-Y.*J
  else
    B = Y*H
    M = eye(nx)-Y*J
  end

  for i=1:nx
    S1 = Interval(0.0)
    S2 = Interval(0.0)
    for j=1:nx
      if (j<i)
        S1 += M[i,j]*(X[j]-x_mid[j])
      elseif (j>i)
        S2 += M[i,j]*(X[j]-x_mid[j])
      end
    end
    if M[i,i].lo*M[i,i].hi > 0.0
      N[i] = x_mid[i] - (B[i]+S1+S2)/M[i,i]
    else
      Ntemp = N
      eD,N[i],Ntemp[i] = extProcess(N[i],X[i],M[i,i],S1,S2,B[i],rtol)
      if eD == 1
        eDflag = true
        Xtemp = X
        Xtemp[i] = Ntemp[i] ∩ X[i]
        X[i] = N[i] ∩ X[i]
        return X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh
      end
    end
    if Strict_XinY(N[i],X[i])
      inclusion[i] = true
      inclusionHigh[i] = true
      inclusionLow[i] = true
    else
      inclusion[i] = false
      inclusionLow[i] = false
      inclusionHigh[i] = false
      if (N[i].lo>X[i].lo)
        inclusionLow[i] = true
      elseif (N[i].hi<X[i].hi)
        inclusionHigh[i] = true
      end
    end
    if ~isdisjoint(N[i],X[i])
      X[i] = N[i] ∩ X[i]
    else
      exclusion = true
      break
    end
  end
  k += 1

  # checks if all components are included
  if (Iflag == false)
    for i=1:nx
      if (inclusion[i] == true)
        Iflag = true
        continue
      elseif (inclusion[i] == false)
        Iflag = false
        break
      end
    end
  end

  #
  while ((k<kmax) && isEqual(X,Xi,etol) == false)

    Xi = copy(X0)
    x_mid = mid.(X)
    H = h(x_mid,P)
    J = hj(X,P)
    Y = Preconditioner(hj,X,P,jac="User")
    if (nx == 1)
      B = Y.*H
      M = eye(nx)-Y.*J
    else
      B = Y*H
      M = eye(nx)-Y*J
    end

    for i=1:nx
      S1 = Interval(0.0)
      S2 = Interval(0.0)
      for j=1:nx
        if (j<i)
          S1 += M[i,j]*(X[j]-x_mid[j])
        elseif (j>i)
          S2 += M[i,j]*(X[j]-x_mid[j])
        end
      end
      if M[i,i].lo*M[i,i].hi > 0.0
        N[i] = x_mid[i] - (B[i]+S1+S2)/M[i,i]
      else
        Ntemp = N
        eD,N[i],Ntemp[i] = extProcess(N[i],X[i],M[i,i],S1,S2,B[i],rtol)
        if eD == 1
          eDflag = true
          Xtemp = X
          Xtemp[i] = Ntemp[i] ∩ X[i]
          X[i] = N[i] ∩ X[i]
          return X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh
        end
      end
      if Strict_XinY(N[i],X[i])
        inclusion[i] = true
        inclusionHigh[i] = true
        inclusionLow[i] = true
      else
        inclusion[i] = false
        inclusionLow[i] = false
        inclusionHigh[i] = false
        if (N[i].lo>X[i].lo)
          inclusionLow[i] = true
        elseif (N[i].hi<X[i].hi)
          inclusionHigh[i] = true
        end
      end
      if ~isdisjoint(N[i],X[i])
        X[i] = N[i] ∩ X[i]
      else
        exclusion = true
        break
      end
    end
    k += 1

    # checks if all components are included
    if (Iflag == false)
      for i=1:nx
        if (inclusion[i] == true)
          Iflag = true
          continue
        elseif (inclusion[i] == false)
          Iflag = false
          break
        end
      end
    end
  end

  if exclusion
    Eflag = true
  end
  Xtemp = copy(X)
  return X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh
end

"""
    PI_KrawczykCW

Uses the interval arithmetic presented in ValidatedNumerics to bound unique
implicit functions x:P->X defined by a system of equations h(z,p)=0 via
a component-wise Krawczyk method. Inputs:
* `X0::Vector{Interval{T}}`: Bounds for dependent variables
* `P::Vector{Interval{T}}`: Bounds for independent variables
* `hj::Function`: Jacobian of h(x,p) with respect to x
* `h::Function`: Equations defining potential implicit function
* `opt::Vector{Any}`: `opt[1]::Int64` is the number of iterations, `opt[2]::Float64`
                       is the equality tolerance for intervals, `opt[3]::Float64`
                       is the amount added to M[i,i] when processing extended
                       interval division.
* `Eflag::Bool`: True if no implicit function exists in box
* `Iflag::Bool`: True if unique implicit function exists in box
Returns a tuple `(X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh)` where
* `X::Vector{Interval{T}}`: Output interval box
* `Xtemp::Vector{Interval{T}}`: Extra box formed from extended division
* `Eflag::Bool`: True if no implicit function exists in box
* `Iflag::Bool`: True if unique implicit function exists in box
* `inclusionLow::Vector{Bool}`: Each component is true if corresponding variable
                                contracted from lower bound
* inclusionHigh::Vector{Bool}`: Each component is true if corresponding variable
                                contracted from upper bound
"""
function PI_KrawczykCW(X0::Vector{Interval{T}},P::Vector{Interval{T}},
                     hj::Function,h::Function,
                     opt::Vector{Any},Eflag::Bool,Iflag::Bool) where {T<:AbstractFloat}
  # unpacks option file
  kmax::Int64 = opt[1]
  etol::Float64 = opt[2]

  # Initializes variables
  nx::Int64 = length(X0)
  S1::Interval{T} = Interval(0.0)
  S2::Interval{T} = Interval(0.0)
  inclusion::Vector{Bool} = [false for i=1:nx]
  inclusionHigh::Vector{Bool} = [false for i=1:nx]
  inclusionLow::Vector{Bool} = [false for i=1:nx]
  exclusion::Bool = false
  Iflag::Bool = false
  Eflag::Bool = false
  eDflag::Bool = false
  X::Vector{Interval{T}} = copy(X0)
  Xi::Vector{Interval{T}} = copy(X)
  N::Vector{Interval{T}} = copy(X)
  x_mid::Vector{T} = mid.(X)
  k::Int64 = 1

  H::Union{Array{Interval{T},2},Vector{Interval{T}}} = h(x_mid,P)
  J::Union{Array{Interval{T},2},Vector{Interval{T}}} = hj(X,P)
  Y::Union{Array{T,2},Vector{T}} = Preconditioner(hj,X,P,jac="User")
  if (nx == 1)
    println("ran to me 1")
    B::Vector{Interval{T}} = Y.*H
    M::Union{Array{Interval{T},2},Vector{Interval{T}}} = eye(nx)-Y.*J
  else
    println("ran to me 2")
    B = Y*H
    M = eye(nx)-Y*J
  end
  println("ran to me 3")

  for i=1:nx
    S1 = Interval(0.0)
    S2 = Interval(0.0)
    for j=1:nx
      if (j<i)
        S1 += M[i,j]*(X[j]-x_mid[j])
      elseif (j>=i)
        S2 += M[i,j]*(X[j]-x_mid[j])
      end
    end
    N[i] = x_mid[i] - B[i] + S1 + S2

    if Strict_XinY(N[i],X[i])
      inclusion[i] = true
      inclusionHigh[i] = true
      inclusionLow[i] = true
    else
      inclusion[i] = false
      inclusionLow[i] = false
      inclusionHigh[i] = false
      if (N[i].lo>X[i].lo)
        inclusionLow[i] = true
      elseif (N[i].hi<X[i].hi)
        inclusionHigh[i] = true
      end
    end
    if ~isdisjoint(N[i],X[i])
      X[i] = N[i] ∩ X[i]
    else
      exclusion = true
      break
    end
  end
  k += 1

  if (Iflag == false)
    for i=1:nx
      if (inclusion[i] == true)
        Iflag = true
        continue
      elseif (inclusion[i] == false)
        Iflag = false
        break
      end
    end
  end

  #
  while ((k<kmax) && (isEqual(X,Xi,etol) == false))
    Xi = X0
    x_mid = mid.(X)
    H = h(x_mid,P)
    J = hj(X,P)
    Y = Preconditioner(hj,X,P,jac="User")
    if (nx == 1)
      B = Y.*H
      M = eye(nx)-Y.*J
    else
      B = Y*H
      M = eye(nx)-Y*J
    end

    for i=1:nx
      S1 = Interval(0.0)
      S2 = Interval(0.0)
      for j=1:nx
        if (j<i)
          S1 += M[i,j]*(X[j]-x_mid[j])
        elseif (j>=i)
          S2 += M[i,j]*(X[j]-x_mid[j])
        end
      end
      N[i] = x_mid[i] - B[i]+S1+S2
      if Strict_XinY(N[i],X[i])
        inclusion[i] = true
        inclusionHigh[i] = true
        inclusionLow[i] = true
      else
        inclusion[i] = false
        inclusionLow[i] = false
        inclusionHigh[i] = false
        if (N[i].lo>X[i].lo)
          inclusionLow[i] = true
        elseif (N[i].hi<X[i].hi)
          inclusionHigh[i] = true
        end
      end
      if ~isdisjoint(N[i],X[i])
        X[i] = N[i] ∩ X[i]
      else
        exclusion = true
        break
      end
    end
    k += 1

    # checks if all components are included
    if (Iflag == false)
      for i=1:nx
        if (inclusion[i] == true)
          Iflag = true
          continue
        elseif (inclusion[i] == false)
          Iflag = false
          break
        end
      end
    end
  end

  if exclusion
    Eflag = true
  end
  return X,Eflag,Iflag,inclusionLow,inclusionHigh
end


"""
    PIn_NewtonGS(X0::Vector{Interval{T}},P::Vector{Interval{T}},hj!::Function,
                 h!::Function,opt::Vector{Any},Eflag::Bool,Iflag::Bool,eDflag::Bool)

Uses the interval arithmetic presented in ValidatedNumerics to bound unique
implicit functions x:P->X defined by a system of equations h(z,p)=0 via
Gauss-Siedel Newton. The jacobian of h w.r.t. z is calculated in placeInputs:
* `X0::Vector{Interval{T}}`: Bounds for dependent variables
* `P::Vector{Interval{T}}`: Bounds for independent variables
* `hj!::Function`: Jacobian of h(x,p) with respect to x (takes input J,x,p)
* `h!::Function`: Equations defining potential implicit function (takes inputs H,x,p)
* `opt::Vector{Any}`: `opt[1]::Int64` is the number of iterations, `opt[2]::Float64`
                       is the equality tolerance for intervals, `opt[3]::Float64`
                       is the amount added to M[i,i] when processing extended
                       interval division.
* `Eflag::Bool`: True if no implicit function exists in box
* `Iflag::Bool`: True if unique implicit function exists in box
* `eDflag::Bool`: True if extended division occured in iteration
Returns a tuple `(X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh)` where
* `X::Vector{Interval{T}}`: Output interval box
* `Xtemp::Vector{Interval{T}}`: Extra box formed from extended division
* `Eflag::Bool`: True if no implicit function exists in box
* `Iflag::Bool`: True if unique implicit function exists in box
* `inclusionLow::Vector{Bool}`: Each component is true if corresponding variable
                                contracted from lower bound
* inclusionHigh::Vector{Bool}`: Each component is true if corresponding variable
                                contracted from upper bound
"""
function PIn_NewtonGS(X0::Vector{Interval{T}},P::Vector{Interval{T}},
                     hj!::Function,h!::Function,
                     opt::Vector{Any},Eflag::Bool,
                     Iflag::Bool,eDflag::Bool) where {T<:AbstractFloat}
  # unpacks option file
  kmax::Int64 = opt[1]
  etol::Float64 = opt[2]
  rtol::Float64 = opt[3]

  # Initializes variables
  nx::Int64 = length(X0)
  S1::Interval{T} = Interval(0.0)
  S2::Interval{T} = Interval(0.0)
  S3::Interval{T} = Interval(0.0)
  inclusion::Vector{Bool} = [false for i=1:nx]
  inclusionHigh::Vector{Bool} = [false for i=1:nx]
  inclusionLow::Vector{Bool} = [false for i=1:nx]
  eD::Int64 = 0
  exclusion::Bool = false
  Iflag::Bool = false
  Eflag::Bool = false
  eDflag::Bool = false
  k::Int64 = 1
  X::Vector{Interval{T}} = copy(X0)
  Xi::Vector{Interval{T}} = copy(X0)
  N::Vector{Interval{T}} = copy(X)

  # creates storage for h(xm,P), J_x(X,P), and preconditioner Y
  H::Vector{Interval{T}} = [Interval(0.0) for i=1:nx]
  J::SparseMatrixCSC{Interval{T},Int64} = spzeros(Interval{T},nx,nx)

  # initalizes storage object for in-place sparse calculations
  SSto = SparseInSto()
  SSto.Xh = copy(H)
  SSto.Yh = copy(H)
  SSto.Zh = copy(H)
  SSto.Xj = spzeros(Interval{T},nx,nx)
  SSto.Yj = spzeros(Interval{T},nx,nx)
  SSto.Zj = spzeros(Interval{T},nx,nx)
  SSto.nx = copy(nx)

  # calculates extension of h(xm,P) and J_x(X,P)
  x_mid::Vector{T} = mid.(X)
  h!(H,x_mid,P)
  hj!(J,X,P)
  Sparse_Precondition!(H,J,mid.(J),SSto)
  Mt = transpose(J)
  for i=1:nx
    S1 = Interval(0.0)
    S2 = Interval(0.0)
    S3 = Interval(0.0)
    for q=(Mt.colptr[i]):(Mt.colptr[i+1]-1)
      if (i < Mt.rowval[q])
        S1 += Mt.nzval[q]*(X[Mt.rowval[q]]-x_mid[Mt.rowval[q]])
      elseif (i > Mt.rowval[q])
        S2 += Mt.nzval[q]*(X[Mt.rowval[q]]-x_mid[Mt.rowval[q]])
      else
        S3 = Mt.nzval[q]
      end
    end
    if S3.lo*S3.hi > 0.0
      N[i] = x_mid[i] - (H[i]+S1+S2)/S3
    else
      Ntemp = N
      eD,N[i],Ntemp[i] = extProcess(N[i],X[i],S3,S1,S2,H[i],rtol)
      if eD == 1
        eDflag = true
        Xtemp = X
        Xtemp[i] = Ntemp[i] ∩ X[i]
        X[i] = N[i] ∩ X[i]
        return X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh
      end
    end
    if Strict_XinY(N[i],X[i])
      inclusion[i] = true
      inclusionHigh[i] = true
      inclusionLow[i] = true
    else
      inclusion[i] = false
      inclusionLow[i] = false
      inclusionHigh[i] = false
      if (N[i].lo>X[i].lo)
        inclusionLow[i] = true
      elseif (N[i].hi<X[i].hi)
        inclusionHigh[i] = true
      end
    end
    if ~isdisjoint(N[i],X[i])
      X[i] = N[i] ∩ X[i]
    else
      exclusion = true
      break
    end
  end
  k += 1

  # checks if all components are included
  if (Iflag == false)
    for i=1:nx
      if (inclusion[i] == true)
        Iflag = true
        continue
      elseif (inclusion[i] == false)
        Iflag = false
        break
      end
    end
  end

  #
  while ((k<kmax) && isEqual(X,Xi,etol) == false)


    x_mid = mid.(X)
    h!(H,x_mid,P)
    hj!(J,X,P)
    Sparse_Precondition!(H,J,mid.(J),SSto)
    Mt = transpose(J)

    for i=1:nx
      S1 = Interval(0.0)
      S2 = Interval(0.0)
      S3 = Interval(0.0)
      for q=(Mt.colptr[i]):(Mt.colptr[i+1]-1)
        if (i < Mt.rowval[q])
          S1 += Mt.nzval[q]*(X[Mt.rowval[q]]-x_mid[Mt.rowval[q]])
        elseif (i > Mt.rowval[q])
          S2 += Mt.nzval[q]*(X[Mt.rowval[q]]-x_mid[Mt.rowval[q]])
        else
          S3 = Mt.nzval[q]
        end
      end
      if S3.lo*S3.hi > 0.0
        N[i] = x_mid[i] - (H[i]+S1+S2)/S3
      else
        Ntemp = N
        eD,N[i],Ntemp[i] = extProcess(N[i],X[i],S3,S1,S2,H[i],rtol)
        if eD == 1
          eDflag = true
          Xtemp = X
          Xtemp[i] = Ntemp[i] ∩ X[i]
          X[i] = N[i] ∩ X[i]
          return X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh
        end
      end
      if Strict_XinY(N[i],X[i])
        inclusion[i] = true
        inclusionHigh[i] = true
        inclusionLow[i] = true
      else
        inclusion[i] = false
        inclusionLow[i] = false
        inclusionHigh[i] = false
        if (N[i].lo>X[i].lo)
          inclusionLow[i] = true
        elseif (N[i].hi<X[i].hi)
          inclusionHigh[i] = true
        end
      end
      if ~isdisjoint(N[i],X[i])
        X[i] = N[i] ∩ X[i]
      else
        exclusion = true
        break
      end
    end
    k += 1

    # checks if all components are included
    if (Iflag == false)
      for i=1:nx
        if (inclusion[i] == true)
          Iflag = true
          continue
        elseif (inclusion[i] == false)
          Iflag = false
          break
        end
      end
    end
  end

  if exclusion
    Eflag = true
  end
  Xtemp = copy(X)
  return X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh
end

"""
    PIn_KrawczykCW(X0::Vector{Interval{T}},P::Vector{Interval{T}},hj!::Function,
                   h!::Function,opt::Vector{Any},Eflag::Bool,Iflag::Bool)

Uses the interval arithmetic presented in ValidatedNumerics to bound unique
implicit functions x:P->X defined by a system of equations h(z,p)=0 via
a component-wise Krawczyk method. Inputs:
* `X0::Vector{Interval{T}}`: Bounds for dependent variables
* `P::Vector{Interval{T}}`: Bounds for independent variables
* `hj!::Function`: Jacobian of h(x,p) with respect to x (takes input J,x,p)
* `h!::Function`: Equations defining potential implicit function (takes inputs H,x,p)
* `opt::Vector{Any}`: `opt[1]::Int64` is the number of iterations, `opt[2]::Float64`
                       is the equality tolerance for intervals, `opt[3]::Float64`
                       is the amount added to M[i,i] when processing extended
                       interval division.
* `Eflag::Bool`: True if no implicit function exists in box
* `Iflag::Bool`: True if unique implicit function exists in box
Returns a tuple `(X,Xtemp,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh)` where
* `X::Vector{Interval{T}}`: Output interval box
* `Xtemp::Vector{Interval{T}}`: Extra box formed from extended division
* `Eflag::Bool`: True if no implicit function exists in box
* `Iflag::Bool`: True if unique implicit function exists in box
* `inclusionLow::Vector{Bool}`: Each component is true if corresponding variable
                                contracted from lower bound
* inclusionHigh::Vector{Bool}`: Each component is true if corresponding variable
                                contracted from upper bound
"""
function PIn_KrawczykCW(X0::Vector{Interval{T}},P::Vector{Interval{T}},
                     hj!::Function,h!::Function,
                     opt::Vector{Any},Eflag::Bool,Iflag::Bool) where {T<:AbstractFloat}
  # unpacks option file
  kmax::Int64 = opt[1]
  etol::Float64 = opt[2]

  # Initializes variables
  nx::Int64 = length(X0)
  S::Vector{Interval{T}} = [Interval(0.0) for i=1:nx]
  inclusion::Vector{Bool} = [false for i=1:nx]
  inclusionHigh::Vector{Bool} = [false for i=1:nx]
  inclusionLow::Vector{Bool} = [false for i=1:nx]
  exclusion::Bool = false
  Iflag::Bool = false
  Eflag::Bool = false
  eDflag::Bool = false
  X::Vector{Interval{T}} = copy(X0)
  Xi::Vector{Interval{T}} = copy(X)
  N::Vector{Interval{T}} = [Interval(0.0) for i=1:nx]
  x_mid::Vector{T} = mid.(X)
  k::Int64 = 1
  # creates storage for h(xm,P), J_x(X,P), and preconditioner Y
  H::Vector{Interval{T}} = [Interval(0.0) for i=1:nx]
  J::SparseMatrixCSC{Interval{T},Int64} = spzeros(Interval{T},nx,nx)

  # initalizes storage object for in-place sparse calculations
  SSto = SparseInSto()
  SSto.Xh = copy(H)
  SSto.Yh = copy(H)
  SSto.Zh = copy(H)
  SSto.Xj = spzeros(Interval{T},nx,nx)
  SSto.Yj = spzeros(Interval{T},nx,nx)
  SSto.Zj = spzeros(Interval{T},nx,nx)
  SSto.nx = copy(nx)

  # calculates extension of h(xm,P) and J_x(X,P)
  h!(H,x_mid,P)
  hj!(J,X,P)
  Sparse_Precondition!(H,J,mid.(J),SSto)
  Mt = transpose(J)
  for i=1:nx
    for q=(Mt.colptr[i]):(Mt.colptr[i+1]-1)
      if (i == Mt.rowval[q])
        N[i] += (1.0-Mt.nzval[q])*(X[Mt.rowval[q]]-x_mid[Mt.rowval[q]])
      else
        N[i] += -Mt.nzval[q]*(X[Mt.rowval[q]]-x_mid[Mt.rowval[q]])
      end
    end
    N[i] = x_mid[i] - H[i] + N[i]

    if Strict_XinY(N[i],X[i])
      inclusion[i] = true
      inclusionHigh[i] = true
      inclusionLow[i] = true
    else
      inclusion[i] = false
      inclusionLow[i] = false
      inclusionHigh[i] = false
      if (N[i].lo>X[i].lo)
        inclusionLow[i] = true
      elseif (N[i].hi<X[i].hi)
        inclusionHigh[i] = true
      end
    end
    if ~isdisjoint(N[i],X[i])
      X[i] = N[i] ∩ X[i]
    else
      exclusion = true
      break
    end
  end
  k += 1

  if (Iflag == false)
    for i=1:nx
      if (inclusion[i] == true)
        Iflag = true
        continue
      elseif (inclusion[i] == false)
        Iflag = false
        break
      end
    end
  end
  #
  while ((k<kmax) && (isEqual(X,Xi,etol) == false))
    Xi = X0
    x_mid = mid.(X)
    h!(H,x_mid,P)
    hj!(J,X,P)
    Sparse_Precondition!(H,J,mid.(J),SSto)
    Mt = transpose(J)
    N = [Interval(0.0) for i=1:nx]
    for i=1:nx
      for q=(Mt.colptr[i]):(Mt.colptr[i+1]-1)
        if (i == Mt.rowval[q])
          N[i] += (1.0-Mt.nzval[q])*(X[Mt.rowval[q]]-x_mid[Mt.rowval[q]])
        else
          N[i] += -Mt.nzval[q]*(X[Mt.rowval[q]]-x_mid[Mt.rowval[q]])
        end
      end
      N[i] = x_mid[i] - H[i] + N[i]

      if Strict_XinY(N[i],X[i])
        inclusion[i] = true
        inclusionHigh[i] = true
        inclusionLow[i] = true
      else
        inclusion[i] = false
        inclusionLow[i] = false
        inclusionHigh[i] = false
        if (N[i].lo>X[i].lo)
          inclusionLow[i] = true
        elseif (N[i].hi<X[i].hi)
          inclusionHigh[i] = true
        end
      end
      if ~isdisjoint(N[i],X[i])
        X[i] = N[i] ∩ X[i]
      else
        exclusion = true
        break
      end
    end
    k += 1

    # checks if all components are included
    if (Iflag == false)
      for i=1:nx
        if (inclusion[i] == true)
          Iflag = true
          continue
        elseif (inclusion[i] == false)
          Iflag = false
          break
        end
      end
    end
  end

  if exclusion
    Eflag = true
  end
  return X,Eflag,Iflag,inclusionLow,inclusionHigh
end
