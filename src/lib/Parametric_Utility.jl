# Will update interface utilities with upgrade to Generalized parametric bisection

#=
function kdel(i,j)
  i == j ? 1.0 : 0.0
end

function Expr_SubSymbols(expr,symdict)
      stack = []
      exprlist = []
      push!(stack,expr)
      cont_flag = false
      while length(stack)>0
            a = pop!(stack)
            if (isempty(exprlist))
                  push!(exprlist,a)
            end
            if ~(typeof(a)<:Number) && ~(typeof(a)==Int64)
              for i=2:length(a.args)
                  if (typeof(a.args[i])==Symbol)
                        a.args[i] = symdict[a.args[i]]
                  elseif typeof(a.args[i])==Expr
                        push!(stack,a.args[i])
                  end
              end
            end
      end
      return exprlist[1]
end

function GenerateJacobianX(ExprArr,SymXArr,SymPArr)

  # creates symbol dictionary
  symdict = Dict{Any,Any}()
  nc = length(SymXArr)
  ns = length(SymPArr)
  for i=1:nc
      symdict[SymXArr[i]] = :(z[$i])
  end
  for i=1:ns
      symdict[SymPArr[i]] = :(p[$i])
  end

  # creates calculates Jacobian symbolically and subs in symbols
  StringArr = string.(ExprArr)
  dexpr_arr = Expr(:vcat)
  for i = 1:length(StringArr)
    out = Calculus.differentiate(StringArr[i],SymXArr)
    expr_out = copy(out)
    for j = 1:length(out)
      expr_out[j] = Expr_SubSymbols(out[j],symdict)
    end
    temp_expr = Expr(:row)
    for j = 1:length(out)
      push!(temp_expr.args,out[j])
    end
    push!(dexpr_arr.args,temp_expr) # populations Jacobian storage
  end

  # puts expression object into function
  @eval Jacob(z,p) = $dexpr_arr
  return Jacob
end

function GenerateH(ExprArr, SymXArr,SymPArr)
    # creates list of symbols
    symdict = Dict{Any,Any}()
    for i=1:length(SymXArr)
        symdict[SymXArr[i]] = :(x[$i])
    end
    for i=1:length(SymPArr)
        symdict[SymPArr[i]] = :(p[$i])
    end
    # creates array of expressions
    htemp1 = deepcopy(ExprArr)
    htemp2 = []
    for i = 1:length(ExprArr)
        push!(htemp2,Expr_SubSymbols(htemp1[i],symdict))
    end
    #
    hexpr = Expr(:vect)
    for i = 1:length(htemp2)
        temp = htemp2[i]
        push!(hexpr.args,Expr(:call,:eval,temp))
    end
    @eval hnew(x,p) = $hexpr
    return hnew
end

# Unpack function
function Unpack(stack)
  cnode = pop!(stack)
  return cnode[1],cnode[2]
end
=#

"""
    Strict_XinY(X::Vector{Interval{T}},Y::Vector{Interval{T}})

Returns true if X is strictly in Y (X.lo>Y.lo && X.hi<Y.hi).
"""
function Strict_XinY(X::Vector{Interval{T}},Y::Vector{Interval{T}}) where {T<:AbstractFloat}
  k::Bool = true
  for i=1:length(X)
    if ((X[i].lo<=Y[i].lo)||
        (X[i].hi>=Y[i].hi))
      k = false
    end
  end
  return k
end

"""
    Strict_XinY(X::Interval{T},Y::Interval{T})

Returns true if X is strictly in Y (X.lo>Y.lo && X.hi<Y.hi).
"""
function Strict_XinY(X::Interval{T},Y::Interval{T}) where {T<:AbstractFloat}
  (X.lo<=Y.lo) && return false
  (X.hi>=Y.hi) && return false
  return true
end

"""
    isEqual(X1::Vector{Interval{T}},X2::Vector{Interval{T}},atol::Float64)

Returns true if X1 and X2 are equal to within tolerance atol in all dimensions.
"""
function isEqual(X1::Vector{Interval{T}},X2::Vector{Interval{T}},
                 atol::Float64) where {T<:AbstractFloat}
  out::Bool = true
  for i=1:length(X1)
    if (abs(X1[i].lo-X2[i].lo)>=atol ||
        abs(X1[i].hi-X2[i].hi)>=atol )
        out = false
        break
    end
  end
  return out
end

"""
    extDivide(A::Interval{T},B::Interval{T},C::Interval{T})

Subfunction to generate output for extended division.
"""
function extDivide(A::Interval{T}) where {T<:AbstractFloat}
  if ((A.lo == -0.0) && (A.hi == 0.0))
    B::Interval{T} = Interval(-Inf,Inf)
    C::Interval{T} = B
    return 0,B,C
  end
  if (A.lo == 0.0)
    B = Interval(1.0/A.hi,Inf)
    C = Interval(Inf,Inf)
    return 1,B,C
  elseif (A.hi == 0.0)
    B = Interval(-Inf,1.0/A.lo)
    C = Interval(-Inf,-Inf)
    return 2,B,C
  else
    B = Interval(-Inf,1.0/A.lo)
    C = Interval(1.0/A.hi,Inf)
    return 3,B,C
  end
end

"""
    extProcess(N::Interval{T},X::Interval{T},Mii::Interval{T},
                    S1::Interval{T},S2::Interval{T},B::Interval{T},rtol::Float64)

Generates output boxes for extended division and flag.
"""
# Error & Case Handling for Extended IA in NewtonGS (CHECKED C++)
function extProcess(N::Interval{T},X::Interval{T},Mii::Interval{T},
                    S1::Interval{T},S2::Interval{T},B::Interval{T},rtol::Float64) where {T<:AbstractFloat}
  v = 1
  Ntemp::Interval{T} = copy(N)
  M::Interval{T} = (B+S1+S2)+Interval(-rtol,rtol)
  if (M.lo<=0 && M.hi>=0)
    println("branch1")
    return 0, Interval(-Inf,Inf), Ntemp
  end
  if (v == 1)
    println("Mii: $Mii")
    k,IML::Interval{T},IMR::Interval{T} = extDivide(Mii)
    println("k: $k")
    if (k == 1)
      println("branch2")
      return 0, (mid(X)-M*IML), Ntemp
    elseif (k == 2)
      println("branch3")
      return 0, (mid(X)-M*IMR), Ntemp
    elseif (k == 3)
      NR = mid(X)-M*IMR
      NL = mid(X)-M*IML
      print("NR: $NR")
      print("NL: $NL")
      if (~isdisjoint(NL,X) && isdisjoint(NR,X))
        println("branch4")
        return 0, NL, Ntemp
      elseif (~isdisjoint(NR,X) && isdisjoint(NL,X))
        println("branch5")
        return 0, NR, Ntemp
      elseif (~isdisjoint(NL,X) && ~isdisjoint(NR,X))
        println("branch6")
        N = NL
        Ntemp = NR
        print("NR: $NR")
        print("NL: $NL")
        return 1, NL, NR
      else
        println("branch7")
        return -1, N, Ntemp
      end
    end
  end
  println("branch8")
  return 0, N, Ntemp
end
