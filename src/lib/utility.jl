############ Defines utilities functions for future calculations ############
"""
    kdel(i,j)

Kronecker's delta function.

# Examples:

'''jldoctest
julia> kdel(1,3)
  0.00
'''

'''jldoctest
julia> kdel(3,3)
  1.00
'''

"""
function kdel(i,j)
  i == j ? 1.0 : 0.0
end

"""
    Expr_SubSymbols(expr,symdict)

Implements a depth-first search a replace routine in which the symbols in 'expr' are substituted for
new symbols using 'symdict' dictionary.

# Example

'''jldoctest
julia> Expr_SubSymbols(:(b1+3),Dict(:(b1) => :a1))
  :(a1 + 3)
'''

"""
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

############ Jacobian Calculations ############
"""
    GenerateJacobianX(ExprArr,SymXArr,SymPArr)

Generates a function that calculates the jacobian of an array of expressions,
'ExprArr', provided an array of symbols corresponding to X, 'SymXArr', and an
array of symbols corresponding to P, 'SymPArr'.

"""
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

"""
    GenerateH(ExprArr,SymXArr,SymPArr)

Generates a function that calculates the h(x,p) from an array of expressions,
'ExprArr', provided an array of symbols corresponding to X, 'SymXArr', and an
array of symbols corresponding to P, 'SymPArr'.
"""
function GenerateH(ExprArr,SymXArr,SymPArr)
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

"""
    JacobianXFD(h,x,p)

Calculates the jacobian of the h(x,p) with respect to the x variable via
forward mode automatic differentation using the ForwardDiff.jl package.
"""
function JacobianXFD(h,x,p)
  nx = length(x)
  np = length(p)
  yvar = vcat(x,p)
  hy(y) = h(y[1:nx],y[(nx+1):(nx+np)])
  ForwardDiff.jacobian(hy,yvar)[1:nx,1:nx]
end

"""
    JacobianXRD(h,x,p)


Calculates the jacobian of the h(x,p) with respect to the x variable via
reverse mode automatic differentation using the ReverseDiff.jl package. 
"""
function JacobianXRD(h,x,p)
  nx = length(x)
  np = length(p)
  yvar = vcat(x,p)
  hy(y) = h(y[1:nx],y[(nx+1):(nx+np)])
  ReverseDiff.jacobian(hy,yvar)[1:nx,1:nx]
end

################# Preconditioning Routines #####################
"""Preconditioner(h,X,P;jac="ForwardDiff") calculate the inverse midpoint Jacobian preconditioner
"""
function Preconditioner(h,X,P;jac="ForwardDiff")
  J = JacobianXFD(h,X,P)
  if (length(X)>1)
    Y = inv(mid.(J))
  else
    Y = (mid.(J))^(-1)
  end
  return Y
end

"""Preconditioner(h,X,P;jac="ReverseDiff") calculate the inverse midpoint Jacobian preconditioner
"""
function Preconditioner(h,X,P;jac="ReverseDiff")
  J = JacobianXRD(h,X,P)
  if (length(X)>1)
    Y = inv(mid.(J))
  else
    Y = (mid.(J))^(-1)
  end
  return Y
end

"""Preconditioner(h,X,P;jac="User") calculate the inverse midpoint Jacobian preconditioner
"""
function Preconditioner(h,X,P;jac="User")
  J = h(X,P)
  if (length(X)>1)
    Y = inv(mid.(J))
  else
    Y = 1.0/(mid(J[1]))
  end
  return Y
end
