# To Do List:
# Finish Documentation

################# Defines Parametric Interval Tests #########################
"""Exclusion_Test(xk1,xk2,method="Miranda") if xk1 and xk2 are distinct boxes
(aka occuring at different iterations) generated from a paramtric interval
method beginning with X0 and the intersection of xk1 and xk2 is empty then the
initial box X0 doesn't contain solutions of h for all p in P. This test is
sufficient but not necessary.
- returns 1 if the box doesn't contain a function x(p)
- returns 0 otherwise
"""
function Exclusion_Test(xk1,xk2;method="Miranda")
  flag = 0
  for i=1:length(xk1)
    if isempty(xk1[i] ∩ xk2[i])
      flag = 1
      break
    end
  end
  return flag
end
"""Exclusion_Test(h,X,P; method="Stuber") provides a test for a partially
enclosed solution
- returns 1 if the box contains a unique function x(p)
- returns 0 otherwise
"""
function Inclusion_Test(h,hj,P,Y,xk1,xk2; method="Stuber")
  # checks to see if solution branch intersect bounds
  for i=1:length(xk1)
    if (xk2[i].lo == xk1[i].lo)
      xk2temp_arr = []
      for j=1:length(xk2)
        if j==i
          push!(xk2temp_arr,Interval(xk2[j].lo,xk2[j].lo))
        else
          push!(xk2temp_arr,xk2[j])
        end
      end
      xk2temp = IntervalBox(tuple(xk2temp_arr...))
      temp = h(xk2temp,P)
      for j=1:length(xk2)
        if (temp[j].lo*temp[j].hi<=0)
          return 0
        end
      end
    elseif (xk2[i].hi == xk1[i].hi)
      xk2temp_arr = []
      for j=1:length(xk2)
        if j==i
          push!(xk2temp_arr,Interval(xk2[j].hi,xk2[j].hi))
        else
          push!(xk2temp_arr,xk2[j])
        end
      end
      xk2temp = IntervalBox(tuple(xk2temp_arr...))
      temp = h(xk2temp,P)
      for j=1:length(xk2)
        if (temp[j].lo*temp[j].hi<=0)
          return 0
        end
      end
    end
  end
  # calculates eigenvalue
  A = (abs.(Y))*radius.(hj(xk2,P))
  if (length(A)>1)
    eigM = eigmax(A)
  else
    eigM = A[1]
  end
  # performs tests
  if eigM < 1
    return 1
  else
    return 0
  end
end
"""Inclusion_Test(xk1,xk2,method="Miranda") if xk1 and xk2 are distinct boxes
(aka occuring at different iterations) generated from a paramtric interval
method beginning with X0 and the intersection of xk1 and xk2 is internal
to xk1 then x(p) is a unique function contained in X. This test is sufficient
but not necessary.
- returns 1 if the box contains a unique function x(p)
- returns 0 otherwise
"""
function Inclusion_Test(xk1,xk2;method="Miranda")
  flag = 1
  for i=1:length(xk1)
    if (xk2[i].lo <= xk1[i].lo) || (xk2[i].hi >= xk1[i].hi)
      flag = 0
    end
  end
  return flag
end

"""Param_Intv_ContactorTest(k,X,P,h,hj,style)
Inputs:
- k number of iterations to run
- X is state space IntervalBox
- P is control space IntervalBox
- h is
- hj is the jacobian of h
- style is "Newton", "NewtonGS", "Krawczyk", "KrawczykCW"
Outputs:
- returns 1 if a unique implicit function is found in box (X,P)
- returns 0 if inconclusive
- returns -1 if no unique implicit function is in the box (X,P)
"""
function Param_Intv_ContactorTest(k,X,P,h,hj,style)

  Xprev = copy(X)
  Xnext = copy(X)

  for i = 1:k

    # calculates Preconditioner
    z = mid.(Xnext)
    Y = Preconditioner(hj,X,P,jac="User") # throws error if singular

    # advances iteration
    if (style=="Krawczyk")
      Xnext = Param_Intv_KrawczykCW(z,Xprev,P,h,hj,Y)
    elseif (style=="KrawczykCW")
      Xnext = Param_Intv_KrawczykCW(z,Xprev,P,h,hj,Y)
    elseif (style=="Newton")
      Xnext = Param_Intv_Newton(z,Xprev,P,h,hj,Y)
    elseif (style=="NewtonGS")
      Xnext = Param_Intv_NewtonGS(z,Xprev,P,h,hj,Y)
    else
      error("Unsupported style of contractor")
    end

    # runs stepwise tests
    inc_test1 = Inclusion_Test(Xprev,Xnext,method="Miranda")
    inc_test2 = Inclusion_Test(h,hj,P,Y,Xprev,Xnext,method="Stuber")
    exc_test1 = Exclusion_Test(Xprev,Xnext,method="Miranda")

    # interprets test results
    if (inc_test1 == 1) || (inc_test2 == 1)
      return 1
    elseif (exc_test1 == 1)
      return -1
    end

    # sets previous to iteration result
    Xprev = Xnext

  end
  return 0
end
