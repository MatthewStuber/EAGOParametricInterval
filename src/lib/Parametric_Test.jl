function Miranda(h,X,P)
  H = h(X,P)
  #println("H: ", H)
  for i=1:length(X)
    if (H[i].lo*H[i].hi>0.0)
      return i
    end
  end
  i = -1
  return i
end

function MirandaExc(h,X,P,Eflag_in)
   Eflag = Eflag_in
   i = Miranda(h,X,P)
   if (i != -1)
     Eflag = true
   end
   return Eflag
end

function partialIncTop(h,X,P,PIflag,incHigh)
  xU = [X[i].hi for i=1:length(X)]
  XU = copy(X)
  Xprev = copy(X)
  k = -1
  ku = -1
  for i=1:length(X)
    #println(" --- --- --- partial top test iter:", i)
    #println("P:  ",P)
    #XU[i] = Xprev[i] # CHECK ME THIS DOESN'T SEEM RIGHT
    XU[i] = Interval(xU[i])
    #println("XU[i]:  ",xU[i])
    #println("htest top:  ", [(3.25-XU[1])/P[1]-XU[3];
    #           XU[1]/P[2]-XU[3];
    #           XU[2]-(XU[1]^2)/(1+XU[1]^2)])
    ku = Miranda(h,XU,P)
    #println("ku:  ", ku)
    if (ku == -1 && ~incHigh[i])
      #println("assignment top:  ", i)
      k = i
      break
    else
      PIflag = false
    end
    XU[i] = Xprev[i]
  end
  return k,PIflag
end

function partialIncBot(h,X,P,PIflag,incLow)
  xL = [X[i].lo for i=1:length(X)]
  XL = copy(X)
  Xprev = copy(X)
  k = -1
  kl = -1
  for i=1:length(X)
    #println(" --- --- --- partial bottom test iter:", i)
    #println("P:  ",P)
    #XL[i] = Xprev[i] # CHECK ME THIS DOESN'T SEEM RIGHT
    XL[i] = Interval(xL[i])
    #println("XL[i]:  ",XL[i])
    #println("htest bottom:  ", [(3.25-XL[1])/P[1]-XL[3];
    #           XL[1]/P[2]-XL[3];
    #           XL[2]-(XL[1]^2)/(1+XL[1]^2)])
    kl = Miranda(h,XL,P)
    #println("kl:  ", kl)
    if (kl == -1 && ~incLow[i])
      #println("assignment bottom:  ", i)
      k = i
      break
    else
      PIflag = false
    end
    XL[i] = Xprev[i]
  end
  return k,PIflag
end


function spectralR(hj,X,P,rho)
  k = true
  J = hj(X,P)
  Y = Preconditioner(hj,X,P,jac="User")
  YJ = abs.(Y)*diam.(J)/2
  if (length(YJ)>1)
    eigval = eigmax(YJ)
  else
    eigval = YJ[1]
  end
  return k,eigval
end

function BoundaryTest(h,hj,X0,X,P,opt,PIcert,PIflag,Iflag,Eflag,inclusionLow,inclusionHigh)
  #println("inclusionLow: ", inclusionLow)
  #println("inclusionHigh: ", inclusionHigh)
  PIflagT = PIflag
  PIflagB = PIflag
  Xtemp = copy(X)
  Xin = copy(X)
  Pin = copy(P)
  ku,PIflagT = partialIncTop(h,Xin,Pin,PIflagT,inclusionHigh)
  kl,PIflagB = partialIncBot(h,Xin,Pin,PIflagB,inclusionLow)
  #println("ku: ", ku)
  #println("kl: ", kl)
  rho = 1.01
  eDflag = false
  #println("X0: ", X0)
  #println("X: ", X)
  #println("Strict_XinY(X,X0): ", Strict_XinY(X,X0))
  #println("((ku==-1) && (kl==-1)): ", ((ku==-1) && (kl==-1)))
  if (Strict_XinY(X,X0) || ((ku==-1) && (kl==-1)))
    blank,rho = spectralR(hj,Xin,Pin,rho)
    PIcert = true
    PIflag = false
    #println("rho: ", rho)
    if (rho<1.0)
      Pmid = mid.(Pin)
      #println("Pmid: ",Pmid)
      if (opt[4]=="NewtonGS")
        #println("Xin: ",Xin)
        Xnew,Xtemp,Eflag,Iflag,eDflag = MC_NewtonGS(Xin,Pmid,hj,h,opt,Eflag,Iflag,eDflag)
        #println("Xnew: ",Xnew)
      elseif (opt[4]=="KrawczykCW")
        #println("Xin: ",Xin)
        Xnew,Eflag,Iflag,blank1,blank2 =  MC_KrawczykCW(Xin,Pmid,hj,h,opt,Eflag,Iflag)
        #println("Xnew: ",Xnew)
      end
      #println("Iflag tight: ",Iflag)
    else
      Iflag = false
    end
  end
  return PIcert,PIflag,Iflag,Eflag
end
#=
function NewBoundaryTest(Xtemp,X,P,Iflag)
  nx = length(X)
  kmax = 1000
  k = 1
  l = 0
  ztemp = 0.0
  eps = 1.0E-12
  Xf = copy(X)
  Xnew = copy(X)
  Pf = copy(P)
  l = SingleStep(Xnew,Pf)
  if (l == 0)
    k = kmax
  end
  for i = 1:nx
    r[2*i-1] = Xf[i].lo - Xnew[i].lo + eps
    r[2*i] = Xnew[i].hi - Xf[i].hi + eps
    z[2*i-1] = Xf[i].lo
    z[2*i] = Xf[i].hi
  end
  # construct Jacobian
  Jr =

end

function SingleStep()
end
=#
