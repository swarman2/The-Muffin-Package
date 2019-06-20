#SV takes m,s as input and returns 
#V, sᵥ, and sᵥ₋₁
#using the V-Conjecture (pg 66) and
#a linear function from pg 68 
#gotten by solving the system:
#   Vsᵥ + (V-1)(sᵥ₋₁) = pieces = 2m
#   sᵥ + sᵥ₋₁ = s
function SV(m,s)
  V=Int64(ceil(2*m/s))
  sᵥ₋₁=V*s-2*m
  sᵥ = 2*m - s*(V-1)
  return V, sᵥ,sᵥ₋₁
end

#FINDEND input m,s alpha and V
#this returns x and y where all V shares are in 
#the interval (𝛂, x) and all V-1 shares are in the 
#interval (y, 1-𝛂)
#This is adapted from psudocode on page 69
function FINDEND(m,s,alpha,V)
  y=m//s -(1-alpha)*(V-2)
  if y>= (1-alpha)
    y=1-alpha
  end
  if y <= alpha
    y=alpha
  end
  x=m//s - alpha*(V-1)
  if x<=alpha
    x=alpha
  end
  if x>= (1-alpha)
    x=1-alpha
  end
  return x,y
end

#VHALF verifies if f(m,s)≦ 𝛂 for some m,s and alpha
#using the half method (i.e. finds all the info and 
#checks to make sure there is a contradiction in "case 5")
#from pg 70
function VHALF(m,s,𝛂)
  if 𝛂<1//3
    println("Bad input")
    return false
  end
  if 𝛂 >1//2
    println("𝛂 > 1//2 is a bad guess, please rethink")
    exit(0)
  end
  (V,sᵥ, sᵥ₋₁)=SV(m,s)
  (x,y)=FINDEND(m,s,𝛂, V)

  #check to see if there are too many pieces smaller than a half
  #or too many pieces larger than a half
  if x<=.5 && V*sᵥ>m
    return true
  end
  if y>=.5 && (V-1)*sᵥ₋₁>m
     return true
  end
return false
end

#HALF 
#input :m,s, output: alpha, 1, or prints not a useful answer
#adapted from pg 71-72
function HALF(m,s)
  if m%s == 0
    return 1
  end  
  V,sᵥ,sᵥ₋₁=SV(m,s)
#  println("V=",V,"  Sv =",sᵥ,"  Sv-1 = ", sᵥ₋₁)
  if (V-1)*sᵥ₋₁>V*sᵥ
    alpha = 1-((m//s-1//2)//(V-2))
    if(alpha<1//3)
      alpha=1//3
    end
    if(VHALF(m,s,alpha))
      return alpha
    end
  elseif (V-1)*sᵥ₋₁<V*sᵥ
    alpha = (m//s-1//2)//(V-1)
    if(alpha<1//3)
      alpha=1//3
    end
    if(VHALF(m,s,alpha))
      return alpha
    end
  else
    println("Not a useful answer")
    return 0
  end
end

#println(VHALF(11,5,13//30))
#println(VHALF(45,26,32/78))
