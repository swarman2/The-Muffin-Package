include("helper_functions.jl")
#for information on this method see chapter 6 of "The Mathematics of Muffins"
#VHALF verifies if f(m,s)≦ alpha for some m,s and alpha
#using the half method (i.e. finds all the info and
#checks to make sure there is a contradiction in "case 5")
function VHALF(m,s,alpha)
  if alpha<1//3
    return false
  end
  if alpha >1//2
    return false
  end
  (V,sᵥ, sᵥ₋₁)=SV(m,s)
  (x,y)=FINDEND(m,s,alpha, V)

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

function HALF(m,s)
  if m%s == 0
    return 1
  end
  V,sᵥ,sᵥ₋₁=SV(m,s)
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
    return 1
  end
end

#Half proof prints out a proof
function Half_proof(m,s,alpha)
  if alpha >1//2
    println("alpha > 1//2 is a bad guess, please rethink")
    exit(0)
  end

  V =Int64( ceil(2*m/s))
  x,y=FINDEND(m,s,alpha,V)
  if x == 1-alpha || y == alpha
    println("WARNING: ENDPOINTS SHIFTED HALF METHOD MAY FAIL")
  end
  pieces=2*m
  println(m," muffins, ",pieces," pieces") #cut every muffin in exactly two pieces
  println("At most ", m, " pieces larger than 1/2")#or smaller than b/c buddy (two lines!)
  println("At most ", m, " pieces smaller than 1/2")
  #check case more than five shares a person
  #if Alice has more than five shares then one of them has to be less than (or equal to?) the amount of muffins
  #Alice has divided by five (amt of muffins Alice has is m/s)

  share= m//(s *(V+1))
  println()
  println("Case 1: ")
  println("if Alice has ",V+1, " or more shares then there exists a share ≤ ", numerator(share),"/",denominator(share))
  #if the share that Alice has is less than our hypthoisized alpha than we know
  #Alice has less than five shares
  less_than_upper = share<alpha #boolean that is true if Alice has less than five shares
  println("We will only consider the case where Alice has less than ", V+1," shares if share is <  alpha")
  println("share < alpha: ", less_than_upper)
  println()
  println("Case 2:")
  #if Bob has less than two shares then one of them has to be greater than
  #1/2 times amt Bob has and its buddy is less than 1-(1/2*amt Bob has)
  share = m//(s *(V-2))
  buddy=(1-share)

  more_than_lower = buddy<alpha
  print("if Bob has ", V-2," or less shares then there exists a share > ", numerator(share),"/",denominator(share))
  println(" whose buddy is < ", numerator(buddy),"/",denominator(buddy))
  println("We will only consider the case where Bob has more than ",V-2," shares if buddy < alpha")
  println("buddy < alpha: ",more_than_lower)
  println()
  if(less_than_upper && more_than_lower)

    println(V-1,"sᵥ₋₁ + ",V,"sᵥ = ", pieces)
    println("sᵥ₋₁ + sᵥ = ", s)
    #These equations were found by manually solving the above matrix equation
    s3=s*V-pieces
    s4=pieces-s*V+s
    println("There are ", s3, " ",V-1,"-students, ", s4," ",V,"-students, ",(V-1)*s3, " ",V-1,"-shares, and ", V*s4, " ",V,"-shares.")
    println("Endpoints are x = ",numerator(x),"/",denominator(x)," and y = ",numerator(y),"/",denominator(y))
    #Alice has a 4-share > x (x=m/s - 3alpha)
    println()

    println("Case 3:")
    println("Alice has a ", V,"-share > ",numerator(x),"/",denominator(x))
    temp=(m-x*s)//s
    println("Alice's other ", V-1," ",V, "-shares add up to < ", m,"/",s, " - ", numerator(x),"/",denominator(x)," = ", numerator(temp),"/",denominator(temp))
    print("so one of Alice's shares must be < ",numerator(temp),"/",denominator(temp), " * 1/",V-1," = ")
    temp=numerator(temp)//(denominator(temp)*(V-1))
    println(numerator(temp),"/",denominator(temp))
    if(temp == alpha)
      println("so one of Alice's shares must be less than alpha")
    elseif(temp<alpha)
      println(numerator(temp),"/",denominator(temp)," < ",numerator(alpha),"/",denominator(alpha),"  so one of Alice's shares must be less than alpha")
    end

    println()
    println("Case 4:")
    #Bob has a 3-share < y (y=m/s - 2*(1-alpha))
   temp=(m-y*s)//s
    println("Bob has a ", V-1,"-share < ",y)
    println("Bob's other ",V-2," ",V-1,"-shares add up to > ",m,"/",s, " - ", y, " = ", numerator(temp),"/",denominator(temp))
    print("so one of Bob's shares must be > ",numerator(temp),"/",denominator(temp)," * 1/",V-2," =")
    temp=temp//(V-2)
    println(numerator(temp),"/",denominator(temp))
    println("whose buddy would be < 1-",temp," = ",numerator(Rational(1-temp)),"/",denominator(Rational(1-temp)))
    if(1-temp == alpha)
      println("so one of Bob's shares must be less than alpha")
    elseif((1-temp)<alpha)
      println(numerator(1-temp),"/",denominator(1-temp)," < ",numerator(alpha),"/",denominator(alpha),"  so one of Bob's shares must be less than alpha")
    end

    println()
    if x<=.5 &&( s4*V)>m #if x is a half and there are more than m shares less than x
      println("Case 5:")
      println("This case is the negation of the other cases")
      println("There are ",s4*V," shares less than a half, this is a contradiction so this case will never happen")
    elseif y>=.5 && (s3*(V-1))>m #if y is greater than a half and there are more than m shares greater than y
      println("Case 5:")
      println("This case is the negation of the other cases")
      println("There are ",s3*(V-1)," shares greater than a half, this is a contradiction so this case will never happen")
    else
      println("We do not know if there are too many shares more/less than a half: PROOF FAILED")
    end
  end
end
