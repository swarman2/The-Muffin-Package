include("helper_functions.jl")
#For information on this method see chapter 8 of "The Mathematics of Muffins"
"""
VINT verifies

INPUTS:
 * m -> number of muffins
 * s -> number of students
 * alphaa -> α
 * proof -> boolean:
   (0 - No Proof)
   (1 - Print Proof)

OUTPUTS:
 * returns true if it INT can prove the upper-bound
 * Print proof (Depending on input (proof))
 """
function VINT(m,s,alpha, proof = false)
  if m%s == 0
    if proof
      println("s divides m, give all students whole muffin(s)")
    end
    return 1
  end
 V,sᵥ,sᵥ₋₁=SV(m,s)
 Vshares=V*sᵥ
 V₋₁shares=(V-1)*sᵥ₋₁
 x,y=FINDEND(m,s,alpha,V)
 if x > y
     if proof
         println("Intervals not disjoint")
     end
     return false, 0
 end

 xbuddy = 1-x
 ybuddy = 1-y

 smallShare = Array{Int64}(undef, V+1)#declare a 1 d array
 for i=1:V+1
   smallShare[i]=-1
 end

 if proof
   print_Intervals(m,s,alpha,false)
 end

  if(V₋₁shares<Vshares)

   #___(_______)________(_____)___|___(_____)____
   #alpha   y-buddy   xbuddy   x  |   y    1-alpha
   #smallShare is an array that holds the possiblities for number of small shares

   num_small_shares = V₋₁shares
   num_large_shares = Vshares - V₋₁shares

   j=1
   smallShare[1]=0
   for numSm = 0:V
    numLg = V-numSm
    if((ybuddy*numSm+x*numLg >m//s) && (alpha*numSm + xbuddy *numLg<m//s))
      smallShare[j]=numSm
      j = j+1
    end
   end
   minSm=Inf
   minLg=Inf
   if proof
      println("Possible distribution of shares: ")
    end
   for i=1:V+1
     if(smallShare[i]!=-1)
       if(smallShare[i]<minSm)
         minSm=smallShare[i]
       end
       if(minLg>V-smallShare[i])
         minLg = V-smallShare[i]
       end
       if proof
         println("    ",smallShare[i]," small share(s) and ",V-smallShare[i]," large share(s)")
       end
     end
   end
   if proof
       if(minSm!=0)
         println("Need at least ", minSm, " small shares")
       end
       if(minLg!=0)
         println("Need at least ",minLg, " large shares")
       end
       if(minSm*sᵥ>num_small_shares)
         println(sᵥ," students need at least ", minSm, " small shares, but there are only ", num_small_shares, " small shares, so f(",m,", ",s,")  ≤ ", numerator(alpha),"/",denominator(alpha))
       end
       if(minLg*sᵥ>num_large_shares)
         println(sᵥ," students need at least ", minLg, " large shares, but there are only ", num_large_shares, " large shares, so f(",m,", ",s,")  ≤ ", numerator(alpha),"/",denominator(alpha))
       end
     end
     if(minSm*sᵥ>num_small_shares)
       return true
     end
     if(minLg*sᵥ>num_large_shares)
       return true
    end
    return false
   elseif(V₋₁shares>Vshares)

   #___(_______)_|__(_________)__________(________)____
   #alpha      x |  y        ybuddy   xbuddy    1-alpha

   num_large_shares = Vshares
   num_small_shares = V₋₁shares - Vshares

   j=1
   smallShare[1]=0
   for numSm = 0:V-1
     numLg = V-numSm -1
     if ((ybuddy*numSm + (1-alpha)*numLg > m//s)&&(y * numSm + xbuddy*numLg<m//s))
      smallShare[j]=numSm
      j=j+1
     end
   end
   minSm=Inf
   minLg=Inf
   for i=1:V
     if(smallShare[i]!=-1)
       if(smallShare[i]<minSm)
         minSm=smallShare[i]
       end
       if(minLg>V-smallShare[i]-1)
         minLg = V-smallShare[i]-1
       end
       if proof
         println(smallShare[i]," small shares and ",V-smallShare[i]-1," large shares works")
       end
     end
   end
   if proof
     if(minSm!=0)
       println("Need at least ", minSm, " small shares")
     end
     if(minLg!=0)
       println("Need at least ",minLg, " large shares")
     end
     if(minSm*sᵥ₋₁>num_small_shares)
       println(sᵥ₋₁," students need at least ", minSm, " small shares, but there are only ", num_small_shares, " small shares, so f(",m,", ",s,") ≤ ", numerator(alpha),"/",denominator(alpha))
     end
     if(minLg*sᵥ₋₁>num_large_shares)
      println(sᵥ₋₁," students need at least ", minLg, " large shares, but there are only ", num_large_shares, " large shares, so f(",m,", ",s,")  ≤ ", numerator(alpha),"/",denominator(alpha))
     end
   end
   if(minSm*sᵥ₋₁>num_small_shares)
    return true
   end
   if(minLg*sᵥ₋₁>num_large_shares)
    return true
   end
   return false
 else
   return false
 end
end
"""
INT is the Interval Method

INPUTS:
 * m -> number of muffins
 * s -> number of students

OUTPUTS:
 * α

"""
function INT(m,s)
  if m%s == 0
    return 1
  end
 V,sᵥ,sᵥ₋₁=SV(m,s)
 Vshares=V*sᵥ
 V₋₁shares=(V-1)*sᵥ₋₁
 alphaPoss= Array{Rational}(undef, V+1)#declare a 1 d array
 for i =1:V+1
   alphaPoss[i]=1
 end

 if(V₋₁shares<Vshares)
   #split the V shares
   num_small_shares= V₋₁shares
   num_large_shares = Vshares-num_small_shares
   #(    )      (    )       (    )
   #𝛂   1-y    1-x    x      y    1-𝛂
   #x = m//s -alpha*(v-1)
   #y = m//s - (1-alpha)*(V-2)
   min_large = Int64(floor(num_large_shares/sᵥ))
   #m//s >= (1-y)*(V-min_large) + x*min_large
  j=1
  for i=0:V
    min_large=V -i
    alpha=(m//s + V + (m*V)//s - V*V -min_large - (2*(min_large)*m)//s+V*min_large)//(-V*V + 2*V - min_large)
    alphaPrime = (m//s + (m*min_large)//s - min_large)//(V-2*min_large+V*min_large)
    alpha = min(alpha,alphaPrime)
    if(VINT(m,s,alpha)==true)
      alphaPoss[j]=alpha
      j=j+1
    end
  end

  alpha = minimum(alphaPoss)
  if(alpha<1//3)
   alpha=1//3
  end
  return (alpha)


 elseif(V₋₁shares>Vshares)
   #(    )      (    )       (    )
   #𝛂   x      y    1-y      1-x    1-𝛂
   #x = m//s -alpha*(v-1)
   #y = m//s - (1-alpha)*(V-2)
   #m//s >= (1-y)*(V-1-min_large) + (1-alpha)*min_large

   num_large_shares = Vshares
   min_large=Int64(floor(num_large_shares/sᵥ₋₁))
   min_large_og = min_large

   j=1
   for i=0:V-1
     #split the V-1 shares
     min_large=V-1-i;
     alpha = ((m*V)//s -V*V+2*V-2*min_large - (m*min_large)//s + V *min_large-1)//(-V*V+3*V+V*min_large-3*min_large-2)
     alphaPrime = ((m*2)//s-(m*V)//s+V*V-3*V+(2*m*min_large)//s-V*min_large + min_large +2)//(V*V - 3*V +min_large +2)
     alpha = min(alpha,alphaPrime)

     if VINT(m,s,alpha) == true
       alphaPoss[j]=alpha
       j=j+1
     end
   end

   alpha = minimum(alphaPoss)
   if(alpha<1//3)
     alpha=1//3
   end
   return (alpha)
 else
      return 1
 end
 #gap is within the maj-shares
end
