include("HALF.jl")
println("Enter: m,s,𝛂")
str=readline()#read in a comma seperated list
strArr=split(str,",")#split the list on the commas (removes the commas)
m=parse(Int64,strArr[1])#assign m to the first thing in the list
s=parse(Int64,strArr[2])#s to the second
alphaArr=split(strArr[3],"/")
num=parse(Int64, alphaArr[1])
denom=parse(Int64,alphaArr[2])
alpha = num//denom
#alpha=parse(Rational,strArr[3])  #attempt to read in a rational number for alpha
#
#treat rationals as pair of integers (reduced)
#
#alpha=32//78#unfortunatly hardcoded in (how do you read in a rational number??)
if alpha >1//2
  println("𝛂 > 1//2 is a bad guess, please rethink")
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
println("if Alice has ",V+1, " or more shares then there exists a share ≦ ", share)
#if the share that Alice has is less than our hypthoisized alpha than we know 
#Alice has less than five shares
less_than_upper = share<alpha #boolean that is true if Alice has less than five shares
println("We will only consider the case where Alice has less than ", V+1," shares if share is <  𝛂")
println("share < 𝛂: ", less_than_upper)
println()
println("Case 2:")
#if Bob has less than two shares then one of them has to be greater than 
#1/2 times amt Bob has and its buddy is less than 1-(1/2*amt Bob has)
share = m//(s *(V-2))
buddy=(1-share)

more_than_lower = buddy<alpha
print("if Bob has ", V-2," or less shares then there exists a share > ", share)
println(" whose buddy is < ", buddy)
println("We will only consider the case where Bob has more than ",V-2," shares if buddy < 𝛂")
println("buddy < 𝛂: ",more_than_lower)
println()
if(less_than_upper && more_than_lower)
  #solve for s3 and s4
  #
  #(V-1)s3+Vs4=pieces
  #s3+s4=students
  #
  #may want to change subscripts on s
  println(V-1,"sᵥ₋₁ + ",V,"sᵥ = ", pieces)
  println("sᵥ₋₁ + sᵥ = ", s)
#  A=[V-1 V; 1 1]; B=[pieces; s];
#  X=A\B #solve the system AX=B
  #These equations were found by manually solving the above matrix equation 
  s3=s*V-pieces
  s4=pieces-s*V+s
#  println("sᵥ₋₁: ", s3," sᵥ: ", s4)
  println("There are ", s3, " ",V-1,"-students, ", s4," ",V,"-students, ",(V-1)*s3, " ",V-1,"-shares, and ", V*s4, " ",V,"-shares.")
  
  println("Endpoints are x = ",x," and y = ",y) 
  #Alice has a 4-share > x (x=m/s - 3𝛂)
  println()
  println("Case 3:")
  #x=Rational(m/s-(3*alpha))
  #x=(m-(V-1)*alpha*s)//s
  

  println("Alice has a ", V,"-share > ",x)
 # temp = Rational((m//s)-x)#this is an intermediate value used frequently
  temp=(m-x*s)//s 
 println("Alice's other ", V-1," ",V, "-shares add up to < ", m//s, " - ", x," = ", temp)
  print("so one of Alice's shares must be < ", temp, " * 1/",V-1," = ")
  temp=numerator(temp)//(denominator(temp)*(V-1))
  println(temp)
  if(temp == alpha)
    println("so one of Alice's shares must be less than alpha")
  elseif(temp<alpha)
    println(temp," < ",alpha,"  so one of Alice's shares must be less than alpha")
  end 

  println()
  println("Case 4:")
  #Bob has a 3-share < y (y=m/s - 2*(1-alpha))
   #y=Rational(m//s - (V-2)*(1-alpha))
  #y=(m-(V-2)*(1-alpha)*s)//s
 #temp = Rational((m//s)-y)
 temp=(m-y*s)//s
  println("Bob has a ", V-1,"-share < ",y)
  println("Bob's other ",V-2," ",V-1,"-shares add up to > ",m//s, " - ", y, " = ", temp)
  print("so one of Bob's shares must be > ",temp,"* 1/",V-2," =")
  temp=numerator(temp)//(denominator(temp)*(V-2))
  println(temp)
  println("whose buddy would be < 1-",temp," = ", Rational(1-temp))
  if(1-temp == alpha)
    println("so one of Bob's shares must be less than alpha")
  elseif((1-temp)<alpha)
    println(1-temp," < ",alpha,"  so one of Bob's shares must be less than alpha")
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