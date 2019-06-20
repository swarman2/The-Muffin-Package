while true
  println("Enter: m,s")
  str=readline()
  strArr=split(str,",")
  if str == "quit"
    break
  end
  m=parse(Int64,strArr[1])
  s=parse(Int64,strArr[2])
  println("f(",m,",",s,") = ", rationalize(max(1/3, min(m/(s*ceil(2*m/s)),1-(m/(s*floor(2*m/s)))))))
end
#println(rationalize(max(1/3,min(5/(3*ceil(2*5/3)),1-(5/(3*floor(2*5/3)))))))
