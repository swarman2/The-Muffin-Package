function FC(m,s)
    if m%s == 0
      return 1
    end
    min_val = rationalize(min(m/(s*ceil(2*m/s)),1-(m/(s*floor(2*m/s)))))
  #  if min_val ==1//3
  #    println("m: ",m,"  s: ",s)
  #  end
return  rationalize(max(1/3, min_val))#min(m/(s*ceil(2*m/s)),1-(m/(s*floor(2*m/s))))))
end

#FC(23,13)
