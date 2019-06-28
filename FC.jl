function FC(m,s)
    if m%s == 0
      return 1
    end
return  rationalize(max(1/3, min(m/(s*ceil(2*m/s)),1-(m/(s*floor(2*m/s))))))
end

#FC(23,13)
