include("Methods/PKG.jl")
using Printf
function menu()
    @printf("\n\n\t MENU: enter 1, 2, 3 or 4")
    @printf("\n1. Given m and s what is the largest possible smallest piece")
    @printf("\n2. Given a range of m and s find all largest possible smallest pieces")
    @printf("\n3. See upper-bound information for a given m, s, and optionally alpha")
    @printf("\n4. Quit\n")
    answer = readline()
    while !(answer == "1" || answer == "2" || answer =="3" || answer == "4")
        println("Enter 1, 2,3 or 4")
        answer = readline()
    end
    return answer
end
function main()
    @printf("\nWelcome to the Muffin Package ")
    @printf("\nThis package is a tool to solve The Muffin Problem: ")
    @printf("\n\n\t\tYou have m muffins and s students. You want to\n
    \t\tdivide the muffins evenly, but no student wants a\n
    \t\ttiny sliver. What division of muffins maximizes\n
    \t\tthe smallest piece?")
    @printf("\n\nNotation used during the program: ")
    @printf("\n\t m : number of muffins")
    @printf("\n\t s : number of students")
    @printf("\n\nAt any user input enter 'q' to go back to the menu")
    println()
    @label start
    answer = menu()
    while answer != "4"
        if answer == "1"
            dual = false
            @printf("Enter m: ")
            @label try_again
            m = readline()
            if m == "q" || m == "Q"
                @goto start
            end
            while true
                try
                    m = parse(Int64,m)
                    break;
                catch
                    @printf("Enter an integer for m: ")
                    m=readline()
                    if m == "q" || m =="Q"
                        @goto start
                    end
                end
            end
            if m<=0
                println("Enter m > 0")
                @goto try_again
            end
            if m>=400
                println("Enter m < 400")
                @goto try_again
            end
            #m=parse(Int64, m)
            @printf("Enter s: ")
            @label try_s
            s = readline()
            if s == "q" || s =="Q"
                @goto start
            end
            while true
                 try
                    s = parse(Int64,s)
                    break;
                catch
                    @printf("Enter an integer for s: ")
                    s=readline()
                    if s == "q" || s =="Q"
                        @goto start
                    end
                end
            end
            if s<=0
                print("Enter s > 0: ")
                @goto try_s
            end
            if s >= 400
                print("Enter s < 400: ")
                @goto try_s
            end
            if m%s == 0
                println()
                println(" m is divisible by s, you can probably do that yourself, try again")
                print("Enter m: ")
                @goto try_again
            end
            if s%m == 0 && s>m
                println()
                println(" s is divisible by m, you can probably do that yourself, try again")
                println(" hint: find f(",s,",",m,") and use the dualtiy thereom: f(s,m) = (s/m)*f(m,s)")
                print("Enter m: ")
                @goto try_again
            end
            if m<=s
                dual = true
                println()
                println("First we solve m = ",s,"  s = ",m,"  then we use the duality thereom: ")
                println("f(s,m) = (s/m)*f(m,s)")
                temp=s
                s=m
                m = temp
            end
            alpha,str,str_eq,time, t_multi, t_solv = FIND_ALPHA(m,s,Inf,Inf,0,0,0,false)
            if str_eq == " TIME OUT "
                str_eq = " ≤ "
            end
            str = str *" NONE "
            valid_methods = split(str)
            leng = length(valid_methods)
            first_letter = Array{String}(undef,0)
            str_first_letter = valid_methods[1][1]
            for i =2:leng
            #    println(valid_methods[i])
                str_first_letter = str_first_letter *" "* valid_methods[i][1]
            end
            first_letter = split(str_first_letter)
            #@printf("\nf(%i,%i) %s %i/%i | Methods: %s  ", m,s,str_eq, numerator(alpha), denominator(alpha), str)
            @label get_method
            if dual
                println()
                println(" \nBy the duality thereom f(",s,",",m,") = ",numerator(s//m*(alpha)),"/",denominator(s//m*(alpha)))
                println(" Enter method for proof [",str,"] (note it is a proof of f(",m,",",s,"))")
            else
                println()
                @printf("\nEnter which method for proof [%s]: ",str)
            end
            method = readline()
            method = uppercase(method)
            if method == "Q"
                @goto start
            end

            while !(method in valid_methods) && !(method in first_letter)
                println("Enter a method from this list: ",str)
                method = readline()
                method = uppercase(method)
                if method == "Q"
                    @goto start
                end
            end
            if method == "H" && "HALF" in valid_methods && "HBM" in valid_methods
                println("HALF or HBM? or enter 1 to go back to methods")
                method = readline()
                if method == "Q" || method == "q"
                    @goto start
                end
                if method == "1"
                    @goto get_method
                end
                method = uppercase(method)
                while !(method == "HALF" || method == "HBM")
                    println("Enter HALF or HBM")
                    method = readline()
                    if method == "1"
                        @goto get_method
                    end
                    if method == "q" || method == "Q"
                        @goto start
                    end
                    method = uppercase(method)
                end
            elseif method == "H" && "HALF" in valid_methods
                method = "HALF"
            elseif method == "H" && "HBM" in valid_methods
                method = "HBM"
            end
            if method != "NONE"
                print_proof(method,m,s,alpha)
                #PROC
                println()
                println("Enter 1 to see procedure or 2 to go back to menu")
                @label  get_see_proc
                see_proc = readline()
                if see_proc == "q" || see_proc == "Q"
                    @goto start
                end
                while see_proc != "1" && see_proc !="2"
                    println("Enter 1 or 2")
                    @goto get_see_proc
                end
                if see_proc == "1"
                    if alpha != 1//3
                        if "GAP" in valid_methods
                            GAP_alpha, GAP_endpoints = GAP(m,s,alpha, true)
                            unique!(GAP_endpoints)
                            proc_bool, err_mess, t_multi, t_solv, length_B = VProc(m,s,alpha, 1000,1000, GAP_endpoints,1)
                        elseif "MID" in valid_methods
                            MID_alpha, MID_endpoints = MID(m,s, alpha, true)
                            unique!(MID_endpoints)
                            proc_bool, err_mess, t_multi, t_solv, length_B = VProc(m,s,alpha, 1000,1000, MID_endpoints,1)

                        else
                            VProc(m,s,alpha,900,900,0,1)
                        end
                    else
                        one_thrd(m,s,1)
                    end
                end
            end
        end
        if answer == "2"
            fileKey = rand(1000:9999)

            println()
            @printf("Enter max m: ")
            @label enter_m
            m = readline()
            if m == "q" || m =="Q"
                @goto start
            end
            while true
                try
                    m = parse(Int64,m)
                    break;
                catch
                    println()
                    @printf("Enter an integer for m: ")
                    m=readline()
                    if m == "q" || m =="Q"
                        @goto start
                    end
                end
            end
            if m <=3
                println()
                println("Enter m > 3")
                @goto enter_m
            end
            if m> 1999
                println()
                println("Enter m < 2000")
                @goto enter_m
            end
            #m=parse(Int64, m)
            @label enter_s
            @printf("Enter min s: ")

            min_s = readline()
            if min_s == "q"  || min_s =="Q"
                @goto start
            end
            while true
                 try
                    min_s = parse(Int64,min_s)
                    break;
                catch
                    println()
                    @printf("Enter an integer for min s: ")
                    min_s=readline()
                    if min_s == "q"  || min_s =="Q"
                        @goto start
                    end
                end
            end
            if min_s<=2
                println()
                @printf("Enter min s > 2: ")
                @goto enter_s
            end
            if min_s >= m
                println()
                @printf("Enter min s < m: ")
                println()
                @goto enter_s
            end
            if min_s >200
                println("Enter min s <= 200")
                @goto enter_s
            end

            @printf("Enter max s: ")
            @label enter_max_s
            max_s = readline()
            if max_s == "q"  || max_s =="Q"
                @goto start
            end
            while true
                 try
                    max_s = parse(Int64,max_s)
                    break;
                catch
                    println()
                    @printf("Enter an integer for max s: ")
                    max_s=readline()
                    if max_s == "q" || max_s =="Q"
                        @goto start
                    end
                end
            end
            if max_s<min_s
                println()
                @printf("Enter max s > min s, enter both again: ")
                println()
                @goto enter_s
            end
            if max_s > 200
                println("Enter s < 200 ")
                @goto enter_max_s
            end
            println()
            println()
            println("Enter the number next to the option to use it")
            println("Options being used are in the right column, press their number to move them back")
            println("You can enter the options in a comma sperated list [ex: 1,4,5]")
            response = "1"
            bool_csv = false
            bool_txt = false
            scott_data = false
            bool_graph = false
            m_le_s2 = false
            stop_at_first = false
            time_limit = 1000
            while response != "0"
                println()
                @label o_menu
                println("Enter '0' to lock in your settings")
                println("-----------Options------------|--------In Use ---------------")
                @printf("%-31s [1] %3d ≤ s ≤ %-3d   s < m ≤ %-4d  "," ",min_s, max_s, m)
                println()
            #    if bool_csv
            #        @printf("%-32s"," ")
            #    end
            #    @printf("[2] Store in csv file")
            #    if bool_csv
            #        @printf(" (run key = %i)",fileKey)
            #    end
            #    println()
                if bool_txt
                    @printf("%-32s"," ")
                end
                @printf("[2] Store in txt file")
                if bool_txt
                    @printf(" (run key = %i)",fileKey)
                end
                println()

                if scott_data
                    @printf("%-31s %s","[3] Use PROC to verify","[3] Use SCOTT to verify")
                else
                    @printf("%-31s %s","[3] Use SCOTT to verify","[3] Use PROC to verify")
                end
                println()
                #if bool_graph
                #    @printf("%-32s"," ")
                #end
                #@printf("[5] Have graphs")
                #if bool_graph
                #    @printf(" (run key = %i)",fileKey)
                #end
                #println()

                @printf("%-31s %s %i sec"," ","[4] Time Limit = ",time_limit)
                println()
                if m_le_s2
                    @printf("%-32s"," ")
                end
                @printf("[5] m ≤ s squared")
                println()
                if stop_at_first
                    @printf("%-32s"," ")
                end
                @printf("[6] stop when one method works")
                println()

                response = readline()

                if response == "q" || response == "Q"
                    @goto start
                end
                r_arr = split(response,",")
                for i=1:length(r_arr)
                    r_arr[i]=strip(r_arr[i])
                end
                if "1" in r_arr
                    println()
                    @printf("Enter max m: ")
                    @label enter_m2
                    m = readline()
                    if m == "q" || m =="Q"
                        @goto start
                    end
                    while true
                        try
                            m = parse(Int64,m)
                            break;
                        catch
                            println()
                            @printf("Enter an integer for m: ")
                            m=readline()
                            if m == "q" || m =="Q"
                                @goto start
                            end
                        end
                    end
                    if m <=3
                        println()
                        println("Enter m > 3")
                        @goto enter_m2
                    end
                    if m> 1999
                        println()
                        println("Enter m < 2000")
                        @goto enter_m2
                    end
                    #m=parse(Int64, m)
                    @label enter_s2
                    @printf("Enter min s: ")

                    min_s = readline()
                    if min_s == "q"  || min_s =="Q"
                        @goto start
                    end
                    while true
                         try
                            min_s = parse(Int64,min_s)
                            break;
                        catch
                            println()
                            @printf("Enter an integer for min s: ")
                            min_s=readline()
                            if min_s == "q"  || min_s =="Q"
                                @goto start
                            end
                        end
                    end
                    if min_s<=2
                        println()
                        @printf("Enter min s > 2: ")
                        @goto enter_s2
                    end
                    if min_s >= m
                        println()
                        @printf("Enter min s < m: ")
                        println()
                        @goto enter_s2
                    end
                    if min_s >200
                        println("Enter min s <= 200")
                        @goto enter_s2
                    end

                    @printf("Enter max s: ")
                    @label enter_max_s2
                    max_s = readline()
                    if max_s == "q"  || max_s =="Q"
                        @goto start
                    end
                    while true
                         try
                            max_s = parse(Int64,max_s)
                            break;
                        catch
                            println()
                            @printf("Enter an integer for max s: ")
                            max_s=readline()
                            if max_s == "q" || max_s =="Q"
                                @goto start
                            end
                        end
                    end
                    if max_s<min_s
                        println()
                        @printf("Enter max s > min s, enter both again: ")
                        println()
                        @goto enter_s2
                    end
                    if max_s > 200
                        println("Enter s < 200 ")
                        @goto enter_max_s2
                    end
                    println()
                end
                if  "2" in r_arr
                    bool_txt = !bool_txt
                end
                if "3" in r_arr
                    scott_data = !scott_data
                end
                if "4" in r_arr

                    println("Enter a new time limit in seconds: ")
                    @label enter_time
                    t_limit = readline()
                     while true
                            try
                                time_limit = parse(Int64,t_limit)
                                break;
                            catch
                                println()
                                @printf("Enter an integer for time limit: ")
                                t_limit=readline()
                                if t_limit == "q" || t_limit =="Q"
                                    @goto start
                                end
                            end
                    end
                    if time_limit < 10
                        println("Enter time limit ≥ 10 sec")
                        @goto enter_time
                    end
                end
                if "5" in r_arr
                    m_le_s2 = !m_le_s2
                end
                if "6" in r_arr
                    stop_at_first = !stop_at_first
                end
                println("\n\n\n\n")
            end

            if false
                println()
                println("Would you like to store the results in a csv file? (yes, no)")
                csv_file = readline()
                csv_file=lowercase(csv_file)
                if csv_file == "q" || csv_file == "Q"
                    @goto start
                end
                while(csv_file != "yes" && csv_file!="no" && csv_file !="y" && csv_file !="n")
                    println()
                    println("Enter yes or no")
                    csv_file = readline()
                    csv_file=lowercase(csv_file)
                    if csv_file == "q" || csv_file == "Q"
                        @goto start
                    end
                end
                println()
                println("Would you like to store the results in a text file? (yes, no)")
                txt_file = readline()
                txt_file=lowercase(txt_file)
                if txt_file == "q" || txt_file == "Q"
                    @goto start
                end
                while(txt_file != "yes" && txt_file!="no" && txt_file!="y"&& txt_file!="n")
                    println()
                    println("Enter yes or no")
                    txt_file = readline()
                    txt_file=lowercase(txt_file)
                    if txt_file == "q" || txt_file == "Q"
                        @goto start
                    end
                end
                println()
                println("Would you like to use SCOTT or PROC to verify? (s, p)")
                verify_type = readline()
                verify_type=lowercase(verify_type)
                if verify_type == "q" || verify_type == "Q"
                    @goto start
                end
                while(verify_type != "scott" && verify_type!="proc" && verify_type != "s" && verify_type !="p")
                    println("Enter SCOTT or PROC")
                    verify_type = readline()
                    verify_type=lowercase(verify_type)
                    if verify_type == "q" || verify_type == "Q"
                        @goto start
                    end
                end
                println()
                bool_graph = false
                graphs = "no"
                #if verify_type != "scott" && verify_type !="s"
                    println("Would you like some graphs? (yes,no)")
                    graphs = readline()
                    graphs=lowercase(graphs)
                    if graphs == "q" || graphs == "Q"
                        @goto start
                    end
                    while(graphs != "yes" && graphs!="no" && graphs!="y"&& graphs!="n")
                        println()
                        println("Enter yes or no")
                        graphs = readline()
                        graphs=lowercase(graphs)
                        if graphs == "q" || graphs == "Q"
                            @goto start
                        end
                    end
                    println()
                #end

                bool_csv = (csv_file == "yes" || csv_file == "y")
                bool_txt = (txt_file == "yes" || txt_file == "y")
                bool_graph = (graphs == "yes" || graphs == "y")
                scott_data = (verify_type == "scott" || verify_type == "s")
            end


                println()
                FIND_ALL(m, min_s, max_s,fileKey, bool_txt, scott_data, bool_csv, bool_graph, time_limit, m_le_s2, stop_at_first )
        end
        if answer == "3"
            @printf("Enter m: ")
            @label tryAgain
            m = readline()
            if m == "q" || m == "Q"
                @goto start
            end
            while true
                try
                    m = parse(Int64,m)
                    break;
                catch
                    @printf("Enter an integer for m: ")
                    m=readline()
                    if m == "q" || m =="Q"
                        @goto start
                    end
                end
            end
            if m<=1
                println("Enter m > 1")
                @goto tryAgain
            end
            if m>=400
                println("Enter m < 400")
                @goto tryAgain
            end
            #m=parse(Int64, m)
            @printf("Enter s: ")
            @label tryS
            s = readline()
            if s == "q" || s =="Q"
                @goto start
            end
            while true
                 try
                    s = parse(Int64,s)
                    break;
                catch
                    @printf("Enter an integer for s: ")
                    s=readline()
                    if s == "q" || s =="Q"
                        @goto start
                    end
                end
            end
            if s<=0
                print("Enter s > 0: ")
                @goto tryS
            end
            if s >= 400
                print("Enter s < 400: ")
                @goto tryS
            end
            if m <= s
                print("Enter m > s: ")
                @goto tryS
            end
            if m%s == 0
                println()
                println(" m is divisible by s, you can probably do that yourself, try again")
                print("Enter m: ")
                @goto tryAgain
            end
            @printf("Enter alpha [x/y] or 0: ")
            @label get_al
            al = readline()
            if al == "q" || al =="Q"
                @goto start
            end
            alpha = 0
            if al != "0"
                while true
                    al_ar = split(al,"/")
                     try
                        n = parse(Int64,al_ar[1])
                        d = parse(Int64,al_ar[2])
                        alpha = n//d
                        break;
                    catch
                        @printf("Enter an integer/integer: ")
                        al=readline()
                        if al == "q" || al =="Q"
                            @goto start
                        end
                    end
                end

                if alpha <= 0 || alpha >=1//2
                    println("Enter alpha between 0 and 1/2: ")
                    @goto get_al
                end
            end

            V,sᵥ,sᵥ₋₁=SV(m,s)
            Vshares=V*sᵥ
            V₋₁shares=(V-1)*sᵥ₋₁
            println()
            println()
            println("    There are ",V,"-students and ",V-1,"-students")
            println("    There are ",sᵥ," ",V,"-students and ",sᵥ₋₁," ",V-1,"-students")
            println()
            println("----------- Please wait while the methods run -----------")
            if al != "0"
               FC_alpha = FC(m,s)
               INT_bool = VINT(m,s,alpha)
               INT_alpha = INT(m,s)
               HALF_bool = VHALF(m,s,alpha)
               HALF_alpha = HALF(m,s)
               EBM_alpha = EBM(m,s)
               HBM_alpha = HBM(m,s)
               MID_bool = VMID(m,s,alpha)
               MID_alpha = MID(m,s)
               GAP_bool = VGAP(m,s,alpha)
               GAP_alpha = GAP(m,s)
               TRAIN_bool = VTRAIN(m,s,alpha)
               TRAIN_alpha = TRAIN(m,s)
               println("  -----  | --- f(",m,", ",s,") ≤ ",al,"----  | --- best upper-bound for f(",m,", ",s,") ---")
               @printf("\n    FC   |           %5s          |          %i/%i       ",FC_alpha == alpha,numerator(FC_alpha),denominator(FC_alpha))
               @printf("\n   INT   |           %5s          |          %i/%i      ",INT_bool,numerator(INT_alpha),denominator(INT_alpha))
               @printf("\n  HALF   |           %5s          |          %i/%i        ",HALF_bool,numerator(HALF_alpha),denominator(HALF_alpha))
               @printf("\n   EBM   |           %5s          |          %i/%i       ",EBM_alpha == alpha,numerator(EBM_alpha),denominator(EBM_alpha))
               @printf("\n   HBM   |           %5s          |          %i/%i       ",HBM_alpha == alpha,numerator(HBM_alpha),denominator(HBM_alpha))
               @printf("\n   MID   |           %5s          |          %i/%i       ",MID_bool,numerator(MID_alpha),denominator(MID_alpha))
               @printf("\n   GAP   |           %5s          |          %i/%i       ",GAP_bool,numerator(GAP_alpha),denominator(GAP_alpha))
               @printf("\n TRAIN   |           %5s          |          %i/%i       ",TRAIN_bool,numerator(TRAIN_alpha),denominator(TRAIN_alpha))
               println()
                @label valid_al
                print_Intervals(m,s,alpha,true)
                while(true)
                    @label get_method2
                    println("  Enter a method to see it's proof or 'none' to go back to menu [FC INT HALF EBM HBM MID GAP TRAIN NONE] ")
                    method = readline()
                    method = uppercase(method)
                    if method == "Q"
                        @goto start
                    end
                    first_letter = ["F","I","H","E","M","G","T","N"]
                    valid_methods = ["FC", "INT", "HALF", "EBM", "HBM", "MID", "GAP", "TRAIN", "NONE"]
                    while !(method in valid_methods)&& !(method in first_letter)
                        println("Enter a method from this list: FC INT HALF EBM HBM MID GAP TRAIN NONE")
                        method = readline()
                        method = uppercase(method)
                        if method == "Q"
                            @goto start
                        end
                    end
                    if method == "H" && "HALF" in valid_methods && "HBM" in valid_methods
                        println("HALF or HBM? or enter 1 to go back to methods")
                        method = readline()
                        if method == "Q"
                            @goto start
                        end
                        if method == "1"
                            @goto get_method2
                        end
                        method = uppercase(method)
                        while !(method == "HALF" || method == "HBM")
                            println("Enter HALF, HBM or enter 1 to go back to methods")
                            method = readline()
                            if method == "Q"
                                @goto start
                            end
                            if method == "1"
                                @goto get_method2
                            end
                            if method == "q" || method == "Q"
                                @goto start
                            end
                            method = uppercase(method)
                        end
                    elseif method == "H" && "HALF" in valid_methods
                        method = "HALF"
                    elseif method == "H" && "HBM" in valid_methods
                        method = "HBM"
                    end
                    if method == "NONE" || method == "N"
                        break
                    end
                    print_proof(method,m,s,alpha)
                end
            else
                FC_alpha = FC(m,s)
                INT_alpha = INT(m,s)
                HALF_alpha = HALF(m,s)
                EBM_alpha = EBM(m,s)
                HBM_alpha = HBM(m,s)
                MID_alpha = MID(m,s)
                GAP_alpha = GAP(m,s)
                TRAIN_alpha = TRAIN(m,s)
                println("  ----- | --- best upper-bound for f(",m,", ",s,") ---")
                @printf("    FC  |          %i/%i       ",numerator(FC_alpha),denominator(FC_alpha))
                @printf("\n   INT  |          %i/%i      ",numerator(INT_alpha),denominator(INT_alpha))
                @printf("\n  HALF  |          %i/%i        ",numerator(HALF_alpha),denominator(HALF_alpha))
                @printf("\n   EBM  |          %i/%i       ",numerator(EBM_alpha),denominator(EBM_alpha))
                @printf("\n   HBM  |          %i/%i       ",numerator(HBM_alpha),denominator(HBM_alpha))
                @printf("\n   MID  |          %i/%i       ",numerator(MID_alpha),denominator(MID_alpha))
                @printf("\n   GAP  |          %i/%i       ",numerator(GAP_alpha),denominator(GAP_alpha))
                @printf("\n TRAIN  |          %i/%i       ",numerator(TRAIN_alpha),denominator(TRAIN_alpha))
                println()
                println()
                while(true)
                    @label get_method3
                    println("  Enter a method to see it's proof or 'none' to go back to menu [FC INT HALF EBM HBM MID GAP TRAIN NONE] ")
                    method = readline()
                    method = uppercase(method)
                    if method == "Q"
                        @goto start
                    end
                    first_letter = ["F","I","H","E","M","G","T","N"]
                    valid_methods = ["FC", "INT", "HALF", "EBM", "HBM", "MID", "GAP", "TRAIN", "NONE"]
                    while !(method in valid_methods)&& !(method in first_letter)
                        println("Enter a method from this list: FC INT HALF EBM HBM MID GAP TRAIN NONE")
                        method = readline()
                        method = uppercase(method)
                        if method == "Q"
                            @goto start
                        end
                    end
                    if method == "H" && "HALF" in valid_methods && "HBM" in valid_methods
                        println("HALF or HBM? or enter 1 to go back to methods")
                        method = readline()
                        if method == "Q"
                            @goto start
                        end
                        if method == "1"
                            @goto get_method3
                        end
                        method = uppercase(method)
                        while !(method == "HALF" || method == "HBM")
                            println("Enter HALF or HBM")
                            method = readline()
                            if method == "Q"
                                @goto start
                            end
                            if method == "1"
                                @goto get_method3
                            end
                            if method == "q" || method == "Q"
                                @goto start
                            end
                            method = uppercase(method)
                        end
                    elseif method == "H" && "HALF" in valid_methods
                        method = "HALF"
                    elseif method == "H" && "HBM" in valid_methods
                        method = "HBM"
                    end
                    if method == "NONE" || method == "N"
                        break
                    end

                    if method == "FC" || method == "F"
                        alpha = FC_alpha
                    elseif method == "INT" || method == "I"
                         alpha = INT_alpha
                    elseif method == "HALF"
                         alpha = HALF_alpha
                    elseif method ==  "MID" || method == "M"
                         alpha =MID_alpha
                    elseif method ==  "EBM" ||method ==  "E"
                         alpha = EBM_alpha
                    elseif method ==  "HBM"
                        alpha = HBM_alpha
                    elseif method ==  "TRAIN" ||method ==  "T"
                         alpha = TRAIN_alpha
                    elseif method ==  "GAP" ||method ==  "G"
                        alpha = GAP_alpha
                    end
                    print_proof(method,m,s,alpha)
                end
            end
        end
            println()

        answer=menu()
    end
end

function print_proof(method, m, s, alpha)
    if method == "FC" || method == "F"
        println("\nFC proof of upperbound")
        print("***********************************************")
        FC(m,s,true)
        println("\n***********************************************")
    elseif method == "INT"|| method == "I"
        println("\nINT proof of upperbound")
        print("***********************************************")
        VINT(m,s,alpha,true)
        println("***********************************************")
    elseif method == "HALF"|| method == "H"
        println("\nHALF proof of upperbound")
        println("***********************************************")
        Half_proof(m,s,alpha)
        println("***********************************************")
    elseif method == "MID"|| method == "M"
        println("\nMID proof of upperbound")
        println("***********************************************")
        VMID(m,s,alpha,true)
        println("***********************************************")
    elseif method == "EBM"|| method == "E"
        println("\nEBM proof of upperbound")
        println("***********************************************")
        EBM(m,s,true)
        println("***********************************************")
    elseif method == "HBM"|| method == "H"
        println("\nHBM proof of upperbound")
        println("***********************************************")
        HBM(m,s,true)
        println("***********************************************")
    elseif method == "GAP"|| method == "G"
        println("Would you like to see the matrices? (yes or no)")
        mat = readline()
        mat = lowercase(mat)
        while  mat != "yes" && mat!="no" && mat!="y"&& mat!="n"
            println("Enter yes or no")
            mat=readline()
            mat=lowercase(mat)
        end
        if mat == "yes" || mat =="y"
            gap_proof = 2
        else
            gap_proof = 1
        end
        println("Would you like to save the proof to a file? (y,n)")
        f = readline()
        f = lowercase(f)
        while  f != "yes" && f!="no" && f!="y"&& f!="n"
            println("Enter yes or no")
            f=readline()
            f=lowercase(f)
        end
        if f == "yes" || f =="y"
            if gap_proof == 2
                gap_proof = 4
            else
                gap_proof = 3
            end
        end
        println("\nGAP proof of upperbound")
        println("***********************************************")
        VGAP(m,s,alpha,gap_proof)
        println("***********************************************")
    elseif method == "TRAIN"||method == "T"
        println("Would you like to save the proof to a file? (y,n)")
        f = readline()
        f = lowercase(f)
        while  f != "yes" && f!="no" && f!="y"&& f!="n"
            println("Enter yes or no")
            f=readline()
            f=lowercase(f)
        end
        if f == "yes" || f =="y"
            file = 3
        else
            file= 1
        end
        println("\nTrain proof of upperbound")
        println("***********************************************")
        VTRAIN(m,s,alpha,file)
        println("***********************************************")
    end
end
main()
