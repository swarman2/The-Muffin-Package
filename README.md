# The Muffin Package
**The Muffin Puzzle**, first proposed by Alan Frank, is as follows.
*You have 5 muffins and 3 students. You want to divide the muffins evenly, but no student wants a tiny sliver. What division of muffins maximizes the smalles piece?*

This package goes with the book *The Mathematics of Muffins* by William Gasarch. Using the methods described by Gasarch and his team of researchs we creatd a package to find the largest possible smallist piece given any number of muffins and students.  


# Requirements
To run this package you need [Julia](https://julialang.org/downloads/) and some version of [Python](https://www.python.org/downloads/).

# Installation
Download and unzip all the files.To run a Julia file open the Julia REPL and type **include(*file_path*)**. Run *setup.jl* to download all the needed packages. To use the program run *user_interface.jl*.

# How to use the package
The menu will look like this

                      MENU: enter 1, 2, 3 or 4
            1. Given m and s what is the largest possible smallest piece
            2. Given a range of m and s find all largest possible smallest pieces
            3. See the proof for a specified m,s, and alpha
            4. Quit

#### Basic functionality:
   * Any user input is not case sensitive
   * To get back to the menu enter "q" at an input prompt
#### Option 1: Given m and s what is the largest possible smallest piece
The following prompts will get you to enter an amount of muffins (m) and an amount of students (s)

* Constraints on m and s: 
* m cannot be divisible by s
* m has to be less than 400
* s has to be less than 400  

Once m and s have been answered the program will calculated the largest possible smallest piece and output something like this:

                  5     |  3     |   FC  INT  HALF  MID  GAP   |  =      5/12    | 0.00 sec
                  
The format is:

* m     |  s     |   methods that provided an upper-bound of alpha   | = largest possible smallest piece    | runtime

The equal will be a less than or equal if the largest possible smallest piece was not found. 

Then you will see this prompt:
                                   
                          Enter which method for proof [ FC  INT  HALF  MID  GAP  NONE ]: 
                          
Type in which method you would like to see the proof for, you can type in the whole method or just the first letter. 

Then this prompt will show up:

                        Enter 1 to see procedure or 2 to go back to menu      
                        
If you press 1 a procedure will be printed 

*warning: this will take roughly the same amount of time as previous step*   

If you enter m less than s, it will calculate f(s,m) and then use the duality thereom. 
<!---
**Example run   m = 11   s = 5**

            Enter m: 11
            Enter s: 5
              11    |  5     |   INT  HALF  MID  GAP   |  =     13/30    | 0.01 sec

            Enter which method for proof [ INT  HALF  MID  GAP  NONE ]: m

            MID proof of upperbound
            ***********************************************

            m  = 11  s = 5
            2 5-students    3 4-students    10 5-shares     12 4-shares
            Numbers assumed to have denominator: 30

            SPLIT THE 4 SHARES

                 (    10   5-shares   )    |    (  2   small  4-shares )         (  10  large  4-shares  )
                13                   14        15                     15        16                      17

            SPLIT THE 4 SHARES AGAIN

                 (   1    4-shares    |    1    4-shares   )         (  10  large  4-shares  )
                15                   15                   15        16                      17

            Possible muffin distributions
            ( 3 3 3 3 )
            System of equations = [1, 1, 10, 3]
            4×1 Array{Int64,2}:
             0
             0
             4
             1
            No solution on the Naturals
            alpha ≤ 13/30
            ***********************************************

            Enter 1 to see procedure or 2 to go back to menu
            1

            Procedure for f(11, 5) = 13/30
            ***********************************************
            All numbers assumed to have denominator: 30
            Cut 8 muffins {  13  17  }
            Cut 2 muffins {  14  16  }
            Cut 1 muffins {  15  15  }

            Give 2 students {  13  13  13  13  14  }
            Give 2 students {  15  17  17  17  }
            Give 1 students {  16  16  17  17  }
            ***********************************************
            --->
#### Option 2: Given a range of m and s find all largest possible smallest pieces
  The goal of this option is to output a range of largest smallest possible pieces with parameters:
  
  
  *s < m ≤ max m         min s ≤ s ≤ max s*    
  
  The following prompts will ask you to enter max_m, min_s and max_s. 
   
              Enter max m: 
              Enter min s: 
              Enter max s: 

  
  Then a settings menu will appear that will allow you to set the settings for your run.
  
            Enter '0' to lock in your settings
            -----------Options------------|--------In Use ---------------
                                            [1]  30 ≤ s ≤ 40    s < m ≤ 50
            [2] Store in txt file
            [3] Use SCOTT to verify         [3] Use PROC to verify
                                            [4] Time Limit =  1000 sec
            [5] m ≤ s squared
            [6] stop when one method works

  
  * You can also enter a comma seperated list of the settings you want (ex: 1,2,4)
 
  * Any files created will be saved in the DATA folder.
  
  * A new folder specific to the run will be created. It's name will end in a "run key" which is a random three digit number. 
  
  Once you have locked in your settings walk aways have a cup of coffee.
  
  When the run finishes some stats are printed in addition to the data
  
* runtime
* the amount of times each method was used
* amount of times none of the methods worked (and if that was due to it timing out or not)
* the amount of times the methods overlapped with each other
* the percent of data solved after adding in a method
  
 *stats are only are saved if you chose to have a txt file*
 <!---
 **Example run  20 ≤ s ≤ 25    s < m ≤ 30**
 
 
             Enter '0' to lock in your settings
            -----------Options------------|--------In Use ---------------
                                            [1]  20 ≤ s ≤ 25    s < m ≤ 30
            [2] Store in txt file
            [3] Use SCOTT to verify         [3] Use PROC to verify
                                            [4] Time Limit =  1000 sec
            [5] m ≤ s squared
            [6] stop when one method works
            0


               m    |   s    |                   Method(s)                  |        α        |       runtime
            ------------------------------------------------------------------------------------------------------
              21    |  20    |   FC  INT                  MID  GAP          |  =      7/20    | 0.23 sec
              23    |  20    |       INT  HALF  EBM  HBM  MID  GAP          |  =      7/20    | 0.04 sec
              27    |  20    |   FC             EBM                         |  =      1/3     | 0.02 sec
              29    |  20    |   FC             EBM                         |  =      1/3     | 0.02 sec
              22    |  21    |                  EBM                         |  =      1/3     | 0.01 sec
              23    |  21    |                  EBM  HBM  MID  GAP          |  =     29/84    | 0.11 sec
              25    |  21    |       INT  HALF  EBM                         |  =      1/3     | 0.01 sec
              26    |  21    |       INT             HBM  MID  GAP          |  =     22/63    | 0.03 sec
              29    |  21    |   FC             EBM                         |  =      1/3     | 0.01 sec
              23    |  22    |                  EBM  HBM       GAP          |  =     15/44    | 0.23 sec
              25    |  22    |                       HBM       GAP          |  =     23/66    | 0.07 sec
              27    |  22    |       INT  HALF  EBM  HBM  MID  GAP          |  =      4/11    | 0.33 sec
              29    |  22    |   FC  INT        EBM       MID  GAP          |  =     15/44    | 0.27 sec
              24    |  23    |   FC  INT                  MID  GAP          |  =      8/23    | 0.11 sec
              25    |  23    |                  EBM                         |  =      1/3     | 0.01 sec
              26    |  23    |       INT        EBM       MID               |  =      8/23    | 0.04 sec
              27    |  23    |       INT  HALF  EBM                         |  =      1/3     | 0.01 sec
              28    |  23    |       INT  HALF  EBM       MID  GAP          |  =     33/92    | 0.07 sec
              29    |  23    |                                 GAP          |  =     49/138   | 0.05 sec
              30    |  23    |   FC  INT        EBM       MID  GAP          |  =      8/23    | 0.09 sec
              25    |  24    |                  EBM                         |  =      1/3     | 0.01 sec
              29    |  24    |       INT  HALF  EBM       MID  GAP          |  =     17/48    | 0.04 sec
              26    |  25    |                  EBM  HBM       GAP          |  =     17/50    | 0.32 sec
              27    |  25    |                  EBM       MID  GAP          |  =     17/50    | 0.15 sec
              28    |  25    |       INT        EBM                         |  =      1/3     | 0.01 sec
              29    |  25    |       INT  HALF  EBM       MID  GAP          |  =     17/50    | 0.15 sec

            ---------------- STATS ----------------
            Total time:  00:00:2.47

            Amount of times each method produced the correct alpha
            FC: 7/26 ==> 26.92 %
            HALF: 7/26 ==> 26.92 %
            INT: 14/26 == >53.85 %
            MID: 13/26 ==> 50.00 %
            EBM: 21/26 ==> 80.77 %
            HBM: 7/26 ==> 26.92 %
            GAP:16/26 ==>  61.54 %
            TRAIN: 0/26 ==> 0.00 %

            Amount of times the correct alpha was not found
            TIME OUT [1000 sec]: 0/26 ==> 0.00 %
            Incorrect upper-bound: 0/26 ==> 0.00 %


            Amount of times methods overlapped:
             FC  INT                  MID  GAP         ==> 2/26.0 = 7.69%
                            EBM  HBM       GAP         ==> 2/26.0 = 7.69%
                            EBM  HBM  MID  GAP         ==> 1/26.0 = 3.85%
                 INT        EBM                        ==> 1/26.0 = 3.85%
             FC             EBM                        ==> 3/26.0 = 11.54%
                 INT  HALF  EBM                        ==> 2/26.0 = 7.69%
                 INT  HALF  EBM  HBM  MID  GAP         ==> 2/26.0 = 7.69%
                 INT        EBM       MID              ==> 1/26.0 = 3.85%
                                 HBM       GAP         ==> 1/26.0 = 3.85%
             FC  INT        EBM       MID  GAP         ==> 2/26.0 = 7.69%
                            EBM       MID  GAP         ==> 1/26.0 = 3.85%
                 INT             HBM  MID  GAP         ==> 1/26.0 = 3.85%
                            EBM                        ==> 3/26.0 = 11.54%
                                           GAP         ==> 1/26.0 = 3.85%
                 INT  HALF  EBM       MID  GAP         ==> 3/26.0 = 11.54%

            Percentage of correct alphas found adding in methods one by one:
            FC    =>    7/26.0 = 26.92%

            FC HALF    =>    14/26.0 = 53.85%

            FC HALF INT    =>    17/26.0 = 65.38%

            FC HALF INT MID    =>    19/26.0 = 73.08%

            FC HALF INT MID EBM    =>    24/26.0 = 92.31%

            FC HALF INT MID EBM HBM    =>    25/26.0 = 96.15%

            FC HALF INT MID EBM HBM GAP    =>    26/26.0 = 100.0%

            FC HALF INT MID EBM HBM GAP TRAIN    =>    26/26.0 = 100.0%

--->
 
#### Option 3: See the proof for a specified m,s, and alpha

You will be prompted to enter m, s and alpha (**enter alpha as an integer/integer**). The program will then print information regarding that m, s and alpha (ex: what types of students there are). Then you can select a method to see how it approaches the problem. You can keep entering methods or enter "none" or "q" to get back to the menu. 
<!---
**Example Run   m = 5   s = 3   alpha = 5/12**

            Enter m: 5
            Enter s: 3
            Enter alpha [x/y] or 0: 5/12


                There are 4-students and 3-students
                There are 1 4-students and 2 3-students

            m  = 5  s = 3
            1 4-students    2 3-students    4 4-shares      6 3-shares
            Numbers assumed to have denominator: 12

            SPLIT THE 3 SHARES

                 (    4    4-shares   )    |    (  2   small  3-shares )         (  4   large  3-shares  )
                5                    5         6                      6         7                       7

            SPLIT THE 3 SHARES AGAIN

                 (   1    3-shares    |    1    3-shares   )         (  4   large  3-shares  )
                6                    6                    6         7                       7

              Enter a method to see it's proof or 'none' to go back to menu [FC INT HALF EBM HBM MID GAP TRAIN NONE]
            f

            FC proof of upperbound
            ***********************************************
            max{ 1/3, min{5/3 * 1/4 , 1- 5/3 * 1/3}
            = max{1/3, 5/12}
            = 5/12
            ***********************************************
              Enter a method to see it's proof or 'none' to go back to menu [FC INT HALF EBM HBM MID GAP TRAIN NONE]
            i

            INT proof of upperbound
            ***********************************************
            m  = 5  s = 3
            1 4-students    2 3-students    4 4-shares      6 3-shares
            Numbers assumed to have denominator: 12

            SPLIT THE 3 SHARES

                 (    4    4-shares   )    |    (  2   small  3-shares )         (  4   large  3-shares  )
                5                    5         6                      6         7                       7

            0 small shares and 3 large shares works
            Need at least 3 large shares
            2 students need at least 3 large shares, but there are only 4 large shares, so alpha ≤ 5/12
--->

# Outside Packages Used
Solvers:
<ul>
  <li> JuMP - a modeling language used to support the solvers</li>
  <li>Cbc</li>
  <li>GLPK</li>
</ul>
Data formatting and file I/O: 
<ul>
  <li>Printf</li>
  <li>Dates</li> 
  <li>Plots</li>
  <li>CSV</li>
  <li>DataFrames</li>
 </ul>
