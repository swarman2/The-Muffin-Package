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

            MENU: enter 1, 2 or 3
    1. Given m and s what is the largest possible smallest piece
    2. Given a range of m and s find all largest possible smallest pieces
    3. Quit

**Basic functionality:**

   * Any user input is not case sensitive
   * To get back to the menu enter "q" at an input prompt

**Option 1 : Given m and s what is the largest possible smallest piece**
    The following prompts will get you to enter an amount of muffins (m) and an amount of students (s)
    * Constraints on m and s: 
        * m cannot be divisible by s
        * m has to be less than 400
        * s has to be less than 400  
     Once m and s have been answered the program will calculated the largest possible smallest piece and output something like this: 
        5     |  3     |   FC  INT  HALF  MID  GAP   |  =      5/12    | 0.00 sec
      The format is 
        m     |  s     |   methods that provided an upper-bound of alpha   | = largest possible smallest piece    | runtime
      The equal will be a less than or equal if the largest possible smallest piece was not found. 
      Then you will see this prompt:
          Enter which method for proof [ FC  INT  HALF  MID  GAP  NONE ]: 
      Type in which method you would like to see the proof for, you can type in the whole method or just the first letter. 
      Then this prompt will show up:
              Enter 1 to see procedure or 2 to go back to menu          
      If you press 1 a procedure will be printed *warning: this will take roughly the same amount of time as previous step*   
      If you enter m less than s, it will calculate f(s,m) and then use the duality thereom. 
**Option 2: Given a range of m and s find all largest possible smallest pieces**
  The goal of this option is to output a range of largest smallest possible pieces with parameters:
                    *s < m $\le$ max m         min s $\le$ s $\le$ max s*    
  The following prompts will ask you to enter max_m, min_s and max_s. 
  
  Then a settings menu will appear that will allow you to set the settings for your run.
    Enter 1 to change the range
    Enter 2 to store in a csv file
    Enter 3 to store in a txt file
    Enter 4 to switch between PROC being used to verify and SCOTT being used to verify
    Enter 5 to change the time limit
    Enter 6 to only go up to m = s squared for each s
  Enter 0 to lock in your settings
  
  You can also enter a comma seperated list of the settings you want (ex: 1,2,4)
 
  Any files created will be saved in the DATA folder.
  A new folder specific to the run will be created. It's name will end in a "run key" which is a random three digit number. 
  
  Once you have locked in your settings walk aways have a cup of coffee.\\
  When the run finishes some stats are printed in addition to the data\\
  <ul>
    <li>runtime<\li>
    <li>the amount of times each method was used\<\li>
    <li>amount of times none of the methods worked (and if that was due to it timing out or not)<\li>
    <li>the amount of times the methods overlapped with each other<\li>
    <li>the percent of data solved after adding in a method<\li>
  <\ul>
 **stats are only are saved if you chose to have a txt file*
    
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
