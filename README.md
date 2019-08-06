# The Muffin Package
**The Muffin Puzzle**, first proposed by Alan Frank, is as follows.
*You have 5 muffins and 3 students. You want to divide the muffins evenly, but no student wants a tiny sliver. What division of muffins maximizes the smalles piece?*

This package goes with the book *The Mathematics of Muffins* by William Gasarch. Using the methods described by Gasarch and his team of researchs we creatd a package to find the largest possible smallist piece given any number of muffins and students.  

The methods folder holds all upper-bound methods as well as FIND_PROC, one third lower bound, and a file name "PKG" which incorporates the upper and lower bound methods together to produce the largest possible smallest piece. The Data folder is empty, it is where the user interface will save its files. 

# Requirements
To run this package you need [Julia](https://julialang.org/downloads/).
To run a Julia file open the Julia REPL and type **include(*file_path*)**.
# Installation
Download and unzip all the files. Run *setup.jl* to download all the needed packages. To use the program run *user_interface.jl*.

# Outside Packages Used

JuMP - a modeling language 
@article{DunningHuchetteLubin2017,
author = {Iain Dunning and Joey Huchette and Miles Lubin},
title = {JuMP: A Modeling Language for Mathematical Optimization},
journal = {SIAM Review},
volume = {59},
number = {2},
pages = {295-320},
year = {2017},
doi = {10.1137/15M1020575},
}
Solvers: [Cbc](https://github.com/JuliaOpt/Cbc.jl) and [GLPK](https://github.com/JuliaOpt/GLPK.jl)

Data formatting and file I/O: Printf, Dates, Plots, CSV, DataFrames 
