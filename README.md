# The Muffin Package
**The Muffin Puzzle**, first proposed by Alan Frank, is as follows.
*You have 5 muffins and 3 students. You want to divide the muffins evenly, but no student wants a tiny sliver. What division of muffins maximizes the smalles piece?*

This package goes with the book *The Mathematics of Muffins* by William Gasarch. Using the methods described by Gasarch and his team of researchs we creatd a package to find the largest possible smallist piece given any number of muffins and students.  


# Requirements
To run this package you need [Julia](https://julialang.org/downloads/).
To run a Julia file open the Julia REPL and type **include(*file_path*)**.
# Installation
Download and unzip all the files. Run *setup.jl* to download all the needed packages. To use the program run *user_interface.jl*.

# Outside Packages Used
Solvers:

  <p>[JuMP](http://www.juliaopt.org/JuMP.jl/v0.17/installation.html) - a modeling language used to support the solvers</p>

  <p>[Cbc](https://github.com/JuliaOpt/Cbc.jl) and [GLPK](https://github.com/JuliaOpt/GLPK.jl)</p>

Data formatting and file I/O: 

  <p>Printf, Dates, Plots, CSV, DataFrames </p>
