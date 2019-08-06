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
<ul>
  <li>[JuMP](http://www.juliaopt.org/JuMP.jl/v0.17/installation.html) - a modeling language used to support the solvers</li>
  <li>[Cbc](https://github.com/JuliaOpt/Cbc.jl)</li>
  <li>[GLPK](https://github.com/JuliaOpt/GLPK.jl)</li>
</ul>
Data formatting and file I/O: 
<ul>
  <li>Printf</li>
  <li>Dates</li> 
  <li>Plots</li>
  <li>CSV</li>
  <li>DataFrames</li>
 </ul>
