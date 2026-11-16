# ripgaitest.jl : reference implementation of
# Projective Geometric Algebra image test
#
# ripga is a Julia port of bivector.net's C++ reference
# implementation of projective geometric algebra
# available at https://bivector.net/tools.html
#
# ripgaitest is an integration test of ripga and Makie.
#
using CairoMakie

# image test
#
function itest()
   f = Figure();          # initialize the Figure (the trailing ; avoids display of blank Figure)
   Axis(f[1,1])           # specify the display should have only one row and one column of plots

   # determine dimension of Euclidean space (i.e., not including ideal dimension)
   nD = trailing_zeros(size(basis,1)) - 1 # -1 discounts ideal dimension
   
   # if 2D Euclidean space
   if nD == 2
      bez2poly(Footprint, 10)
	  
   # else if 3D Euclidean space
   elseif nD == 3
      println("no 3D itest() available yet")
	  
   # else if 4D Euclidean space
   elseif nD == 4
      println("no 4D itest() available yet")
	  
   end
#=
   pts = Float32[0f0 2f0 3f0 1f0;0f0 0f0 1f0 1f0]
   poly(                  # poly() is for plotting closed polygons
      pts,                # define the coordinates of the points in the closed polygon
      color=:red,         # specify the fill color of the closed polygon
      strokecolor=:black, # specify the border color
      strokewidth=1)      # specify the border width
=#
end # ripga2d itest()