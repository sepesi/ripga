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
   fig = Figure(size=(1000,400)); # initialize the Figure (the trailing ; avoids display of blank Figure)
   ax1 = Axis(fig[1,1],           # specify the display in left column of plots
      autolimitaspect = 1,        # 1:1 is x:y screen length aspect ratio
	  title = "footprint from Bezier curves")
   ax2 = Axis(fig[1,2],           # specify the display in right column of plots
      autolimitaspect = 1,        # 1:1 is x:y screen length aspect ratio
	  title = "infinity symbol from Bezier curves")

   # determine dimension of Euclidean space (i.e., not including ideal dimension)
   nD = trailing_zeros(size(basis,1)) - 1 # -1 discounts ideal dimension
   
   # scale down size of Footprint by 200
   Footprint200 = copy(Footprint)
   Footprint200[1:end-1,:] ./= 200
   
   # if 2D Euclidean space
   if nD == 2
      VRP = bez2poly(Footprint200, 16)
      nVRP = length(VRP)
      for iVRP = 1:nVRP
         P = reshape(VRP[iVRP],:,2)'
         # println("polynomial $iVRP> $P")
         poly!(ax1,             # poly() is for plotting closed polygons
            P,                  # define the coordinates of the points in the closed polygon
            color=:black,       # specify the fill color of the closed polygon
            strokecolor=:black, # specify the border color
            strokewidth=1)      # specify the border width
	  end
      
      VRP = bez2poly(Infinity, 16)
      nVRP = length(VRP)
      for iVRP = 1:nVRP
         P = reshape(VRP[iVRP],:,2)'
         # println("polynomial $iVRP> $P")
         poly!(ax2,             # poly() is for plotting closed polygons
            P,                  # define the coordinates of the points in the closed polygon
            color=:white,       # specify the fill color of the closed polygon
            strokecolor=:black, # specify the border color
            strokewidth=1,      # specify the dotted line outline
            linestyle = :dot)
	  end
        
	  save("./res/Bezier.png", fig)
	  
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