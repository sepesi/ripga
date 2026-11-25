# ripgaplot.jl : plotting functions for
#   reference implementation of
#   Projective Geometric Algebra image test
#
# ripga is a Julia port of bivector.net's C++ reference
# implementation of projective geometric algebra
# available at https://bivector.net/tools.html
#

# bezlerp : Bezier Linear Interpolation
#
# arguments:
#  BCP1: Bezier Control Point matrix (each column is Euclidean dimension of Bezier curve Control Points)
#  TSM: Time Sample Matrix (stored in a Matrix{Float32}) for linear interpolation
#
# return value:
#  PC: Polynomial Coordinates (stored in a Matrix{Float32})
#
function bezlerp(BCP1::Matrix{Float32},TSM::Matrix{Float32})

   # Allocate buffer for return of Polynomial Coordinates in a Matrix{Float32}.
   nCoord = size(BCP1,2)
   nTS = size(TSM,2)
   PC = Matrix{Float32}(undef,nTS,nCoord)

   # for each coordinate in polynomial
   for iCoord = 1:nCoord
      PC[:,iCoord] = sum(TSM .* BCP1[:,iCoord],dims=1)'
   end

   return PC # each column of matrix holds a Euclidean dimension of polynomial coordinates
end # bezlerp()

# bez2poly : convert cubic Bezier curves into
#            rolled vectors of polynomial coordinates
#
# arguments:
#  BCP: Bezier Control Point matrix where each column contains all control points of Bezier curve
#  nT: number of sampling periods in each Bezier curve segment
#
# return value:
#  Vector{Vector{Float32}} where each "column" is a set of polynomial
#    coordinates rolled into a single vector
#
function bez2poly(BCP::Matrix{Float32}, nT::Int)

   # find start point and length of each contiguous Bezier curve
   Status = trunc.(Int16, BCP[end,:])
   Istart = findall(x -> (x & 1) == 1, Status)
   push!(Istart, length(Status)+1) # append data sentinel
   Ilength = diff(Istart)
   pop!(Istart) # remove data sentinel

   # initialize data
   TS = range(0f0, step=1f0/nT, length=nT) # define Time Samples
   TSM = Matrix{Float32}([(1 .- TS).^3 (3*TS) .* (1 .- TS).^2 3*TS.^2 .* (1 .- TS) TS.^3]') # TS Matrix
   VRP = Vector{Vector{Float32}}()   # initialize Vector of Rolled Polynomials
   nCoord = (size(BCP,1) & ~1) >> 2  # discount appended status field, divide by # of control points (4)
   # println("TSM: $TSM")

   # for each set of closed polynomials in BCP
   nPoly = length(Istart)
   for iPoly = 1:nPoly
     # println("iPoly $iPoly> start at $(Istart[iPoly]) with length $(Ilength[iPoly])")
     PC = Matrix{Float32}(undef, 0, 2) # initialize Polynomial Coordinates

     # for each cubic Bezier curve within the closed polynomial
     nCol = Ilength[iPoly]
     for iCol = 1:nCol
        # linear interpolate between control points of single Bezier curve (i.e., single column)
        iBez1 = Istart[iPoly] + iCol - 1
        BCP1 = reshape(BCP[1:end-1,iBez1],:,nCoord)
        # println("iBez1 $iBez1; iCol $iCol> $BCP1")
        PC1 = bezlerp(BCP1, TSM)
        # println("Bez curve $iBez1 matrix> $PC1")
        PC = vcat(PC, PC1)
     end

     # roll up PC (Polynomial Coordinates) into a single vector and push it onto VRP
     push!(VRP, PC[:])
   end

   return VRP
end # bez2poly()

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
end # itest()