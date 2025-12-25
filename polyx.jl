# polyx.jl
# polygon intersection testing
# The Separating Axis Theorem (SAT) applies to only
# convex polygons. However, any non-convex polygon
# can be partitioned into convex polygons, as done
# in this example.
#
# NOTE:
# In this example, the vertices of each polygon should 
# be indexed in the counter clockwise direction (i.e.,
# according to the right hand rule).
#
# NOTE2:
# This implementation of sat() is efficient in that it
# does not require the intermediate storage of all the 
# projected distances and then a comparison of the 
# minimum and maximum values of those projected distances.
# Instead, because of the choice of the origin and the
# ordering of the polygon vertices, the projected 
# distances of the vertices within the polygon containing
# the origin are guaranteed to be less than or equal to
# zero and therefore don't need to be calculated at all. 
# A separating axis is detected only when all the sides
# of the other polygon (not including the chosen origin)
# have projection distances greater than zero.
#
using GLMakie
using Images # for RGBA
include("ripga2d.jl") # needed for satpga() not sat()

# separating axis test
function sat(A::Matrix{Float32}, B::Matrix{Float32})
 nA = size(A, 2)
 nB = size(B, 2)
 TOPERP = [0 1; -1 0] # rotate to outward perpendicular
 for iA = 1:nA-1
  P0 = A[:, iA]
  DA = (TOPERP * (A[:, iA+1] - P0))'
  count = 0
  for iB = 1:nB-1
   d = DA * (B[:, iB] - P0)
   count = (d > 0) ? count + 1 : break
  end
  if (count == nB-1) return iA; end
 end
 for iB = 1:nB-1
  P0 = B[:, iB]
  DB = (TOPERP * (B[:, iB+1] - P0))'
  count = 0
  for iA = 1:nA-1
   d = DB * (A[:, iA] - P0)
   count = (d > 0) ? count + 1 : break
  end
  if (count == nA-1) return nA + iB; end
 end
 return 0
end

# separating axis test ... PGA version
#
# NOTE:
# The separating axis test is fundamentally a bunch
# of 2D dot products. Employing PGA to implement 2D
# dot products is overkill but it can be done and is
# only a little slower. For example, my laptop computer
# takes 3.1 seconds to run the PGA version to generate 
# the 15 second video and 3.0 seconds to run the nonPGA 
# version to generate the same 15 second video.
#
function satpga(A::Matrix{Float32}, B::Matrix{Float32})
 nA = size(A, 2)
 nB = size(B, 2)
 AX = point(A)
 BX = point(B)
 for iA = 1:nA-1
#  lineProj = Float32.(rotor(pi/2, AX[:,iA])) >>> 
#   ga"(AX[:,iA]∨AX[:,iA+1])"
  lineProj = Float32.(rotor(pi/2, AX[:,iA])) >>> 
   (AX[:,iA] & AX[:,iA+1])
  count = 0
  for iB = 1:nB-1
#   d = ga"((AX[:,iA]∨BX[:,iB])·lineProj)[1]"
   d = ((AX[:,iA] & BX[:,iB]) | lineProj)[1]
   count = (d > 0) ? count + 1 : break
  end
  if (count == nB-1) return iA; end
 end
 for iB = 1:nB-1
#  lineProj = Float32.(rotor(pi/2, BX[:,iB])) >>> 
#   ga"(BX[:,iB]∨BX[:,iB+1])"
  lineProj = Float32.(rotor(pi/2, BX[:,iB])) >>> 
   (BX[:,iB] & BX[:,iB+1])
  count = 0
  for iA = 1:nA-1
#   d = ga"((BX[:,iB]∨AX[:,iA])·lineProj)[1]"
   d = ((BX[:,iB] & AX[:,iA]) | lineProj)[1]
   count = (d > 0) ? count + 1 : break
  end
  if (count == nA-1) return nA + iB; end
 end
 return 0
end

# polygon intersection testing
function polyx()

 # initialize figure
 fig = Figure(resolution = (800, 800))
 LIM = (-8,8, -8,8)
 ax2d = GLMakie.Axis(fig[1,1], limits=LIM, aspect=1,
  title = "polygon intersection testing")
 
 # define shapes and paths
 iFrame = 1
 nFrame = 360
 nPath = nFrame + 1 # +1 to avoid duplicate frame
      # in repeating animated gif
 rMajor = 5
 rMinor = 3
 THETAP = LinRange(0, 2*pi, nPath)'
 PATH = [rMinor.*cos.(THETAP); rMajor.*sin.(THETAP)]
 THETAS = LinRange(0, 6*pi, nPath)'
 
 nA = 3
 rA = 2
 xA = 0
 yA = 0
 THETAA = LinRange(0, 2*pi, nA+1)'
 A = Float32.(
  [xA .+ rA.*cos.(THETAA); # polygon A
   yA .+ rA.*sin.(THETAA)])
 
 nB = 6
 rB = 4
 xB = rMinor
 yB = 0
 THETAB = LinRange(0, 2*pi, nB+1)'
 BV = [rB.*cos.(THETAB); # polygon B
   rB.*sin.(THETAB)]
 BV[:,div(nB,2,RoundDown)+1] .= 0
 B = Matrix{Float32}[]
 push!(B, hcat(BV[:,1:4], BV[:,1]))
 push!(B, hcat(BV[:,1], BV[:,4], BV[:,5:end]))
 
 # initialize plot containing observable data
 MB = Matrix{Float32}[] # mobile convex polygon partitions
 push!(MB, Float32.(B[1] .+ PATH[:,iFrame])) 
 push!(MB, Float32.(B[2] .+ PATH[:,iFrame]))
 OBS_MB = Observable[]
 push!(OBS_MB, Observable(MB[1]))
 push!(OBS_MB, Observable(MB[2]))
 OBS_C = Observable[]
 push!(OBS_C, Observable(RGBA(0, 1, 0, 0.5)))
 push!(OBS_C, Observable(RGBA(0, 1, 0, 0.5)))
 poly!(ax2d, A, color = RGBA(1, 0, 0, 0.5))
 poly!(ax2d, OBS_MB[1], color = OBS_C[1])
 poly!(ax2d, OBS_MB[2], color = OBS_C[2])
 lines!(ax2d, PATH, linestyle = :dot, color = :gray)
 iSep = sat(A,MB[1])
 fig
 
 # record spinning polynomial following path
 record(fig, "polyx.mp4", 1:nFrame) do iFrame
  ROT = [cos(THETAS[iFrame]) -sin(THETAS[iFrame]);
      sin(THETAS[iFrame])  cos(THETAS[iFrame])]
  nPartition = length(B)
  for iPartition = 1:nPartition
   MB[iPartition] = Float32.((ROT * B[iPartition]) 
    .+ PATH[:,iFrame])
   OBS_C[iPartition][] = (satpga(A, MB[iPartition]) == 0) ?
    RGBA(1, 0, 0, 0.5) : RGBA(0, 1, 0, 0.5)
   OBS_MB[iPartition][] = MB[iPartition]
  end # for each partition
 end # for each video frame
end # polyx()

# quick and dirty ffmpeg conversion from polyx.mp4 to polyx.gif:
# ffmpeg -i polyx.mp4 -r 24 -s 480x480 -loop 0 polyx.gif