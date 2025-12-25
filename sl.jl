# sl.jl
# 3D slicing application example
#
# The test tetrahedron has three faces (the base face
# is empty). All edges are of length 2.0. The peak is 
# directly above (i.e., z offset) the origin. All the
# test tetrahedron's vertices and faces are defined,
# according to the wavefront 3D object file format,
# by seven lines of text (within the following 
# multiline comment): 
#=
v 1.0 -0.5773502692 0.0
v 0.0 1.15470053838 0.0
v -1.0 -0.5773502692 0.0
v 0.0 0.0 1.632993162
f 1 2 4
f 2 3 4
f 3 1 4
=#

using GLMakie, GLMakie.FileIO
using GeometryBasics # for normal_mesh()
include("ripga3d.jl")

function tetratest()
 zmax = 1.5
 zdel = 0.1
 rmax = 1.0
 rdel = 0.15
 nTheta = 7
 THETA = LinRange(-pi/6, 23*pi/6, nTheta)
 open("tetratest.obj", "w") do io
  for iTheta = 1:nTheta
   r = rmax - (iTheta-1) * rdel
   xTri = r*cos(THETA[iTheta])
   yTri = r*sin(THETA[iTheta])
   zTri = (iTheta < nTheta) ?
    (iTheta - 1) * zdel : 0.15
   println(io, "v $xTri $yTri $zTri")
  end
  println(io, "v 0 0 $zmax") # vertex 8
  println(io, "v 0 0 1.48")  # vertex 9
  println(io, "v 0 0 1.46")  # vertex 10
  println(io, "v 0 0 1.44")  # vertex 11
  println(io, "v 0 0 1.42")  # vertex 12
  println(io, "v 0 0 1")     # vertex 13
  println(io, "f 1 2 8")
  println(io, "f 2 3 9")
  println(io, "f 3 4 10")
  println(io, "f 4 5 11")
  println(io, "f 5 6 12")
  println(io, "f 6 7 13")
 end 
end

function zslice(zCut::Float32, F::GeometryBasics.Mesh,
 TZ::Matrix{Float32}, TI::Matrix{Int32},
 SEG::Matrix{Float32}, MEASURE::Vector{Float32})
 
 # initialize
 nCol = 0
 MEASURE[:] .= 0
 SUM1 = zeros(Float32, 16)
# iT0 = findfirst(x -> x>=zCut, TZ[:,1])
 
 # for each triangle (face sorted by height) in 3D object
 nF = length(F)
 for iT = 1:nF
  dz2 = TZ[iT,2] - TZ[iT,1]
  dzc = zCut - TZ[iT,1]
  
  # if triangle intersects with cut plane
  if (dzc >= 0) && (dzc <= dz2)
   iLo = TI[iT,1]
   iHi = TI[iT,2]
   iMid = xor(iLo, iHi)
   iF = TI[iT,3]
   P0 = F[iF][iLo][1:2]
   P1 = F[iF][iMid][1:2]
   P2 = F[iF][iHi][1:2]
   dzm = F[iF][iMid][3] - TZ[iT,1]
   PM = P0 + (dzm/dz2) .* (P2 - P0)
   
   # if cut in lower section of triangle
   if dzc <= dzm
    if dzm == 0
     nCol += 1
     SEG[:,nCol] = F[iF][iLo][1:3]
     nCol += 1
     SEG[:,nCol] = F[iF][iMid][1:3]
     MEASURE[1] += sqrt(
      (SEG[1,nCol] - SEG[1,nCol-1])^2 +
      (SEG[2,nCol] - SEG[2,nCol-1])^2)
     MEASURE[2] += norm(
      point(SEG[:,nCol-1]) &
      point(SEG[:,nCol]))
     P0X = point(SEG[:,nCol-1])
     P1X = point(SEG[:,nCol])
     X = TI[iT,5] * (P0X & P1X)
     SUM1 += X
     nCol += 1 # skip NaN32 column
    else
     s = dzc / dzm
     nCol += 1
     SEG[:,nCol] = [P0+s.*(PM-P0); zCut]
     nCol += 1
     SEG[:,nCol] = [P0+s.*(P1-P0); zCut]
     MEASURE[1] += sqrt(
      (SEG[1,nCol] - SEG[1,nCol-1])^2 +
      (SEG[2,nCol] - SEG[2,nCol-1])^2)
     MEASURE[2] += norm(
      point(SEG[:,nCol-1]) &
      point(SEG[:,nCol]))
     P0X = point(SEG[:,nCol-1])
     P1X = point(SEG[:,nCol])
     X = TI[iT,5] * (P0X & P1X)
     SUM1 += X
     nCol += 1 # skip NaN32 column
    end
     
   # else cut in upper section of triangle
   else
    if dz2 - dzm == 0
     nCol += 1
     SEG[:,nCol] = F[iF][iHi][1:3]
     nCol += 1
     SEG[:,nCol] = F[iF][iMid][1:3]
     MEASURE[1] += sqrt(
      (SEG[1,nCol] - SEG[1,nCol-1])^2 +
      (SEG[2,nCol] - SEG[2,nCol-1])^2)
     MEASURE[2] += norm(
      point(SEG[:,nCol-1]) &
      point(SEG[:,nCol]))
     P0X = point(SEG[:,nCol-1])
     P1X = point(SEG[:,nCol])
     X = TI[iT,5] * (P0X & P1X)
     SUM1 += X
     nCol += 1 # skip NaN32 column
    else
     s = (dzc - dzm) / (dz2 - dzm)
     nCol += 1
     SEG[:,nCol] = [PM+s.*(P2-PM); zCut]
     nCol += 1
     SEG[:,nCol] = [P1+s.*(P2-P1); zCut]
     MEASURE[1] += sqrt(
      (SEG[1,nCol] - SEG[1,nCol-1])^2 +
      (SEG[2,nCol] - SEG[2,nCol-1])^2)
     MEASURE[2] += norm(
      point(SEG[:,nCol-1]) &
      point(SEG[:,nCol]))
     P0X = point(SEG[:,nCol-1])
     P1X = point(SEG[:,nCol])
     X = TI[iT,5] * (P0X & P1X)
     SUM1 += X
     nCol += 1 # skip NaN32 column
    end
   end # else cut in upper section of triangle
  end # if cut slice intersects a triangle
 end # for each triangle
 MEASURE[3] = normIdeal(SUM1) / 2
 return nCol
end

function sl()
 # load and scale (by 13) faces of bunny object
 #F = load("tetrahedron.obj")
 #F = load("tetratest.obj")
 F = load("xbunny.obj")
 nF = length(F)
 
 # initialize volume and surface area measurements
 surface_area = 0
 VOL = zeros(Float32, 16)

 # sort triangles by height within object
 # TZ (i.e., Triangle Z-values) has 2 columns:
 # 1: min z value of triangle (in sorted order)
 # 2: max z value of triangle (not necessarily in order)
 TZ = Matrix{Float32}(undef, (nF,2))
 # TI (i.e., Triangle Indices) has 4 columns:
 # 1: in-triangle index of min z value
 # 2: in-triangle index of max z value
 # 3: index of face corresponding to sorted triangle
 # 4: rank order of triangle's max z value
 # 5: orientation of sorted triangle
 TI = Matrix{Int32}(undef, (nF,5))
 for iF = 1:nF
  FACE = 
   point(F[iF][1][1],F[iF][1][2],F[iF][1][3]) &
   point(F[iF][2][1],F[iF][2][2],F[iF][2][3]) &
   point(F[iF][3][1],F[iF][3][2],F[iF][3][3])
  surface_area += norm(FACE)
  VOL += FACE
   
  if F[iF][1][3] <= F[iF][2][3]
   if F[iF][1][3] <= F[iF][3][3]
    TZ[iF,1] = F[iF][1][3]
    TI[iF,1] = 1
    if F[iF][2][3] <= F[iF][3][3]
     TZ[iF,2] = F[iF][3][3]
     TI[iF,2] = 3
     TI[iF,5] = 1
    else
     TZ[iF,2] = F[iF][2][3]
     TI[iF,2] = 2
     TI[iF,5] = -1
    end
   else
    TZ[iF,1] = F[iF][3][3]
    TI[iF,1] = 3
    TZ[iF,2] = F[iF][2][3]
    TI[iF,2] = 2
    TI[iF,5] = 1
   end
  else
   if F[iF][3][3] <= F[iF][2][3]
    TZ[iF,1] = F[iF][3][3]
    TI[iF,1] = 3
    TZ[iF,2] = F[iF][1][3]
    TI[iF,2] = 1
    TI[iF,5] = -1
   else
    TZ[iF,1] = F[iF][2][3]
    TI[iF,1] = 2
    if F[iF][1][3] <= F[iF][3][3]
     TZ[iF,2] = F[iF][3][3]
     TI[iF,2] = 3
     TI[iF,5] = -1
    else
     TZ[iF,2] = F[iF][1][3]
     TI[iF,2] = 1
     TI[iF,5] = 1
    end
   end
  end
 end
 TI[:,4] = sortperm(TZ[:,2])
 P = sortperm(TZ[:,1])
 TZ = TZ[P,:]
 TI = TI[P,:]
 TI[:,3] = P
 surface_area /= 2
 volume = normIdeal(VOL, 3) / 3
 
 # initialize figures
 fig = Figure(resolution = (1000, 500))
 ax3d = Axis3(fig[1,1],
  elevation = pi/16,
  azimuth = -5*pi/8,
  viewmode = :fit,
  zlabel = "z (cm)",
  aspect = (1,1,1))
 m = mesh!(ax3d, normal_mesh(F),
  color = :chocolate4,
  shading = true)
 ax2d = Axis(fig[1,2],
  limits = (-1,1, -1,1),
  xlabel = "x (cm)",
  ylabel = "y (cm)")

 # plot initial observable data
 SEG = fill(NaN32, 3, 3*nF)  # line SEGment buffer
 MEASURE = zeros(Float32, 3) # slice MEASUREments
 zCut::Float32 = 1.395f0
 nCol = zslice(zCut, F, TZ, TI, SEG, MEASURE)
 SEG_obs = Observable(SEG)
 slice2d = @lift @view $SEG_obs[1:2,1:nCol]
 slice3d = @lift @view $SEG_obs[:,1:nCol]
 strHeight = @sprintf("""
  slice height = %.3f cm
  circumference = %.2f cm
  circumference = %.2f cm
  area = %.3f sq cm
  """,
  zCut, MEASURE[1], MEASURE[2], MEASURE[3])
 str_obs = Observable(strHeight)
 strWholeObject = @sprintf("""
  surface area = %.2f sq cm
  volume = %.2f cc
  """,
  surface_area, volume)
 lines!(ax3d, slice3d,
  linewidth = 5,
  color = :black)
 lines!(ax2d, slice2d,
  color = :black)
 text!(ax2d, str_obs,
  position = (0.05, 0.60),
  space = :data)
 text!(ax2d, strWholeObject,
  position = (0.05, -0.98),
  space = :data)
 fig
 
 # generate video of slicing
 zMax = 2.0
 nFrame = 800
 record(fig, "sl.mp4", 1:nFrame) do iFrame
  ax3d.azimuth[] = -3pi/4 - 2pi*(iFrame-1)/nFrame
  zCut = zMax - abs(zMax - 2*zMax*(iFrame-1)/nFrame)
  nCol = zslice(zCut, F, TZ, TI, SEG, MEASURE)
  notify(SEG_obs)
  str_obs[] = @sprintf("""
   slice height = %.3f cm
   circumference = %.2f cm
   circumference = %.2f cm
   area = %.3f sq cm
   """,
   zCut, MEASURE[1], MEASURE[2], MEASURE[3])
 end # for each video frame
end # sl()

# quick and dirty ffmpeg conversion from sl.mp4 to sl.gif:
# ffmpeg -i sl.mp4 -r 24 -s 460x240 -loop 0 sl.gif
# ffmpeg -i sl.mp4 -r 24 -s 920x480 -loop 0 sl.gif