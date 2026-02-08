# polyx3.jl
# polygon 3D intersection testing
# The Separating Axis Theorem (SAT) applies to only
# convex polygons. However, any non-convex polygon
# can be partitioned into convex polygons.
#
# NOTE:
# In this example, the vertices of each face should 
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
# The tetrahedron has three faces (the base face is
# empty), all edges of length 2.0, and a peak that
# is directly above (i.e., z offset) the origin. All
# the tetrahedron's vertices and faces are defined,
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
f 1 3 2
=#
#
# The unit cube has four sides (the bottom and top
# sides are empty). All the unit cube's vertices and 
# faces are defined, according to the wavefront 3D 
# object file format, by 16 lines of text (within the following 
# multiline comment): 
#=
v  0.5  0.5 0
v -0.5  0.5 0
v -0.5 -0.5 0
v  0.5 -0.5 0
v  0.5  0.5 1
v -0.5  0.5 1
v -0.5 -0.5 1
v  0.5 -0.5 1
f 1 2 6
f 6 5 1
f 2 3 7
f 7 6 2
f 3 4 8
f 8 7 3
f 4 1 5
f 5 8 4
f 1 3 2
f 3 1 4
f 5 6 7
f 7 8 5
=#
using GLMakie, GLMakie.FileIO
using GeometryBasics # for normal_mesh()
using Images # for RGBA
include("ripga3d.jl") # needed for satpga() not sat()

function point(A::Point{3, Float32})
 return A[1]*e032 + A[2]*e013 + A[3]*e021 + e123
end

# separating axis test
# arguments: 
#  vertices, faces, face normals of 2 3D objects
function sat(
 AV::Vector{Point{3, Float32}},
 AF::Vector{NgonFace{3, OffsetInteger{-1, UInt32}}},
 AN::Matrix{Float32},
 BV::Vector{Point{3, Float32}},
 BF::Vector{NgonFace{3, OffsetInteger{-1, UInt32}}},
 BN::Matrix{Float32})

 # determine number of vertices and faces in 3D objects
 nAV = length(AV)
 nAF = length(AF)
 nBV = length(BV)
 nBF = length(BF)
 
 # for each face of first 3D object
 for iAF = 1:nAF
  P0 = AV[AF[iAF]][1]
  P1 = P0 + AN[:,iAF]
  R01 = (P1 - P0)'
  count = 0
  for iBV = 1:nBV
   P2 = BV[iBV]
   d = R01 * (P2 - P0)
   count = (d > 0) ? count + 1 : break
  end
  if (count == nBV) 
   return iAF; 
  end
 end
 
 # for each face of second 3D object
 for iBF = 1:nBF
  P0 = BV[BF[iBF]][1]
  P1 = P0 + BN[:,iBF]
  R01 = (P1 - P0)'
  count = 0
  for iAV = 1:nAV
   d = R01 * (AV[iAV] - P0)
   count = (d > 0) ? count + 1 : break
  end
  if (count == nAV)
   return nAF + iBF; 
  end
 end
 return 0 # no gaps detected -> an intersection
end

# separating axis test using PGA
# arguments: 
#  vertices, faces, face normals of 2 3D objects
function satpga(
 AV::Vector{Point{3, Float32}},
 AF::Vector{NgonFace{3, OffsetInteger{-1, UInt32}}},
 AN::Matrix{Float32},
 BV::Vector{Point{3, Float32}},
 BF::Vector{NgonFace{3, OffsetInteger{-1, UInt32}}},
 BN::Matrix{Float32})

 # determine number of vertices and faces in 3D objects
 nAV = length(AV)
 nAF = length(AF)
 nBV = length(BV)
 nBF = length(BF)
 
 # for each face of first 3D object
 for iAF = 1:nAF
  P0 = point(AV[AF[iAF]][1])
  P1 = point(AV[AF[iAF]][1] + AN[:,iAF])
  count = 0
  for iBV = 1:nBV
   P2 = point(BV[iBV])
#   Plane = ga"P0 · P1 · P2"
#   proj = ga"P0 ∨ P1 · Plane"
#   d = ga"(proj · (P0 ∨ P2 · Plane))[1]"
   Plane = P0 & P1 & P2
   proj = P0 & P1 | Plane
   d = (proj | (P0 & P2 | Plane))[1]
   count = (d > 0) ? count + 1 : break
  end
  if (count == nBV)
   return iAF; 
  end
 end
 
 # for each face of second 3D object
 for iBF = 1:nBF
  P0 = point(BV[BF[iBF]][1])
  P1 = point(BV[BF[iBF]][1] + BN[:,iBF])
  count = 0
  for iAV = 1:nAV
   P2 = point(AV[iAV])
#   Plane = ga"P0 · P1 · P2"
#   proj = ga"P0 ∨ P1 · Plane"
#   d = ga"(proj · (P0 ∨ P2 · Plane))[1]"
   Plane = P0 & P1 & P2
   proj = P0 & P1 | Plane
   d = (proj | (P0 & P2 | Plane))[1]
   count = (d > 0) ? count + 1 : break
  end
  if (count == nAV)
   return nAF + iBF; 
  end
 end
 return 0 # no gaps detected -> an intersection
end

# calculate unit vectors normal to 3D object's faces
function calc_normals(X)
 C,F = coordinates(X),faces(X)
 nF = length(F)
 res = Matrix{Float32}(undef,3,nF)
 for iF = 1:nF
  P = C[F[iF]]
  V1 = P[2] - P[1]
  V2 = P[3] - P[1]
  N = [
   V1[2]*V2[3] - V1[3]*V2[2];
    -(V1[1]*V2[3] - V1[3]*V2[1]);
   V1[1]*V2[2] - V1[2]*V2[1]]
  mag = sqrt(N[1]^2 + N[2]^2 + N[3]^2)
  res[:,iF] = N ./ mag
 end
 return res
end

# polygon 3D intersection testing
function polyx3()

 # define path
 nFrame = 360
 nPath = nFrame + 1 # +1 to avoid duplicate frame
      # in repeating animated gif
 rMajor = 1.2
 rMinor = 0.6
 zMax = 1.1
 THETAP = LinRange(0, 2*pi, nPath)'
 PATH = Float32.([
  rMajor.*cos.(THETAP); 
  rMinor.*sin.(THETAP);
  zMax.*sin.(THETAP)])
 
 # define spin of object and elevation of camera
 THETAS = Float32.(LinRange(0, 6*pi, nPath)')
 nFrame4 = div(nFrame,4,RoundDown)
 THETAEL = Vector{Float32}(undef, nPath)
 THETAEL[1:nFrame4] = 
  LinRange(pi/8, 0.95*pi/2, nFrame4)
 THETAEL[nFrame4+1:3*nFrame4] =
  LinRange(0.95*pi/2,0.05,2*nFrame4)
 THETAEL[3*nFrame4+1:end] = 
  LinRange(0.05, pi/8, nFrame4+1)
 
 # initialize figure
 fig = Figure(resolution = (800, 800))
 ax3d = GLMakie.Axis3(fig[1,1],
  elevation = pi/8,
  azimuth = -pi/4,
  viewmode = :fit,
  aspect = (1,1,1),
  limits = (-2.,2., -2.,2., -2.,2.),
  title = "polygon 3D intersection testing")
 
 # load meshes of two convex 3D objects
 C = load("ucube.obj")  # Cube (stationary)
 CV,CF,CN = coordinates(C),faces(C),calc_normals(C)
# T = load("tetrahedron.obj") # Tetrahedron (on path)
 T = load("ucube.obj")
# T = deepcopy(C)
 TV,TF = coordinates(T),faces(T)

 # initialize plot containing observable data
 iFrame = 45
# str = @sprintf("iFrame %d> ROT=%f; PATH=[%f %f %f]",
#  iFrame, THETAS[iFrame], 
#  PATH[1,iFrame], PATH[2,iFrame], PATH[3,iFrame])
# println(str)
 ROT = Float32.([
  cos(THETAS[iFrame]) -sin(THETAS[iFrame]) 0;
  sin(THETAS[iFrame])  cos(THETAS[iFrame]) 0;
  0      0      1])
 TV2 = Vector{Point{3, Float32}}(
  map(c2 -> PATH[:,iFrame] + c2, # translate second
  map(c -> ROT * c, TV))) # rotate first
 T2 = GeometryBasics.Mesh(TV2, TF)
 T_obs = Observable(T2)
 TN2 = calc_normals(T2)
# iSepFace = sat(CV,CF,CN, TV2,TF,TN2)
 iSepFace = satpga(CV,CF,CN, TV2,TF,TN2)
 T_color_obs = iSepFace==0 ?
  Observable(RGBA(1, 0, 0, 0.5)) :
  Observable(RGBA(0, 1, 0, 0.5))
 lines!(ax3d, PATH, color = :gold)
 mesh!(ax3d, C, 
  color = RGBA(1, 0, 0, 0.5))
 mesh!(ax3d, T_obs, 
  color = T_color_obs)
 fig
 
 # generate video
 record(fig, "polyx3.mp4", 1:nFrame) do iFrame
  ax3d.elevation[] = THETAEL[iFrame]
  ROT = Float32.([
   cos(THETAS[iFrame]) -sin(THETAS[iFrame]) 0;
   sin(THETAS[iFrame])  cos(THETAS[iFrame]) 0;
   0      0      1])
  TV2 = Vector{Point{3, Float32}}(
   map(c2 -> PATH[:,iFrame] + c2, # translate second
   map(c -> ROT * c, TV))) # rotate first 
  T2 = GeometryBasics.Mesh(TV2, TF)
  T_obs[] = T2
  TN2 = calc_normals(T2)
#  T_color_obs[] = sat(CV,CF,CN, TV2,TF,TN2)==0 ?
  T_color_obs[] = satpga(CV,CF,CN, TV2,TF,TN2)==0 ?
   RGBA(1, 0, 0, 0.5) :
   RGBA(0, 1, 0, 0.5)
 end # for each video frame
end # polyx()

# quick and dirty ffmpeg conversion from polyx3.mp4 to polyx3.gif:
# ffmpeg -i polyx3.mp4 -r 24 -s 480x480 -loop 0 polyx3.gif