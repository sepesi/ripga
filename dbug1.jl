# dbug1.jl : debug ga macro error
#
# usage:
# > include("ripga3d.jl")
# > include("ripgand.jl")
# > include("dbug1.jl")
# > dbug1()

function dbug1(flgSimplify::Bool=false)
 # allocate geometric objects
 nField = length(basis)
 axis_z = Vector{Float32}(undef,nField)
 origin = Vector{Float32}(undef,nField)
 px = Vector{Float32}(undef,nField)
 line = Vector{Float32}(undef,nField)
 p = Vector{Float32}(undef,nField)
 rot = Vector{Float32}(undef,nField)
 rot_point = Vector{Float32}(undef,nField)

 # calculate geometric objects with programming syntax
#  axis_z = e1 ^ e2
#  origin = axis_z ^ e3
#  px = point(1, 0, 0)
#  line = origin & px
#  p = plane(2,0,1,-3)
#  rot = rotor(pi/2, e1*e2)
#  rot_point = rot * px * ~rot

 # calculate geometric objects with math syntax
 ga"axis_z = e1 ∧ e2"
 ga"origin = axis_z ∧ e3"
 px = point(1, 0, 0)
 ga"line = origin ∨ px"
 p = plane(2, 0, 1,-3)
 ga"rot = rotor(pi/2, e1 e2)"
 ga"rot_point = rot px ~rot"

 # define tests
 toStr2 = flgSimplify ? toStr : toStr1
 S = Matrix{String}(undef,5,3) # 3 columns:
 S[1,1] = " point          : "  #  1) label
 S[1,2] = toStr2(px)            #  2) toStr() or toStr1()
 S[1,3] = flgSimplify ?         #  3) expected string
  "e032 + e123" :
  "1e032 + 1e123"
 
 S[2,1] = " line           : "
 S[2,2] = toStr2(line)
 S[2,3] = flgSimplify ?
  "-e23" :
  "-1e23"
 
 S[3,1] = " plane          : "
 S[3,2] = toStr2(p)
 S[3,3] = flgSimplify ?
  "-3e0 + 2e1 + e3" :
  "-3e0 + 2e1 + 1e3"
 
 S[4,1] = " rotor          : "
 S[4,2] = toStr2(rot)
 S[4,3] = flgSimplify ?
  "0.7071068 + 0.7071068e12" :
  "0.7071068 + 0.7071068e12"
 
 S[5,1] = " rotated point  : "
 S[5,2] = toStr2(rot_point)
 S[5,3] = flgSimplify ?
  "-0.9999999e013 + 0.9999999e123" :
  "-0.9999999e013 + 0.9999999e123"

 # run tests
 nError = 0
 nTest = size(S,1)
 for iTest = 1:nTest
  isError = S[iTest,2] != S[iTest,3]
  xChar = isError ? 'x' : ' '
  println(xChar * S[iTest,1] * S[iTest,2])
  if isError
   println(' ' * S[iTest,1] * S[iTest,3])
   nError += 1
  end
 end # for iTest
end # dbug1()