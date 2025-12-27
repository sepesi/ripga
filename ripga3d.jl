# ripga3d.jl : reference implementation of 
# Projective Geometric Algebra for 3D
#
# This is a Julia port of pga3d.cpp, bivector.net's C++ 
# reference implementation of projective geometric algebra
# available at https://bivector.net/tools.html
#
using Printf
 
# define multivector basis names 
# 0 denotes projective dimension (e.g., e0 * e0 = 0)
basis = [ # iField
 "1"  #  1 scalar
 "e0" #  2 grade 1 vectors
 "e1" #  3
 "e2" #  4
 "e3" #  5
 "e01" #  6 grade 2 vectors (bivectors)
 "e02" #  7
 "e03" #  8
 "e12" #  9
 "e31" # 10
 "e23" # 11
 "e021" # 12 grade 3 vectors (trivectors)
 "e013" # 13
 "e032" # 14
 "e123" # 15
 "e0123"]# 16 pseudoscalar

# define basis multivectors
nField = 2^4 # 4 = 3D + 1 dimensions 
e0 =    zeros(Float32, nField); e0[2] = 1
e1 =    zeros(Float32, nField); e1[3] = 1
e2 =    zeros(Float32, nField); e2[4] = 1
e3 =    zeros(Float32, nField); e3[5] = 1
e01 =   zeros(Float32, nField); e01[6] = 1
e02 =   zeros(Float32, nField); e02[7] = 1
e03 =   zeros(Float32, nField); e03[8] = 1
e12 =   zeros(Float32, nField); e12[9] = 1
e31 =   zeros(Float32, nField); e31[10] = 1
e23 =   zeros(Float32, nField); e23[11] = 1
e021 =  zeros(Float32, nField); e021[12] = 1
e013 =  zeros(Float32, nField); e013[13] = 1
e032 =  zeros(Float32, nField); e032[14] = 1
e123 =  zeros(Float32, nField); e123[15] = 1
e0123 = zeros(Float32, nField); e0123[16] = 1

# geometric product
function Base.:*(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=b[1]*a[1]+b[3]*a[3]+b[4]*a[4]+b[5]*a[5]-b[9]*a[9]-b[10]*a[10]-b[11]*a[11]-b[15]*a[15]
 res[2]=b[2]*a[1]+b[1]*a[2]-b[6]*a[3]-b[7]*a[4]-b[8]*a[5]+b[3]*a[6]+b[4]*a[7]+b[5]*a[8]+b[12]*a[9]+b[13]*a[10]+b[14]*a[11]+b[9]*a[12]+b[10]*a[13]+b[11]*a[14]+b[16]*a[15]-b[15]*a[16]
 res[3]=b[3]*a[1]+b[1]*a[3]-b[9]*a[4]+b[10]*a[5]+b[4]*a[9]-b[5]*a[10]-b[15]*a[11]-b[11]*a[15]
 res[4]=b[4]*a[1]+b[9]*a[3]+b[1]*a[4]-b[11]*a[5]-b[3]*a[9]-b[15]*a[10]+b[5]*a[11]-b[10]*a[15]
 res[5]=b[5]*a[1]-b[10]*a[3]+b[11]*a[4]+b[1]*a[5]-b[15]*a[9]+b[3]*a[10]-b[4]*a[11]-b[9]*a[15]
 res[6]=b[6]*a[1]+b[3]*a[2]-b[2]*a[3]-b[12]*a[4]+b[13]*a[5]+b[1]*a[6]-b[9]*a[7]+b[10]*a[8]+b[7]*a[9]-b[8]*a[10]-b[16]*a[11]-b[4]*a[12]+b[5]*a[13]+b[15]*a[14]-b[14]*a[15]-b[11]*a[16]
 res[7]=b[7]*a[1]+b[4]*a[2]+b[12]*a[3]-b[2]*a[4]-b[14]*a[5]+b[9]*a[6]+b[1]*a[7]-b[11]*a[8]-b[6]*a[9]-b[16]*a[10]+b[8]*a[11]+b[3]*a[12]+b[15]*a[13]-b[5]*a[14]-b[13]*a[15]-b[10]*a[16]
 res[8]=b[8]*a[1]+b[5]*a[2]-b[13]*a[3]+b[14]*a[4]-b[2]*a[5]-b[10]*a[6]+b[11]*a[7]+b[1]*a[8]-b[16]*a[9]+b[6]*a[10]-b[7]*a[11]+b[15]*a[12]-b[3]*a[13]+b[4]*a[14]-b[12]*a[15]-b[9]*a[16]
 res[9]=b[9]*a[1]+b[4]*a[3]-b[3]*a[4]+b[15]*a[5]+b[1]*a[9]+b[11]*a[10]-b[10]*a[11]+b[5]*a[15]
 res[10]=b[10]*a[1]-b[5]*a[3]+b[15]*a[4]+b[3]*a[5]-b[11]*a[9]+b[1]*a[10]+b[9]*a[11]+b[4]*a[15]
 res[11]=b[11]*a[1]+b[15]*a[3]+b[5]*a[4]-b[4]*a[5]+b[10]*a[9]-b[9]*a[10]+b[1]*a[11]+b[3]*a[15]
 res[12]=b[12]*a[1]-b[9]*a[2]+b[7]*a[3]-b[6]*a[4]+b[16]*a[5]-b[4]*a[6]+b[3]*a[7]-b[15]*a[8]-b[2]*a[9]+b[14]*a[10]-b[13]*a[11]+b[1]*a[12]+b[11]*a[13]-b[10]*a[14]+b[8]*a[15]-b[5]*a[16]
 res[13]=b[13]*a[1]-b[10]*a[2]-b[8]*a[3]+b[16]*a[4]+b[6]*a[5]+b[5]*a[6]-b[15]*a[7]-b[3]*a[8]-b[14]*a[9]-b[2]*a[10]+b[12]*a[11]-b[11]*a[12]+b[1]*a[13]+b[9]*a[14]+b[7]*a[15]-b[4]*a[16]
 res[14]=b[14]*a[1]-b[11]*a[2]+b[16]*a[3]+b[8]*a[4]-b[7]*a[5]-b[15]*a[6]-b[5]*a[7]+b[4]*a[8]+b[13]*a[9]-b[12]*a[10]-b[2]*a[11]+b[10]*a[12]-b[9]*a[13]+b[1]*a[14]+b[6]*a[15]-b[3]*a[16]
 res[15]=b[15]*a[1]+b[11]*a[3]+b[10]*a[4]+b[9]*a[5]+b[5]*a[9]+b[4]*a[10]+b[3]*a[11]+b[1]*a[15]
 res[16]=b[16]*a[1]+b[15]*a[2]+b[14]*a[3]+b[13]*a[4]+b[12]*a[5]+b[11]*a[6]+b[10]*a[7]+b[9]*a[8]+b[8]*a[9]+b[7]*a[10]+b[6]*a[11]-b[5]*a[12]-b[4]*a[13]-b[3]*a[14]-b[2]*a[15]+b[1]*a[16]
 return res
end # geometric product (*)

# regressive product: vee operator (&, \vee)
function Base.:&(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[16]=a[16]*b[16]
 res[15]=a[15]*b[16]+a[16]*b[15]
 res[14]=a[14]*b[16]+a[16]*b[14]
 res[13]=a[13]*b[16]+a[16]*b[13]
 res[12]=a[12]*b[16]+a[16]*b[12]
 res[11]=a[11]*b[16]+a[14]*b[15]-a[15]*b[14]+a[16]*b[11]
 res[10]=a[10]*b[16]+a[13]*b[15]-a[15]*b[13]+a[16]*b[10]
 res[9]=a[9]*b[16]+a[12]*b[15]-a[15]*b[12]+a[16]*b[9]
 res[8]=a[8]*b[16]+a[13]*b[14]-a[14]*b[13]+a[16]*b[8]
 res[7]=a[7]*b[16]-a[12]*b[14]+a[14]*b[12]+a[16]*b[7]
 res[6]=a[6]*b[16]+a[12]*b[13]-a[13]*b[12]+a[16]*b[6]
 res[5]=a[5]*b[16]+a[8]*b[15]-a[10]*b[14]+a[11]*b[13]+a[13]*b[11]-a[14]*b[10]+a[15]*b[8]+a[16]*b[5]
 res[4]=a[4]*b[16]+a[7]*b[15]+a[9]*b[14] -a[11]*b[12]-a[12]*b[11]+a[14]*b[9] +a[15]*b[7]+a[16]*b[4]
 res[3]=a[3]*b[16]+a[6]*b[15]-a[9]*b[13]+a[10]*b[12]+a[12]*b[10]-a[13]*b[9]+a[15]*b[6]+a[16]*b[3]
 res[2]=a[2]*b[16]-a[6]*b[14]-a[7]*b[13]-a[8]*b[12] -a[12]*b[8] -a[13]*b[7]-a[14]*b[6]+a[16]*b[2]
 res[1]=a[1]*b[16]-a[2]*b[15]-a[3]*b[14]-a[4]*b[13] -a[5]*b[12] +a[6]*b[11]+a[7]*b[10]+a[8]*b[9]+a[9]*b[8]+a[10]*b[7]+a[11]*b[6]+a[12]*b[5]+a[13]*b[4]+a[14]*b[3]+a[15]*b[2]+a[16]*b[1]
 return res
end # regressive product; vee operator (&, \vee)

# inner product (|)
function Base.:|(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=b[1]*a[1]+b[3]*a[3]+b[4]*a[4]+b[5]*a[5]-b[9]*a[9]-b[10]*a[10]-b[11]*a[11]-b[15]*a[15]
 res[2]=b[2]*a[1]+b[1]*a[2]-b[6]*a[3]-b[7]*a[4]-b[8]*a[5]+b[3]*a[6]+b[4]*a[7]+b[5]*a[8]+b[12]*a[9]+b[13]*a[10]+b[14]*a[11]+b[9]*a[12]+b[10]*a[13]+b[11]*a[14]+b[16]*a[15]-b[15]*a[16]
 res[3]=b[3]*a[1]+b[1]*a[3]-b[9]*a[4]+b[10]*a[5]+b[4]*a[9]-b[5]*a[10]-b[15]*a[11]-b[11]*a[15]
 res[4]=b[4]*a[1]+b[9]*a[3]+b[1]*a[4]-b[11]*a[5]-b[3]*a[9]-b[15]*a[10]+b[5]*a[11]-b[10]*a[15]
 res[5]=b[5]*a[1]-b[10]*a[3]+b[11]*a[4]+b[1]*a[5]-b[15]*a[9]+b[3]*a[10]-b[4]*a[11]-b[9]*a[15]
 res[6]=b[6]*a[1]-b[12]*a[4]+b[13]*a[5]+b[1]*a[6]-b[16]*a[11]-b[4]*a[12]+b[5]*a[13]-b[11]*a[16]
 res[7]=b[7]*a[1]+b[12]*a[3]-b[14]*a[5]+b[1]*a[7]-b[16]*a[10]+b[3]*a[12]-b[5]*a[14]-b[10]*a[16]
 res[8]=b[8]*a[1]-b[13]*a[3]+b[14]*a[4]+b[1]*a[8]-b[16]*a[9]-b[3]*a[13]+b[4]*a[14]-b[9]*a[16]
 res[9]=b[9]*a[1]+b[15]*a[5]+b[1]*a[9]+b[5]*a[15]
 res[10]=b[10]*a[1]+b[15]*a[4]+b[1]*a[10]+b[4]*a[15]
 res[11]=b[11]*a[1]+b[15]*a[3]+b[1]*a[11]+b[3]*a[15]
 res[12]=b[12]*a[1]+b[16]*a[5]+b[1]*a[12]-b[5]*a[16]
 res[13]=b[13]*a[1]+b[16]*a[4]+b[1]*a[13]-b[4]*a[16]
 res[14]=b[14]*a[1]+b[16]*a[3]+b[1]*a[14]-b[3]*a[16]
 res[15]=b[15]*a[1]+b[1]*a[15]
 res[16]=b[16]*a[1]+b[1]*a[16]
 return res
end # inner product (|)

# outer product; wedge operator (^)
function Base.:^(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=b[1]*a[1]
 res[2]=b[2]*a[1]+b[1]*a[2]
 res[3]=b[3]*a[1]+b[1]*a[3]
 res[4]=b[4]*a[1]+b[1]*a[4]
 res[5]=b[5]*a[1]+b[1]*a[5]
 res[6]=b[6]*a[1]+b[3]*a[2]-b[2]*a[3]+b[1]*a[6]
 res[7]=b[7]*a[1]+b[4]*a[2]-b[2]*a[4]+b[1]*a[7]
 res[8]=b[8]*a[1]+b[5]*a[2]-b[2]*a[5]+b[1]*a[8]
 res[9]=b[9]*a[1]+b[4]*a[3]-b[3]*a[4]+b[1]*a[9]
 res[10]=b[10]*a[1]-b[5]*a[3]+b[3]*a[5]+b[1]*a[10]
 res[11]=b[11]*a[1]+b[5]*a[4]-b[4]*a[5]+b[1]*a[11]
 res[12]=b[12]*a[1]-b[9]*a[2]+b[7]*a[3]-b[6]*a[4]-b[4]*a[6]+b[3]*a[7]-b[2]*a[9]+b[1]*a[12]
 res[13]=b[13]*a[1]-b[10]*a[2]-b[8]*a[3]+b[6]*a[5]+b[5]*a[6]-b[3]*a[8]-b[2]*a[10]+b[1]*a[13]
 res[14]=b[14]*a[1]-b[11]*a[2]+b[8]*a[4]-b[7]*a[5]-b[5]*a[7]+b[4]*a[8]-b[2]*a[11]+b[1]*a[14]
 res[15]=b[15]*a[1]+b[11]*a[3]+b[10]*a[4]+b[9]*a[5]+b[5]*a[9]+b[4]*a[10]+b[3]*a[11]+b[1]*a[15]
 res[16]=b[16]*a[1]+b[15]*a[2]+b[14]*a[3]+b[13]*a[4]+b[12]*a[5]+b[11]*a[6]+b[10]*a[7]+b[9]*a[8]+b[8]*a[9]+b[7]*a[10]+b[6]*a[11]-b[5]*a[12]-b[4]*a[13]-b[3]*a[14]-b[2]*a[15]+b[1]*a[16]
 return res
end # outer product; wedge operator (^)

# reverse operator (~)
function Base.:~(a::Vector{Float32})
 res = copy(a)
 res[6:15] .*= -1
 return res
end

# conjugate
function conjugate(a::Vector{Float32})::Vector{Float32}
 res = copy(a)
 res[2:11] .*= -1
 return res
end

function plane(a::Number, b::Number, c::Number, d::Number)::Vector{Float32}
 return a*e1 + b*e2 + c*e3 + d*e0
end

function circle(t::Number, radius::Number, line::Vector{Float32})
 return rotor(t*2*pi,line) * translator(radius,e1*e0)
end

function torus(s::Number, t::Number,  
 r1::Number, l1::Vector{Float32},
 r2::Number, l2::Vector{Float32})
 return circle(s,r2,l2) * circle(t,r1,l1)
end
 
function normIdeal(a::Vector{Float32},nd::Int64=2)
 if nd == 2 # for area
  return sqrt(a[6]^2 + a[7]^2 + a[8]^2)
 else # for volume
  return sqrt(a[2]^2)
 end
end

# unit test
# arguments:
# - nLoop repeats a section of the PGA calculations for benchmarking
# - flgSimplify set to true diverges from the text output by pga3d.cpp
# usage notes:
# - @time utest(1) checks on whether the unit test output of ripga.jl
# exactly matches the unit test output of pga3d.cpp. The comparison
# ends with the printing of the number of tests in the unit test that
# don't match. 0 indicates success.
# - @time utest(1,true) outputs a slightly simplified version of the 
# unit test output that does not match the unit test output by pga3d.cpp.
# - @btime utest() is a test for execution speed of ripga3d.jl.
#   (NOTE: requires using BenchmarkTools)
function utest(nLoop=100, flgSimplify::Bool=false)
 nField = length(basis)
 axis_z = Vector{Float32}(undef,nField)
 origin = Vector{Float32}(undef,nField)
 px = Vector{Float32}(undef,nField)
 line = Vector{Float32}(undef,nField)
 p = Vector{Float32}(undef,nField)
 rot = Vector{Float32}(undef,nField)
 rot_point = Vector{Float32}(undef,nField)
 rot_line = Vector{Float32}(undef,nField)
 rot_plane = Vector{Float32}(undef,nField)
 point_on_plane = Vector{Float32}(undef,nField)
 to = Vector{Float32}(undef,nField)
 point_on_torus = Vector{Float32}(undef,nField)
 tst1 = Vector{Float32}(undef,nField)
 tst2 = Vector{Float32}(undef,nField)
 
 for iLoop = 1:nLoop
  # geometric algebra equations in programming syntax
  axis_z = e1 ^ e2
  origin = axis_z ^ e3
  
  px = point(1, 0, 0)
  line = origin & px
  p = plane(2,0,1,-3)
  rot = rotor(pi/2, e1*e2)
  rot_point = rot * px * ~rot
  rot_line = rot * line * ~rot
  rot_plane = rot * p * ~rot
  point_on_plane = (p | px) * p
  to = torus(0,0, 0.25,e1*e2, 0.6,e1*e3)
  point_on_torus = to * e123 * ~to
  
  tst1 = e0 - 1
  tst2 = 1 - e0
  
#  # geometric algebra equations in math syntax
#  ga"axis_z = e1 ∧ e2"
#  ga"origin = axis_z ∧ e3"
#  
#  px = point(1, 0, 0)
#  ga"line = origin ∨ px"
#  p = plane(2, 0, 1,-3)
#  ga"rot = rotor(pi/2, e1 e2)"
#  ga"rot_point = rot px ~rot"
#  ga"rot_line = rot line ~rot"
#  ga"rot_plane = rot p ~rot"
#  ga"point_on_plane = (p·px) p"
#  ga"to = torus(0,0, 0.25,e1 e2, 0.6,e1 e3)"
#  ga"point_on_torus = to e123 ~to"
#  
#  tst1 = e0 - 1
#  tst2 = 1 - e0
 end
 
 # if verbose/slow output of unit test results wanted
 if nLoop == 1
  nError = 0
  toStr2 = flgSimplify ? toStr : toStr1
  
  S = Matrix{String}(undef,11,3) # 3 columns:
  S[1,1] = " point          : "  #  1) label
  S[1,2] = toStr2(px)            #  2) toStr1()
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
  
  S[6,1] = " rotated line   : "
  S[6,2] = toStr2(rot_line)
  S[6,3] = flgSimplify ?
   "0.9999999e31" :
   "0.9999999e31"
  
  S[7,1] = " rotated plane  : "
  S[7,2] = toStr2(rot_plane)
  S[7,3] = flgSimplify ?
   "-3e0 - 2e2 + 0.9999999e3" :
   "-3e0 + -2e2 + 0.9999999e3"
  
  S[8,1] = " point on plane : "
  S[8,2] = toStr2(normalize(point_on_plane))
  S[8,3] = flgSimplify ?
   "0.2e021 + 1.4e032 + e123" :
   "0.2e021 + 1.4e032 + 1e123"
  
  S[9,1] = " point on torus : "
  S[9,2] = toStr2(point_on_torus)
  S[9,3] = flgSimplify ?
   "0.85e032 + e123" :
   "0.85e032 + 1e123"
  
  S[10,1] = " toStr1 test 1  : "
  S[10,2] = toStr2(tst1)
  S[10,3] = flgSimplify ?
   "-1 + e0" :
   "-1 + 1e0"
  
  S[11,1] = " toStr1 test 2  : "
  S[11,2] = toStr2(tst2)
  S[11,3] = flgSimplify ?
   "1 - e0" :
   "1 + -1e0"
  
  # print unit test results;
  # print 'x' in first column in tests with errors
  nTest = size(S,1)
  for iTest = 1:nTest
   isError = S[iTest,2] != S[iTest,3]
   xChar = isError ? 'x' : ' '
   println(xChar * S[iTest,1] * S[iTest,2])
   if isError
    println(' ' * S[iTest,1] * S[iTest,3])
    nError += 1
   end
  end
  
  return nError # return unit test results
 end
end # rip3d utest()