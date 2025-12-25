# ripga2d.jl : reference implementation of 
# Projective Geometric Algebra for 2D
#
# This is a Julia port of bivector.net's C++ reference
# implementation of projective geometric algebra
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
 "e01" #  5 grade 2 vectors (bivectors)
 "e20" #  6
 "e12" #  7
 "e012"] #  8 pseudoscalar

# define the basis elements
nField = 2^3 # 3 = 2D + 1 dimensions
e0 =   zeros(Float32, nField); e0[2] = 1
e1 =   zeros(Float32, nField); e1[3] = 1
e2 =   zeros(Float32, nField); e2[4] = 1
e01 =  zeros(Float32, nField); e01[5] = 1
e20 =  zeros(Float32, nField); e20[6] = 1
e12 =  zeros(Float32, nField); e12[7] = 1
e012 = zeros(Float32, nField); e012[8] = 1

function Base.:*(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=a[1]*b[1]+a[3]*b[3]+a[4]*b[4]-a[7]*b[7]
 res[2]=a[1]*b[2]+a[2]*b[1]+a[5]*b[3]-a[3]*b[5]+a[4]*b[6]-a[6]*b[4]-a[7]*b[8]-a[8]*b[7]
 res[3]=a[1]*b[3]+a[3]*b[1]-a[4]*b[7]+a[7]*b[4]
 res[4]=a[1]*b[4]+a[4]*b[1]+a[3]*b[7]-a[7]*b[3]
 res[5]=a[1]*b[5]+a[5]*b[1]+a[2]*b[3]-a[3]*b[2]+a[4]*b[8]+a[8]*b[4]+a[6]*b[7]-a[7]*b[6]
 res[6]=a[1]*b[6]+a[6]*b[1]-a[2]*b[4]+a[4]*b[2]+a[3]*b[8]+a[8]*b[3]-a[5]*b[7]+a[7]*b[5]
 res[7]=a[1]*b[7]+a[7]*b[1]+a[3]*b[4]-a[4]*b[3]
 res[8]=a[1]*b[8]+a[8]*b[1]+a[2]*b[7]+a[7]*b[2]+a[3]*b[6]+a[6]*b[3]+a[4]*b[5]+a[5]*b[4]
 return res
end

function Base.:&(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=a[1]*b[8]+a[2]*b[7]+a[3]*b[6]+a[4]*b[5]+a[5]*b[4]+a[6]*b[3]+a[7]*b[2]+a[8]*b[1]
 res[2]=a[2]*b[8]+a[8]*b[2]+a[5]*b[6]-a[6]*b[5]
 res[3]=a[3]*b[8]+a[8]*b[3]-a[5]*b[7]+a[7]*b[5]
 res[4]=a[4]*b[8]+a[8]*b[4]+a[6]*b[7]-a[7]*b[6]
 res[5]=a[5]*b[8]+a[8]*b[5]
 res[6]=a[6]*b[8]+a[8]*b[6]
 res[7]=a[7]*b[8]+a[8]*b[7]
 res[8]=a[8]*b[8]
 return res
end

function Base.:|(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=a[1]*b[1]+a[3]*b[3]+a[4]*b[4]-a[7]*b[7]
 res[2]=a[1]*b[2]+a[2]*b[1]+a[5]*b[3]-a[6]*b[4]-a[3]*b[5]+a[4]*b[6]-a[7]*b[8]-a[8]*b[7]
 res[3]=a[1]*b[3]+a[3]*b[1]-a[4]*b[7]+a[7]*b[4]
 res[4]=a[1]*b[4]+a[4]*b[1]+a[3]*b[7]-a[7]*b[3]
 res[5]=a[1]*b[5]+a[5]*b[1]+a[4]*b[8]+a[8]*b[4]
 res[6]=a[1]*b[6]+a[6]*b[1]+a[3]*b[8]+a[8]*b[3]
 res[7]=a[1]*b[7]+a[7]*b[1]
 res[8]=a[1]*b[8]+a[8]*b[1]
 return res
end

function Base.:^(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=a[1]*b[1]
 res[2]=a[1]*b[2]+a[2]*b[1]
 res[3]=a[1]*b[3]+a[3]*b[1]
 res[4]=a[1]*b[4]+a[4]*b[1]
 res[5]=a[1]*b[5]+a[5]*b[1]+a[2]*b[3]-a[3]*b[2]
 res[6]=a[1]*b[6]+a[6]*b[1]-a[2]*b[4]+a[4]*b[2]
 res[7]=a[1]*b[7]+a[7]*b[1]+a[3]*b[4]-a[4]*b[3]
 res[8]=a[1]*b[8]+a[2]*b[7]+a[3]*b[6]+a[4]*b[5]+a[5]*b[4]+a[6]*b[3]+a[7]*b[2]+a[8]*b[1]
 return res
end

function Base.:~(a::Vector{Float32}) # reverse operator
 res = copy(a)
 res[5:8] .*= -1
 return res
end

function conjugate(a::Vector{Float32})::Vector{Float32}
 res = copy(a)
 res[2:7] .*= -1
 return res
end

function normIdeal(a::Vector{Float32})
 return sqrt(a[2]^2)
end

# unit test
# arguments:
# - nLoop repeats a section of the PGA calculations for benchmarking
# usage notes:
# - @time utest(1) checks on whether the unit test output of ripga.jl
# exactly matches the unit test output of pga3d.cpp. The comparison
# ends with the printing of the number of tests in the unit test that
# don't match. 0 indicates success.
# - @time utest(1,true) outputs a slightly simplified version of the 
# unit test output that does not match the unit test output by pga3d.cpp.
# - @btime utest() is a test for execution speed of ripga.jl.
#   (NOTE: requires using BenchmarkTools)
function utest(nLoop=100, flgSimplify::Bool=false)
 nField = length(basis)
 P0 = Vector{Float32}(undef,nField)
 P1 = Vector{Float32}(undef,nField)
 P2 = Vector{Float32}(undef,nField)
 P3 = Vector{Float32}(undef,nField)
 line0 = Vector{Float32}(undef,nField)
 line1 = Vector{Float32}(undef,nField)
 x = Vector{Float32}(undef,nField)
 tst1 = Vector{Float32}(undef,nField)
 tst2 = Vector{Float32}(undef,nField)

 for iLoop = 1:nLoop
  # define some points
  P0 = point(0,0)
  P1 = point(1,0)
  P2 = point(0,1)
  P3 = point(1,1)
  
  # calculate intersection of parallel lines
  line0 = P0 & P1
  line1 = P2 & P3
  x = line0 ^ line1
  
  tst1 = e0 - 1f0
  tst2 = 1f0 - e0
 
#  # geometric algebra equations in math syntax
#  ga"axis_z = e1 ∧ e2"
#  ga"origin = axis_z ∧ e3"
#  
#  px = point(1f0, 0f0, 0f0)
#  ga"line = origin ∨ px"
#  p = plane(2f0,0f0,1f0,-3f0)
#  ga"rot = rotor(Float32(pi/2), e1 e2)"
#  ga"rot_point = rot px ~rot"
#  ga"rot_line = rot line ~rot"
#  ga"rot_plane = rot p ~rot"
#  
#  tst1 = e0 - 1f0
#  tst2 = 1f0 - e0
 end

 # if (slow) output of unit test results wanted
 if nLoop == 1
  nError = 0

  S = Matrix{String}(undef,9,3)  # 3 columns:
  S[1,1] = " P0             : "  # 1) label
  S[1,2] = toStr(P0)             # 2) toStr()
  S[1,3] = "e12"                 # 3) expected string
  
  S[2,1] = " P1             : "
  S[2,2] = toStr(P1)
  S[2,3] = "e20 + e12"
  
  S[3,1] = " P2             : "
  S[3,2] = toStr(P2)
  S[3,3] = "e01 + e12"
  
  S[4,1] = " P3             : "
  S[4,2] = toStr(P3)
  S[4,3] = "e01 + e20 + e12"
  
  S[5,1] = " line0          : "
  S[5,2] = toStr(line0)
  S[5,3] = "-e2"
  
  S[6,1] = " line1          : "
  S[6,2] = toStr(line1)
  S[6,3] = "e0 - e2"
  
  S[7,1] = " intersection   : "
  S[7,2] = toStr(x)
  S[7,3] = "-e20"
  
  S[8,1] = " toStr test 1   : "
  S[8,2] = toStr(tst1)
  S[8,3] = "-1 + e0"
  
  S[9,1] = " toStr test 2   : "
  S[9,2] = toStr(tst2)
  S[9,3] = "1 - e0"
  
  # print unit test results;
  nTest = size(S,1)
  for iTest = 1:nTest
   if S[iTest,2] == S[iTest,3]
    mark = " "
   else
    mark = "x"
    nError += 1
   end
   println("$mark" * S[iTest,1] * S[iTest,2])
  end

  return nError # return unit test results
 end
end # ripga2d utest()