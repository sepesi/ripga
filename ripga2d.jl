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

# define basis multivectors
nField = 2^3 # 3 = 2D + 1 dimensions
e0 =   zeros(Float32, nField); e0[2] = 1
e1 =   zeros(Float32, nField); e1[3] = 1
e2 =   zeros(Float32, nField); e2[4] = 1
e01 =  zeros(Float32, nField); e01[5] = 1
e20 =  zeros(Float32, nField); e20[6] = 1
e12 =  zeros(Float32, nField); e12[7] = 1
e012 = zeros(Float32, nField); e012[8] = 1

# geometric product
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
end # geometric product (*)

# regressive product: vee operator (&, \vee)
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
end # regressive product; vee operator (&, \vee)

# inner product: dot operator (|)
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
end # inner product (|)

# outer product; wedge operator (^)
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
end # outer product; wedge operator (^)

# reverse operator (~)
function Base.:~(a::Vector{Float32}) # reverse operator
 res = copy(a)
 res[5:8] .*= -1
 return res
end

# conjugate
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
# - flgSimplify set to true diverges from the text output by pga2d.cpp
# usage notes:
# - @time utest(1) checks on whether the unit test output of ripga2d.jl
# exactly matches the unit test output of pga2d.cpp. The comparison
# ends with the printing of the number of tests in the unit test that
# don't match. 0 indicates success of the unit test.
# - @time utest(1,true) outputs a slightly simplified version of the 
# unit test output that does not match the unit test output by pga2d.cpp.
# - @time utest(1,true,true) same as the above case except that
# the geometric objects are calculated with math syntax.
# - @btime utest() is a test for execution speed of ripga2d.jl.
#   (NOTE: requires using BenchmarkTools)
function utest(nLoop=100,
  flgSimplify::Bool=false,
  flgMathSyntax::Bool=false)

 # allocate some multivectors
 nField = length(basis)
 P0 = Vector{Float32}(undef,nField)
 P1 = Vector{Float32}(undef,nField)
 P2 = Vector{Float32}(undef,nField)
 P3 = Vector{Float32}(undef,nField)
 line0 = Vector{Float32}(undef,nField)
 line1 = Vector{Float32}(undef,nField)
 tst1 = Vector{Float32}(undef,nField)
 tst2 = Vector{Float32}(undef,nField)
 x = Vector{Float32}(undef,nField)

 for iLoop = 1:nLoop
  P0 = point(0,0)
  P1 = point(1,0)
  P2 = point(0,1)
  P3 = point(1,1)

  if flgMathSyntax == false
   # geometric objects in programming syntax
   (nLoop == 1) && println("  # calculated with programming syntax")
   line0 = P0 & P1
   line1 = P2 & P3
   x = line0 ^ line1
   
   tst1 = e0 - 1
   tst2 = 1 - e0
 
  else # flgMathSyntax == true
   # geometric objects in math syntax
   (nLoop == 1) && println("  # calculated with math syntax")
   line0 = ga"P0 ∨ P1"
   line1 = ga"P2 ∨ P3"
   x = ga"line0 ^ line1"
   
   tst1 = e0 - 1
   tst2 = 1 - e0
  end # flgMathSyntax
 end # iLoop

 # if verbose/slow output of unit test results wanted
 if nLoop == 1
  nError = 0
  toStr2 = flgSimplify ? toStr : toStr1

  S = Matrix{String}(undef,9,3) # 3 columns:
  S[1,1] = " P0           : "   #  1) label
  S[1,2] = toStr2(P0)           #  2) toStr() or toStr1()
  S[1,3] = flgSimplify ?        #  3) expected string
   "e12" :
   "1e12"

  S[2,1] = " P1           : "
  S[2,2] = toStr2(P1)
  S[2,3] = flgSimplify ?
   "e20 + e12" :
   "1e20 + 1e12"
  
  S[3,1] = " P2           : "
  S[3,2] = toStr2(P2)
  S[3,3] = flgSimplify ?
   "e01 + e12" :
   "1e01 + 1e12"
  
  S[4,1] = " P3           : "
  S[4,2] = toStr2(P3)
  S[4,3] = flgSimplify ?
   "e01 + e20 + e12" :
   "1e01 + 1e20 + 1e12"
  
  S[5,1] = " line0        : "
  S[5,2] = toStr2(line0)
  S[5,3] = flgSimplify ?
   "-e2" :
   "-1e2"
  
  S[6,1] = " line1        : "
  S[6,2] = toStr2(line1)
  S[6,3] = flgSimplify ?
   "e0 - e2" :
   "1e0 + -1e2"
  
  S[7,1] = " intersection : "
  S[7,2] = toStr2(x)
  S[7,3] = flgSimplify ?
   "-e20" :
   "-1e20"
  
  S[8,1] = " toStr test 1 : "
  S[8,2] = toStr2(tst1)
  S[8,3] = flgSimplify ?
   "-1 + e0" :
   "-1 + 1e0"
  
  S[9,1] = " toStr test 2 : "
  S[9,2] = toStr2(tst2)
  S[9,3] = flgSimplify ?
   "1 - e0" :
   "1 + -1e0"

  # print unit test results
  #  'x' in first column denotes tests with errors
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
end # ripga2d utest()
