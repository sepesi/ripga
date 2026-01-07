# ripga1d.jl : reference implementation of
# Projective Geometric Algebra for 1D
#
# This is a Julia port of bivector.net's C++ reference
# implementation of projective geometric algebra
# available at https://bivector.net/tools.html
#
using Printf

# define multivector basis names
# 0 denotes projective dimension (e.g., e0 * e0 = 0)
basis = [ # iField
 "1"  #  1 scalar vector (i.e., vector with only nonzero scalar term)
 "e0" #  2 grade 1 vectors
 "e1" #  3
 "e01"]# 4 pseudoscalar

# define basis multivectors
nField = 2^2 # 2 = 1D + 1 dimensions
e0 =  zeros(Float32, nField); e0[2] = 1
e1 =  zeros(Float32, nField); e1[3] = 1
e01 = zeros(Float32, nField); e01[4] = 1

# geometric product
function Base.:*(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=a[1]*b[1]+a[3]*b[3]
 res[2]=a[1]*b[2]+a[2]*b[1]-a[3]*b[4]+a[4]*b[3]
 res[3]=a[1]*b[3]+a[3]*b[1]
 res[4]=a[1]*b[4]+a[2]*b[3]-a[3]*b[2]+a[4]*b[1]
 return res
end # geometric product (*)

# regressive product: vee operator (&, \vee)
function Base.:&(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=a[1]*b[4]-a[2]*b[3]+a[3]*b[2]+a[4]*b[1]
 res[2]=a[2]*b[4]+a[4]*b[2]
 res[3]=a[3]*b[4]+a[4]*b[3]
 res[4]=a[4]*b[4]
 return res
end # regressive product; vee operator (&, \vee)

# inner product: dot operator (|)
function Base.:|(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=a[1]*b[1]+a[3]*b[3]
 res[2]=a[1]*b[2]+a[2]*b[1]-a[3]*b[4]+a[4]*b[3]
 res[3]=a[1]*b[3]+a[3]*b[1]
 res[4]=a[1]*b[4]+a[4]*b[1]
 return res
end # inner product (|)

# outer product; wedge operator (^)
function Base.:^(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 res = similar(a)
 res[1]=a[1]*b[1]
 res[2]=a[1]*b[2]+a[2]*b[1]
 res[3]=a[1]*b[3]+a[3]*b[1]
 res[4]=a[1]*b[4]+a[2]*b[3]-a[3]*b[2]+a[4]*b[1]
 return res
end # outer product; wedge operator (^)

# reverse operator (~)
function Base.:~(a::Vector{Float32}) # reverse operator
 res = copy(a)
 res[4] .*= -1
 return res
end

# conjugate
function conjugate(a::Vector{Float32})::Vector{Float32}
 res = copy(a)
 res[2:4] .*= -1
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
 a    = Vector{Float32}(undef,nField)
 b    = Vector{Float32}(undef,nField)
 res1 = Vector{Float32}(undef,nField)
 res2 = Vector{Float32}(undef,nField)
 res3 = Vector{Float32}(undef,nField)
 res4 = Vector{Float32}(undef,nField)
 res5 = Vector{Float32}(undef,nField)
 res6 = Vector{Float32}(undef,nField)
 res7 = Vector{Float32}(undef,nField)
 res8 = Vector{Float32}(undef,nField)
 res9 = Vector{Float32}(undef,nField)
 res10= Vector{Float32}(undef,nField)
 res11= Vector{Float32}(undef,nField)
 res12= Vector{Float32}(undef,nField)
 res13= Vector{Float32}(undef,nField)
 res14= Vector{Float32}(undef,nField)

 tst1 = Vector{Float32}(undef,nField)
 tst2 = Vector{Float32}(undef,nField)

 for iLoop = 1:nLoop

  if flgMathSyntax == false
   # geometric objects in programming syntax
   (nLoop == 1) && println("  # calculated with programming syntax")

   # note: multiplicationn operator needed to specify vectors
   a = 1 + 2*e0 + 3*e1 + 4*e01 # not "1 + 2e0 + 3e1 + 4e01"
   b = 2 + 4*e0 + 6*e1 + 8*e01 # not "2 + 4e0 + 6e1 + 8e01"

   # calculate results to be tested
   res1 =  (1*e0) * (1*e0) # = 0
   res2 =  (1*e1) * (1*e1) # = 1
   res3 =  (1*e0) ^ (1*e1) # = e01
   res4 =  !(1*e01)        # = 1
   res5 =  (1*e0) & (1*e1) # = -1
   res6 =  (1*e0) | (1*e1) # = 0
   res7 =  (1*e1) | (1*e1) # = 1
   res8 =  a + b   # = 3 + 6*e0 + 9*e1 + 12*e01
   res9 =  a - b   # = -1 - 2*e0 - 3*e1 - 4*e01
   res10 = a * b   # = 20 + 8*e0 + 12*e1 + 16*e01
   res11 = a ^ b   # = 2 + 8*e0 + 12*e1 + 16*e01
   res12 = a & b   # = 16 + 32*e0 + 48*e1 + 32*e01
   res13 = a | b   # = 20 + 8*e0 + 12*e1 + 16*e01
     
   tst1 = e0 - 1
   tst2 = 1 - e0
 
  else # flgMathSyntax == true
   # geometric objects in math syntax
   (nLoop == 1) && println("  # calculated with math syntax")

   # note: multiplicationn operator needed to specify vectors
   a = ga"1 + 2 e0 + 3 e1 + 4 e01" # not "1 + 2e0 + 3e1 + 4e01"
   b = ga"2 + 4 e0 + 6 e1 + 8 e01" # not "2 + 4e0 + 6e1 + 8e01"

   # calculate results to be tested
   res1 =  ga"(1 e0)   (1 e0)" # = 0
   res2 =  ga"(1 e1)   (1 e1)" # = 1
   res3 =  ga"(1 e0) ^ (1 e1)" # = e01
   res4 =  ga"(1 e01)∗"        # = 1
   res5 =  ga"(1 e0) ∨ (1 e1)" # = -1
   res6 =  ga"(1 e0) ⋅ (1 e1)" # = 0
   res7 =  ga"(1 e1) ⋅ (1 e1)" # = 1
   res8 =  ga"a + b"   # = 3 + 6*e0 + 9*e1 + 12*e01
   res9 =  ga"a - b"   # = -1 - 2*e0 - 3*e1 - 4*e01
   res10 = ga"a   b"   # = 20 + 8*e0 + 12*e1 + 16*e01
   res11 = ga"a ^ b"   # = 2 + 8*e0 + 12*e1 + 16*e01
   res12 = ga"a ∨ b"   # = 16 + 32*e0 + 48*e1 + 32*e01
   res13 = ga"a ⋅ b"   # = 20 + 8*e0 + 12*e1 + 16*e01

   tst1 = e0 - 1
   tst2 = 1 - e0
  end # flgMathSyntax
 end # iLoop

 # if verbose/slow output of unit test results wanted
 if nLoop == 1
  nError = 0
  toStr2 = flgSimplify ? toStr : toStr1

  S = Matrix{String}(undef,15,3) # 3 columns:
  S[1,1] = " res1         : "    #  1) label
  S[1,2] = toStr2(res1)          #  2) toStr() or toStr1()
  S[1,3] = flgSimplify ?         #  3) expected string
   "0" :
   "0"

  S[2,1] = " res2         : "
  S[2,2] = toStr2(res2)
  S[2,3] = flgSimplify ?
   "1" :
   "1"
  
  S[3,1] = " res3         : "
  S[3,2] = toStr2(res3)
  S[3,3] = flgSimplify ?
   "e01" :
   "1e01"
  
  S[4,1] = " res4         : "
  S[4,2] = toStr2(res4)
  S[4,3] = flgSimplify ?
   "1" :
   "1"
  
  S[5,1] = " res5         : "
  S[5,2] = toStr2(res5)
  S[5,3] = flgSimplify ?
   "-1" :
   "-1"
  
  S[6,1] = " res6         : "
  S[6,2] = toStr2(res6)
  S[6,3] = flgSimplify ?
   "0" :
   "0"
  
  S[7,1] = " res7         : "
  S[7,2] = toStr2(res7)
  S[7,3] = flgSimplify ?
   "1" :
   "1"
  
  S[8,1] = " res8         : "
  S[8,2] = toStr2(res8)
  S[8,3] = flgSimplify ?
   "3 + 6e0 + 9e1 + 12e01" :
   "3 + 6e0 + 9e1 + 12e01"
  
  S[9,1] = " res9         : "
  S[9,2] = toStr2(res9)
  S[9,3] = flgSimplify ?
   "-1 - 2e0 - 3e1 - 4e01" :
   "-1 + -2e0 + -3e1 + -4e01"
  
  S[10,1]= " res10        : "
  S[10,2]= toStr2(res10)
  S[10,3]= flgSimplify ?
   "20 + 8e0 + 12e1 + 16e01" :
   "20 + 8e0 + 12e1 + 16e01"
  
  S[11,1]= " res11        : "
  S[11,2]= toStr2(res11)
  S[11,3]= flgSimplify ?
   "2 + 8e0 + 12e1 + 16e01" :
   "2 + 8e0 + 12e1 + 16e01"
  
  S[12,1]= " res12        : "
  S[12,2]= toStr2(res12)
  S[12,3]= flgSimplify ?
   "16 + 32e0 + 48e1 + 32e01" :
   "16 + 32e0 + 48e1 + 32e01"
  
  S[13,1]= " res13        : "
  S[13,2]= toStr2(res13)
  S[13,3]= flgSimplify ?
   "20 + 8e0 + 12e1 + 16e01" :
   "20 + 8e0 + 12e1 + 16e01"
  
  S[14,1]= " toStr test 1 : "
  S[14,2]= toStr2(tst1)
  S[14,3]= flgSimplify ?
   "-1 + e0" :
   "-1 + 1e0"
  
  S[15,1]= " toStr test 2 : "
  S[15,2]= toStr2(tst2)
  S[15,3]= flgSimplify ?
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
end # ripga1d utest()
