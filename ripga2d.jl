# ripga2d.jl : reference implementation of 
# Projective Geometric Algebra for 2D
#
# This is a Julia port of bivector.net's C++ reference
# implementation of projective geometric algebra
# available at https://bivector.net/tools.html
#
using Printf

# define multivector basis names
basis = [ 															# iField
 "1"    "#1: scalar (specified as eu in vector form)" 				# 1
 "e0"   "#2: ideal line (line at infinity, encloses the 2D space)"	# 2
 "e1"   "#3: y-axis line (i.e., the x=0 line)"						# 3
 "e2"   "#4: x-axis line (i.e., the y=0 line)"						# 4
 "e01"  "#5: ideal point in y-direction"							# 5
 "e20"  "#6: ideal point in x-direction"							# 6
 "e12"  "#7: Euclidean point at origin (x=0,y=0)"					# 7
 "e012" "#8: pseudoscalar (the entire 2D space)"]					# 8

# define basis multivectors
nField = 2^3+1 # 3 = 2 dimensions + extra dimension; trailing +1 is a status field 
eu =   zeros(Float32, nField); eu[1] = 1
e0 =   zeros(Float32, nField); e0[2] = 1
e1 =   zeros(Float32, nField); e1[3] = 1
e2 =   zeros(Float32, nField); e2[4] = 1
e01 =  zeros(Float32, nField); e01[5] = 1
e20 =  zeros(Float32, nField); e20[6] = 1
e12 =  zeros(Float32, nField); e12[7] = 1
e012 = zeros(Float32, nField); e012[8] = 1

# alternative basis index ordering
e10 = -e01
e02 = -e20
e21 = -e12

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

# dual operator (!)
function Base.:!(a::Vector{Float32})::Vector{Float32}
 [reverse(a[1:end-1]); a[end]] # keep status field at end
end

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

# index to grade lookup table
#
# Construction of 2D PGA index to grade lookup table (I2G):
#  Z = zeros(Int,8) # 8 is nBasis for 2D PGA
#  B = binomial.(3,0:3) # 3 is nD + 1, 0:3 is range of grades
#  BC = cumsum(B)[1:end-1] .+ 1 # B Cumulated
#  Z[BC] .= 1;
#  I2G = cumsum(Z)
#
function grade(i::Int)
 [0;1;1;1;2;2;2;3][i]
end

# convert Euclidean coordinates to PGA expression
function point(
 x::Number,
 y::Number)::Vector{Float32}
 return x*e20 + y*e01 + e12
end
function point(M::Matrix{Float32})::Matrix{Float32}
 nPoint = size(M,2) # M is coordinate Matrix where each column is a point
 res = Matrix{Float32}(undef, 8+1, nPoint) # nBasis is 8, +1 for appended status
 for iPoint=1:nPoint
  res[:,iPoint] =
   M[1,iPoint]*e20 +
   M[2,iPoint]*e01 + e12
 end
 return res # each column of result matrix is a PGA expression of a point
end

# convert PGA expression to Euclidean coordinates
function toCoord(V::Vector{Float32})
 res = Vector{Float32}(undef, 2) # nD = 2
 res[1] = V[6] # e20 element is x component
 res[2] = V[5] # e01 element is y component
 return res
end
function toCoord(M::Matrix{Float32})::Matrix{Float32}
 nPoint = size(M,2) # M is PGA Matrix, each column is PGA expression of point
 res = Matrix{Float32}(undef, 2, nPoint) # nD = 2
 for iPoint=1:nPoint
  res[1,iPoint] = M[6,iPoint] # e20 element is x component
  res[2,iPoint] = M[5,iPoint] # e01 element is y component
 end
 return res # each column of result matrix is a coordinate of a point
end

# unit test
# arguments:
# - nLoop repeats a section of the PGA calculations for benchmarking
# - flgMathSyntax set to true tests ga macro
# usage notes:
# - @time utest(1) checks on whether the unit test
# exactly matches the expected output. The
# comparison ends with the printing of the
# number of tests in the unit test that don't match.
# 0 indicates success of the unit test.
# - @time utest(1,true) calculates the geometric
# objects using math syntax.
# - @btime utest() is a test for execution speed of ripga2d.jl.
#   (NOTE: requires using BenchmarkTools)
function utest(nLoop=100,
  flgMathSyntax::Bool=false)

 # allocate some multivectors
 nField = size(basis,1)+1 # +1 is for status field
 P0    = Vector{Float32}(undef,nField)
 P1    = Vector{Float32}(undef,nField)
 P2    = Vector{Float32}(undef,nField)
 P3    = Vector{Float32}(undef,nField)
 line0 = Vector{Float32}(undef,nField)
 line1 = Vector{Float32}(undef,nField)
 tst1  = Vector{Float32}(undef,nField)
 tst2  = Vector{Float32}(undef,nField)
 x     = Vector{Float32}(undef,nField)

 B =  [eu e0 e1 e2 e01 e20 e12 e012]
 BR = [e012 e12 e20 e01 e2 e1 e0 eu]
 BBR = geoprodset(B[1:end-1,:],BR[1:end-1,:])

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

  S = Matrix{String}(undef,17,3) # 3 columns:
  S[1,1] = " P0           : "   #  1) label
  S[1,2] = toStr(P0)            #  2) toStr() or toStr1()
  S[1,3] = "e12"		        #  3) expected string

  S[2,1] = " P1           : "
  S[2,2] = toStr(P1)
  S[2,3] = "e20 + e12"
  
  S[3,1] = " P2           : "
  S[3,2] = toStr(P2)
  S[3,3] = "e01 + e12"
  
  S[4,1] = " P3           : "
  S[4,2] = toStr(P3)
  S[4,3] = "e01 + e20 + e12"
  
  S[5,1] = " line0        : "
  S[5,2] = toStr(line0)
  S[5,3] = "-e2"
  
  S[6,1] = " line1        : "
  S[6,2] = toStr(line1)
  S[6,3] = "e0 - e2"
  
  S[7,1] = " intersection : "
  S[7,2] = toStr(x)
  S[7,3] = "-e20"
  
  S[8,1] = " !!e0         : "
  S[8,2] = toStr(!!e0)
  S[8,3] = "e0"
  
  S[9,1] = " !!!!e0       : "
  S[9,2] = toStr(!!!!e0)
  S[9,3] = "e0"
  
  S[10,1]= " toStr test 1 : "
  S[10,2]= toStr(tst1)
  S[10,3]= "-1 + e0"
  
  S[11,1]= " toStr test 2 : "
  S[11,2]= toStr(tst2)
  S[11,3]= "1 - e0"

  S[12,1]= " point test   : "
  S[12,2]= toStr(point(5,6))
  S[12,3]= "6e01 + 5e20 + e12"
  
  S[13,1]= " toCoord      : "
  S[13,2]= string(toCoord(point(5,6)))
  S[13,3]= "Float32[5.0, 6.0]"
  
  S[14,1]= " point test 2 : "
  S[14,2]= string(toCoord(point([5f0 10f0; 6f0 11f0])))
  S[14,3]= "Float32[5.0 10.0; 6.0 11.0]"
  
  S[15,1]= " BBR[end,:]   : "
  S[15,2]= toStr(BBR[end,:])
  S[15,3]= "1 + e0 + e1 + e2 + e01 + e20 + e12 + e012"

  S[16,1]= " min(ZBBR)    : "
  S[16,2]= string(minimum(BBR[1:end-1,:][:]))
  S[16,3]= "0.0"

  S[17,1]= " max(ZBBR)    : "
  S[17,2]= string(maximum(BBR[1:end-1,:][:]))
  S[17,3]= "0.0"

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
