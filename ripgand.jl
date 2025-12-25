# ripgand.jl : reference implementation of
# Projective Geometric Algebra common for nD
#
# This is a Julia port of pga3d.cpp, bivector.net's C++
# reference implementation of projective geometric algebra
# available at https://bivector.net/tools.html
#
using Printf

function Base.:*(a::Number, b::Vector{Float32})::Vector{Float32}
 res = copy(b)
 res .*= Float32(a)
 return res
end

function cumgprod(a::Matrix{Float32})::Matrix{Float32}
 res = similar(a)
 res[:,1] = a[:,1]
 nCol = size(a,2)
 for iCol = 2:nCol
  res[:,iCol] = res[:,iCol-1] * a[:,iCol]
 end
 return res  # cumulative geometric product
end    # similar to cumulative product cumprod()

function Base.:&(a::Matrix{Float32},b::Matrix{Float32})::Matrix{Float32}
 nCol = size(a,2)
 if nCol >= size(b,2)
  res = similar(a)
 else
  res = similar(b)
  nCol = size(b,2)
 end
 for iCol = 1:nCol
  res[:,iCol] = a[:,iCol] & b[:,iCol]
 end
 return res
end

function Base.:^(a::Matrix{Float32},b::Vector{Float32})::Matrix{Float32}
 res = similar(a)
 nCol = size(res,2)
 for iCol = 1:nCol
  res[:,iCol] = a[:,iCol] ^ b
 end
 return res
end

function Base.:^(a::Vector{Vector{Float32}},b::Vector{Float32})::Matrix{Float32}
 nCol = length(a)
 res = Matrix{Float32}(undef,16,nCol)
 for iCol = 1:nCol
  res[:,iCol] = a[iCol] ^ b
 end
 return res
end

function Base.:^(a::Vector{Float32},b::Matrix{Float32})::Matrix{Float32}
 res = similar(b)
 nCol = size(res,2)
 for iCol = 1:nCol
  res[:,iCol] = a ^ b[:,iCol]
 end
 return res
end

function Base.:^(a::Vector{Float32},b::Vector{Vector{Float32}})::Matrix{Float32}
 nCol = length(b)
 res = Matrix{Float32}(undef,16,nCol)
 for iCol = 1:nCol
  res[:,iCol] = a ^ b[iCol]
 end
 return res
end

function Base.:+(a::Vector{Float32},b::Number)::Vector{Float32}
 res = copy(a)
 res[1] += Float32(b)
 return res
end

function Base.:+(a::Number,b::Vector{Float32})::Vector{Float32}
 res = copy(b)
 res[1] += Float32(a)
 return res
end

function Base.:-(a::Vector{Float32},b::Number)::Vector{Float32}
 res = copy(a)
 res[1] -= Float32(b)
 return res
end

function Base.:-(a::Number,b::Vector{Float32})::Vector{Float32}
 res = copy(.-b)
 res[1] += Float32(a)
 return res
end

function Base.:>>>(a::Vector{Float32},b::Vector{Float32})::Vector{Float32}
 return a * b * ~a
end

function Base.:>>>(a::Vector{Float32},b::Matrix{Float32})
 res = similar(b)
 nCol = size(b,2)
 for iCol = 1:nCol
  res[:,iCol] = a * b[:,iCol] * ~a
 end
 return res
end

function Base.:!(a::Vector{Float32})::Vector{Float32}
 return reverse(a)
end

function Base.:!(a::Matrix{Float32})
 return mapslices(!, a, dims=1)
end

# mask off all vectors except those with specified grade (k)
function grade(a::Vector{Float32},k::Int64)
 # if specified grade out of range
 nBasis = length(a)
 nD = Int(log2(nBasis)) - 1
 nGrade = nD + 2
 res = zeros(Float32,nBasis)
 if (k < 0) || (k >= nGrade)
  return res
 end
 
 # generate CN: ending index of each grade
 N = zeros(Int32,nGrade) # basis vectors with grade
 for iGrade = 1:nGrade
  N[iGrade] = binomial(nD+1,iGrade-1)
 end
 CN = cumsum(N)
 
 # copy vectors with specified grade
 i2 = CN[k+1]
 i0 = i2 - N[k+1] + 1
 res[i0:i2] = a[i0:i2]
 return res
end

# convert Euclidean coordinates to PGA expression
function point(x::Number, y::Number)::Vector{Float32}
 return x*e20 + y*e01 + e12
end
function point(x::Number, y::Number, z::Number)::Vector{Float32}
 return x*e032 + y*e013 + z*e021 + e123
end
function point(
 x::Number,
 y::Number,
 z::Number,
 w::Number)::Vector{Float32}
 return x*e0234 +
  y*e0134 +
  z*e0124 +
  w*e0123 +
  e1234
end
function point(V::Vector{Float32},
 nBasis::Int64=0)::Vector{Float32}
 
 if nBasis == 0
  nBasis = length(basis)
 end
 if nBasis == 8 # 2D
  return V[1]*e20 +
   V[2]*e01 +
   e12
 elseif nBasis == 16 # 3D
  return V[1]*e032 +
   V[2]*e013 +
   V[3]*e021 +
   e123
 elseif nBasis == 32 # 4D
  return V[1]*e0234 +
   V[2]*e0134 +
   V[3]*e0124 +
   V[4]*e0123 +
   e1234
 end
end
function point(M::Matrix{Float32})::Matrix{Float32}
 nCol = size(M,2)
 res = Matrix{Float32}(undef, length(basis), nCol)
 for iCol = 1:nCol
  res[:,iCol] = point(M[:,iCol])
 end
 return res
end

# convert PGA expressions to Euclidean coordinates
function toCoord(M::Matrix{Float32},
 keepIdeal::Bool=false,
 nBasis::Int64=0)
 
 if nBasis == 0
  nBasis = length(basis)
 end
 
 nB1 = nBasis - 1
 nD = Int(log2(nBasis)) - 1
 MC = Matrix{Float32}(undef, (nD,size(M,2)))
 B = M[nB1,:] .!= 0
 S = sign.(M[nB1,:])
 S[.!B] .= 1.0
 for iRow = 1:nD
  MC[iRow,:] = M[nB1-iRow,:] .* S # account for orientation
 end
 return keepIdeal ? MC : MC[:,B]
end
function toCoord(V::Vector{Float32})
 nBasis = size(V,1)
 nB1 = nBasis - 1
 nD = Int(log2(nBasis)) - 1
 res = Vector{Float32}(undef, nD)
 for iRow = 1:nD
  res[iRow] = V[nB1-iRow]
 end
 return res
end

function rotor(angle::Number, line::Vector{Float32})::Vector{Float32}
 return Float32(cos(angle/2)) +
  Float32(sin(angle/2))*normalize(line)
end

function translator(dist::Number, line::Vector{Float32})::Vector{Float32}
 return 1 + Float32(dist/2) * line
end

function norm(a::Vector{Float32})
 return sqrt(abs((a * conjugate(a))[1]))
end

function norm(a::Matrix{Float32})
 return mapslices(norm, a, dims=1)
end

function normalize(a::Vector{Float32})::Vector{Float32}
 return a ./ norm(a)
end

# exponential, restricted to case B ^ B = 0
function E(alpha::Number, B::Vector{Float32})
 s = (B * B)[1]
 if s == 0
  return 1 + alpha*B
 elseif s < 0
  return cos(alpha) + sin(alpha)*B
 else
  return cosh(alpha) + sinh(alpha)*B
 end
end

# convert multivector fields to string
# usage notes:
# - Use toStr1 to generate text that exactly
# matches the text generated by pga3d.cpp
# for example to unit test ripga.jl.
# - Use toStr to generate text that simplifies
# signs (e.g., a + -b simplifies to a - b)
# and avoids overloading the parsing of
# constants with exponentials to represent
# basis vectors (e.g., 1e1 with overloading
# simplifies to e1 without overloading).
function toStr1(V::Vector{Float32})
 nField = size(V,1)
 nNZField = 0
 S = String[]
 for iField = 1:nField
  if V[iField] != 0
   if nNZField != 0
    push!(S, @sprintf(" + %0.7g%s",
     V[iField], basis[iField]))
   else
    push!(S, @sprintf("%0.7g%s",
     V[iField],
     (iField==1) ? "" : basis[iField]))
   end
   nNZField += 1
  end
 end
 return string(S...)
end

function toStr(V::Vector{Float32})
 nField = size(V,1)
 nNZField = 0
 S = String[]
 for iField = 1:nField
  if V[iField] != 0
   if nNZField != 0
    if V[iField] < 0
     if V[iField] == -1
      push!(S, @sprintf(" - %s",
       basis[iField]))
     else
      push!(S, @sprintf(" - %0.7g%s",
       abs(V[iField]), basis[iField]))
     end
    else
     if V[iField] == 1
      push!(S, @sprintf(" + %s",
       basis[iField]))
     else
      push!(S, @sprintf(" + %0.7g%s",
       V[iField], basis[iField]))
     end
    end
   else
    if V[iField] == 1
     push!(S, @sprintf("%s",
      (iField==1) ? "1" : basis[iField]))
    elseif V[iField] == -1
     push!(S, @sprintf("-%s",
      (iField==1) ? "1" : basis[iField]))
    elseif iField == 1
     push!(S, @sprintf("%0.7g",
      V[iField]))
    else
     push!(S, @sprintf("%0.7g%s",
      V[iField], basis[iField]))
    end
   end
   nNZField += 1
  end
 end
 if length(S) == 0
  push!(S, "0")
 end
 return string(S...)
end # toStr()

# convert GA math syntax string to GA programming syntax expression
macro ga_str(str)
 C = collect(str)
 n = length(C)
 for i = 1:n
  if C[i] == ' '  # \thinspace for geometric product
   C[i] = '*'
  elseif C[i] == '∧' # \wedge for outer product
   C[i] = '^'
  elseif C[i] == '∨' # \vee for regressive product
   C[i] = '&'
  elseif C[i] == '·' # \cdotp for inner product
   C[i] = '|'
  elseif C[i] == '∗' # \ast for dual (suffix)
   j = i-1
   if C[j] == ')'
    nDepth = 1
    C[j+1] = C[j]
    j -= 1
    # shift from suffix to prefix of parentheses
    while (j > 0) && (nDepth > 0)
     if C[j] == ')'
      nDepth += 1
     elseif C[j] == '('
      nDepth -= 1
     end
     C[j+1] = C[j]
     j -= 1
    end
   else
    # shift from suffix to prefix of variable
    while (j > 0) && (isletter(C[j]) || isnumeric(C[j]) || C[j]=='_')
     C[j+1] = C[j]
     j -= 1
    end
   end
   C[j+1] = '!' # prefix '!'
  end
 end
 return esc(Meta.parse(String(C)))
end # ga_str()

# switch dimension of geometry workspace
function xdimension(nD::Int64, isVerbose::Bool=false)
 # if specified dimension is current dimension
 if 2^(nD+1) == length(basis)
  if isVerbose
   println("dimension is already $nD")
  end
  
 # else if specified dimension is implemented 
 elseif nD in [2, 3, 4]
  strFile = @sprintf("ripga%dd.jl", nD)
  include(strFile)
 
 # else specified dimension not yet implemented
 else
  if isVerbose
   println("xDimension error: dimension $nD not yet implemented")
  end
 end
end # xdimension()

# rtest: random test
function rtest()
 nBasis = length(basis)
 M = Float32.(rand(1:9,nBasis,3))
 M[:,3] = M[:,1] * M[:,2]
# M[:,3] = normalize(M[:,1])
 display([1:nBasis basis M])
 println("expression to copy to bivector.net evaluator:")
 println("(" * toStr(M[:,1]) * ") * (" *
  toStr(M[:,2]) * ")")
end # rtest()