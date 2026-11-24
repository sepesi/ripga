# ripgand.jl : reference implementation of
# Projective Geometric Algebra common for nD (i.e., 2D,3D,4D)
#
# This is a Julia port of bivector.net's C++ reference
# implementation of projective geometric algebra
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

function geoprodset(a::Matrix{Float32},b::Matrix{Float32})::Matrix{Float32}
 n = min(size(a,2),size(b,2))
 res = Matrix{Float32}(undef, size(a,1), n)
 for i=1:n
  res[:,i] = a[:,i] * b[:,i]
 end
 return res  # geometric product set (1 geometric product per column)
end

function outprodset(a::Matrix{Float32},b::Matrix{Float32})::Matrix{Float32}
 n = min(size(a,2),size(b,2))
 res = Matrix{Float32}(undef, size(a,1), n)
 for i=1:n
  res[:,i] = a[:,i] ^ b[:,i]
 end
 return res  # outer product set (1 outer product per column)
end

function regprodset(a::Matrix{Float32},b::Matrix{Float32})::Matrix{Float32}
 n = min(size(a,2),size(b,2))
 res = Matrix{Float32}(undef, size(a,1), n)
 for i=1:n
  res[:,i] = a[:,i] & b[:,i]
 end
 return res  # regressive product set (1 regressive product per column)
end

function inprodset(a::Matrix{Float32},b::Matrix{Float32})::Matrix{Float32}
 n = min(size(a,2),size(b,2))
 res = Matrix{Float32}(undef, size(a,1), n)
 for i=1:n
  res[:,i] = a[:,i] | b[:,i]
 end
 return res  # inner product set (1 inner product per column)
end

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
function epga(alpha::Number, B::Vector{Float32})
 s = (B * B)[1]
 if s == 0
  return 1 - alpha*B
 elseif s < 0
  return cos(alpha) - sin(alpha)*B
 else
  return cosh(alpha) - sinh(alpha)*B
 end
end

# bezlerp : Bezier Linear Interpolation
#
# arguments:
#  BCP1: Bezier Control Point matrix (each column is Euclidean dimension of Bezier curve Control Points)
#  TSM: Time Sample Matrix (stored in a Matrix{Float32}) for linear interpolation
#
# return value:
#  PC: Polynomial Coordinates (stored in a Matrix{Float32})
#
function bezlerp(BCP1::Matrix{Float32},TSM::Matrix{Float32})
  
  # Allocate buffer for return of Polynomial Coordinates in a Matrix{Float32}.
  nCoord = size(BCP1,2)
  nTS = size(TSM,2)
  PC = Matrix{Float32}(undef,nTS,nCoord)
  
  # for each coordinate in polynomial
  for iCoord = 1:nCoord
    PC[:,iCoord] = sum(TSM .* BCP1[:,iCoord],dims=1)'
  end
  
  return PC # each column of matrix holds a Euclidean dimension of polynomial coordinates
end

# bez2poly : convert cubic Bezier curves into 
#            rolled vectors of polynomial coordinates
#
# arguments:
#  BCP: Bezier Control Point matrix where each column contains all control points of Bezier curve
#  nT: number of sampling periods in each Bezier curve segment
#
# return value:
#  Vector{Vector{Float32}} where each "column" is a set of polynomial
#    coordinates rolled into a single vector
# 
function bez2poly(BCP::Matrix{Float32}, nT::Int)

  # find start point and length of each contiguous Bezier curve
  Status = trunc.(Int16, BCP[end,:])
  Istart = findall(x -> (x & 1) == 1, Status)
  push!(Istart, length(Status)+1) # append data sentinel
  Ilength = diff(Istart)
  pop!(Istart) # remove data sentinel
  
  # initialize data
  TS = range(0f0, step=1f0/nT, length=nT) # define Time Samples
  TSM = Matrix{Float32}([(1 .- TS).^3 (3*TS) .* (1 .- TS).^2 3*TS.^2 .* (1 .- TS) TS.^3]') # TS Matrix
  VRP = Vector{Vector{Float32}}()   # initialize Vector of Rolled Polynomials
  nCoord = (size(BCP,1) & ~1) >> 2  # discount appended status field, divide by # of control points (4)
  # println("TSM: $TSM")

  # for each set of closed polynomials in BCP
  nPoly = length(Istart)
  for iPoly = 1:nPoly
    # println("iPoly $iPoly> start at $(Istart[iPoly]) with length $(Ilength[iPoly])")
    PC = Matrix{Float32}(undef, 0, 2) # initialize Polynomial Coordinates

    # for each cubic Bezier curve within the closed polynomial
    nCol = Ilength[iPoly]
    for iCol = 1:nCol
       # linear interpolate between control points of single Bezier curve (i.e., single column)
       iBez1 = Istart[iPoly] + iCol - 1
       BCP1 = reshape(BCP[1:end-1,iBez1],:,nCoord)
       # println("iBez1 $iBez1; iCol $iCol> $BCP1")
       PC1 = bezlerp(BCP1, TSM)
       # println("Bez curve $iBez1 matrix> $PC1")
       PC = vcat(PC, PC1)
	end
	
	# roll up PC (Polynomial Coordinates) into a single vector and push it onto VRP
	push!(VRP, PC[:])
  end
  
  return VRP
end

# calculate dual of basis
# arg:
# 1) V::Vector{String} is the first column of basis
#
# usage:
# basis_dual(basis[:,1])
#
# returns matrix with 2 columns of Strings:
# 1) V (the original basis)
# 2) the dual (!) of the basis
#
function basis_dual(V::Vector{String})::Matrix{String}
  # expand the basis so each element covers the entire space
  nBasis = length(V)
  nChar = length(V[end]) - 1 # -1 avoids counting the leading 'e'
  BC = chop.(V, head=1, tail=0) # Basis Chopped (to remove leading 'e')
  IS = BC .* reverse(BC) # Index Strings

  # convert IS from Vector{String} to Matrix{Int} using Julia comprehension
  M = [parse(Int, str[j]) for str in IS, j in 1:nChar] # indices Matrix
  M = [-ones(Int, nBasis) M] # prepend data sentinel (for simpler sorting)
  nChar += 1 # account for added data sentinel
  
  # count left shifts for each expanded element in basis
  #  while performing insertion sort
  for iBasis = 1:nBasis
    lshift = 0
    for iChar = 3:nChar
      v = M[iBasis,iChar]
      jChar = iChar - 1
      while v < M[iBasis,jChar]
        M[iBasis,jChar+1] = M[iBasis,jChar]
        jChar -= 1
        lshift += 1
      end
      M[iBasis,jChar+1] = v
    end
    M[iBasis,1] = lshift # replace data sentinel with left shift count
  end
  
  # use Julia comprehension to generate sign strings based upon number of left shifts
  S = [v&0x1 > 0 ? "-" : "" for v in M[:,1]]
  
  return [V S.*reverse(V)] # return basis (col 1) and its dual (col 2)
end

# calculate reverse of basis
# arg:
# 1) V::Vector{String} is the first column of basis
#
# usage:
# basis_reverse(basis[:,1])
#
# returns matrix with 2 columns of Strings:
# 1) V (the original basis)
# 2) the reverse (~) of the basis
#
function basis_reverse(V::Vector{String})::Matrix{String}
  # generate reverse position matrix
  nBasis = length(V)
  nChar = length(V[end]) # -1 removes 'e', +1 for data sentinel
  R = zeros(Int,1,nChar) # empty Row
  copyto!(R,nChar:-1:1)
  R[1] = -1 # put data sentinel (for simpler sort) in leftmost position
  maxLength = nChar - 1
  M = repeat(R, maxLength)
  
  # sort reverse position matrix for each possible length
  for iLength = 1:maxLength
    lshift = 0
    for iChar = 3:iLength+1
      v = M[iLength,iChar]
      jChar = iChar - 1
      while v < M[iLength,jChar]
        M[iLength,jChar+1] = M[iLength,jChar]
        jChar -= 1
        lshift += 1
      end
      M[iLength,jChar+1] = v
    end
    M[iLength,1] = lshift # replace data sentinel with left shift count
  end
  
  # define prepending sign string based upon left shift count
  S = fill("",nBasis)
  for iBasis = 2:nBasis
    n = length(V[iBasis]) - 1 # -1 removes leading each
	(M[n,1]&0x1 > 0) && (S[iBasis] = "-")
  end
  
  return [V S.*V] # return basis (col 1) and its reverse (col 2)
end

# convert multivector fields to string
function toStr(V::Vector{Float32})
 nField = size(V,1)&~1 # clearing LSB ignores status field
 nNZField = 0
 S = String[]
 for iField = 1:nField
  # if field/component is not empty
  if V[iField] != 0
   # if not the first nonzero field
   if nNZField != 0
    # if nonzero field is negative
    if V[iField] < 0
     if V[iField] == -1
      push!(S, @sprintf(" - %s",
       basis[iField,1]))
     else
      push!(S, @sprintf(" - %0.7g%s",
       abs(V[iField]), basis[iField,1]))
     end
    # else nonzero field is positive
    else
     if V[iField] == 1
      push!(S, @sprintf(" + %s",
       basis[iField,1]))
     else
      push!(S, @sprintf(" + %0.7g%s",
       V[iField], basis[iField,1]))
     end
    end
   # else the first nonzero field
   else
    if V[iField] == 1
     push!(S, @sprintf("%s",
      (iField==1) ? "1" : basis[iField,1]))
    elseif V[iField] == -1
     push!(S, @sprintf("-%s",
      (iField==1) ? "1" : basis[iField,1]))
    elseif iField == 1
     push!(S, @sprintf("%0.7g",
      V[iField]))
    else
     push!(S, @sprintf("%0.7g%s",
      V[iField], basis[iField,1]))
    end
   end
   nNZField += 1
  end
 end
 (length(S) == 0) && push!(S, "0")
 return string(S...)
end # toStr()

# convert GA math syntax string to GA programming syntax expression
macro ga_str(str)
 C = collect(str)
 n = length(C)
 for i = 1:n
  if C[i] == ' '  # \thinspace for geometric product
   C[i] = '*'
  elseif C[i] == '∧' # \wedge for outer product
   C[i] = '^'
  elseif C[i] == '∨' # \vee for regressive product
   C[i] = '&'
  elseif C[i] == '⋅' # \cdot for inner product
   C[i] = '|'
#= simplify by avoiding postfix operator
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
=#
  end # ifs for conversion of math operators
 end # for each character in string
 return esc(Meta.parse(String(C)))
end # ga_str()

# dump expression
# args:
#  x: Expr::expression
#  d: Integer::depth
function dumpexpr(x::Expr, d::Integer=0)
 n = 0
 d += 2
 for a in x.args
  n += 1
  t = typeof(x.args[n])
  println("$(" "^d)$n> $a ($t)")
  (t == Expr) && m2p(x.args[n], d)
 end # for each arg
end # dumpexpr()

# switch dimension of geometry workspace
function xdimension(nD::Int64, isVerbose::Bool=false)
 # if specified dimension is current dimension
 if 2^(nD+1) == size(basis,1)
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
 nBasis = size(basis,1)
 M = Float32.(rand(1:9,nBasis,3))
 M[:,3] = M[:,1] * M[:,2]
# M[:,3] = normalize(M[:,1])
 display([1:nBasis basis M])
 println("expression to copy to bivector.net evaluator:")
 println("(" * toStr(M[:,1]) * ") * (" *
  toStr(M[:,2]) * ")")
end # rtest()
