# ripgand.jl : reference implementation of
# Projective Geometric Algebra common for nD (i.e., 2D,3D,4D)
#
# This is a Julia port of bivector.net's C++ reference
# implementation of projective geometric algebra
# available at https://bivector.net/tools.html
#
using Printf

# 2D cubic Bezier curves of left footprint
# NOTES:
# 1) There is one set of Bezier curve control points per column
# 2) Appended to the end of each column is a status field of flags:
#    bit 0: start of a group of columns (e.g., a group of adjacent curves)
#    bit 1: not currently used
#    bit 2: the group of columns makes a closed curve (the end connects to the beginning)
# 3) Each Bezier curve (i.e., each column) has 4 Bezier curve control points:
#    P0: start point
#    P1: mid control point
#    P2: mid control point
#    P3: end point
# 4) However, P3 is not stored because the endpoint of one segment is identical to the
#    start point of the next segment (i.e., Footprint2D[6,k] == Footprint2D[1,k+1] and
#    Footprint2D[7,k] == Footprint2D[2,k+1]) because the adjacent segments in a group
#    are connected.
# 5) Also, because adjacent curves in a group have G1 geometric continuity (i.e., the slope
#    stays the same across the boundary between two adjacent Bezier curves in a group) which
#    means that the previous curve's P2 is colinear with the new curve's P0 and P1.
# 6) If a Bezier curve is not closed, then there needs to be a final column for defining the
#    curve's endpoint.
# 7) The first 10 columns in the following matrix specify the Bezier control points of the
#    outline of the toes, starting with the big toe (columns 1 and 2) and ending with the little
#    toe (columns 9 and 10).
# 8) The remaining 7 columns (columns 11 through 17) specify the Bezier control points of the
#    outline of the sole of the foot.
#=
Footprint2D = Float32[  # uncentered 2D cubic Bezier control points of footprint
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  # column numbers
 80  82  59  60  43  45  29  32  17  21  48  22  16  28  78  74  64; # P0 x
212 250 212 242 205 233 195 220 181 202  12  44 128 180 201 132  46; # P0 y
 96  66  69  50  52  36  37  24  24  14  32  22  14  42  95  50  72; # P1 x
210 252 211 243 203 235 193 222 178 205  12  68 154 191 199 104  28; # P1 y
 98  64  70  49  54  34  40  21  28  10  22  18  15  62  96  52  64; # P2 x
248 214 241 213 231 207 218 197 199 184  24  98 170 203 157  72  12; # P2 y
0x5   0 0x5   0 0x5   0 0x5   0 0x5   0 0x5   0   0   0   0   0   0; # status flags
]
=#
Footprint2D = Float32[  # centered 2D cubic Bezier control points of footprint
#  1    2    3    4    5    6    7    8    9   10    11   12   13   14  15   16    17  # column numbers
  26   28    5    6  -11   -9  -25  -22  -37  -33    -6  -32  -38  -26  24   20    10  # P0 x
  80  118   80  110   73  101   63   88   49   70  -120  -88   -4   48  69    0   -86  # P0 y
  42   12   15   -4   -2  -18  -17  -30  -30  -40   -22  -32  -40  -12  41   -4    18  # P1 x
  78  120   79  111   71  103   61   90   46   73  -120  -64   22   59  67  -28  -104  # P1 y
  44   10   16   -5    0  -20  -14  -33  -26  -44   -32  -36  -39    8  42   -2    10  # P2 x
 116   82  109   81   99   75   86   65   67   52  -108  -34   38   71  25  -60  -120  # P2 y
 0x5    0  0x5    0  0x5    0  0x5    0  0x5    0   0x5    0    0    0   0    0     0  # status flags
]

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
  return 1 + alpha*B
 elseif s < 0
  return cos(alpha) + sin(alpha)*B
 else
  return cosh(alpha) + sinh(alpha)*B
 end
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
