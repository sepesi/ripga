# origami.jl : animate folding of kabuto (Samurai helmet)
using GLMakie
using GeometryBasics
using Images # for RGBA
include("ripga3d.jl")

# convert coordinate Point{3, Float32} to/from PGA point
function point(A::Point{3, Float32})::Vector{Float32}
 return A[1]*e032 + A[2]*e013 + A[3]*e021 + e123
end
function point(A::Vector{Point{3, Float32}})
 nCol = length(A)
 res = Matrix{Float32}(undef, (16,nCol))
 for iCol = 1:nCol
  res[:,iCol] = A[iCol][1]*e032 + A[iCol][2]*e013 +
   A[iCol][3]*e021 + e123
 end
 return res
end
function toPoint(A::Vector{Float32})
 return Point{3, Float32}(A[14],A[13],A[12])
end
function toPoint(A::Matrix{Float32})
 res = Point{3, Float32}[]
 nCol = size(A,2)
 for iCol = 1:nCol
  push!(res, Point{3, Float32}(A[14,iCol],
   A[13,iCol], A[12,iCol]))
 end
 return res
end

# fold
#
# NOTE:
# The mesh's vertices (SV) do not change, but the
# mesh's faces (SF) do change: new folds modify some 
# previous faces and generate new face(s). During a
# fold, some faces don't rotate and some do. To help
# with the bookkeeping, the faces that rotate during 
# the fold are typically put at the end of the face 
# vector. 
#
function fold(iFold::Int64,
 SV::Vector{Point{3, Float32}},   # Sheet Vertices (i/o)
 SF::Vector{TriangleFace{Int64}}, # Sheet Faces (i/o)
 PL::Vector{Float32},  # Pivot Line (output)
 MVI::Vector{Int64})   # Moving Vertex Indices (output)
 
 # initialize data
 nSweep = 1
 
 # fold 1: diagonal fold that puts 
 # vertex 1 on top of vertex 3
 if iFold == 1
  # calculate Pivot Line
  PL[:] = point(SV[2]) & point(SV[4])
  
  # generate new faces in mesh
  push!(SF, [2,3,4]) # 1.
  push!(SF, [1,2,4]) # 2. rotating
  
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 1)
  push!(MVI, 7)
  push!(MVI, 8)
  push!(MVI, 10)
  push!(MVI, 11)
  push!(MVI, 16)
  push!(MVI, 17)
  push!(MVI, 20)
  push!(MVI, 21)
  push!(MVI, 23)
  push!(MVI, 24)
  push!(MVI, 25)
  push!(MVI, 26)
  push!(MVI, 27)
  push!(MVI, 28)
 
 # fold 2: horizontal fold that puts
 # vertex 4 on top of vertex 3
 elseif iFold == 2
  # calculate Pivot Line
  PL[:] = point(SV[5]) & point(SV[6])
  
  # modify existing faces in mesh
  deleteat!(SF, 1); insert!(SF, 1, [2,3,5]) # 1.
  deleteat!(SF, 2); insert!(SF, 2, [1,2,5]) # 2.
  
  # generate new faces in mesh
  push!(SF, [3,6,5]) # 3.
  push!(SF, [1,5,8]) # 4.
  push!(SF, [4,5,6]) # 5. rotating
  push!(SF, [4,8,5]) # 6. rotating
  
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 4)
  push!(MVI, 14)
  push!(MVI, 15)
  push!(MVI, 16)
  push!(MVI, 19)
  push!(MVI, 20)
 
 # fold 3: vertical fold that puts
 # vertex 2 on top of vertex 3
 elseif iFold == 3
  # calculate Pivot Line
  PL[:] = point(SV[9]) & point(SV[5])
  
  # modify existing faces in mesh
  deleteat!(SF, 1); insert!(SF, 1, [9,3,5]) # 1.
  deleteat!(SF, 2); insert!(SF, 2, [1,5,7]) # 2.
  
  # generate new faces in mesh
  push!(SF, [2,9,5]) # 7. rotating
  push!(SF, [5,7,2]) # 8. rotating
  
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 2)
  push!(MVI, 12)
  push!(MVI, 17)
  push!(MVI, 18)
  push!(MVI, 21)
  push!(MVI, 22)
 
 # fold 4: diagonal fold that puts 
 # vertex 4 on top of vertex 5
 elseif iFold == 4
  # calculate Pivot Line
  PL[:] = point(SV[6]) & point(SV[9])
  
  # modify existing faces in mesh
  deleteat!(SF, 5); insert!(SF, 5, [5,6,14]) # 5.
  deleteat!(SF, 6); insert!(SF, 6, [5,8,14]) # 6.
  
  # generate new faces in mesh
  push!(SF, [4,14,6]) # 9. rotating
  push!(SF, [4,14,8]) #10. rotating
  
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 4)
  push!(MVI, 15)
  push!(MVI, 16)
  push!(MVI, 19)
  push!(MVI, 20)
 
 # fold 5: diagonal fold that puts 
 # vertex 2 on top of vertex 5
 elseif iFold == 5
  # calculate Pivot Line
  PL[:] = point(SV[6]) & point(SV[9])
  
  # modify existing faces in mesh
  deleteat!(SF, 7); insert!(SF, 7, [5,9,12]) # 7.
  deleteat!(SF, 8); insert!(SF, 8, [5,7,12]) # 8.
  
  # generate new faces in mesh
  push!(SF, [9,12,2]) #11. rotating
  push!(SF, [7,12,2]) #12. rotating
  
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 2)
  push!(MVI, 17)
  push!(MVI, 18)
  push!(MVI, 21)
  push!(MVI, 22)
 
 # fold 6: vertical crease
 elseif iFold == 6
  # calculate Pivot Line
  PL[:] = point(SV[14]) & point(SV[16])
  
  # modify existing faces in mesh
  deleteat!(SF, 9); insert!(SF, 9, [6,14,15]) # 9.
  deleteat!(SF, 10); insert!(SF, 10, [8,14,16]) # 10.
  
  # generate new faces in mesh
  push!(SF, [4,14,16]) #13. rotating
  push!(SF, [4,14,15]) #14. rotating
  
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 4)
  push!(MVI, 19)
  push!(MVI, 20)
  
  # sweep back and forth
  nSweep = 2
 
 # fold 7: horizontal crease
 elseif iFold == 7
  # calculate Pivot Line
  PL[:] = point(SV[18]) & point(SV[12])
  
  # modify existing faces in mesh
  deleteat!(SF, 11); insert!(SF, 11, [9,18,12]) # 11.
  deleteat!(SF, 12); insert!(SF, 12, [7,17,12]) # 12.
  
  # generate new faces in mesh
  push!(SF, [12,2,18]) #15. rotating
  push!(SF, [12,2,17]) #16. rotating
 
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 2)
  push!(MVI, 21)
  push!(MVI, 22)
  
  # sweep back and forth
  nSweep = 2
 
 # fold 8: near diagonal fold that aligns 
 # vertex 4 with vertices 14 and 16
 elseif iFold == 8
  # calculate Pivot Line
  PL[:] = point(SV[14]) & point(SV[20])
  
  # modify existing faces in mesh
  deleteat!(SF, 13); insert!(SF, 13, [20,14,16]) # 13.
  deleteat!(SF, 14); insert!(SF, 14, [19,14,15]) # 14.
  
  # generate new faces in mesh
  push!(SF, [4,14,20]) #17. rotating
  push!(SF, [4,14,19]) #18. rotating
 
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 4)
 
 # fold 9: near diagonal fold that aligns 
 # vertex 2 with vertices 12 and 17
 elseif iFold == 9
  # calculate Pivot Line
  PL[:] = point(SV[22]) & point(SV[12])
  
  # modify existing faces in mesh
  deleteat!(SF, 15); insert!(SF, 15, [22,12,18]) # 15.
  deleteat!(SF, 16); insert!(SF, 16, [21,12,17]) # 16.
  
  # generate new faces in mesh
  push!(SF, [2,12,22]) #19. rotating
  push!(SF, [2,12,21]) #20. rotating
 
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 2)
 
 # fold 10: diagonal fold that puts 
 # vertex 1 on top of vertex 10
 elseif iFold == 10
  # calculate Pivot Line
  PL[:] = point(SV[24]) & point(SV[23])
  
  # modify existing faces in mesh
  deleteat!(SF, 1); insert!(SF, 1, [9,13,5]) # 1.
  deleteat!(SF, 2); insert!(SF, 2, [7,11,5]) # 2.
  deleteat!(SF, 3); insert!(SF, 3, [13,6,5]) # 3.
  deleteat!(SF, 4); insert!(SF, 4, [11,5,8]) # 4.
  
  # generate new faces in mesh
  push!(SF, [9,3,13])  # 21.
  push!(SF, [3,6,13])  # 22.
  push!(SF, [7,23,11]) # 23.
  push!(SF, [23,25,11]) # 24.
  push!(SF, [23,1,25]) # 25.
  push!(SF, [8,24,11]) # 26.
  push!(SF, [24,25,11]) # 27.
  push!(SF, [24,1,25]) # 28.
  
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 1)
  push!(MVI, 26)
  push!(MVI, 27)
  push!(MVI, 28)
 
 # fold 11: diagonal fold to finish front brim
 elseif iFold == 11
  # calculate Pivot Line
  PL[:] = point(SV[8]) & point(SV[7])
  
  # modify existing faces in mesh
  deleteat!(SF, 25); insert!(SF, 25, [26,1,28]) # 25.
  deleteat!(SF, 28); insert!(SF, 28, [27,1,28]) # 28.
  
  # generate new faces in mesh
  push!(SF, [23,26,28]) # 29.
  push!(SF, [28,25,23]) # 30.
  push!(SF, [24,27,28]) # 31.
  push!(SF, [28,25,24]) # 32.

  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 23)
  push!(MVI, 24)
  push!(MVI, 25)
 
 # fold 12: diagonal fold to finish back brim
 elseif iFold == 12
  # calculate Pivot Line
  PL[:] = point(SV[7]) & point(SV[8])
  
  # identify the rotating vertices
  empty!(MVI)
  push!(MVI, 3)
 end
 return nSweep
end

function origami()
 # initialize data about sheet
 nFold = 12
 x = 1 + tan(pi/8)
 SV = Point{3, Float32}[ # Sheet Vertices
  [-2, 0, 2],  # 1. upper left corner of unfolded sheet
  [-2, 0, -2], # 2. lower left corner of unfolded sheet
  [2, 0, -2],  # 3. lower right corner of unfolded sheet
  [2, 0, 2],  # 4. upper right corner of unfolded sheet
  [0, 0, 0],  # 5. center of unfolded sheet
  [2, 0, 0],  # 6. center of right edge of unfolded sheet
  [-2, 0, 0],  # 7. center of left edge of unfolded sheet
  [0, 0, 2],  # 8. center of top edge of unfolded sheet
  [0, 0, -2],  # 9. center of bottom edge of unfolded sheet
  [-1/2, 0, 1/2], #10. 
  [-1, 0, 1],  #11. upper left center of folded square
  [-1, 0, -1], #12. lower left center of folded square
  [1, 0, -1],  #13. lower right center of folded square
  [1, 0, 1],  #14. upper right center of folded square
  [2, 0, 1],  #15. right edge quarter length
  [1, 0, 2],  #16. top edge quarter length
  [-2, 0, -1], #17. left edge quarter length
  [-1, 0, -2], #18. bottom edge quarter length
  [2, 0, x],  #19. right edge eighth length
  [x, 0, 2],  #20. top edge eighth length
  [-2, 0, -x], #21. left edge eighth length
  [-x, 0, -2], #22. bottom edge eighth length
  [-2, 0, 1/2], #23. 
  [-1/2, 0, 2], #24. 
  [-5/4, 0, 5/4], #25. 
  [-2, 0, 1],  #26. 
  [-1, 0, 2],  #27.
  [-3/2, 0, 3/2]] #28.
 SF = NgonFace{3, Int64}[] # Sheet Faces (no folds -> no faces)
 S = GeometryBasics.Mesh(SV,SF) # mesh of Sheet
 S_obs = Observable(S)
 
 # precalculate some rotation angles
 nFrameFold = 100 # number of video frames per fold
 nFrameShow = 400 # number of video frames to show kabuto
 THETAFOLD = LinRange(0, pi, nFrameFold)
 THETAFOLD2 = pi .* (1 .-
  abs.(2 .* ((1:nFrameFold) .- 1) ./ nFrameFold .- 1))
 phi0 = -pi/4
 PHISHOW = LinRange(phi0, phi0+2*pi, nFrameShow)
 
 # initialize figure
 strTitle = @sprintf("fold %d of %d", 0, nFold)
 str_obs = Observable(strTitle)
 fig = Figure(resolution = (600, 650))
 ax3d = Axis3(fig[1,1],
  elevation = pi/16,
  azimuth = phi0,
  viewmode = :fit,
  limits = (-2,2, -2,2, -2,2),
  aspect = (1,1,1),
  title = str_obs)
 poly!(ax3d, S_obs,
  color = RGBA(1, 1, 0.9, 0.8),
  strokewidth = 1)
 fig
 
 # generate video of folding
 iFold = 0
 nFrame = nFold * nFrameFold + nFrameShow
 nSweep = 0
 SP = Point{3, Float32}[]# Starting Point(s)
 PL = zeros(Float32, 16) # Pivot Line PGA expression
 MVI = zeros(Int64, 1) # Moving Vertex Indices
 record(fig, "origami.mp4", 1:nFrame) do iFrame
  iFrameMod = mod(iFrame-1, nFrameFold)
  
  # if time for a new fold
  if iFrameMod == 0
   iFold += 1
   nSweep = fold(iFold, SV, SF, PL, MVI)
   SP = point(SV[MVI]) # Starting Point
   str_obs[] = (iFold <= nFold) ?
    @sprintf("fold %d of %d", iFold, nFold) :
    @sprintf("kabuto (Samurai helmet)")
  end
  
  # if folding origami
  if iFold <= nFold
   angle = nSweep == 1 ?
    THETAFOLD[iFrameMod+1] :
    THETAFOLD2[iFrameMod+1]
   R = rotor(angle, PL)
   MP = R >>> SP # Moving Point
   SV[MVI] = toPoint(MP)
   S_obs[] = GeometryBasics.Mesh(SV, SF)
  
  # else showing finished origami
  else
   ax3d.azimuth[] = PHISHOW[iFrame - nFold*nFrameFold]
  end
 end # for each video frame
end # origami()

# quick and dirty ffmpeg conversion from origami.mp4 to origami.gif:
# ffmpeg -i origami.mp4 -r 24 -s 460x480 -loop 0 origami.gif
# ffmpeg -i origami.mp4 -r 24 -s 600x650 -loop 0 origami.gif