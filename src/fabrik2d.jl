# fabrik2d.jl
# interactive graphical demonstration of Inverse Kinematics
# iterative solver algorithm called FABRIK defined by
# Andreas Aristidou and Joan Lasenby in their paper at
# http://www.andreasaristidou.com/publications/papers/FABRIK.pdf
using GLMakie
using Printf
include("ripga2d.jl")

# translate distance along line
function xlator(line::Vector{Float32},dist::Number)
# return ga"1 - dist/2 (e0 normalize(line) e0∗)"
 return 1 - dist/2*(e0*normalize(line)*!e0)
end

# inverse kinematics algorithm; plot of convergence rate
# arguments:
# - coordinates of target
# - number of links in robot arm
function ik(target::Vector{Float64}=[2.5,2.,0.], nLink::Int64=6)
 nEP = nLink + 1 # number of link endpoints
 armLength = 3 # max reach of robot arm
 
 # initialize figure
 fig = Figure(size = (1800, 800))
 C = ["red"; "green"; "blue"] # line colors
 LIM = (-2,3, -2.5,2.5)
 YTIC = (-2:1:2)
 AX = [
  Axis(fig[1,1], limits=LIM, yticks=YTIC, aspect=1,
   title = "1. initial endpoints,\n" *
    "the target is separate from robot arm");
  Axis(fig[1,2], limits=LIM, yticks=YTIC, aspect=1,
   title = "2. endpoints after 1st pass\n" *
    "of backward relaxation");
  Axis(fig[1,3], limits=LIM, yticks=YTIC, aspect=1,
   title = "3. endpoints after 1st pass\n" *
    "of forward relaxation");
  Axis(fig[1,4], limits=LIM, yticks=YTIC, aspect=1,
   title = "4. endpoints after 2nd pass\n" *
    "of backward relaxation");
  Axis(fig[1,5], limits=LIM, yticks=YTIC, aspect=1,
   title = "5. endpoints after 2nd pass\n" *
    "of forward relaxation");
  Axis(fig[2,2], limits=LIM, yticks=YTIC, aspect=1,
   title = "6. endpoints after 3rd pass\n" *
    "of backward relaxation");
  Axis(fig[2,3], limits=LIM, yticks=YTIC, aspect=1,
   title = "7. endpoints after 3rd pass\n" *
    "of forward relaxation");
  Axis(fig[2,4], limits=LIM, yticks=YTIC, aspect=1,
   title = "8. endpoints after 4th pass\n" *
    "of backward relaxation");
  Axis(fig[2,5], limits=LIM, yticks=YTIC, aspect=1,
   title = "9. endpoints after 4th pass\n" *
    "of forward relaxation");
 ]
 
 # allocate endpoint PGA expressions
 # (appended endpoint in last column is target)
 PX = Matrix{Float32}(undef, (length(e0),nEP+1))
 
 # define link endpoints
 linkLength::Float32 = armLength / nLink
 for iEP = 1:nEP
#  PX[:,iEP] = ga"(e0 + (iEP linkLength - 1.5) e1)∗"
  PX[:,iEP] = !(e0 + (iEP*linkLength - 1.5)*e1)
 end
 PX[:,nEP+1] = point(target[1], target[2], target[3])
 
 # plot link endpoints and target point
 P = toCoord(PX)
 iAx = 1
 scatterlines!(AX[iAx], P[1:3,1:nEP], color="black")
 scatterlines!(AX[iAx], # target point is at end
  [P[1,end]], [P[2,end]], [P[3,end]],
  color="black")
 
 # plot results of each relaxation loop
 for iRelax = 1:4
  # set tip to target, changing length of last link
  PX[:,nEP] = PX[:,nEP+1]
  P = toCoord(PX)
  
  # restore link lengths from back to front
  iAx += 1
  scatterlines!(AX[iAx],
   P[1:3,1:nEP], color = "light gray")
  for jLink = 1:nEP-2
   i = nEP - jLink
   iColor = mod(jLink-1,3) + 1
#   XL = xlator(ga"PX[:,i+1] ∨ PX[:,i]", linkLength)
#   PX[:,i] = ga"XL PX[:,i+1] ~XL" # perform translation
   XL = xlator( # define translation along line
    PX[:,i+1] & PX[:,i], linkLength)
   PX[:,i] = XL >>> PX[:,i+1] # perform translation
   P = toCoord(PX)
   if i > 2
    scatterlines!(AX[iAx],
     P[1:3,i:i+1], color = C[iColor])
   else
    scatterlines!(AX[iAx],
     P[1:3,i-1:i+1], color = C[iColor])
   end
  end
  
  # restore link lengths from front to back
  iAx += 1
  scatterlines!(AX[iAx],
   P[1:3,1:nEP], color = "light gray")
  for i = 2:nEP
   iColor = mod(i-2,3) + 1
#   XL = xlator(ga"PX[:,i-1] ∨ PX[:,i]", linkLength)
#   PX[:,i] = ga"XL PX[:,i-1] ~XL" # perform translation
   XL = xlator( # define translation along line
    PX[:,i-1] & PX[:,i], linkLength)
   PX[:,i] = XL >>> PX[:,i-1] # perform translation
   P = toCoord(PX)
   scatterlines!(AX[iAx],
    P[1:3,i-1:i], color = C[iColor])
  end
 end
 fig
end

function ik_solver(PX::Matrix{Float32},
 linkLength::Float32,
 SRC::Matrix{Float32},
 nPass::Int64,
 alpha::Float64,
 pressure::Float64 = 0)
 nEP = size(PX,2) - 1 # -1 because last column is target
 nSRC = size(SRC,2)
 
 # for each relaxation pass
 for iPass = 1:nPass
  # set tip to target, changing length of last link
  PX[:,nEP] = PX[:,nEP+1]
  
  # initially nudge inner pivot points away from point sources
  if alpha > 0.15 && iPass == 1
   NUDGE = zeros(Float32, size(PX,1))
   for i = 2:nEP-1
    for iSRC = 1:nSRC
     L = SRC[:,iSRC] & PX[:,i]
     d = sqrt(L[6]^2 + L[5]^2)
     XL = xlator(L, d+pressure/(1+d^2))
     PN = XL >>> PX[:,i]
     NUDGE += (PN - PX[:,i])
    end
    PX[:,i] += NUDGE
   end # for
  end # if
  
  # restore link lengths from back to front
  for jLink = 1:nEP-2
   i = nEP - jLink
#   XL = xlator(ga"PX[:,i+1] ∨ PX[:,i]", linkLength)
#   PX[:,i] = ga"XL PX[:,i+1] ~XL" # perform translation
   XL = xlator( # define translation along line
    PX[:,i+1] & PX[:,i], linkLength)
   PX[:,i] = XL >>> PX[:,i+1] # perform translation
  end
  
  # restore link lengths from front to back
  for i = 2:nEP
#   XL = xlator(ga"PX[:,i-1] ∨ PX[:,i]", linkLength)
#   PX[:,i] = ga"XL PX[:,i-1] ~XL" # perform translation
   XL = xlator( # define translation along line
    PX[:,i-1] & PX[:,i], linkLength)
   PX[:,i] = XL >>> PX[:,i-1] # perform translation
  end
 end
end

# interactive inverse kinematics
# arguments:
# - coordinates of target of robot arm
# - number of links in robot arm
function iik(target::Vector{Float64}=[2.5,2.,0.], nLink::Int64=6)
 nEP = nLink + 1 # number of link endpoints
 armLength = 3 # max reach of robot arm

 # initialize figure
 fig = Figure(size = (800, 800))
 LIM = (-2,3, -2.5,2.5)
 YTIC = (-2:1:2)
 ax1 = Axis(fig[1,1], limits=LIM, yticks=YTIC, aspect=1,
  title = "Interactive demonstration of inverse kinematics algorithm.\n" *
   "(The target point is initially separated from the robot arm.\n" *
   "Drag that target point around to see how the robot arm reacts.)")

 # allocate and define endpoint PGA expressions
 # (appended endpoint in last column is target)
 PX = Matrix{Float32}(undef, (length(e0),nEP+1))
 linkLength::Float32 = armLength / nLink
 ANCHOR = [-1; 0; 0]
 for iEP = 1:nEP
#  PX[:,iEP] = ga"(e0 + ((iEP-1) linkLength + ANCHOR[1]) e1)∗"
  PX[:,iEP] = !(e0 + ((iEP-1)*linkLength + ANCHOR[1])*e1)
 end
 PX[:,nEP+1] = point(target[1], target[2], target[3])

 # calculate inverse kinematics
 ik_solver(PX, linkLength)
 P = toCoord(PX) # convert PGA expressions to Euclidean coordinates

 # define observables for plotting
 RCOORD = Observable(P[1:3,1:nEP]) # robot coordinates
 TCOORD = Observable([ANCHOR P[1:3,end]]) # target coordinates

 # plot robot and target coordinates
 scatterlines!(ax1, RCOORD, color="black")
 scatter!(ax1, TCOORD, color="red")

 deregister_interaction!(ax1, :rectanglezoom)
 register_interaction!(ax1, :my_mouse_interaction) do event::MouseEvent, axis
  if Makie.is_mouseinside(ax1.scene)
   if event.type === MouseEventTypes.leftdrag
    PX[:,end] = point(event.data[1], event.data[2], 0)
    ik_solver(PX, linkLength)
    P = toCoord(PX)
    RCOORD[] = P[1:3,1:end-1] # update the plotted observables to 
    TCOORD[] = [ANCHOR P[1:3,end]] # automatically update the plot
   end
  end
 end
 fig
end

# collides: detects collisions of non-adjacent links
# arguments:
# 1. PX: robot arm pivot point expressions
# 2. LX: line segment expressions
function collides(PX::Matrix{Float32},LX::Matrix{Float32})::Bool
 nLink = size(LX,2)
 for iLink = 1:nLink-2
  for jLink = iLink+2:nLink
   # calculate areas that have negative product
   # when the two points are on opposite side of line
   A0 = (LX[:,iLink] & PX[:,jLink])[1]
   A1 = (LX[:,iLink] & PX[:,jLink+1])[1]
   if A0 * A1 < 0
    # calculate areas that have negative product
    # when the two points are on opposite side of line
    A2 = (LX[:,jLink] & PX[:,iLink])[1]
    A3 = (LX[:,jLink] & PX[:,iLink+1])[1]
    if A2 * A3 < 0
     return true
    end
   end
  end
 end
 return false
end

# inverse kinematics video
# arguments:
# - number of relaxation passes in ik solver
# - number of links in robot arm
function ikv(nPass::Int64=8, nLink::Int64=12)
 nEP = nLink + 1 # number of link endpoints
 armLength = 3 # max reach of robot arm
 
 # initialize figure
 fig = Figure(size = (800, 800))
 LIM = (-3,3, -3,3)
 ax2d = GLMakie.Axis(fig[1,1],
  limits=LIM,
  aspect=1,
  titlesize = 22,
  title = "inverse kinematics using " *
   "fast iterative algorithm called FABRIK\n" *
   "(see http://www.andreasaristidou.com/publications/papers/FABRIK.pdf)")
 
 # plot elliptical path of target
 nFrame = 360 # sample points from ellipse
 a = 2.5   # cos coefficient (1/2 major axis)
 b = 1   # sin coefficient (1/2 minor axis)
 x0 = 1   # ellipse x offset
 y0 = 0   # ellipse y offset
 phi = pi/4   # tilt of ellipse (major axis tilt)
 PRESSURE = [ # collision avoidance pressure
  0.0; 0.001]
 THETA = LinRange(0, 2*pi,
  nFrame+1) # +1 to avoid adjacent duplicate
     # frames in repeating animated gif
 T = zeros(Float32,2,nFrame+1) # Target (destination)
 T[1,:] = (a*cos(phi)) .* cos.(THETA) -
    (b*sin(phi)) .* sin.(THETA) .+ x0
 T[2,:] = (a*sin(phi)) .* cos.(THETA) +
    (b*cos(phi)) .* sin.(THETA) .+ y0
 lines!(ax2d, T,
  linestyle = :dot,
  color = :gray)
 
 # plot initial robot links
 ANCHOR = [-1.5; 0]
# SRC = [point(x0, y0) point(ANCHOR[1], 0)]
 SRC = Matrix{Float32}(undef, (length(e0),1))
 SRC[:,1] = point(x0, y0)
 println("SRC: $SRC")
 UX = point(-1,0) & point(0,0)
 UY = point(-1,0) & point(-1,1)
 TCOORD = Observable([ANCHOR T[:,1]])
 P = zeros(Float32, 2, nEP+1) # robot Pivots + 1 for target
 P[1,1:nEP] = LinRange(ANCHOR[1],ANCHOR[1]+armLength,nEP)
 PCOORD = Observable(P[:,1:nEP])
 scatterlines!(ax2d, PCOORD, color="black")
 scatter!(ax2d, TCOORD, color="red")
 PRESSURELABEL = Observable("")
 text!(-2.75, 2.8,
  text = PRESSURELABEL,
  align = (:left, :center),
  fontsize = 25)
 
 # pre-allocate for animation calculations
 PX = Matrix{Float32}(undef, (length(e0),nEP+1))
 linkLength::Float32 = armLength / nLink
 for iEP = 1:nEP
#  PX[:,iEP] = ga"(e0+((iEP-1) linkLength+ANCHOR[1]) e1)∗"
  PX[:,iEP] = !(e0 + ((iEP-1)*linkLength + ANCHOR[1])*e1)
 end
 PX[:,end] = point(T[1,1], T[2,1])
 fig
 LX = Matrix{Float32}(undef, (length(e0),nLink))
 nRev = 2
 LINKANGLE = zeros(nFrame, nEP, nRev)
 
 # record video of robot following target path
 iRev = 0
 fn = @sprintf("fabrik2dv%02d.mp4", nPass)
 record(fig, fn, 1:2*nFrame) do jFrame
  # check for beginning of ellipse (traversed twice)
  iFrame = mod(jFrame-1, nFrame) + 1
  alpha = iFrame / nFrame
  if iFrame == 1
   iRev += 1
   PRESSURELABEL[] = @sprintf(
    "collision avoidance nudge \"pressure\" = %.4f",
    PRESSURE[iRev])
  end
  
  # move target T along ellipse
  PX[:,end] = point(T[1,iFrame], T[2,iFrame])
  ik_solver(PX, linkLength, SRC, nPass,
   alpha, PRESSURE[iRev])
  P = toCoord(PX)
  PCOORD[] = P[:,1:end-1]
  TCOORD[] = [ANCHOR P[:,end]] # automatically update the plot
  
  # precalculate link lines for link angle calculations
  for iLink = 1:nLink
   LX[:,iLink] = PX[:,iLink] & PX[:,iLink+1]
  end
  
  # calculate link angles (for adjacent link collisions)
  for iLink = 1:nLink
   LINKANGLE[iFrame,iLink+1,iRev] = # +1 for column 1 of zeros
    180/pi * atan(
     (UY|LX[:,iLink])[1],
     (UX|LX[:,iLink])[1])
  end
 end
 
 # for each animated revolution
 for iRev = 1:nRev
  # unwrap wrapped link angles
  isAdjacentLinkCollision = false
  LAD = diff(LINKANGLE[:,:,iRev], dims=2) # Link Angle Diff
  DELTA = diff(LAD, dims=1)
  for iFrame = 1:nFrame-1
   for iLink = 1:nLink
    if DELTA[iFrame,iLink] < -0.9*360
     LAD[iFrame+1:end,iLink] .+= 360
     isAdjacentLinkCollision = true
    elseif DELTA[iFrame,iLink] > 0.9*360
     LAD[iFrame+1:end,iLink] .-= 360
     isAdjacentLinkCollision = true
    end
   end
  end
  
  # plot link angles
  fig2 = Figure(size = (800, 800))
  strTitle = @sprintf("""link angles (unit: degrees)
   top to bottom plots correspond to anchor to end links
   %d relaxation passes
   collision avoidance nudge \"pressure\" = %0.4f""", 
   nPass, PRESSURE[iRev])
  for iLink = 1:nLink
   if iLink == 1
    ax = GLMakie.Axis(fig2[iLink,1],
     xticks = (0:90:360),
     titlesize = 22,
     title = strTitle)
    hidexdecorations!(ax,
     grid = false,
     ticks = false)
   elseif iLink < nLink
    ax = GLMakie.Axis(fig2[iLink,1],
     xticks = (0:90:360))
    hidexdecorations!(ax,
     grid = false,
     ticks = false)
   else
    strXlabel = "position along elliptical path " *
     "(unit: degrees)"
    ax = GLMakie.Axis(fig2[iLink,1],
     xticks = (0:90:360),
     xlabel = strXlabel)
   end
   lines!(ax, LAD[:,iLink], color=:black)
  end
  fn = @sprintf("fabrik2dang%02d_%d.png", 
   nPass, iRev)
  save(fn, fig2)
  fig2
 end
end

# quick and dirty ffmpeg conversion from .mp4 to .gif:
# /tool/ffmpeg-2022-07/bin/ffmpeg.exe -i fabrik2dv.mp4 -r 24 -s 480x480 -loop 0 fabrik2dv.gif