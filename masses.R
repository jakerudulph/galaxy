# A SIMPLE GALAXY - Produced 2011 and converted to R 1/26/2024

# SET UP INITIAL CONDITIONS AND CONSTANTS
rm(list=ls()) # Clear environment in R

n <- as.numeric(readline(prompt = "Enter the number of stars in your galaxy. A wise number is between 50 and 500. "))
sm <- 1.9891E30 # KILOGRAMS 
G <- 6.672E-11 # Nm^2/kg^2 (GRAVITATIONAL CONSTANT)
sol <- 299792458 # SPEED OF LIGHT, IN M/S 
secondsyear <- 365*24*60*60
ly <- sol * secondsyear # ONE LIGHTYEAR
parsec <- 3.26 * ly
gr <- 15000 # GALAXY RADIUS IN parsecs
bhm <- 1.73E11 # BLACK HOLE MASS
bhm <- bhm * sm # SET IN UNITS OF SOLAR MASS
origin <- 0 # GALACTIC CENTER ORIGIN FOR BOTH X AND Y
tstep <- 2000000 # OUR STEP TIME IN MILLIONS OF YEARS
tinc <- tstep * secondsyear # THIS TIME IN MILLIONS OF YEARS*SECONDS/YR (SPEED OF LIGHT IN M/S)
by <- 5 # NUMBER OF BILLIONS OF YEARS TO RUN THE MODEL
ntime <- by * 10^9 / tstep # NUMBER OF TIME STEPS IN OUR BILLIONS OF YEARS

rO <- gr * parsec # THE OUTER BOUNDARY FOR GENERATING STARS IN PARSECS
rIn <- (gr / 6) * parsec # THRESHOLD

# LET'S SET UP THE GALAXY CENTER
masses <- matrix(0, n, 6)
masses[1, 1] <- bhm
masses[1, 2:3] <- origin * parsec
masses[1, 4:6] <- 0

# FOR ALL MASSES WHICH ARE NOT THE BLACK HOLE (FIRST ROW)
for (i in 2:n) {
  mi <- runif(1)
  m <- (mi + 0.8)^3 * sm # NORMALIZE SO SOME STARS LARGER, SOME SMALLER THAN OUR SUN
  masses[i, 1] <- m
  
  ri <- runif(1)
  r <- ri * (rO - rIn) + rIn
  
  ang0 <- runif(1)
  ang <- ang0 * 2 * pi
  x <- r * cos(ang)
  y <- r * sin(ang)
  dx <- cos(ang + pi / 2)
  dy <- sin(ang + pi / 2)
  
  v <- sqrt(G * bhm / r)
  masses[i, 2] <- x + origin * parsec
  masses[i, 3] <- y + origin * parsec
  masses[i, 4] <- dx * v
  masses[i, 5] <- dy * v
  masses[i, 6] <- r
}

# PLOTTING
par(mfrow=c(2,2))
plot(masses[,2], masses[,3], col="blue", pch=20, main="Initial Positions, Galactic Center in Red", xlab="X Plane, Meters", ylab="Y plane, Meters", asp=1)
points(masses[1,2], masses[1,3], col="red", pch=20)

hist(masses[-1,1] / sm, main="Histogram of Star Masses", xlab="Stellar Mass Equivalent", ylab="Frequency")

# SIMULATION OVER TIME
for (j in 1:tinc) {
  for (k in 1:n) {
    xd <- masses[k, 2]
    yd <- masses[k, 3]
    acc <- c(0, 0)
    dr2 <- xd^2 + yd^2
    dstars <- matrix(1, n, 2)

	massk <- matrix(0, n, 7)
	massk[,1:3] <- masses[,1:3]

    for (q in 2:n) {
      massk[q,4] <- massk[q,2] - xd
      massk[q,5] <- massk[q,3] - yd
      massk[q,6] <- (massk[q,4] / sqrt(dr2)) * (G * (massk[q,1])) / dr2
      massk[q,7] <- (massk[q,5] / sqrt(dr2)) * (G * (massk[q,1])) / dr2
    }
    starG <- c(sum(massk[,6]), sum(massk[,7]))
    dcenter <- c(origin, origin) - c(xd, yd)
    acc <- acc + starG + (dcenter / sqrt(dr2)) * ((G * (bhm + masses[k,1])) / dr2)
    masses[k,4:5] <- masses[k,4:5] + acc * tinc
    masses[k,2:3] <- masses[k,2:3] + masses[k,4:5] * tinc
    masses[k,6] <- sqrt(dcenter[1]^2 + dcenter[2]^2)
  }
  
	# After the simulation loop
	par(mfrow=c(2,2)) # Prepare a 2x2 plotting area
	
	# Plot initial positions with the galactic center in red
	plot(masses[,2], masses[,3], col="blue", pch=20, main="Initial Positions, Galactic Center in Red", xlab="X Plane, Meters", ylab="Y Plane, Meters", asp=1)
	points(masses[1,2], masses[1,3], col="red", pch=20)
	
	# Histogram of star masses
	hist(masses[-1,1] / sm, main="Histogram of Star Masses", xlab="Stellar Mass Equivalent", ylab="Frequency", breaks=20)
	
	# Simulate and plot positions over time inside the loop
	# (The real-time plotting inside the loop would require an interactive R session. For long simulations, consider plotting only key frames.)
	
	# Galaxy Rotation (Velocity vs. Distance from Center) Plot
	velocity <- sqrt(masses[,4]^2 + masses[,5]^2)
	plot(masses[,6]/ly, velocity, col="blue", pch=20, main="Galaxy Rotation Curve", xlab="Distance from Center, Lightyears", ylab="Velocity, m/s", xlim=c(0.5e4, 6e4), ylim=c(0, 6e5))
	
	# Galactic rotation visualization at the final state
	plot(masses[,2], masses[,3], col="white", bg="black", pch=20, main="Galactic Rotation", xlab="X Plane, Meters", ylab="Y Plane, Meters", asp=1)
	points(masses[1,2], masses[1,3], col="red", pch=20) # Mark the galactic center
	
	# Set the plotting area limits based on the galaxy radius and ensure aspect ratio is 1
	xlim <- c(-2*gr*parsec, 2*gr*parsec)
	ylim <- c(-2*gr*parsec, 2*gr*parsec)
	  
	  Sys.sleep(0.05) # Pause for animation effect, equivalent to MATLAB's pause(0.05)
}

# The R code for plotting the galaxy's rotation and the rotation curve will be similar to the MATLAB plotting code,
# using plot(), points(), and hist() functions for visualization.
