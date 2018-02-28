#Load libraries
library(circular) #For correlated random walk
library(sp) #For spatial plotting (and nest coordinates)
library(rgdal) #For geospatial data
library(RColorBrewer) #For plotting
library(rgeos) #For placing buffer around nests

#Test:
# attack radius of 10, 50, 100m
# movement speed: 10, 50, 100m/min
# turn angle: 0.1, 0.5, 0.9

#Predator attributes
predator_turn_angle = 0.95 # 0.001 = diffusion, 0.999 = straight lines
predator_speed = 33 #meters per min (2km/h)
attack_radius = 50 #meters


#LPB nest locations (use 2013)
nestloc = read.csv("C:/Users/koonslt/Documents/Iles/University/PhD_Project/Data/COEI/Database/spatial_data/gps_data/coei_2013_gps_REVISED.csv")

#Select eider nests
nestloc = nestloc[1:404,]

#Convert to UTM
coordinates(nestloc) <- c("lon", "lat")
cord.dec = SpatialPoints(nestloc, proj4string=CRS("+proj=longlat"))
res <- spTransform(cord.dec, CRS("+proj=utm +zone=15 ellps=WGS84"))

#Initial nest dataframe
nest_init = as.data.frame(res)
nest_init$ID = paste("Nest",1:nrow(nest_init))
nest_init$min_failed = NA

#Nest phenology
nest_init$init = as.integer(rnorm(nrow(nest_init),5,3)) #Dispersed phenology
#nest_init$init = 1 #Perfectly synchronous phenology

nest_init$hatch = nest_init$init + 30 #30 day nesting period

days = 60
min_season = 60*24*days #convert days to minutes (length of simulation)

camera_proportion = 0.1 #proportion of nests with cameras on them
cameras = round(nrow(nest_init) * camera_proportion)
cam_nests = sample(1:nrow(nest_init),cameras, replace=FALSE) #Randomly select nests to place cameras on
nest_init$cam = 0
nest_init$cam[cam_nests] = 1

########################################################
# Place reflexive boundary around colony
########################################################

boundary_x = c(min(res$lon)-20000, max(res$lon)+20000)
boundary_y = c(min(res$lat)-50000, max(res$lat)+2000)

#Initially place predator on coast, about 5 km away
x = boundary_x[2]
y = runif(1,min(boundary_y),max(boundary_y))
y_dir = abs((y - min(boundary_y))/(max(boundary_y)-min(boundary_y)) + 0.5)
Phi = y_dir*pi

#Random initial location and direction of predator (choose either this or conditions above)
x = runif(1,boundary_x[1],boundary_x[2])
y = runif(1,boundary_y[1],boundary_y[2])
y_dir = runif(1,0,2*pi)

x_vec = x
y_vec = y

#Dynamic dataframe that will update at each timestep
nest_dyn = nest_init
nest_dyn = subset(nest_dyn, init <= 1) #Only select nests that have initiated so far

day = 1

# Start the clock to record processing time for this simulation
ptm <- proc.time()

#Loop through minutes of study
for (i in 2:min_season){

    day_update = round(i/60/24) #convert i to days (i is recorded in minutes)

    #Check if it is a new day
    #If so, add new nests, remove hatched ones, and update the day counter
    if (day_update > day){

        day = day+1
        #Add nests that have initiated today
        if (length(which(nest_dyn$init == day)>0)) nest_dyn = rbind(nest_dyn, subset(nest_init, init == day))

        #Remove nests that have hatched today
        if (length(which(nest_dyn$hatch == day)>0))nest_dyn = nest_dyn[-which(nest_dyn$hatch == day),]
    }

    #If all nests have hatched or failed, stop the simulation
    if (day > max(subset(nest_init, is.na(min_failed))$hatch)) break

    #If nest_dyn is empty, create large dist (to keep predator moving)
    if (nrow(nest_dyn) < 1 | is.null(nrow(nest_dyn))) dist = 10000 else dists = sqrt( (nest_dyn$lon - x)^2 + (nest_dyn$lat - y)^2)

    #If predator is within attack radius of any nests
    if (sum(dists<attack_radius) > 0){

        #"Teleport" to closest nest within attack radius
        close_nest = which(dists == min(dists))

        nest_id = nest_dyn$ID[close_nest]

        nest_init$min_failed[nest_init$ID == nest_id] = i

        x = nest_dyn$lon[close_nest]
        y = nest_dyn$lat[close_nest]
        nest_dyn = nest_dyn[-close_nest,]

        x_vec = c(x_vec,x)
        y_vec = c(y_vec,y)

        #After eating nest, head in new random direction
        Phi <- rwrappedcauchy(1,circular(0),0.001)

        print(i)

        #If otherwise, continue normal movement
    } else{

        # make weibull distributed steps
        steps <- 150 #rweibull(1,5,160)

        # check out their distribution
        #hist(steps)

        # make clustered turning angles
        theta <- rwrappedcauchy(1,circular(0),predator_turn_angle)

        # cumulative angle (absolute orientation)
        Phi <- Phi + theta

        # step length components
        dX <- steps*cos(Phi)
        dY <- steps*sin(Phi)

        # actual X-Y values
        x = x+dX
        y = y+dY

        # Make sure predator stays inside boundary box (home range)
        while (x < min(boundary_x) | x > max(boundary_x) | y < min(boundary_y) | y > max(boundary_y)){
            x = x-dX
            y = y-dY

            # If outside boundary, make random turn angles until re-entering bounding box
            theta <- rwrappedcauchy(1,circular(0),0.001)

            # cumulative angle (absolute orientation)
            Phi <- Phi + theta

            # step length components
            dX <- steps*cos(Phi)
            dY <- steps*sin(Phi)

            # actual X-Y values
            x = x+dX
            y = y+dY

        }

        #Store predator steps
        x_vec = c(x_vec,x)
        y_vec = c(y_vec,y)

        print(i)
    }

}
# Stop the clock
sim.time = proc.time() - ptm
print(sim.time) #Print processing time

#Determine failed nests
nest_init$dep = 0
nest_init$dep[nest_init$min_failed > 0] = 1

#Buffer around nests - to delineate study area
buf1 <- gBuffer(res, width=attack_radius, byid=TRUE)
buf2 <- gUnaryUnion(buf1)

#Predator in colony
pred_pts = data.frame(lat = y_vec, lon = x_vec, time = 1:length(x_vec))
coordinates(pred_pts) <- c("lon", "lat")
pred_pts = SpatialPoints(pred_pts, proj4string=CRS("+proj=utm +zone=15 ellps=WGS84"))

occupied = over( pred_pts , buf2 , fn = NULL)
times_occupied = which(!is.na(occupied))
days_occupied = floor(times_occupied/60/24 + 1)

#Predation dynamics
nest_init$day_failed = floor(nest_init$min_failed/60/24 + 1)

nest_init$camera_captures = nest_init$cam * nest_init$day_failed

if (sum(nest_init$camera_captures == 0, na.rm=TRUE)>0) nest_init$camera_captures[nest_init$camera_captures == 0] = NA

par(mfrow=c(1,1))
par(mar=c(4,4,2,1))

predator_walk_col = colorRampPalette(colors=c("gray25","gray95"))(length(x_vec))

#Predator walk
par(fig = c(0,0.5, 0.5, 1))
plot(y_vec~x_vec, xlim = boundary_x, ylim = boundary_y, type="l",
     xlab = "UTM x", ylab = "UTM y", col = "gray85", main = "Track of Predator", cex = 0.01)
points(lat~lon, col = "blue", pch = 19, cex = 0.2, data = nest_init)
points(lat~lon, col = "red", pch = 19, data = na.omit(nest_init[,c("lon","lat","day_failed")]), cex = 0.2)
#points(y_vec[1]~x_vec[1], col = "black", pch = 4, cex = 0.5)

text(x = boundary_x[1]+diff(boundary_x),y = boundary_y[1]+diff(boundary_y)*0.9, labels = paste("speed=",predator_speed,"m/min",sep=""), adj = 1)
text(x = boundary_x[1]+diff(boundary_x),y = boundary_y[1]+diff(boundary_y)*0.8, labels = paste("angle=",predator_turn_angle,sep=""), adj = 1)
text(x = boundary_x[1]+diff(boundary_x),y = boundary_y[1]+diff(boundary_y)*0.7, labels = paste("radius=",attack_radius,"m",sep=""), adj = 1)

par(mar=c(2,2,2,1))
par(fig = c(0.05,0.3, 0.7, 1), new = T)
rose.diag(rwrappedcauchy(10000,circular(0),predator_turn_angle),bins=24, col = "white", cex = 0.8)

par(fig = c(0.5,1, 0.5, 1), new = T)
par(mar=c(5,5,1,1))
#hist(days_occupied, breaks = seq(0,days), xlab = "Day of Season", ylab = "Steps inside Colony by Predator", main = "Presence in Colony Per Day", col = "gray75")
hist(nest_init$day_failed, breaks = seq(0,days), xlab="Day of Season", ylab = "Number of Nests Eaten", main = "Nests Eaten per Day", col = "red")

par(fig = c(0,0.5, 0, 0.5), new = T)
par(mar=c(4,4,2,2))
plot(buf2,xlab = "", ylab = "", main = "Close-up of Colony", col = "white")
#plot(lat~lon, data = nest_init, pch = 19, col = "blue",xlab = "UTM x", ylab = "UTM y", main = "Close-up of Colony")

points(lat~lon, col = "blue", pch = 19, data = nest_init)
points(lat~lon, col = "red", pch = 19, data = na.omit(nest_init[,c("lon","lat","day_failed")]))
points(lat~lon, col = "black", pch = 0, cex = 1.5, data = subset(nest_init, cam == 1))
lines(y_vec~x_vec, lty = 2, col = "gray75")
points(y_vec~x_vec, pch=19, col = predator_walk_col, cex = 0.5)
points(y_vec[which(y_vec %in% nest_init$lat)]~x_vec[which(x_vec %in% nest_init$lon)], pch=19, col = "red", cex = 0.8)

par(fig = c(0.5,1, 0, 0.5), new = T)
par(mar=c(5,5,2,1))
#plota = hist(nest_init$day_failed, breaks = seq(0,days), xlab="Day of Season", ylab = "Number of Nests Eaten", main = "Nests Eaten per Day", col = "gray75", plot = FALSE)
hist(nest_init$camera_captures, breaks = seq(0,days), xlab="Day of Season", ylab = "Number of Times Caught on Camera", main = "Camera Captures", col = "black")

par(mfrow=c(1,1))
