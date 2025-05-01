## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install------------------------------------------------------------------
# install.packages("ecoregime")
# devtools::install_github(repo = "MSPinillos/ecoregime", dependencies = T, build_vignettes = T)

## ----setup--------------------------------------------------------------------
library(ecoregime)

## ----citation-----------------------------------------------------------------
citation("ecoregime")

## ----data---------------------------------------------------------------------
# ID of the sampling units and observations
ID_sampling <- LETTERS[1:10]
ID_obs <- 1:3

# Define initial species abundances
set.seed(123)
initial <- data.frame(sampling_units = ID_sampling, 
                   sp1 = round(runif(10), 2), 
                   sp2 = round(runif(10), 2))

# Define the parameters for the Lotka-Volterra model
parms <- c(r1 = 1, r2 = 0.1, a11 = 0.02, a21 = 0.03, a22 = 0.02, a12 = 0.01)

# We can use primer and deSolve to run the simulations
library(primer)
simulated_abun <- lapply(1:nrow(initial), function(isampling){
  # Initial abundance in the sampling unit i
  initialN <- c(initial$sp1[isampling], initial$sp2[isampling])
  # Simulate community dynamics
  simulated_abun <- data.frame(ode(y = initialN, times = ID_obs, func = lvcomp2, parms = parms))
  # Add the names of the sampling units
  simulated_abun$sampling_unit <- ID_sampling[isampling]
  # Calculate relative abundances
  simulated_abun$sp1 <- simulated_abun$X1/rowSums(simulated_abun[, 2:3])
  simulated_abun$sp2 <- simulated_abun$X2/rowSums(simulated_abun[, 2:3])
  
  return(simulated_abun[, c("sampling_unit", "time", "sp1", "sp2")])
})

# Compile species abundances of all sampling units in the same data frame
abundance <- do.call(rbind, simulated_abun)


## ----abundance----------------------------------------------------------------
head(abundance)

## ----state_dissim-------------------------------------------------------------
# Generate a matrix containing dissimilarities between states
state_dissim <- vegan::vegdist(abundance[, c("sp1", "sp2")], method = "bray")
as.matrix(state_dissim)[1:6, 1:6]


## ----state_mds, fig = TRUE, fig.height = 5, fig.width = 5, fig.align = "center"----
# Multidimensional scaling
state_mds <- smacof::smacofSym(state_dissim, ndim = 2)
state_mds <- data.frame(state_mds$conf)

# Define different colors for each trajectory
traj.colors <- RColorBrewer::brewer.pal(10, "Paired")

# Plot the distribution of the states in the state space
plot(state_mds$D1, state_mds$D2,
     col = rep(traj.colors, each = 3), # use different colors for each sampling unit
     pch = rep(ID_sampling, each = 3), # use different symbols for each sampling unit
     xlab = "Axis 1", ylab = "Axis 2",
     main = "State space")


## ----traj_state_mds, fig = TRUE, fig.height = 5, fig.width = 5, fig.align = "center"----
## -->> This code shows you the process step by step. You could directly use 
## -->> ecotraj::trajectoryPlot()

# Plot the distribution of the states in the state space
plot(state_mds$D1, state_mds$D2,
     col = rep(traj.colors, each = 3), # use different colors for each sampling unit
     pch = rep(ID_sampling, each = 3), # use different symbols for each sampling unit
     xlab = "Axis 1", ylab = "Axis 2",
     main = "State space")

# Link trajectory states in chronological order
for (isampling in seq(1, 30, 3)) {
  # From observation 1 to observation 2
  shape::Arrows(state_mds[isampling, 1], state_mds[isampling, 2],
         state_mds[isampling + 1, 1], state_mds[isampling + 1, 2],
         col = traj.colors[1+(isampling-1)/3], arr.adj = 1)
  # From observation 2 to observation 3
  shape::Arrows(state_mds[isampling + 1, 1], state_mds[isampling + 1, 2],
         state_mds[isampling + 2, 1], state_mds[isampling + 2, 2],
         col = traj.colors[1+(isampling-1)/3], arr.adj = 1)
}



## ----traj_dissim, fig = TRUE, fig.height = 5, fig.width = 5, fig.align = "center"----
# Generate a matrix containing dissimilarities between trajectories
traj_dissim <- ecotraj::trajectoryDistances(
  ecotraj::defineTrajectories(state_dissim, 
                              sites = rep(ID_sampling, each = 3),
                              surveys = rep(ID_obs, 10)),
  distance.type = "DSPD"
)

as.matrix(traj_dissim)[1:6, 1:6]


## ----traj_mds, fig = TRUE, fig.height = 5, fig.width = 5, fig.align = "center"----
# Multidimensional scaling
traj_mds <- smacof::smacofSym(traj_dissim, ndim = 2)
traj_mds <- data.frame(traj_mds$conf)

# Plot the distribution of the trajectories in the trajectory space
plot(traj_mds$D1, traj_mds$D2,
     col = traj.colors, # use different colors for each sampling unit
     pch = ID_sampling, # use different symbols for each sampling unit
     xlab = "Axis 1", ylab = "Axis 2",
     main = "Trajectory space")



## ----sp_variation, fig = TRUE, fig.height = 5, fig.width = 5, fig.align = "center"----
# Plot species abundances over time
plot(abundance$time, abundance$sp1, type = "n", ylim = c(-0.01, 1),
     xlab = "Observation", ylab = "Species abundance",
     main = "Temporal variation in species composition")

for (i in seq_along(ID_sampling)) {
  points(abundance[which(abundance$sampling_unit == ID_sampling[i]), ]$time,
         abundance[which(abundance$sampling_unit == ID_sampling[i]), ]$sp1,
         col = traj.colors[i], pch = ID_sampling[i])
  lines(abundance[which(abundance$sampling_unit == ID_sampling[i]), ]$time,
         abundance[which(abundance$sampling_unit == ID_sampling[i]), ]$sp1,
         col = "black", pch = ID_sampling[i])
  
  points(abundance[which(abundance$sampling_unit == ID_sampling[i]), ]$time,
         abundance[which(abundance$sampling_unit == ID_sampling[i]), ]$sp2,
         col = traj.colors[i], pch = ID_sampling[i])
  lines(abundance[which(abundance$sampling_unit == ID_sampling[i]), ]$time,
         abundance[which(abundance$sampling_unit == ID_sampling[i]), ]$sp2,
         col = "red", pch = ID_sampling[i], lty = 2)
}

legend("bottomleft", legend = c("sp1", "sp2"), col = 1:2, lty = 1:2, bg = "white")



## -----------------------------------------------------------------------------
# Correlation of the first axis of the state space with the abundance of sp1
cor(state_mds$D1, abundance$sp1)

# Correlation of the first axis of the state space with the abundance of sp2
cor(state_mds$D1, abundance$sp2)


## ----RETRA-EDR----------------------------------------------------------------
# Use set.seed to obtain reproducible results of the segment space in RETRA-EDR
set.seed(123)

# Apply RETRA-EDR
repr_traj <- retra_edr(d = state_dissim, trajectories = rep(ID_sampling, each = 3),
                   states = rep(ID_obs, 10), minSegs = 2)



## ----retra_summary------------------------------------------------------------
summary(repr_traj)

## ----retra_segments-----------------------------------------------------------
lapply(repr_traj, "[[", "Segments")

## ----plot_retra, fig = TRUE, fig.height = 5, fig.width = 5, fig.align = "center"----
# Plot the representative trajectories of an EDR
plot(repr_traj, # <-- This is a RETRA object returned by retra_edr()
     # data to generate the state space
     d = state_mds, trajectories = rep(ID_sampling, each = 3), states = rep(ID_obs, 10),
     # use the colors previously used for individual trajectories. 
     # (make them more transparent to highlight representative trajectories)
     traj.colors = alpha(traj.colors, 0.3),
     # display representative trajectories in blue    
     RT.colors = "blue",
     # select T2 to be displayed with a different color (black)
     select_RT = "T2", sel.color = "black",
     # Identify artificial links using dashed lines (default) and a different color (red)
     link.lty = 2, link.color = "red",
     # We can use other arguments in plot()
     main = "Representative trajectories")

# Add a legend
legend("bottomright", c("T1", "T2", "Link"), 
       col = c("blue", "black", "red"), lty = c(1, 1, 2))


## ----define_retra_df----------------------------------------------------------
# Generate a data frame indicating the states forming the new trajectories
new_traj_df <- data.frame(
  # name of the new trajectories (as many times as the number of states)
  RT = c(rep("T1.1", 4), rep("T2.1", 5), rep("T1.2", 2)),
  # name of the trajectories (sampling units)
  RT_traj = c(rep("B", 2), rep("I", 2), # for the first trajectory (T1.1)
              rep("C", 3), rep("I", 2), # for the second trajectory (T2.1)
              rep("E", 2)), # for the third trajectory (T1.2)
  # states in each sampling unit 
  RT_states = c(1:2, 2:3, # for the first trajectory (T1.1)
                1:3, 2:3, # for the second trajectory (T2.1)
                2:3), # for the third trajectory (T1.2)
  # representative trajectories obtained in retra_edr()
  RT_retra = c(rep("T1", 4), rep("T2", 5), 
               rep("T1", 2)) # The last segment belong to both (T1, T2), choose one
)

new_traj_df


## ----define_retra_ls----------------------------------------------------------
# List including the sequence of segments for each new trajectory
new_traj_ls <- list(
  # First part of T1 excluding the last segment
  repr_traj$T1$Segments[1:(length(repr_traj$T1$Segments) - 1)],
  # First part of T2 excluding the last segment
  repr_traj$T2$Segments[1:(length(repr_traj$T2$Segments) - 1)],
  # Last segment of T1 and T2: segment composed of states 2 and 3 of the sampling unit E ("E[2-3]")
  "E[2-3]"
)

new_traj_ls


## ----define_retra-------------------------------------------------------------
# Define representative trajectories using a data frame
new_repr_traj <- define_retra(data = new_traj_df, 
             # Information of the state space
             d = state_dissim, trajectories = rep(ID_sampling, each = 3), 
             states = rep(ID_obs, 10), 
             # RETRA object returned by retra_edr()
             retra = repr_traj)

# Define representative trajectories using a list with sequences of segments
new_repr_traj_ls <- define_retra(data = new_traj_ls, 
             # Information of the state space
             d = state_dissim, trajectories = rep(ID_sampling, each = 3), 
             states = rep(ID_obs, 10), 
             # RETRA object returned by retra_edr()
             retra = repr_traj)

if (all.equal(new_repr_traj, new_repr_traj_ls)) {
  print("Yes, both are equal!")
}



## ----plot_newretra, fig = TRUE, fig.height = 5, fig.width = 5, fig.align = "center"----
plot(new_repr_traj, # <-- This is the RETRA object returned by define_retra()
     # data to generate the state space
     d = state_mds, trajectories = rep(ID_sampling, each = 3), states = rep(ID_obs, 10),
     # display individual trajectories in light blue
     traj.colors = "lightblue",
     # display representative trajectories in dark blue    
     RT.colors = "darkblue",
     # select T1.2 to be displayed in a different color (red)
     select_RT = "T1.2", sel.color = "red",
     # Identify artificial links using dashed lines (default), but use the same 
     # color than the representative trajectories (default)
     link.lty = 2, link.color = NULL,
     # We can use other arguments in plot()
     main = "Defined representative trajectories")

# Add a legend 
legend("bottomright", c("T1.1, T2.1", "T1.2", "Link"), 
       col = c("darkblue", "red", "darkblue"), lty = c(1, 1, 2))


## ----sp_diversity-------------------------------------------------------------
# Set an ID in the abundance matrix
abundance$ID <- paste0(abundance$sampling_unit, abundance$time)

# Calculate the Shannon index in all trajectory states
abundance$Shannon <- vegan::diversity(abundance[, c("sp1", "sp2")], index = "shannon")

# Identify the states forming both representative trajectories
traj_states <- lapply(repr_traj, function(iRT){
  segments <- iRT$Segments
  seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", segments)), "-")
  traj_states <- vapply(seg_components, function(iseg){
    c(paste0(iseg[1], iseg[2]), paste0(iseg[1], iseg[3]))
  }, character(2))
  traj_states <- unique(as.character(traj_states))
  traj_states <- data.frame(ID = traj_states, RT_states = 1:length(traj_states))
})

# Associate the states of the representative trajectories with their corresponding
# value of the Shannnon index
RT_diversity <- lapply(traj_states, function(iRT){
  data <- merge(iRT, abundance, by = "ID", all.x = T)
  data <- data[order(data$RT_states), ]
})

RT_diversity$T2


## ----plot_sp_diversity, fig = TRUE, fig.height = 5, fig.width = 5, fig.align = "center"----
# Plot the variation of species diversity in T1
plot(x = RT_diversity$T1$RT_states + 1, y = RT_diversity$T1$Shannon,
     type = "l", col = "blue", xlim = c(1, 7), ylim = c(0, 0.7), 
     xlab = "RT state", ylab = "Shannon index",
     main = "Variation of species diversity")
# Add the variation of species diversity in T2
lines(x = RT_diversity$T2$RT_states, y = RT_diversity$T2$Shannon,
      col = "black")
# Add a legend
legend("topright", c("T1", "T2"), col = c("blue", "black"), lty = 1)



## ----dDis_data----------------------------------------------------------------
# Change the trajectory identifier and the number of the state
abundance_T2 <- RT_diversity$T2
abundance_T2$sampling_unit <- "T2"
abundance_T2$time <- abundance_T2$RT_states

# Add representative trajectories' data to the abundance matrix
abundance_T2 <- rbind(abundance, abundance_T2[, names(abundance)])

# Calculate state dissimilarities including the representative trajectory
state_dissim_T2 <- vegan::vegdist(abundance_T2[, c("sp1", "sp2")], method = "bray")


## ----dDis_T2------------------------------------------------------------------
# Compute dDis taking T2 as reference
dDis_T2 <- dDis(d = state_dissim_T2, d.type = "dStates", 
         trajectories = abundance_T2$sampling_unit, states = abundance_T2$time, 
         reference = "T2")
dDis_T2


## ----dDis_I_F-----------------------------------------------------------------
# dDis: reference C
dDis_I <- dDis(d = traj_dissim, d.type = "dTraj", 
         trajectories = ID_sampling, 
         reference = "I")
# dDis: reference F
dDis_F <- dDis(d = traj_dissim, d.type = "dTraj", 
         trajectories = ID_sampling, 
         reference = "F")

dDis_I
dDis_F


## ----dDis_w-------------------------------------------------------------------
# Define w values
initial_sp2 <- abundance[which(abundance$time == 1), ]$sp2

# Identify the index of the reference trajectories
ind_I <- which(ID_sampling == "I")
ind_F <- which(ID_sampling == "F")

# Compute dDis with weighted trajectories:
# Considering I as reference
dDis_I_w <- dDis(d = traj_dissim, d.type = "dTraj", 
         trajectories = ID_sampling, reference = "I", 
         w.type = "precomputed", w.values = initial_sp2[-ind_I])
# Considering F as reference
dDis_F_w <- dDis(d = traj_dissim, d.type = "dTraj", 
         trajectories = ID_sampling, reference = "F", 
         w.type = "precomputed", w.values = initial_sp2[-ind_F])
dDis_I_w
dDis_F_w


## ----dBD----------------------------------------------------------------------
# Calculate dBD
dBD(d = state_dissim, d.type = "dStates", 
    trajectories = rep(ID_sampling, each = 3), states = rep(ID_obs, 10))


## ----dEve---------------------------------------------------------------------
# Calculate dEve
dEve(d = traj_dissim, d.type = "dTraj", trajectories = ID_sampling)


## ----dEvew--------------------------------------------------------------------
# Calculate dEve weighting trajectories by the initial abundance of sp2
dEve(d = traj_dissim, d.type = "dTraj", trajectories = ID_sampling,
         w.type = "precomputed", w.values = initial_sp2)


## ----data2--------------------------------------------------------------------
# ID of the sampling units for EDR2
ID_sampling2 <- paste0("inv_", LETTERS[1:10])

# Define initial species abundances for sp3
set.seed(321)
initial$sp3 <- round(runif(10, 0, 0.5), 2)

# Define the parameters for the Lotka-Volterra model
parms2 <- c(r1 = 1, r2 = 0.1, a11 = 0.02, a21 = 0.03, a22 = 0.02, a12 = 0.01,
           r3 = 1.5, a33 = 0.02, a31 = 0.01, a32 = 0.01, a13 = 0.1, a23 = 0.1)

# We can use primer to run the simulations
simulated_abun2 <- lapply(1:nrow(initial), function(isampling){
  # Initial abundance in the sampling unit i
  initialN <- c(initial$sp1[isampling], initial$sp2[isampling], initial$sp3[isampling])
  # Simulate community dynamics
  simulated_abun <- data.frame(ode(y = initialN, times = ID_obs, func = lvcomp3, parms = parms2))
  # Add the names of the sampling units
  simulated_abun$sampling_unit <- ID_sampling2[isampling]
  # Calculate relative abundances
  simulated_abun$sp1 <- simulated_abun$X1/rowSums(simulated_abun[, c("X1", "X2", "X3")])
  simulated_abun$sp2 <- simulated_abun$X2/rowSums(simulated_abun[, c("X1", "X2", "X3")])
  simulated_abun$sp3 <- simulated_abun$X3/rowSums(simulated_abun[, c("X1", "X2", "X3")])
  
  return(simulated_abun[, c("sampling_unit", "time", "sp1", "sp2", "sp3")])
})

# Compile species abundances of all sampling units in the same data frame
abundance2 <- do.call(rbind, simulated_abun2)


## ----data3--------------------------------------------------------------------
# ID of the sampling units for EDR3
ID_sampling3 <- LETTERS[11:20]

# Define initial species abundances
set.seed(3)
initial3 <- data.frame(sampling_units = ID_sampling3, 
                   sp4 = round(runif(10), 2))

# Define the parameters for the Lotka-Volterra model
parms3 <- c(r1 = 1, r2 = 0.1, a11 = 0.2, a21 = 0.1, a22 = 0.02, a12 = 0.01)

# We can use primer to run the simulations
simulated_abun3 <- lapply(1:nrow(initial), function(isampling){
  # Initial abundance in the sampling unit i
  initialN <- c(initial3$sp4[isampling], initial$sp2[isampling])
  # Simulate community dynamics
  simulated_abun <- data.frame(ode(y = initialN, times = ID_obs, func = lvcomp2, parms = parms3))
  # Add the names of the sampling units
  simulated_abun$sampling_unit <- ID_sampling3[isampling]
  # Calculate relative abundances
  simulated_abun$sp4 <- simulated_abun$X1/rowSums(simulated_abun[, c("X1", "X2")])
  simulated_abun$sp2 <- simulated_abun$X2/rowSums(simulated_abun[, c("X1", "X2")])

  return(simulated_abun[, c("sampling_unit", "time", "sp4", "sp2")])
})

# Compile species abundances of all sampling units in the same data frame
abundance3 <- do.call(rbind, simulated_abun3)


## ----traj_dissim_allEDR-------------------------------------------------------
# Bind all abundance matrices
abundance_allEDR <- data.table::rbindlist(list(abundance, abundance2, abundance3), fill = T)
abundance_allEDR[is.na(abundance_allEDR)] <- 0

# Calculate state dissimilarities including states in the three EDRs
state_dissim_allEDR <- vegan::vegdist(abundance_allEDR[, paste0("sp", 1:4)], 
                                      method = "bray")

# Calculate trajectory dissimilarities including trajectories in the three EDRs
traj_dissim_allEDR <- ecotraj::trajectoryDistances(
  ecotraj::defineTrajectories(state_dissim_allEDR, 
                              sites = abundance_allEDR$sampling_unit,
                              surveys = abundance_allEDR$time))


## ----state_mds_allEDR, fig = TRUE, fig.height = 5, fig.width = 5, fig.align = "center"----
# Multidimensional scaling
st_mds_all <- smacof::smacofSym(state_dissim_allEDR, 
                                  ndim = nrow(as.matrix(state_dissim_allEDR))-1)
st_mds_all <- data.frame(st_mds_all$conf)

# Plot ecological trajectories in the state space
state.colorsEDRs <- rep(RColorBrewer::brewer.pal(3, "Set1"), each = 30)
# Set an empty plot
plot(st_mds_all$D1, st_mds_all$D2, type = "n",
     xlab = "Axis 1", ylab = "Axis 2",
     main = "EDRs in the state space")
# Add arrows
for (isampling in seq(1, 90, 3)) {
  # From observation 1 to observation 2
  shape::Arrows(st_mds_all[isampling, 1], st_mds_all[isampling, 2],
         st_mds_all[isampling + 1, 1], st_mds_all[isampling + 1, 2],
         col = state.colorsEDRs[isampling], arr.adj = 1)
  # From observation 2 to observation 3
  shape::Arrows(st_mds_all[isampling + 1, 1], st_mds_all[isampling + 1, 2],
         st_mds_all[isampling + 2, 1], st_mds_all[isampling + 2, 2],
         col = state.colorsEDRs[isampling], arr.adj = 1)
}

# Add a legend
legend("bottomleft", paste0("EDR", 1:3), col = unique(state.colorsEDRs), lty = 1)


## ----EDR_dissim---------------------------------------------------------------
# Compute the dissimilarities between EDRs
EDR_dissim <- dist_edr(d = traj_dissim_allEDR, d.type = "dTraj",
                       edr = rep(c("EDR1", "EDR2", "EDR3"), each = 10), 
                       metric = "dDR")
round(EDR_dissim, 2)


