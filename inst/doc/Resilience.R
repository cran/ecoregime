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

## ----dDis---------------------------------------------------------------------
# Species abundances in undisturbed states of the disturbed trajectories
# (Undisturbed states are identified by disturbed_states = 0)
abun_undist <- EDR_data$EDR3_disturbed$abundance[disturbed_states == 0] 
selcols <- names(EDR_data$EDR1$abundance)

## EDR1 ------------------------------------------------------------------------

# Species abundances in EDR1 and the undisturbed states of disturbed trajectories
abun1_undist <- rbind(EDR_data$EDR1$abundance, abun_undist[, ..selcols])

# State dissimilarities in EDR1 and the undisturbed states of disturbed trajectories
d1_undist <- vegan::vegdist(x = abun1_undist[, paste0("sp", 1:12)], method = "bray")

# dDis of the disturbed trajectories in relation to EDR1
dDis1 <- sapply(unique(abun_undist$traj), function(iundist){
  dDis(d = d1_undist, d.type = "dStates", 
       trajectories = abun1_undist$traj, states = abun1_undist$state, 
       reference = as.character(iundist))
})

## EDR2 ------------------------------------------------------------------------

# Species abundances in EDR2 and the undisturbed states of disturbed trajectories
abun2_undist <- rbind(EDR_data$EDR2$abundance, abun_undist[, ..selcols])

# State dissimilarities in EDR2 and the undisturbed states of disturbed trajectories
d2_undist <- vegan::vegdist(x = abun2_undist[, paste0("sp", 1:12)], method = "bray")

# dDis of the disturbed trajectories in relation to EDR2
dDis2 <- sapply(unique(abun_undist$traj), function(iundist){
  dDis(d = d2_undist, d.type = "dStates", 
       trajectories = abun2_undist$traj, states = abun2_undist$state, 
       reference = as.character(iundist))
})

## EDR3 ------------------------------------------------------------------------

# Species abundances in EDR3 and the undisturbed states of disturbed trajectories
abun3_undist <- rbind(EDR_data$EDR3$abundance, abun_undist[, ..selcols])

# State dissimilarities in EDR3 and the undisturbed states of disturbed trajectories
d3_undist <- vegan::vegdist(x = abun3_undist[, paste0("sp", 1:12)], method = "bray")

# dDis of the disturbed trajectories in relation to EDR3
dDis3 <- sapply(unique(abun_undist$traj), function(iundist){
  dDis(d = d3_undist, d.type = "dStates", 
       trajectories = abun3_undist$traj, states = abun3_undist$state, 
       reference = as.character(iundist))
})

## Compare dynamic dispersion --------------------------------------------------

# Compare dDis values for the three EDRs
dDis_df <- data.frame(EDR1 = dDis1, EDR2 = dDis2, EDR3 = dDis3)


## ----dDis_comparison, echo=FALSE----------------------------------------------
knitr::kable(dDis_df, digits = 3)


## ----retra--------------------------------------------------------------------
# State dissimilarities for EDR3 (considering only the undisturbed trajectories)
d_EDR3 <- vegan::vegdist(EDR_data$EDR3$abundance[, paste0("sp", 1:12)])

# Representative trajectories
retra <- retra_edr(d = d_EDR3, 
                   trajectories = EDR_data$EDR3$abundance$traj,
                   states = EDR_data$EDR3$abundance$state, minSegs = 5)


## ----summary_retra------------------------------------------------------------
# Summarize retra
summary(retra)

# Define T4 as the unique representative trajectory and generate an object of class 'RETRA'
retra_ref <- define_retra(data = retra$T4$Segments, d = d_EDR3, 
                          trajectories = EDR_data$EDR3$abundance$traj, 
                          states = EDR_data$EDR3$abundance$state,
                          retra = retra)


## ----plot_retra, fig.dim = c(6,5)---------------------------------------------
# Plot EDR3 and its representative trajectories
plot(retra, d = d_EDR3, 
     trajectories = EDR_data$EDR3$abundance$traj, 
     states = EDR_data$EDR3$abundance$state, select_RT = "T4",
     main = "Representative trajectories in EDR3")
legend("topleft", c("Representative trajectory 'T4'", 
                    "Other representative trajectories", 
                    "Individual trajectories in EDR3"),
       lty = 1, col = c("red", "black", "grey"), bty = "n")

## ----graphical_indices, echo=FALSE, fig.dim=c(6,5)----------------------------
par(mar = c(0,0,0,0))

# Representative trajectory
plot(x = 0:10, y = c(0,0,0:8), type = "n", axes = F, xlab = "", ylab = "")
for (i in 0:9) {
  shape::Arrows(x0 = i, y0 = 1, x1 = i + 1, y1 = 1, lwd = 2,  arr.adj = 1)
}
text(x = 0.2, y = 1.3, "RT")

# Reference
lines(x = c(0, 10), y = c(2.5, 2.5), col = RColorBrewer::brewer.pal(3, "Set1")[1], lty = 3, lwd = 2)

# Disturbed trajectory
for (i in 0.5:1.5) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 2.5, y0 = 2.5, x1 = 4, y1 = 7, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
shape::Arrows(x0 = 4, y0 = 7, x1 = 5, y1 = 6, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
shape::Arrows(x0 = 5, y0 = 6, x1 = 5.5, y1 = 4, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
shape::Arrows(x0 = 5.5, y0 = 4, x1 = 6.5, y1 = 3.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
shape::Arrows(x0 = 6.5, y0 = 3.5, x1 = 7.5, y1 = 3.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 0.7, y = 2.8, "DT", col = RColorBrewer::brewer.pal(3, "Set1")[1])

# Resistance
lines(x = c(2.5, 4)-0.2, y = c(2.5, 7), col = RColorBrewer::brewer.pal(3, "Set1")[2])
text(x = 2.8, y = 5, "Rt", col = RColorBrewer::brewer.pal(3, "Set1")[2])

# Amplitude
plotrix::draw.arc(x = 2.5, y = 2.5, radius = 0.7, angle1 = 0, deg2 = 70,
                  col = RColorBrewer::brewer.pal(3, "Set1")[3])
text(x = 3.3, y = 3.1, "A", col = RColorBrewer::brewer.pal(3, "Set1")[3])

# Recovery
lines(x = c(4, 5), y = c(7, 7), lty = 3, lwd = 2, 
      col = RColorBrewer::brewer.pal(4, "Set1")[4])
plotrix::draw.arc(x = 4, y = 7, radius = 0.6, deg1 = 0, deg2 = -45,
                  col = RColorBrewer::brewer.pal(4, "Set1")[4])
text(x = 4.9, y = 6.7, "Rc", col = RColorBrewer::brewer.pal(4, "Set1")[4])

# Net change
lines(x = c(2.5, 6.5), y = c(2.5, 3.5), lty = 3, lwd = 2,
      col = RColorBrewer::brewer.pal(9, "Set1")[9])
plotrix::draw.arc(x = 2.5, y = 2.5, radius = 1.5, angle1 = 0, deg2 = 15,
                  col = RColorBrewer::brewer.pal(9, "Set1")[9])
text(x = 4.3, y = 2.7, "NC", col = RColorBrewer::brewer.pal(9, "Set1")[9])

# Legend
legend("topright",
       c("RT: Representative trajectory",
         "DT: Disturbed trajectory",
         "Rt: Resistance",
         "A: Amplitude",
         "Rc: Recovery",
         "NC: Net change"),
       col = c("black", RColorBrewer::brewer.pal(9, "Set1")[c(1:4, 9)]),
       lty = 1, text.col = c("black", RColorBrewer::brewer.pal(9, "Set1")[c(1:4, 9)]), bty = "n")

# States
text(x = c(2.5, 4, 5.2, 5.3, 6.5, 7.5), y = c(2.2, 7.4, 6, 4, 3.8, 3.8),
     0:5, cex = 0.9, col = RColorBrewer::brewer.pal(5, "Set1")[1])


## ----graphical_resistance, echo=FALSE, fig.dim=c(6,5)-------------------------

par(mar = c(0,0,0,0))

# Representative trajectory
plot(x = c(0, 10), y = c(0, 4.5), type = "n", axes = F, xlab = "", ylab = "")
for (i in 0:9) {
  shape::Arrows(x0 = i, y0 = 0, x1 = i + 1, y1 = 0, lwd = 2,  arr.adj = 1)
}
text(x = 0.2, y = 0.2, "RT")

# Reference
lines(x = c(0, 10), y = c(2.5, 2.5), col = RColorBrewer::brewer.pal(3, "Set1")[1], lty = 3, lwd = 2)

# Disturbed trajectories
for (i in 5.5:6.5) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 7.5, y0 = 2.5, x1 = 8.5, y1 = 4, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 5.7, y = 2.7, "DT2", col = RColorBrewer::brewer.pal(3, "Set1")[1])

for (i in 0.5:1.5) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 2.5, y0 = 2.5, x1 = 8, y1 = 1, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 0.7, y = 2.7, "DT1", col = RColorBrewer::brewer.pal(3, "Set1")[1])


# Resistance
lines(x = c(7.5, 8.5)-0.2, y = c(2.5, 4)+0.05, lwd = 2, col = RColorBrewer::brewer.pal(3, "Set1")[2])
text(x = 7.2, y = 3.4, "High Rt", col = RColorBrewer::brewer.pal(3, "Set1")[2])
lines(x = c(2.5, 8)-0.1, y = c(2.5, 1)-0.2, lwd = 2, col = RColorBrewer::brewer.pal(3, "Set1")[2])
text(x = 5, y = 1.3, "Low Rt", col = RColorBrewer::brewer.pal(3, "Set1")[2])


# Legend
legend("topleft",
       c("RT: Representative trajectory",
         "DT: Disturbed trajectories",
         "Rt: Resistance"),
       col = c("black", RColorBrewer::brewer.pal(5, "Set1")),
       lty = 1, text.col = c("black", RColorBrewer::brewer.pal(5, "Set1")), bty = "n")

# States
text(x = c(2.5, 7.5, 8.7, 8.3), y = c(2.7, 2.3, 4.2, 1),
     c(0, 0, 1, 1), cex = 0.9, col = RColorBrewer::brewer.pal(5, "Set1")[1])


## ----resistance---------------------------------------------------------------
# To calculate resistance, we need a state dissimilarity matrix for the disturbed trajectories
d_disturbed <- vegan::vegdist(EDR_data$EDR3_disturbed$abundance[, paste0("sp", 1:12)], 
                              method = "bray")

# Compute resistance
# Note that the disturbed states are identified by disturbed_states = 1
Rt <- resistance(d = d_disturbed, 
                 trajectories = EDR_data$EDR3_disturbed$abundance$traj, 
                 states = EDR_data$EDR3_disturbed$abundance$state, 
                 disturbed_trajectories = unique(EDR_data$EDR3_disturbed$abundance$traj), 
                 disturbed_states = EDR_data$EDR3_disturbed$abundance[disturbed_states == 1]$state)


## ----Rt_results, echo=FALSE---------------------------------------------------
knitr::kable(Rt, row.names = F, col.names = c("Disturbed trajectories", "Rt"), digits = 3)

## ----graphical_amplitude, echo=FALSE, fig.dim=c(6,5)--------------------------

par(mar = c(0,0,0,0))

# Representative trajectory
plot(x = c(0, 16), y = c(0, 5.5), type = "n", axes = F, xlab = "", ylab = "")
for (i in 0:14) {
  shape::Arrows(x0 = i, y0 = 0, x1 = i + 1, y1 = 0, lwd = 2, arr.adj = 1)
}
text(x = 0.2, y = 0.2, "RT")

# Reference
lines(x = c(0, 15), y = c(2.5, 2.5), col = RColorBrewer::brewer.pal(3, "Set1")[1], lty = 3, lwd = 2)

# Disturbed trajectory 1
for (i in 0.5:1.5) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 2.5, y0 = 2.5, x1 = 4, y1 = 0.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 0.7, y = 2.7, "DT1", col = RColorBrewer::brewer.pal(3, "Set1")[1])

# Disturbed trajectory 2
for (i in 5:6) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 7, y0 = 2.5, x1 = 14, y1 = 4.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 5, y = 2.7, "DT2", col = RColorBrewer::brewer.pal(3, "Set1")[1])



# # Amplitude
lines(x = c(4, 4), y = c(0.5, 2.5), lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[3])
text(x = 5, y = 1.5, expression(A[abs1] ~ "< 0"), col = RColorBrewer::brewer.pal(5, "Set1")[3])
plotrix::draw.arc(x = 2.5, y = 2.5, radius = 0.5, deg2 = -70, lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[3])
text(x = 3.8, y = 2.2, expression(A[rel1] ~ "< 0"), col = RColorBrewer::brewer.pal(5, "Set1")[3])

lines(x = c(14, 14), y = c(2.5, 4.5), lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[3])
text(x = 15, y = 3.5, expression(A[abs2] ~ "> 0"), col = RColorBrewer::brewer.pal(5, "Set1")[3])
plotrix::draw.arc(x = 7, y = 2.5, radius = 1.5, deg2 = 32, lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[3])
text(x = 9.5, y = 2.75, expression(A[rel2] ~ "> 0"), col = RColorBrewer::brewer.pal(5, "Set1")[3])

text(x = c(0), y = c(4.5), adj = 0, expression("|"~A[abs1]~"| = |"~A[abs2]~"|"), cex = 0.9, 
     col = RColorBrewer::brewer.pal(5, "Set1")[3])
text(x = c(0), y = c(4.2), adj = 0, expression("|"~A[rel1]~"| > |"~A[rel2]~"|"), cex = 0.9, 
     col = RColorBrewer::brewer.pal(5, "Set1")[3])

# Legend
legend("topleft",
       c("RT: Representative trajectory",
         "DT: Disturbed trajectories",
         "A: Amplitude"),
       col = c("black", RColorBrewer::brewer.pal(5, "Set1")[c(1, 3)]),
       lty = 1, text.col = c("black", RColorBrewer::brewer.pal(5, "Set1")[c(1, 3)]), bty = "n")

# States
text(x = c(2.5, 7, 4, 14), y = c(2.7, 2.3, 0.3, 4.7),
     c(0, 0, 1, 1), cex = 0.9, col = RColorBrewer::brewer.pal(5, "Set1")[1])


## ----amplitude----------------------------------------------------------------
# We need a state dissimilarity matrix containing the states of the disturbed 
# trajectories and the representative trajectory taken as the reference
abun <- rbind(EDR_data$EDR3$abundance, EDR_data$EDR3_disturbed$abundance, fill = T)
d <- vegan::vegdist(abun[, paste0("sp", 1:12)], method = "bray")

# Compute amplitude
A <- amplitude(d = d,
               trajectories = abun$traj,
               states = abun$state,
               disturbed_trajectories = abun[disturbed_states == 1]$traj,
               disturbed_states = abun[disturbed_states == 1]$state,
               reference = retra_ref, method = "nearest_state")


## ----A_results, echo=FALSE----------------------------------------------------
knitr::kable(A, col.names = c("Disturbed trajectories", "Reference", "A~abs~", "A~rel~"), digits = 3)

## ----graphical_recovery, echo=FALSE, fig.dim=c(6,5)---------------------------

par(mar = c(0,0,0,0))

# Representative trajectory
plot(x = c(0, 16), y = c(0, 5.5), type = "n", axes = F, xlab = "", ylab = "")
for (i in 0:14) {
  shape::Arrows(x0 = i, y0 = 1, x1 = i + 1, y1 = 1, lwd = 2, arr.adj = 1)
}
text(x = 0.2, y = 1.2, "RT")

# Reference
lines(x = c(0, 15), y = c(2.5, 2.5), col = RColorBrewer::brewer.pal(3, "Set1")[1], lty = 3, lwd = 2)

# Disturbed trajectory 1
for (i in 0.5:1.5) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 2.5, y0 = 2.5, x1 = 4, y1 = 4.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
shape::Arrows(x0 = 4, y0 = 4.5, x1 = 5, y1 = 3.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 0.7, y = 2.7, "DT1", col = RColorBrewer::brewer.pal(3, "Set1")[1])

# Disturbed trajectory 2
for (i in 6:7) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 8, y0 = 2.5, x1 = 9.5, y1 = 4.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
shape::Arrows(x0 = 9.5, y0 = 4.5, x1 = 14, y1 = 3.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 6.2, y = 2.7, "DT2", col = RColorBrewer::brewer.pal(3, "Set1")[1])

# States
text(x = c(2.5, 3.5, 5, 8, 9.5, 14), y = c(2.3, 4.5, 3.3, 2.3, 4.7, 3.3),
     c(0:2, 0:2), cex = 0.9, col = RColorBrewer::brewer.pal(5, "Set1")[1])


# Recovery DT1
lines(x = c(4, 5), y = c(4.5, 4.5), lty = 3, col = RColorBrewer::brewer.pal(5, "Set1")[4])
lines(x = c(5, 5), y = c(3.5, 4.5), lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[4])
text(x = 5.2, y = 4, adj = 0, expression(Rc[abs1] ~ "> 0"), col = RColorBrewer::brewer.pal(5, "Set1")[4])

# lines(x = c(4, 4), y = c(3.5, 4.5), lty = 3, col = RColorBrewer::brewer.pal(5, "Set1")[4])
plotrix::draw.arc(x = 4, y = 4.5, radius = 0.5, deg1 = 0, deg2 = -65, lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[4])
text(x = 4, y = 4.7, adj = 0, expression(Rc[rel1] ~ "> 0"), col = RColorBrewer::brewer.pal(5, "Set1")[4])

# Recovery DT2
lines(x = c(9.5, 14), y = c(4.5, 4.5), lty = 3, col = RColorBrewer::brewer.pal(5, "Set1")[4])
lines(x = c(14, 14), y = c(3.5, 4.5), lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[4])
text(x = 15, y = 4, expression(Rc[abs2] ~ "> 0"), col = RColorBrewer::brewer.pal(5, "Set1")[4])

# lines(x = c(9.5, 9.5), y = c(3.5, 4.5), lty = 3, col = RColorBrewer::brewer.pal(5, "Set1")[4])
plotrix::draw.arc(x = 9.5, y = 4.5, radius = 1, deg1 = 0, deg2 = -25, lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[4])
text(x = 11, y = 4.3, adj = 0, expression(Rc[rel2] ~ "> 0"), col = RColorBrewer::brewer.pal(5, "Set1")[4])

text(x = c(0), y = c(4.5), adj = 0, expression(Rc[abs1]~" = "~Rc[abs2]), cex = 0.9, 
     col = RColorBrewer::brewer.pal(5, "Set1")[4])
text(x = c(0), y = c(4.2), adj = 0, expression(Rc[rel1]~" > "~Rc[rel2]), cex = 0.9, 
     col = RColorBrewer::brewer.pal(5, "Set1")[4])

# Legend
legend("topleft",
       c("RT: Representative trajectory",
         "DT: Disturbed trajectories",
         "Rc: Recovery"),
       col = c("black", RColorBrewer::brewer.pal(5, "Set1")[c(1, 4)]),
       lty = 1, text.col = c("black", RColorBrewer::brewer.pal(5, "Set1")[c(1, 4)]), bty = "n")



## ----graphical_recovery2, echo=FALSE, fig.dim=c(6,5)--------------------------

par(mar = c(0,0,0,0))

# Representative trajectory
plot(x = c(0, 16), y = c(0, 5.5), type = "n", axes = F, xlab = "", ylab = "")
for (i in 0:14) {
  shape::Arrows(x0 = i, y0 = 0, x1 = i + 1, y1 = 0, lwd = 2, arr.adj = 1)
}
text(x = 0.2, y = 0.2, "RT")

# Reference
lines(x = c(0, 15), y = c(2.5, 2.5), col = RColorBrewer::brewer.pal(3, "Set1")[1], lty = 3, lwd = 2)

# Disturbed trajectory 3
for (i in 0.5:1.5) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 2.5, y0 = 2.5, x1 = 4, y1 = 0.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
shape::Arrows(x0 = 4, y0 = 0.5, x1 = 5, y1 = 1.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 0.7, y = 2.7, "DT3", col = RColorBrewer::brewer.pal(3, "Set1")[1])

# Disturbed trajectory 4
for (i in 6:7) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 8, y0 = 2.5, x1 = 9.5, y1 = 4.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
shape::Arrows(x0 = 9.5, y0 = 4.5, x1 = 14, y1 = 5.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 6.2, y = 2.7, "DT4", col = RColorBrewer::brewer.pal(3, "Set1")[1])

# States
text(x = c(2.5, 3.4, 5.5, 8, 9.2, 14.2), y = c(2.7, 0.6, 1.5, 2.3, 4.5, 5.6),
     c(0:2, 0:2), cex = 0.9, col = RColorBrewer::brewer.pal(5, "Set1")[1])


# Recovery DT3
lines(x = c(4, 5), y = c(0.5, 0.5), lty = 3, col = RColorBrewer::brewer.pal(5, "Set1")[4])
lines(x = c(5, 5), y = c(0.5, 1.5), lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[4])
text(x = 5.2, y = 1, adj = 0, expression(Rc[abs3] ~ "< 0"), col = RColorBrewer::brewer.pal(5, "Set1")[4])

# lines(x = c(4, 4), y = c(0.5, 1.5), lty = 3, col = RColorBrewer::brewer.pal(5, "Set1")[4])
plotrix::draw.arc(x = 4, y = 0.5, radius = 0.5, deg1 = 0, deg2 = 65, lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[4])
text(x = 3.8, y = 0.3, adj = 0, expression(Rc[rel3] ~ "< 0"), col = RColorBrewer::brewer.pal(5, "Set1")[4])

# Recovery DT4
lines(x = c(9.5, 14), y = c(4.5, 4.5), lty = 3, col = RColorBrewer::brewer.pal(5, "Set1")[4])
lines(x = c(14, 14), y = c(5.5, 4.5), lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[4])
text(x = 15, y = 5, expression(Rc[abs4] ~ "< 0"), col = RColorBrewer::brewer.pal(5, "Set1")[4])

# lines(x = c(9.5, 9.5), y = c(5.5, 4.5), lty = 3, col = RColorBrewer::brewer.pal(5, "Set1")[4])
plotrix::draw.arc(x = 9.5, y = 4.5, radius = 1.3, deg1 = 0, deg2 = 25, lwd = 2, col = RColorBrewer::brewer.pal(5, "Set1")[4])
text(x = 11.1, y = 4.7, adj = 0, expression(Rc[rel4] ~ "< 0"), col = RColorBrewer::brewer.pal(5, "Set1")[4])

text(x = c(0), y = c(5.5), adj = 0, expression(Rc[abs3]~" = "~Rc[abs4]), cex = 0.9, 
     col = RColorBrewer::brewer.pal(5, "Set1")[4])
text(x = c(0), y = c(5.2), adj = 0, expression(Rc[rel3]~" < "~Rc[rel4]~" < "~Rc[rel2]~" < "~Rc[rel1]), cex = 0.9, 
     col = RColorBrewer::brewer.pal(5, "Set1")[4])



## ----recovery-----------------------------------------------------------------
# Compute recovery using the same data used for amplitude
Rc <- recovery(d = d, 
               trajectories = abun$traj,
               states = abun$state, 
               disturbed_trajectories = abun[disturbed_states == 1]$traj,
               disturbed_states = abun[disturbed_states == 1]$state, 
               reference = retra_ref, method = "nearest_state")


## ----plot_recovery, fig.dim=c(6,5)--------------------------------------------
# Number of states after the disturbed state
Rc <- data.table::data.table(Rc)
Rc[, ID_post := 1:(.N), by = disturbed_trajectories]

# Plot absolute recovery over time
plot(x = range(Rc$ID_post), y = range(Rc$Rc_abs), type = "n",
     xlab = "Nb. states after disturbance", ylab = "Absolute recovery",
     main = "Variation of absolute recovery")
for (i in unique(Rc$disturbed_trajectories)) {
  lines(Rc[disturbed_trajectories == i, c("ID_post", "Rc_abs")],
        col = which(unique(Rc$disturbed_trajectories) %in% i) + 1)
}
legend("bottomleft", legend = unique(Rc$disturbed_trajectories), lty = 1, 
       col = seq_along(unique(Rc$disturbed_trajectories)) + 1, bty = "n")

# Plot relative recovery over time
plot(x = range(Rc$ID_post), y = range(Rc$Rc_rel), type = "n",
     xlab = "Nb. states after disturbance", ylab = "Relative recovery",
     main = "Variation of relative recovery")
for (i in unique(Rc$disturbed_trajectories)) {
  lines(Rc[disturbed_trajectories == i, c("ID_post", "Rc_rel")],
        col = which(unique(Rc$disturbed_trajectories) %in% i) + 1)
}
legend("topright", legend = unique(Rc$disturbed_trajectories), lty = 1, 
       col = seq_along(unique(Rc$disturbed_trajectories)) + 1, bty = "n")




## ----graphical_NetChange, echo=FALSE, fig.dim=c(6,5)--------------------------

par(mar = c(0,0,0,0))

# Representative trajectory
plot(x = c(0, 16), y = c(0, 5.5), type = "n", axes = F, xlab = "", ylab = "")
for (i in 0:14) {
  shape::Arrows(x0 = i, y0 = 1, x1 = i + 1, y1 = 1, lwd = 2, arr.adj = 1)
}
text(x = 0.2, y = 1.2, "RT")

# Reference
lines(x = c(0, 15), y = c(2.5, 2.5), col = RColorBrewer::brewer.pal(3, "Set1")[1], lty = 3, lwd = 2)

# Disturbed trajectory 1
for (i in 0.5:1.5) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 2.5, y0 = 2.5, x1 = 4, y1 = 4.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
shape::Arrows(x0 = 4, y0 = 4.5, x1 = 5, y1 = 3.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 0.7, y = 2.7, "DT1", col = RColorBrewer::brewer.pal(3, "Set1")[1])

# Disturbed trajectory 2
for (i in 6:7) {
  shape::Arrows(x0 = i, y0 = 2.5, x1 = i + 1, y1 = 2.5, lwd = 2, arr.adj = 1,
                col = RColorBrewer::brewer.pal(3, "Set1")[1])
}
shape::Arrows(x0 = 8, y0 = 2.5, x1 = 9.5, y1 = 3.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
shape::Arrows(x0 = 9.5, y0 = 3.5, x1 = 14, y1 = 1.5, lwd = 2, arr.adj = 1,
              col = RColorBrewer::brewer.pal(3, "Set1")[1])
text(x = 6.2, y = 2.7, "DT2", col = RColorBrewer::brewer.pal(3, "Set1")[1])

# States
text(x = c(2.5, 4, 5.5, 8, 9.5, 14), y = c(2.3, 4.7, 3.5, 2.3, 3.7, 1.3),
     c(0:2, 0:2), cex = 0.9, col = RColorBrewer::brewer.pal(5, "Set1")[1])


# Net change DT1
lines(x = c(2.5, 5), y = c(2.5, 3.5), lty = 3, col = RColorBrewer::brewer.pal(9, "Set1")[9])
lines(x = c(5, 5), y = c(2.5, 3.5), lwd = 2, col = RColorBrewer::brewer.pal(9, "Set1")[9])
text(x = 5.2, y = 3, adj = 0, expression(NC[abs1] ~ "> 0"), col = RColorBrewer::brewer.pal(9, "Set1")[9])

plotrix::draw.arc(x = 2.5, y = 2.5, radius = 0.6, deg2 = 40, lwd = 2, col = RColorBrewer::brewer.pal(9, "Set1")[9])
text(x = 3.2, y = 2.7, adj = 0, expression(NC[rel1] ~ "> 0"), col = RColorBrewer::brewer.pal(9, "Set1")[9])

# Net change DT2
lines(x = c(8, 14), y = c(2.5, 1.5), lty = 3, col = RColorBrewer::brewer.pal(9, "Set1")[9])
lines(x = c(14, 14), y = c(2.5, 1.5), lwd = 2, col = RColorBrewer::brewer.pal(9, "Set1")[9])
text(x = 15, y = 2, expression(NC[abs2] ~ "< 0"), col = RColorBrewer::brewer.pal(9, "Set1")[9])

plotrix::draw.arc(x = 8, y = 2.5, radius = 1.7, deg2 = -20, lwd = 2, col = RColorBrewer::brewer.pal(9, "Set1")[9])
text(x = 10, y = 2.3, adj = 0, expression(NC[rel2] ~ "< 0"), col = RColorBrewer::brewer.pal(9, "Set1")[9])

text(x = c(0), y = c(4.5), adj = 0, expression("|"~NC[abs1]~"| = |"~NC[abs2]~"|"), cex = 0.9, 
     col = RColorBrewer::brewer.pal(9, "Set1")[9])
text(x = c(0), y = c(4.2), adj = 0, expression("|"~NC[rel1]~"| > |"~NC[rel2]~"|"), cex = 0.9, 
     col = RColorBrewer::brewer.pal(9, "Set1")[9])

# Legend
legend("topleft",
       c("RT: Representative trajectory",
         "DT: Disturbed trajectories",
         "NC: Net change"),
       col = c("black", RColorBrewer::brewer.pal(9, "Set1")[c(1, 9)]),
       lty = 1, text.col = c("black", RColorBrewer::brewer.pal(9, "Set1")[c(1, 9)]), bty = "n")



## ----net_change---------------------------------------------------------------
# Compute net change using the same data used for amplitude
NC <- net_change(d = d, 
                 trajectories = abun$traj,
                 states = abun$state,
                 disturbed_trajectories = abun[disturbed_states == 1]$traj,
                 disturbed_states = abun[disturbed_states == 1]$state,
                 reference = retra_ref, method = "nearest_state")


## ----plot_NetChange, fig.dim=c(6,5)-------------------------------------------
# ID post-disturbance states
NC <- data.table::data.table(NC)
NC[, ID_post := 1:(.N), by = disturbed_trajectories]

# Plot absolute net change over time
plot(x = range(NC$ID_post), y = range(NC$NC_abs), type = "n",
     xlab = "Nb. states after disturbance", ylab = "Absolute net change",
     main = "Variation of absolute net change")
for (i in unique(NC$disturbed_trajectories)) {
  lines(NC[disturbed_trajectories == i, c("ID_post", "NC_abs")],
        col = which(unique(NC$disturbed_trajectories) %in% i) + 1)
}
legend("topleft", legend = unique(NC$disturbed_trajectories), lty = 1, 
       col = seq_along(unique(NC$disturbed_trajectories)) + 1, bty = "n")

# Plot relative net change over time
plot(x = range(NC$ID_post), y = range(NC$NC_rel), type = "n",
     xlab = "Nb. states after disturbance", ylab = "Relative net change",
     main = "Variation of relative net change")
for (i in unique(NC$disturbed_trajectories)) {
  lines(NC[disturbed_trajectories == i, c("ID_post", "NC_rel")],
        col = which(unique(NC$disturbed_trajectories) %in% i) + 1)
}
legend("bottomleft", legend = unique(NC$disturbed_trajectories), lty = 1, 
       col = seq_along(unique(NC$disturbed_trajectories)) + 1, bty = "n")


## ----A_Rc_NC------------------------------------------------------------------
# Merge the results for resistance, amplitude, recovery, and net change
results <- Reduce(function(x, y) merge(x, y, all = T),
                  list(Rt, A, Rc, NC))
results <- results[which(results$ID_post == 14), 
                   c("disturbed_trajectories", "Rt", "A_abs", "Rc_abs", "NC_abs")]

## ----all_results, echo=FALSE--------------------------------------------------
knitr::kable(results, row.names = F, digits = 3, 
             col.names = c("Disturbed trajectories", "Rt", "A~abs~", "Rc~abs~", "NC~abs~"))

## ----dDis_post, echo = T------------------------------------------------------
# Species abundances in post-disturbance states of the disturbed trajectories
# The states after the release of the disturbance are identified by disturbed_states > 1
abun_post <- EDR_data$EDR3_disturbed$abundance[disturbed_states > 1] 
selcols <- names(EDR_data$EDR3$abundance)

# Species abundances in EDR3 and the post-disturbance states of disturbed trajectories
abun3_post <- rbind(EDR_data$EDR3$abundance, abun_post[, ..selcols])

# State dissimilarities in EDR3 and the post-disturbance states of disturbed trajectories
d3_post <- as.matrix(vegan::vegdist(x = abun3_post[, paste0("sp", 1:12)], method = "bray"))

# dDis of each disturbed trajectory
dDis_dist <- sapply(unique(abun_post$traj), function(idist){
  rm_disturbed <- unique(abun_post$traj[which(abun_post$traj != idist)])
  irm <- which(abun3_post$traj %in% rm_disturbed)
  dDis(d = d3_post[-irm, -irm], d.type = "dStates", 
       trajectories = abun3_post$traj[-irm], 
       states = abun3_post$state[-irm], 
       reference = as.character(idist))
})


## ----dDis_post_tb, echo = F---------------------------------------------------
knitr::kable(data.frame(disturbed_trajectory = 31:33, dDis = dDis_dist), digits = 3,
             row.names = F, col.names = c("Disturbed trajectories", "dDis"))

## ----plot_disturbed, fig.dim=c(7,6)-------------------------------------------
# Plot EDR3 and its representative trajectories
plot(retra_ref, d = d, 
     trajectories = abun$traj, 
     states = abun$state, 
     traj.colors = c(rep("grey", 30), 2:4), 
     main = "Comparison of disturbed trajectories and EDR3")
legend("topleft", c("Representative trajectory", 
                    "Individual trajectories in EDR3",
                    "Trajectory 31",
                    "Trajectory 32",
                    "Trajectory 33"),
       lty = 1, col = c("black", "grey", 2:4), bty = "n")


