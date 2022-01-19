library(mvtnorm)
library(sp)
library(lattice)
library(spatstat)
library(data.table)
library(raster)
library(rasterVis)
library(ROCR)
library(gridExtra)
library(grid)
library(latticeExtra)
library(NLMR)
library(viridisLite)
library(abind)
library(TMB)
library(purrr)
library(ppmlasso)

source("Share Functions_original.R") # This function comes from :
# Renner, IW, Louvrier, J, Gimenez, O. Combining multiple data sources in species distribution models while accounting for spatial dependence and overfitting with combined penalized likelihood maximization. Methods Ecol Evol. 2019; 10: 2118– 2128. https://doi-org.ezproxy.newcastle.edu.au/10.1111/2041-210X.13297

source("DynSharfunction_TMB18062020_GPlasso.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("ERROR - no arguments", call.=FALSE)
}


# Simulation of data

## Covariates
### Environment

XY = expand.grid(seq(0, 30, 0.5), seq(0, 30, 0.5))
names(XY) = c("X", "Y")
X = XY[,1]
Y = XY[,2]

quad = expand.grid(seq(0, 30, 0.5), seq(0, 30, 0.5))
names(quad) = c("X", "Y")

win   = owin(xrange = c(0, 30), yrange = c(0, 30))


# Setting up covariates


set.seed(3)
n_centers = 2000
centers = matrix(runif(n_centers, 1, 29), n_centers/2, 2)


set.seed(3)
x1 = makecovar(110)
set.seed(19)
x2 = makecovar(200)
set.seed(61)
x1.2 = makecovar(130)
set.seed(11)
x2.2 = makecovar(90)
set.seed(7)
x3 = makecovar(180)
set.seed(88)
x4 = makecovar(220)


Datav1 = scale(x1)
Datav2 = scale(x2)
Datav1.2 = scale(x1.2)
Datav2.2 = scale(x2.2)
Datav3 = scale(x3)
Datav4 = scale(x4)


### Bias

## Create road network for observer bias

set.seed(3)
network1 = addnetwork(0, 15, angle = 0, p_split = 0.03, split_decay = 0.2,
                      nmult = 0.6, min_angle = 60, sdmult = 1.5, max_levels = 3)
set.seed(12)
network2 = addnetwork(0, 15, angle = 40, p_split = 0.035, split_decay = 0.25,
                      nmult = 0.6, min_angle = 60, sdmult = 1.5, max_levels = 3)
set.seed(5)
network3 = addnetwork(0, 15, angle = -40, p_split = 0.015, split_decay = 0.2,
                      nmult = 0.6, min_angle = 60, sdmult = 1.5, max_levels = 3)
combined = c(network1, network2, network3)
combinednew = fixnodes(combined)
plotroadlist(combinednew, main = "Simulated Road Network", asp = 1)

combinednewnodes = nodelist(combinednew)
summary(combinednewnodes)

PO.data = ppp(x = c(quad$X, combinednewnodes$X),
              y = c(quad$Y, combinednewnodes$Y),
              marks = c(rep("Quad", dim(quad)[1]), rep("Road", dim(combinednewnodes)[1])),
              window = owin(c(0, 30), c(0, 30)))
nndists = nndist(PO.data, k = 1, by = as.factor(marks(PO.data)))

# distance to the nearest road
d_rd1 = nndists[1:dim(quad)[1], 2] 

ldrd = levelplot(d_rd1 ~ quad[,1] + quad[,2], ylab="", xlab="",
                 main = "d_rd: bias road", asp = "iso", colorkey=F)
# Vegetation/Habitat facilitate the detection
set.seed(28)

z3 = makecovar(90) # detection for occupancy
l3 = levelplot(z3 ~ quad[,1] + quad[,2], main = "z3: detection covariate", 
               asp = "iso", ylab="", xlab="", colorkey=F)
l3

### Covariates set up

# from raster to dataframe for the landscape variables
Datav1 = scale(as.data.frame(flip(v1,2)))
Datav2 = scale(as.data.frame(flip(v2,2)))
Datav1.2 = scale(as.data.frame(flip(v1.2,2)))
Datav2.2 = scale(as.data.frame(flip(v2.2,2)))
Datav3 = scale(as.data.frame(flip(v3,2)))
Datav4 = scale(as.data.frame(flip(v4,2)))
d_rd = scale(sqrt(d_rd1))
v8 = scale(z3)

set.seed(50)
d1 = makecovar(40)
set.seed(105)
d2 = makecovar(40)
d1 = scale(d1)
d2 = scale(d2)


covframe = data.frame(Datav1, Datav2, Datav1.2, Datav2.2, 
                      Datav3, Datav4, d_rd, v8, d1, d2)
cor(covframe)

vmat = as.matrix(data.frame(1, Datav1, Datav2, Datav1.2, 
                            Datav2.2, Datav3, Datav4, d_rd, v8))

cor(vmat[,-1])

colnames(vmat)=c("int", "v1", "v2", "v1.1", "v2.2", "v3", "v4","d_rd", "v8")
vmat.df= as.data.frame(vmat)

# Matrice without covariates for bias and the dummy cov
env.mat = vmat[,1:5]


## True species intensity
# Generate true PPM coefficients based on linear and quadratic terms 
# of 2 covariates
sp_coef = c(1.8, 0.5, 0.6, -0.5, -0.4) # If positive integer the intensity is really high..but more points
sp_int = exp(env.mat %*% sp_coef)

# True intensity surface
sp_int_im = as.im(data.frame(x = quad$X, y = quad$Y, z = sp_int))

# Generate points from the true intensity
set.seed(14)
sp_sim = rpoispp(sp_int_im)
sp_sim

#sp_xy = data.frame(X = sp_sim$x, Y = sp_sim$y)
plot(sp_sim)



## Year 1 introduce colonization/extinction using the covar v5 and v6
PO1d_xy = data.frame(X = sp_sim$x, Y = sp_sim$y)

# environmental info at quad points and species locations
quads = data.frame(X = quad$X, Y = quad$Y, V1 = vmat.df$v1, V2 = vmat.df$v2,
                   V1.1 = vmat.df$v1.1, V2.2 = vmat.df$v2.2, D1=d1, D2=d2,
                   V3= vmat.df$v3, V4= vmat.df$v4, d_rd = vmat.df$d_rd,
                   v8=v8)

PO1d_env = newenv.var(sp.xy = PO1d_xy, env.grid = quads,
                      env.scale = 0.5, coord = c("X", "Y"), file.name = NA)

po1d_X = data.frame(Intercept = 1,PO1d_env[,9])
po1ext_beta = c(1.8,0.6) # extinction coef
po1Ext_intensity = 1/(1+exp(-as.matrix(po1d_X) %*% po1ext_beta))


# Species intensity is used to consider the more suitable habitat
prob_death = 0.65
P_retention1 = as.matrix(po1Ext_intensity)

breakloop = 0
i=1
while(breakloop == 0){
  set.seed(as.numeric(args[1])+i)
  notdead1 = sample(1:sp_sim$n, round(rnorm(1, 1-prob_death, 0.02)*sp_sim$n),
                    prob = P_retention1)
  
  PO1_yr1 = sp_sim[notdead1]
  
  ext.ratio = PO1_yr1$n/sp_sim$n
  ext.ratio
  
  i=i+1
  i
  
  if (ext.ratio < (1-0.64) & ext.ratio > (1-0.66))
  {
    breakloop = 1
  }
  breakloop
}



PO1_yr1


# birth-colonization

po1col_X = data.frame(Intercept = 1,PO1d_env[,10])
po1col_beta = c(1.8,0.5)
po1col_int = as.matrix(1/(1+exp(-as.matrix(po1col_X) %*% po1col_beta)))

prob_col = 0.65

breakloop = 0
i=1
while(breakloop == 0){
  set.seed(as.numeric(args[1])+i)
  OffS1 = sample(1:sp_sim$n, round(rnorm(1, prob_col, 0.05)*sp_sim$n), prob=po1col_int)
  Offsprg1 = sp_sim[OffS1]
  Offsprg1
  
  Col.ratio = Offsprg1$n/sp_sim$n
  Col.ratio
  
  i=i+1
  i
  
  if (Col.ratio > 0.64 & Col.ratio < 0.66)
  {
    breakloop = 1
  }
  breakloop
}

Offsprg1
POtrue_year1 = superimpose(PO1_yr1, Offsprg1)



## Species data Year 1

### PO data

int1 = -2.2

# Biased intensity for observer bias
po1_beta = c(int1, 0.5, 0.6, -0.5, -0.4, -1)
po1_intensity = exp(vmat[,c(1:5, 8)] %*% po1_beta)


# PO intensity surface
po1_int_im = as.im(data.frame(x = quad$X, y = quad$Y, z = po1_intensity))

# environmental info at quad points and species locations
quads = data.frame(X = quad$X, Y = quad$Y, V1 = vmat.df$v1, V2 = vmat.df$v2,
                   V1.1 = vmat.df$v1.1, V2.2 = vmat.df$v2.2, D1=d1, D2=d2,
                   V3= vmat.df$v3, V4= vmat.df$v4, d_rd = vmat.df$d_rd,
                   v8=v8)


sp_xy = data.frame(X = POtrue_year1$x, Y = POtrue_year1$y)

sp_env = newenv.var(sp.xy = sp_xy, env.grid = quads,
                    env.scale = 0.5, coord = c("X", "Y"), file.name = NA)

po1_X = data.frame(Intercept = 1, sp_env[,c(3:6, 11)])
po1_beta = c(int1, 0.5, 0.6, -0.5, -0.4, -1)
po1_intensity = exp(as.matrix(po1_X) %*% po1_beta)
PO1_rows = c() #vector indicating the sampled rows
X_add = sample(1:POtrue_year1$n, 300, prob = po1_intensity)
PO1_rows = c(PO1_rows, X_add)
PO1 = POtrue_year1[PO1_rows]
PO1


### Occ data

# *Generate occupancy data from species true intensity*
# Creating the grid sites
XY.occ =expand.grid(seq(1.5, 28.5, 3), seq(1.5, 28.5, 3))
X.occ = XY.occ[,1]
Y.occ = XY.occ[,2]

# Simulation a point pattern
sp_and_occ = ppp(x = c(X.occ, POtrue_year1$x), y = c(Y.occ, POtrue_year1$y),
                 marks = c(rep("Occ", length(X.occ)), rep("Sp", POtrue_year1$n)),
                 window = owin(c(-0.25, 30.25), c(-0.25, 30.25)))
dist_sp_occ = nndist(sp_and_occ, k = 1, by = as.factor(marks(sp_and_occ)))

# compute distance to nearest species
sp_occ_dists = dist_sp_occ[1:length(X.occ), 2]

# species considered present at site if nearest one is within 0.25 units
occ_present = as.numeric(sp_occ_dists <= 0.30)
table(occ_present) # we get about 2/3 of the sites occupied

occ_paste = paste(X.occ, Y.occ)
quad_paste = paste(quad$X, quad$Y)
occ_quadrow1 = match(occ_paste, quad_paste)


p_detect = clogloginv(v8[occ_quadrow1])*occ_present 

#p_detect = expit(v8[occ_quadrow1])*occ_present 

n_visits = 5
sim_history = matrix(as.integer(matrix(rep(p_detect, times = n_visits),
                                       length(X.occ), n_visits) >
                                  matrix(runif(length(X.occ)*n_visits),
                                         length(X.occ), n_visits)),
                     length(X.occ), n_visits)


# Occupancy plot
site_sum = apply(sim_history, 1, sum)
plot_xy =expand.grid(seq(0, 30, 0.5), seq(0, 30, 0.5))
plot_x = plot_xy[,1]
plot_y = plot_xy[,2]
plot_id = paste(plot_x, plot_y)
occ_id = paste(X.occ, Y.occ)
occ_match = match(occ_id, plot_id)
plot_z = rep(NA, length(plot_x))
plot_z[occ_match] = site_sum




## Species data Year 2


### PO data

## Year 1 introduce colonization/extinction using the covar v5 and v6
# Start from pattern from year1
PO2d_xy = data.frame(X = POtrue_year1$x, Y = POtrue_year1$y)

# environmental info at quad points and species locations
quads = data.frame(X = quad$X, Y = quad$Y, V1 = vmat.df$v1, V2 = vmat.df$v2,
                   V1.1 = vmat.df$v1.1, V2.2 = vmat.df$v2.2, D1=d1, D2=d2,
                   V3= vmat.df$v3, V4= vmat.df$v4, d_rd = vmat.df$d_rd,
                   v8=v8)

PO2d_env = newenv.var(sp.xy = PO2d_xy, env.grid = quads,
                      env.scale = 0.5, coord = c("X", "Y"), file.name = NA)

po2d_X = data.frame(Intercept = 1,PO2d_env[,9])
po2ext_beta = c(1.8,0.5) # extinction coef
po2Ext_intensity = 1/(1+exp(-as.matrix(po2d_X) %*% po2ext_beta))


# Species intensity is used to consider the more suitable habitat
prob_death2 = 0.65
P_retention2 = as.matrix(po2Ext_intensity)

breakloop2 = 0
i=1
while(breakloop2 == 0){
  set.seed(as.numeric(args[2])+i)
  notdead2 = sample(1:POtrue_year1$n, round(rnorm(1, 1-prob_death2, 0.02)*POtrue_year1$n),
                    prob = P_retention2)
  
  PO_yr2 = POtrue_year1[notdead2]
  
  ext.ratio2 = PO_yr2$n/POtrue_year1$n
  ext.ratio2
  
  i=i+1
  i
  
  if (ext.ratio2 < (1-0.64) & ext.ratio2 > (1-0.66))
  {
    breakloop2 = 1
  }
  breakloop2
}



PO_yr2


# birth-colonization

po2col_X = data.frame(Intercept = 1,PO2d_env[,10])
po2col_beta = c(1.8,0.6)
po2col_int = as.matrix(1/(1+exp(-as.matrix(po2col_X) %*% po2col_beta)))

prob_col2 = 0.65

breakloop2 = 0
i=1
while(breakloop2 == 0){
  set.seed(as.numeric(args[2])+i)
  OffS2 = sample(1:POtrue_year1$n, round(rnorm(1, prob_col2, 0.05)*POtrue_year1$n), prob=po2col_int)
  Offsprg2 = POtrue_year1[OffS2]
  Offsprg2
  
  Col.ratio2 = Offsprg2$n/POtrue_year1$n
  Col.ratio2
  
  i=i+1
  i
  
  if (Col.ratio2 > 0.64 & Col.ratio2 < 0.66)
  {
    breakloop2 = 1
  }
  breakloop2
}

Offsprg2
POtrue_year2 = superimpose(PO_yr2, Offsprg2)


int2 = -2.2
# Biased intensity for observer bias
po2_beta = c(int2, 0.5, 0.6, -0.5, -0.4, -1)
po2_intensity = exp(vmat[,c(1:5,8)] %*% po2_beta)


# PO intensity surface
po2_int_im = as.im(data.frame(x = quad$X, y = quad$Y, z = po2_intensity))

# environmental info at quad points and species locations
quads2 = data.frame(X = quad$X, Y = quad$Y, V1 = vmat.df$v1, V2 = vmat.df$v2,
                    V1.1 = vmat.df$v1.1, V2.2 = vmat.df$v2.2, D1=d1, D2=d2,
                    V3= vmat.df$v3, V4= vmat.df$v4, d_rd = vmat.df$d_rd,
                    v8=v8)

sp_xy2 = data.frame(X = POtrue_year2$x, Y = POtrue_year2$y)

sp_env2 = newenv.var(sp.xy = sp_xy2, env.grid = quads2,
                     env.scale = 0.5, coord = c("X", "Y"), file.name = NA)

po2_X = data.frame(Intercept = 1, sp_env2[,c(3:6, 11)])
po2_beta = c(int2, 0.5, 0.6, -0.5, -0.4, -1)
po2_intensity = exp(as.matrix(po2_X) %*% po2_beta)
PO2_rows = c() #vector indicating the sampled rows
X_add2 = sample(1:POtrue_year2$n, 300, prob = po2_intensity)
PO2_rows = c(PO2_rows, X_add2)
PO2 = POtrue_year2[PO2_rows]
PO2


### Occ data

# # Creating the grid sites
XY.occ =expand.grid(seq(1.5, 28.5, 3), seq(1.5, 28.5, 3))
X.occ = XY.occ[,1]
Y.occ = XY.occ[,2]

#for occupancy data

sp_and_occ2 = ppp(x = c(X.occ, POtrue_year2$x), y = c(Y.occ, POtrue_year2$y),
                  marks = c(rep("Occ", length(X.occ)), rep("Sp", POtrue_year2$n)),
                  window = owin(c(-0.25, 30.25), c(-0.25, 30.25)))
dist_sp_occ2 = nndist(sp_and_occ2, k = 1, by = as.factor(marks(sp_and_occ2)))

# compute distance to nearest species
sp_occ_dists2 = dist_sp_occ2[1:length(X.occ), 2]

# species considered present at site if nearest one is within 0.25 units
occ_present2 = as.numeric(sp_occ_dists2 <= 0.30)
table(occ_present2) 

occ_paste2 = paste(X.occ, Y.occ)
quad_paste2 = paste(quad$X, quad$Y)
occ_quadrow2 = match(occ_paste, quad_paste)


p_detect2 = clogloginv(v8[occ_quadrow2])*occ_present2
#p_detect2 = expit(v8[occ_quadrow2])*occ_present2

n_visits = 5
sim_history2 = matrix(as.integer(matrix(rep(p_detect2, times = n_visits),
                                        length(X.occ), n_visits) >
                                   matrix(runif(length(X.occ)*n_visits),
                                          length(X.occ), n_visits)),
                      length(X.occ), n_visits)


# Occupancy plot
site_sum2 = apply(sim_history2, 1, sum)
plot_xy =expand.grid(seq(0, 30, 0.5), seq(0, 30, 0.5))
plot_x = plot_xy[,1]
plot_y = plot_xy[,2]
plot_id = paste(plot_x, plot_y)
occ_id = paste(X.occ, Y.occ)
occ_match = match(occ_id, plot_id)
plot_z2 = rep(NA, length(plot_x))
plot_z2[occ_match] = site_sum2



## Species data Year 3


### PO data

## Year 1 introduce colonization/extinction using the covar v5 and v6
# Start from pattern from year1
PO3d_xy = data.frame(X = POtrue_year2$x, Y = POtrue_year2$y)

# environmental info at quad points and species locations
quads = data.frame(X = quad$X, Y = quad$Y, V1 = vmat.df$v1, V2 = vmat.df$v2,
                   V1.1 = vmat.df$v1.1, V2.2 = vmat.df$v2.2, D1=d1, D2=d2,
                   V3= vmat.df$v3, V4= vmat.df$v4, d_rd = vmat.df$d_rd,
                   v8=v8)

PO3d_env = newenv.var(sp.xy = PO3d_xy, env.grid = quads,
                      env.scale = 0.5, coord = c("X", "Y"), file.name = NA)

po3d_X = data.frame(Intercept = 1,PO3d_env[,9])
po3ext_beta = c(1.8,0.5) # extinction coef
po3Ext_intensity = 1/(1+exp(-as.matrix(po3d_X) %*% po3ext_beta))


# Species intensity is used to consider the more suitable habitat
prob_death3 = 0.65
P_retention3 = as.matrix(po3Ext_intensity)

breakloop3 = 0
i=1
while(breakloop3 == 0){
  set.seed(as.numeric(args[3])+i)
  notdead3 = sample(1:POtrue_year2$n, round(rnorm(1, 1-prob_death3, 0.02)*POtrue_year2$n),
                    prob = P_retention3)
  
  PO_yr3 = POtrue_year2[notdead3]
  
  ext.ratio3 = PO_yr3$n/POtrue_year2$n
  ext.ratio3
  
  i=i+1
  i
  
  if (ext.ratio3 < (1-0.64) & ext.ratio3 > (1-0.66))
  {
    breakloop3 = 1
  }
  breakloop3
}



PO_yr3


# birth-colonization

po3col_X = data.frame(Intercept = 1,PO3d_env[,10])
po3col_beta = c(1.8,0.6)
po3col_int = as.matrix(1/(1+exp(-as.matrix(po3col_X) %*% po3col_beta)))

prob_col3 = 0.65

breakloop3 = 0
i=1
while(breakloop3 == 0){
  set.seed(as.numeric(args[3])+i)
  OffS3 = sample(1:POtrue_year2$n, round(rnorm(1, prob_col2, 0.05)*POtrue_year2$n), prob=po3col_int)
  Offsprg3 = POtrue_year2[OffS3]
  Offsprg3
  
  Col.ratio3 = Offsprg3$n/POtrue_year2$n
  Col.ratio3
  
  i=i+1
  i
  
  if (Col.ratio3 > 0.64 & Col.ratio3 < 0.66)
  {
    breakloop3 = 1
  }
  breakloop3
}

Offsprg3
POtrue_year3 = superimpose(PO_yr3, Offsprg3)


int3 = -2.2
# Biased intensity for observer bias
po3_beta = c(int3, 0.5, 0.6, -0.5, -0.4, -1)
po3_intensity = exp(vmat[,c(1:5,8)] %*% po3_beta)


# PO intensity surface
po3_int_im = as.im(data.frame(x = quad$X, y = quad$Y, z = po3_intensity))

# environmental info at quad points and species locations
quads3 = data.frame(X = quad$X, Y = quad$Y, V1 = vmat.df$v1, V2 = vmat.df$v2,
                    V1.1 = vmat.df$v1.1, V2.2 = vmat.df$v2.2, D1=d1, D2=d2,
                    V3= vmat.df$v3, V4= vmat.df$v4, d_rd = vmat.df$d_rd,
                    v8=v8)

sp_xy3 = data.frame(X = POtrue_year3$x, Y = POtrue_year3$y)

sp_env3 = newenv.var(sp.xy = sp_xy3, env.grid = quads3,
                     env.scale = 0.5, coord = c("X", "Y"), file.name = NA)

po3_X = data.frame(Intercept = 1, sp_env3[,c(3:6, 11)])
po3_beta = c(int3, 0.5, 0.6, -0.5, -0.4, -1)
po3_intensity = exp(as.matrix(po3_X) %*% po3_beta)
PO3_rows = c() #vector indicating the sampled rows
X_add3 = sample(1:POtrue_year3$n, 300, prob = po3_intensity)
PO3_rows = c(PO3_rows, X_add3)
PO3 = POtrue_year3[PO3_rows]
PO3


### Occ data

# # Creating the grid sites
XY.occ =expand.grid(seq(1.5, 28.5, 3), seq(1.5, 28.5, 3))
X.occ = XY.occ[,1]
Y.occ = XY.occ[,2]

#for occupancy data

sp_and_occ3 = ppp(x = c(X.occ, POtrue_year3$x), y = c(Y.occ, POtrue_year3$y),
                  marks = c(rep("Occ", length(X.occ)), rep("Sp", POtrue_year3$n)),
                  window = owin(c(-0.25, 30.25), c(-0.25, 30.25)))
dist_sp_occ3 = nndist(sp_and_occ3, k = 1, by = as.factor(marks(sp_and_occ3)))

# compute distance to nearest species
sp_occ_dists3 = dist_sp_occ3[1:length(X.occ), 2]

# species considered present at site if nearest one is within 0.25 units
occ_present3 = as.numeric(sp_occ_dists3 <= 0.30)
table(occ_present3) 

occ_paste3 = paste(X.occ, Y.occ)
quad_paste3 = paste(quad$X, quad$Y)
occ_quadrow3 = match(occ_paste, quad_paste)


p_detect3 = clogloginv(v8[occ_quadrow3])*occ_present3

n_visits = 5
sim_history3 = matrix(as.integer(matrix(rep(p_detect3, times = n_visits),
                                        length(X.occ), n_visits) >
                                   matrix(runif(length(X.occ)*n_visits),
                                          length(X.occ), n_visits)),
                      length(X.occ), n_visits)


# Occupancy plot
site_sum3 = apply(sim_history3, 1, sum)
plot_xy =expand.grid(seq(0, 30, 0.5), seq(0, 30, 0.5))
plot_x = plot_xy[,1]
plot_y = plot_xy[,2]
plot_id = paste(plot_x, plot_y)
occ_id = paste(X.occ, Y.occ)
occ_match = match(occ_id, plot_id)
plot_z3 = rep(NA, length(plot_x))
plot_z3[occ_match] = site_sum3



# Fitting the model

## Set up some arguments

# year 1
po1_env = sp_env[PO1_rows,]
occ_env = quads[occ_quadrow1,]

# year 2
# get the rows
po2_env = sp_env2[PO2_rows,]
occ_env2 = quads2[occ_quadrow2,]

# year 3
# get the rows
po3_env = sp_env3[PO3_rows,]
occ_env3 = quads3[occ_quadrow3,]



## list for the years
env_y1 = list(po1_env, occ_env)
env_y2 = list(po2_env, occ_env2)  
env_y3 = list(po3_env, occ_env3) 

# quad information
is.data.frame(quads)
is.data.frame(quads2)

Quad.dyn = list(quads, quads2, quads3)

# Other elements
env_formula = ~ V1 + V2 + V1.1 + V2.2+ D1 + D2 + V3 + V4

bias_formula = list(list(~ d_rd, ~ v8), list(~ d_rd, ~ v8), list(~ d_rd, ~ v8))

quad_data = Quad.dyn
sp_data = list(env_y1, env_y2, env_y3)

sp_y = list(list(rep(1, nrow(po1_env)),sim_history), 
            list(rep(1, nrow(po2_env)), sim_history2),
            list(rep(1, nrow(po3_env)), sim_history3))
dat.type = c("PO", "Occ")
site.area = pi*0.18^2


n.time = length(sp_data)

# Function test #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

testcomb = comb_lasso(env_formula, bias_formula, intercept_env = list(1, 1,1), 
                      intercept_bias = list(NA, NA, NA), quad_data, sp_data, sp_y, 
                      dat.type, coord = c("X", "Y"), sp_res = 0.5, penalty_vec = NULL, 
                      alpha = 1, gamma = 0, init.coef = list(NA, NA, NA), 
                      standardise = TRUE, criterion = "BIC", family = "poisson", 
                      tol = 1.e-7, b.min = 1.e-6, max.it = 25, n.fits = 15, 
                      noshrink = rep(list(NULL), length(quad_data)), method = "BFGS", 
                      link = "logit", site.area = site.area, obs.scores = NULL, 
                      extinct_comp=c(8), coloni_comp=c(9), area.int = FALSE, 
                      r = NULL, wt.vec = NULL, pen.min = 1.e-6, verbose = TRUE)


testcomb$runtime
testcomb$penalty
testcomb$beta


