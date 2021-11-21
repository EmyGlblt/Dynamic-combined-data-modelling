# functions for generating covariates

trunc.Gauss.density = function(X, Mean, Sigma, trunc.sd = 3) # truncated Gaussian density
{
  X.dist    = spDists(as.matrix(X), t(as.matrix(Mean)), longlat = FALSE)
  stdev     = sqrt(Sigma[1, 1])       
  in.kernel = which(X.dist <= trunc.sd*stdev)
  denom     = pmvnorm(lower = Mean - trunc.sd*stdev, upper = Mean + trunc.sd*stdev, mean = Mean, sigma = Sigma)
  dens.out  = rep(0, length(X.dist))
  dens.out[in.kernel] = dmvnorm(x = X[in.kernel,], mean = Mean, sigma = Sigma)/denom
  dens.out
}

makecovar = function(n_patches, patch_centers = centers, quads = quad) # generate environmental covariate
{
  patches = sample(1:(dim(centers)[1]), n_patches)
  hab.centers = patch_centers[patches,]
  hab.qual = rep(0, dim(quad)[1])
  for (i in 1:n_patches)
  {
    hab.qual = hab.qual + trunc.Gauss.density(X = quads, Mean = hab.centers[i,], Sigma = diag(2)*(runif(1, 1.5, 2)^2), trunc.sd = runif(1, 2, 4))
  }
  levelplot(hab.qual ~ quads[,1] + quads[,2])
  hab.qual
}

# Road network functions

addroad = function(x0, y0, angle, sd = 3, dist = 0.1, n = 1000, xlim = c(0, 30), ylim = c(0, 30), level = 1, plot = FALSE) # add a road
{
  x.out = y.out = segangle = rep(NA, n)
  x.out[1] = x0
  y.out[1] = y0
  segangle[1] = angle
  
  for (seg in 2:n)
  {
    segangle[seg] = rnorm(1, mean = segangle[seg - 1], sd = sd)
    dx = dist*cos(segangle[seg]*pi/180)
    dy = dist*sin(segangle[seg]*pi/180)
    x.out[seg] = x.out[seg - 1] + dx
    y.out[seg] = y.out[seg - 1] + dy
  }
  if (plot != FALSE)
  {
    plot(x.out, y.out, type = "l", xlim = xlim, ylim = ylim)
  }
  outpoints = data.frame(x.out, y.out)
  out_of_window = which(x.out < xlim[1] | x.out > xlim[2] | y.out < ylim[1] | y.out > ylim[2])
  if (length(out_of_window) > 0)
  {
    outpoints = outpoints[1:(min(out_of_window) - 1),]
    segangle = segangle[1:(min(out_of_window) - 1)]
  }
  
  return(list(xy = outpoints, angles = segangle,
              sd = sd, dist = dist, n = length(segangle), xlim = xlim, ylim = ylim, level = level))
}

roadsplits = function(road, max_split, p_split, min_dist = 10) # split road into two roads
{
  proposed_gaps = rgeom(max_split, p_split)
  proposed_gaps[proposed_gaps < min_dist] = min_dist
  proposed_splits = cumsum(proposed_gaps)
  splits = proposed_splits[proposed_splits <= length(road$angles)]
  splits
}

addlevel = function(road, max_split, p_split, sdmult = 1, distmult = 1, nmult = 1, xlim = c(0, 30), ylim = c(0, 30), min_angle = 30) # add level of roads
{
  if (is.null(road$level) == FALSE) # single road
  {
    road = list(road)
  }
  levels = unlist(lapply(road, function(x) x$level))
  roadlist = road
  
  for (r in which(levels == max(levels)))
  {
    splits = roadsplits(road[[r]], max_split = max_split, p_split = p_split)
    if (length(splits) > 0)
    {
      for (i in 1:length(splits))
      {
        roadlist[[length(roadlist) + 1]] = addroad(x0 = road[[r]]$xy$x.out[splits[i]], 
                                                   y0 = road[[r]]$xy$y.out[splits[i]], 
                                                   angle = runif(1, min = road[[r]]$angles[splits[i]] + min_angle, max = road[[r]]$angles[splits[i]] + 180 - min_angle) + rbinom(1, 1, 0.5)*180,
                                                   sd = road[[r]]$sd*sdmult, dist = road[[r]]$dist*distmult, n = ceiling(road[[r]]$n*nmult),
                                                   level = road[[r]]$level + 1, xlim = xlim, ylim = ylim)
      }
    }
  }
  
  return(roadlist)
}

addnetwork = function(x0, y0, angle, sd = 3, dist = 0.1, n = 1000, p_split = 0.01, split_decay = 0.5, max_split = 10, max_levels = 4, sdmult = 1, distmult = 1, nmult = 1, xlim = c(0, 30), ylim = c(0, 30), min_angle = 30, plot = FALSE) # add network of roads
{
  road = addroad(x0, y0, angle = angle, sd = sd, dist = dist, n = n, level = 1)
  for (l in 1:(max_levels - 1))
  {
    road = addlevel(road, max_split = max_split, p_split = p_split*(split_decay)^(l - 1),
                    sdmult = sdmult, distmult = distmult,
                    nmult = nmult, xlim = xlim, ylim = ylim, min_angle = min_angle)
  }
  if (plot != FALSE)
  {
    plotroadlist(road)
  }
  return(road)
}

plotroadlist = function(roadlist, main = "", asp = 1) # plot list of roads
{
  n.roads = length(roadlist)
  plot(roadlist[[1]]$xy$x.out, roadlist[[1]]$xy$y.out, type = "l",
       xlab = "X", ylab = "Y", xlim = roadlist[[1]]$xlim, ylim = roadlist[[1]]$ylim,
       main = main, asp = asp)
  for (r in 2:n.roads)
  {
    points(roadlist[[r]]$xy$x.out, roadlist[[r]]$xy$y.out, type = "l", col = roadlist[[r]]$level)
  }
}

fixnodes = function(roadlist) # aggregate nodes
{
  newroadlist = roadlist
  allx = ally = c()
  for (i in 1:length(roadlist))
  {
    allx = c(allx, roadlist[[i]]$xy$x.out)
    ally = c(ally, roadlist[[i]]$xy$y.out)
  }
  n.nodes = unlist(lapply(roadlist, function(x) {length(x$angles)}))
  roadmarks = rep(1:length(roadlist), times = n.nodes)
  roadwindow = owin(xrange = roadlist[[1]]$xlim, yrange = roadlist[[1]]$ylim)
  roadppp = ppp(allx, ally, marks = roadmarks, window = roadwindow)
  nndists  = nndist(roadppp, by = as.factor(marks(roadppp)))
  dist = roadlist[[1]]$dist
  for (i in 1:length(roadlist))
  {
    droads = apply(nndists[roadmarks == i,], 2, min)
    closeroads = which(droads < dist & droads > 0)
    closeroads = closeroads[closeroads > i]
    if (length(closeroads > 0))
    {
      for (j in closeroads)
      {
        closedist = min(nndists[roadmarks == i, j])
        closenode_i = which(nndists[roadmarks == i, j] == closedist)
        closenode_j = which(nndists[roadmarks == j, i] == closedist)
        newroadlist[[i]]$xy[closenode_i,] = roadlist[[j]]$xy[closenode_j,]
      }
    }
  }
  newroadlist
}

nodelist = function(roadlist) # create list of road network junction nodes
{
  allx = ally = c()
  for (i in 1:length(roadlist))
  {
    allx = c(allx, roadlist[[i]]$xy$x.out)
    ally = c(ally, roadlist[[i]]$xy$y.out)
  }
  allxy = paste(allx, ally)
  n.nodes = unlist(lapply(roadlist, function(x) {length(x$angles)}))
  r.levels = unlist(lapply(roadlist, function(x) {x$level}))
  roadmarks = rep(1:length(roadlist), times = n.nodes)
  nodelevels = rep(r.levels, times = n.nodes)
  max_neighbors = max(table(allxy))*2
  neighbors = matrix(NA, length(roadmarks), max_neighbors)
  
  for (node in 1:length(roadmarks))
  {
    noderoad = roadmarks[node]
    neighbor_candidates = which(roadmarks == noderoad)
    node_neighbors = neighbor_candidates[which(abs(neighbor_candidates - node) == 1)]
    neighbors[node, 1:length(node_neighbors)] = node_neighbors
  }
  dupenodes = duplicated(allxy)
  n_unique_nodes = sum(dupenodes == FALSE)
  uniqueneighbors = data.frame(Node = 1:n_unique_nodes, Level = NA, X = allx[dupenodes == FALSE], Y = ally[dupenodes == FALSE], matrix(NA, n_unique_nodes, max_neighbors))
  allxysub = paste(uniqueneighbors$X, uniqueneighbors$Y)
  
  for (node in 1:n_unique_nodes)
  {
    samenodes = which(allxy == allxysub[node])
    oldneighbors = sort(as.vector(neighbors[samenodes,]))
    uniqueneighbors[node, 5:(length(oldneighbors) + 4)] = match(allxy[oldneighbors], allxysub)
    uniqueneighbors[node, 2] = min(nodelevels[samenodes])
  }
  names(uniqueneighbors)[5:length(uniqueneighbors)] = paste("N", 1:max_neighbors, sep = "")
  uniqueneighbors
}

# combined penalised likelihood functions

comb_lasso = function(env_formula, bias_formula = NULL, intercept_env = NULL, intercept_bias = NULL, 
                          quad_data = NULL, sp_data, sp_y, dat.type = "PO", coord = c("X", "Y"), sp_res = 1, 
                          penalty_vec = NULL, alpha = 1, gamma = 0, init.coef = NA, standardise = TRUE, criterion = "BIC",
                          family = "poisson", tol = 1.e-7, b.min = 1.e-6, 
                          max.it = 25, n.fits = 100, noshrink = NULL, method = "BFGS", 
                          link = "logit", site.area = 1,
                          area.int = FALSE, r = NULL, wt.vec = NULL, pen.min = 1.e-2, verbose = FALSE)
{
  
  # To do:
  # Automatic choice of penalty vector
  # output standardisation means and standard deviations
  # block CV
  
  start.time = Sys.time()
  formula.out = list(env_formula, bias_formula)
  dat.use = dataprep(env_formula = env_formula, bias_formula = bias_formula, 
                     intercept_env = intercept_env, intercept_bias = intercept_bias, 
                     quad_data = quad_data, sp_data = sp_data, 
                     sp_y = sp_y, coord = coord, sp.res = sp_res, dat.type = dat.type, 
                     standardise = standardise, area.int = area.int, r = r)
  
  noshrink   = c("Intercept", noshrink)
  all_names  = c(colnames(dat.use$X_env[[1]]), unlist(lapply(dat.use$X_bias, colnames)))
  
  # Delete below 2 lines!
  #is.intercept = which(all_names == "Intercept")
  #all_names[is.intercept[2:length(is.intercept)]] = "1"
  
  noshrink   = which(all_names %in% noshrink)
  
  n_components = dat.use$n_components
  
  if (is.null(wt.vec))
  {
    wt.vec = rep(1, n_components)
  }
  
  env_p     = dim(dat.use$X_env[[1]])[2]
  bias_p    = unlist(lapply(dat.use$X_bias, function(x) dim(x)[[2]]))
  num_p     = env_p + sum(bias_p)
  bias_ind1 = env_p + 1 + cumsum(bias_p) - bias_p
  bias_ind2 = env_p + cumsum(bias_p)
  env_ix    = rep(list(1:env_p), n_components)
  bias_ix   = list()
  for (i in 1:n_components)
  {
    bias_ix[[i]] = bias_ind1[i]:bias_ind2[i]
  }
  
  if (is.na(init.coef)[1])
  {
    gamma = 0
  }
  adapt.weights = if (gamma == 0) rep(1, num_p) else 1/abs(init.coef)^gamma
  
  if (is.null(penalty_vec))
  {
    intmodel  = score.int_new(dat.use$X_env, dat.use$X_bias, dat.use$y_list, dat.use$ob_wt, dat.type = dat.type, site.area = site.area, lasso = 0, env_ix = env_ix, bias_ix = bias_ix, coord = coord, wt.vec = wt.vec, noshrink = noshrink)
    score0    = intmodel$score_int
    b.init    = intmodel$newpar
    
    new.score  = abs(score0/adapt.weights)
    sub.score  = new.score[-noshrink]
    pen.max    = max(sub.score[is.infinite(sub.score) == FALSE])
    if (is.na(pen.max) == FALSE)
    {
      penalty_vec = c(0, exp(seq(log(pen.min), log(pen.max + 1.e-5), length.out = (n.fits - 1))))
    }
  }
  
  if (is.na(init.coef)[1] == FALSE)
  {
    init.coef[init.coef == 0] = b.min
  }
  
  fit.0 = combsingle.lasso_new(dat.use$y_list, dat.use$X_env, dat.use$X_bias, lamb = penalty_vec[2], ob.wt = dat.use$ob_wt, 
                               alpha = alpha, b.init = b.init, intercept = NA, family = family, 
                               tol = tol, gamma = gamma, init.coef = init.coef, mu.min = 1.e-16, mu.max = 1/mu.min, area.int = area.int, 
                               interactions, max.it = 100, standardise = standardise, 
                               noshrink = noshrink, method = method, b.min = b.min, 
                               link = link, site.area = site.area, dat.type = dat.type, wt.vec = wt.vec, grad = FALSE)
  
  b.init = fit.0$b
  
  fitpath  = matrix(NA, num_p, length(penalty_vec))
  likepath = rep(NA, length(penalty_vec))
  
  for (i in 1:length(penalty_vec))
  {
    it.max = if (i == 1) 100 else max.it
    
    #function(y, Xenv, Xbias = NULL, lamb, ob.wt = rep(1, length(y)), alpha = 1, b.init = NA, intercept = NA, family = "gaussian", tol = 1.e-9, gamma = 0, init.coef = NA, mu.min = 1.e-16, mu.max = 1/mu.min, area.int = FALSE, interactions, max.it = 25, standardise = TRUE, noshrink = c(1), method = "BFGS", b.min = 1.e-4, link = "logit", site.area = 1, dat.type = "PO", wt.vec = NULL)
    
    fit.i = combsingle.lasso_new(dat.use$y_list, dat.use$X_env, dat.use$X_bias, lamb = penalty_vec[i], ob.wt = dat.use$ob_wt, 
                                 alpha = alpha, b.init = b.init, intercept = NA, family = family, 
                                 tol = tol, gamma = gamma, init.coef = init.coef, mu.min = 1.e-16, mu.max = 1/mu.min, area.int = area.int, 
                                 interactions, max.it = it.max, standardise = standardise, 
                                 noshrink = noshrink, method = method, b.min = b.min, 
                                 link = link, site.area = site.area, dat.type = dat.type, wt.vec = wt.vec)
    fitpath[,i] = fit.i$b
    likepath[i] = fit.i$pen.likelihood
    b.init      = fit.i$b
    if (verbose == TRUE)
    {
      cat(paste(i, "\n"))
      flush.console()
    }
  }
  
  rownames(fitpath) = all_names
  
  variable_set = abs(fitpath) >= b.min
  n_variables  = apply(variable_set, 2, sum)
  
  AICs = 2*likepath + 2*n_variables
  BICs = 2*likepath + log(sum(dat.use$y_n))*n_variables
  HQCs = 2*likepath + 2*(n_variables + 1)*log(log(sum(dat.use$y_n)))
  AICcs = AICs + 2*(n_variables + 1)*(n_variables + 2)/(sum(dat.use$y_n) - n_variables - 2)
  
  criterion.matrix        = data.frame(AICs, AICcs, BICs, HQCs)
  names(criterion.matrix) = c("AIC", "AICc", "BIC", "HQC")
  
  meth.id   = paste(criterion, "s", sep = "")
  choice.id = max(which.min(get(meth.id)))
  
  beta.hat   = fitpath[,choice.id]
  like.hat   = likepath[choice.id]
  penalty    = penalty_vec[choice.id]
  
  mu   = list()
  bias = list()
  psi  = list()
  p_detect = list()
  
  for (i in 1:n_components)
  {
    mu[[i]] = exp(dat.use$X_env[[i]] %*% beta.hat[env_ix[[i]]])
    psi[[i]] = 1 - exp(-mu[[i]]*site.area)
    if (dat.type[i] == "PO")
    {
      bias[[i]] = exp(dat.use$X_bias[[i]] %*% beta.hat[bias_ix[[i]]])
      p_detect[[i]] = list(NULL)
    }
    else
    {
      bias[[i]] = list(NULL)
      lin_det  = dat.use$X_bias[[i]] %*% beta.hat[bias_ix[[i]]]
      if (link == "logit")
      {
        p_detect[[i]] = expit(lin_det)
      }
      if (link == "cloglog")
      {
        p_detect[[i]] = clogloginv(lin_det)
        #p_detect[[i]] = 1 - exp(-exp(lin_det))
      }
    }
  }
  end.time = Sys.time()
  runtime = end.time - start.time
  return(list(formula = formula.out, betas = fitpath, pen_likelihoods = likepath, penalty_vec = penalty_vec, 
              criterion.matrix = criterion.matrix, beta = beta.hat, pen_likelihood = like.hat, 
              penalty = penalty, mu = mu, bias = bias, psi = psi, p_detect = p_detect, 
              X_env = dat.use$X_env, X_bias = dat.use$X_bias, ob_wt = dat.use$ob_wt, y = dat.use$y_list, 
              quad_data = quad_data, sp_data = sp_data, dat.type = dat.type,
              quad_xy = dat.use$quad_xy, sp_xy = dat.use$sp_xy, site.area = site.area, runtime = runtime))
}

dataprep = function(env_formula, bias_formula, intercept_env, intercept_bias, quad_data, sp_data, sp_y, coord, sp.res, dat.type, standardise, area.int, r)
{
  ob_wt = list()
  y_list = list()
  if (class(sp_data) != "list")
  {
    sp_data = list(sp_data)
  }
  sp_data_xycol = lapply(lapply(sp_data, names), function(x) match(coord, x))
  n_components  = length(sp_data)
  sp_data_col =  sapply(lapply(sp_data, dim), function(x) x[2]) 
  y_n         = sapply(lapply(sp_data, dim), function(x) x[1])
  if (is.null(quad_data) == FALSE)
  {
    quad_coord = match(coord, names(quad_data))
    quad_vars  = names(quad_data[-quad_coord])
    for (i in 1:n_components)
    {
      need.interp = setdiff(quad_vars, names(sp_data[[i]]))
      if (length(need.interp) > 0) # change this to interpolate only needed variables?
      {
        sp_data[[i]] = newenv.var(sp.xy = sp_data[[i]], env.scale = sp.res, env.grid = quad_data, coord = coord, file.name = "SpEnvData")
      }
    }
  }
  quad_data = sample.quad(quad_data, sp.res)
  for (i in 1:n_components)
  {
    if (dat.type[i] == "PO")
    {
      ob_wt[[i]] = weights(sp_data[[i]], quad_data)
      pres.i     = rep(0, length(ob_wt[[i]]))
      pres.i[1:length(sp_y[[i]])] = sp_y[[i]]
      y_list[[i]] = pres.i/ob_wt[[i]]
    }
    else
    {
      ob_wt[[i]] = list(NULL)
      y_list[[i]] = sp_y[[i]]
    }
  }
  
  X_env = list()
  X_bias = list()
  
  #env_all = rbind(do.call("rbind", sp_data), quad_data)
  env_all = if (any(dat.type == "PO") == TRUE) rbind(do.call("rbind", sp_data), quad_data) else do.call("rbind", sp_data)
  row_start = c(1, cumsum(y_n) + 1)[1:n_components]
  row_end   = cumsum(y_n)
  sp_rows   = if (n_components == 1) list(row_start:row_end) else mapply(function(x, y) x:y, row_start, row_end)
  #quad_rows = (max(row_end) + 1):(dim(env_all)[1])
  quad_rows = if (any(dat.type == "PO") == TRUE) (max(row_end) + 1):(dim(env_all)[1]) else NULL
  for (i in 1:n_components)
  {
    if (dat.type[i] == "PO")
    {
      sp_rows[[i]] = c(sp_rows[[i]], quad_rows)
    }
  }
  
  mf_env     = model.frame(env_formula, env_all)
  mt_env     = attr(mf_env, "terms")
  attr(mt_env, "intercept") = 0
  X_env_all = if (!is.empty.model(mt_env)) model.matrix(mt_env, mf_env, contrasts)else matrix(, NROW(X_env_all), 0L)
  
  if (is.null(bias_formula) == FALSE)
  {
    for (i in 1:n_components)
    {
      if (dat.type[i] == "PO")
      {
        data.i      = rbind(sp_data[[i]], quad_data)
        bf.i = model.frame(bias_formula[[i]], data.i)
        bt.i = attr(bf.i, "terms")
        attr(bt.i, "intercept") = 0
        X_bias[[i]] = if (!is.empty.model(bt.i)) model.matrix(bt.i, bf.i, contrasts)
        else matrix(, NROW(sp_data[[i]]), 0L)
      }
      else
      {
        bf.i = model.frame(bias_formula[[i]], sp_data[[i]])
        bt.i = attr(bf.i, "terms")
        attr(bt.i, "intercept") = 0
        X_bias[[i]] = if (!is.empty.model(bt.i)) model.matrix(bt.i, bf.i, contrasts)
        else matrix(, NROW(sp_data[[i]]), 0L)
      }
    }
  }
  
  s_means = s_sds = NULL
  
  if (standardise == TRUE)
  {
    X_env_all_stand = standardise.X(X_env_all)
    s_means_env = X_env_all_stand$dat.means
    s_sds_env   = X_env_all_stand$dat.sds
    X_env_s     = X_env_all_stand$X
    X_env = lapply(sp_rows, function(x) X_env_s[x,])
    
    s_means_bias = s_sds_bias = c()
    if (is.null(bias_formula) == FALSE)
    {
      for (i in 1:n_components)
      {
        X_bias_stand_i = standardise.X(X_bias[[i]])
        s_means_bias   = c(s_means_bias, X_bias_stand_i$dat.means)
        s_sds_bias     = c(s_sds_bias, X_bias_stand_i$dat.sds)
        X_bias[[i]]    = X_bias_stand_i$X
      }
    }
    s_means = c(s_means_env, s_means_bias)
    s_sds   = c(s_sds_env, s_sds_bias)
  }
  
  
  # add intercepts
  
  #intercept_env = 1
  #intercept_bias = list(NA, NA, 1)
  
  if (is.null(intercept_env) == FALSE)
  {
    X_env = lapply(X_env, function(x) cbind(Intercept = 1, x))
  }
  
  if (is.null(intercept_bias) == FALSE)
  {
    intercept_add = which(is.na(intercept_bias) == FALSE)
    for (add.i in intercept_add)
    {
      X_bias[[add.i]] = cbind(Intercept = 1, X_bias[[add.i]])
    }
  }
  
  X_env = lapply(X_env, polynamesadd)
  X_bias = lapply(X_bias, polynamesadd)
  
  if (area.int == TRUE)
  {
    for (i in 1:n_components)
    {
      if (dat.type[i] == "PO" & is.na(r[i]) == FALSE)
      {
        int_i = point.interactions(env = quad_data, pres = sp_data[[i]], r = r[i], coord = coord)
        if (standardise == TRUE)
        {
          int_i = scale(int_i)
        }
        X_bias[[i]] = cbind(X_bias[[i]], Interaction = int_i)
        dimnames(X_bias[[i]])[[2]][dim(X_bias[[i]])[2]] = "Interaction"
      }
    }
  }
  sp_xy = lapply(sp_data, function(x) {x[,match(coord, names(x))]})
  quad_xy = quad_data[,match(coord, names(quad_data))]
  return(list(quad_data = quad_data, sp_data = sp_data, y_list = y_list, ob_wt = ob_wt, y_n = y_n, n_components = n_components, X_env = X_env, X_bias = X_bias, quad_xy = quad_xy, sp_xy = sp_xy))
}

m_to_km = function(data, coord = c("X", "Y"))
{
  data_km = data
  xy_col  = match(coord, names(data_km))
  data_km[,xy_col] = data[,xy_col]/1000
  data_km
}

score.int_new = function(Xenv, Xbias = NULL, y, ob.wt = NULL, link = "logit", sp_res = 1, site.area = 1, dat.type = "PO", lasso = 0, env_ix, bias_ix, coord = c("x", "y"), wt.vec, noshrink, method = "BFGS", b.min = 1.e-5)
{
  if (class(y) == "list"){
    num_env  = dim(Xenv[[1]])[2]
    num_bias = unlist(lapply(lapply(Xbias, as.matrix), ncol))
    num_p    = num_env + sum(num_bias)
  }else{
    num_p = dim(Xenv)[2] + dim(Xbias)[2]
  }
  
  par = rep(0, num_p)
  beta0 = rep(0, length(y))
  
  for (i in 1:length(y))
  {
    if (dat.type[i] == "PO")
    {
      beta0[i] = sum(y[[i]] > 0)/(sum(y[[i]] == 0)*sp_res^2)
    }
    if (dat.type[i] == "Occ")
    {
      beta0[i] = sum(apply(y[[i]], 1, sum) > 0)/(dim(y[[i]])[1]*site.area)
    }
  }
  
  par[noshrink[which(noshrink %in% env_ix[[1]])]] = log(mean(beta0))
  
  intercept_model = optim(par = par, fn = LL_Lasso_Wt_Comb_new, #gr = scoreeq_new,
                          Xenv = Xenv, 
                          Xbias = Xbias, 
                          y = y, ob.wt = ob.wt, link = link, coord = coord,
                          env_ix = env_ix,
                          bias_ix = bias_ix,
                          site.area = site.area, 
                          lasso = rep(0, num_p), wt.vec = wt.vec, dat.type = dat.type, is.in = 1:num_p %in% noshrink,
                          method = method, control = list(ndeps = rep(b.min, num_p))) # betas
  
  # July 17 2018: behaves strangely with 0 as par for occ lasso with large site area. 
  # Better starting point than 0 (such as log(31/(73*site.area))) seems to work, but how to generalize?
  # Maybe an average intensity from PO plus occ?
  
  score_int = scoreeq_new(Xenv, Xbias = Xbias, y = y, par = intercept_model$par, ob.wt = ob.wt, 
                          link = link, site.area = site.area, dat.type = dat.type, 
                          lasso = 0, env_ix = env_ix, bias_ix = bias_ix, coord = coord, 
                          wt.vec = wt.vec, is.in = rep(TRUE, num_p))
  return(list(score_int = score_int, newpar = intercept_model$par))
}

combsingle.lasso_new = function(y, Xenv, Xbias = NULL, lamb, ob.wt = rep(1, length(y)), alpha = 1, b.init = NA, intercept = NA, family = "gaussian", tol = 1.e-9, gamma = 0, init.coef = NA, mu.min = 1.e-16, mu.max = 1/mu.min, area.int = FALSE, interactions, max.it = 25, standardise = TRUE, noshrink = c(1), method = "BFGS", b.min = 1.e-4, link = "logit", site.area = 1, dat.type = "PO", wt.vec = NULL, grad = TRUE)
{
  # To do:
  # Fix behvaior with AI models
  # Other estimation procedures (IRLS, etc.?)
  
  if (class(y) == "list") # determine individual data sources and corresponding indices
  {
    num_env  = dim(Xenv[[1]])[2]
    num_bias = unlist(lapply(lapply(Xbias, as.matrix), ncol))
    num_p    = num_env + sum(num_bias)
    n.components = length(num_bias)
    bias_ind1 = num_env + 1 + cumsum(num_bias) - num_bias
    bias_ind2 = num_env + cumsum(num_bias)
    ix = list()
    env_ix = rep(list(1:num_env), n.components)
    bias_ix = list()
    for (i in 1:n.components)
    {
      ix[[i]] = c(1:num_env, bias_ind1[i]:bias_ind2[i])
      bias_ix[[i]] = bias_ind1[i]:bias_ind2[i]
    }
    X = mapply(cbind, Xenv, Xbias, SIMPLIFY = FALSE) 
  }
  else
  {
    env_cov  = 1:dim(Xenv)[2]
    bias_cov = if (is.null(Xbias) == FALSE) (max(env_cov) + 1):(max(env_cov) + dim(Xbias)[2])
    else c()
    
    #bias_cov = (max(env_cov) + 1):(max(env_cov) + dim(Xbias)[2])
    
    X = cbind(Xenv, Xbias)
    num_p = length(c(env_cov, bias_cov))
  }
  
  if (is.null(wt.vec))
  {
    wt.vec = rep(1, n.components)
  }
  
  error.flag = FALSE
  
  # eta.min = family$linkfun(mu.min)
  # eta.max = family$linkfun(mu.max)
  
  adapt.weights = (if(gamma == 0) rep(1, num_p) # check area int
                   else 1/abs(init.coef)^gamma)
  
  #b.glm   = rep(NA, (dim(X)[2] + 1))
  b.glm   = rep(NA, num_p)
  b.lasso = b.glm
  
  lambda.start = (if(length(lamb) == 1) rep(lamb, num_p)
                  else rep(0, num_p))
  lambda.start[noshrink] = 0
  lambda = as.array(lambda.start)
  lambda = lambda * abs(adapt.weights[1:length(lambda)])
  
  # if (area.int == TRUE)
  # {
  #   X      = cbind(X, interactions)
  #   is.in  = rep(TRUE, dim(X)[2])
  # }
  
  killed  = is.infinite(lambda)
  keep = which(killed == FALSE)
  
  ix.keep = lapply(ix, function(x) killed[x] == FALSE)
  #env.keep = lapply(env_ix, function(x) x %in% keep)
  #bias.keep = lapply(bias_ix, function(x) x %in% keep)
  
  Xenv  = subsetlist(Xenv, env_ix, keep)
  Xbias = subsetlist(Xbias, bias_ix, keep)
  
  
  #X       = X[,killed == FALSE]
  b.lasso = b.lasso[killed == FALSE]
  lambda  = lambda[killed == FALSE]
  
  # if (length(lamb) == 1 & lamb[1] == 0)
  # {
  #   b.init    = NA
  #   intercept = NA
  # }
  
  if (any(is.na(b.init)) & is.na(intercept))
  {
    
    # b.est.old = optim(par = rep(0, sum(killed == FALSE)), fn = LL_Lasso_Wt_Comb,
    #               Xenv = Xenv, Xbias = Xbias, y = y, ob.wt = ob.wt,
    #               env_ix = env_ix, bias_ix = bias_ix, site.area = site.area,
    #               lasso = lambda[killed == FALSE], wt.vec = wt.vec, dat.type = dat.type, is.in = killed == FALSE,
    #               method = method, control = list(ndeps = rep(b.min, sum(killed == FALSE))))
    
    b.est = optim(par = rep(0, sum(killed == FALSE)), fn = LL_Lasso_Wt_Comb_new, gr = if (grad == TRUE) scoreeq_new else NULL,
                  Xenv = Xenv, 
                  Xbias = Xbias, 
                  y = y, ob.wt = ob.wt, link = link, coord = coord,
                  env_ix = env_ix,
                  bias_ix = bias_ix,
                  site.area = site.area, 
                  lasso = lambda[killed == FALSE], wt.vec = wt.vec, dat.type = dat.type, is.in = killed == FALSE,
                  method = method, control = list(ndeps = rep(b.min, sum(killed == FALSE))))
    
    b.lasso = b.est$par
    like1   = b.est$value
    
    
    # mu.old  = rep(mean(y), dim(X)[1])
    # eta.old = eta.from.mu(mu.old, family = family, mu.min, mu.max, eta.min, eta.max)
    # l.old   = like.calc(X, family, ob.wt, mu = mu.old, y, alpha, lambda, beta = rep(0, dim(X)[2]), penalty = FALSE)
    # b.old   = rep(1, dim(X)[2])
    # diff    = 10
    # 
    # while (diff > tol)
    # {
    #   update   = irls.update(y, X, ob.wt, is.in, signs, eta = eta.old, mu = mu.old, alpha = alpha, lambda = lambda, beta.old = b.old, penalty = FALSE, family = family, mu.min, mu.max, eta.min, eta.max, tol)
    #   l.new    = like.calc(X, family, ob.wt, mu = update$mu, y, alpha, lambda, beta = update$beta, penalty = FALSE)
    #   diff     = abs(l.new - l.old)
    #   mu.old   = update$mu
    #   eta.old  = update$eta
    #   b.old    = update$beta
    #   l.old    = l.new
    # }
    # 
    # b.glm   = update$beta
    # b.lasso = b.glm
  }
  
  if (any(is.na(b.init)) == FALSE)
  {
    b.lasso = b.init[killed == FALSE]
    # like1   = LL_Lasso_Wt_Comb(par = b.lasso, Xenv = Xenv, Xbias = Xbias, y = y, ob.wt = ob.wt,
    #                                env_ix = env_ix, bias_ix = bias_ix, site.area = site.area, 
    #                                lasso = lambda[killed == FALSE], wt.vec = wt.vec, dat.type = dat.type)
    
    like1   = LL_Lasso_Wt_Comb_new(par = b.lasso, Xenv = Xenv, Xbias = Xbias, y = y, ob.wt = ob.wt,
                                   env_ix = env_ix, bias_ix = bias_ix, site.area = site.area, 
                                   lasso = lambda[killed == FALSE], wt.vec = wt.vec, dat.type = dat.type, is.in = killed == FALSE)
    
    # Need to calculate likelihood here too
  }
  
  if (is.na(intercept) == FALSE)
  {
    if (family$family == "poisson")
    {
      #b.lasso = c(log(mean(y)), rep(0, (dim(X)[2] - 1 - area.int)), rep(0, area.int))
      b.lasso = c(log(mean(y)), rep(0, (dim(X)[2])))
    }
    
    if (family$family == "binomial")
    {
      b.lasso = c(log(mean(y)/(1 - mean(y))), rep(0, (dim(X)[2] - 1 - area.int)), rep(1, area.int))
    }
    
    if (family$family == "gaussian")
    {
      b.lasso = c(mean(y), rep(0, (dim(X)[2] - 1 - area.int)), rep(1, area.int))
    }
  }
  
  is.in        = abs(b.lasso) > b.min
  is.in[noshrink] = TRUE
  b.lasso[is.in == FALSE] = 0
  sign.change  = rep(-1, num_p)
  is.different = is.in
  
  diff   = 2*tol
  num.it = 0
  
  signs = sign(b.lasso)
  
  betas   = c(b.lasso)
  scores  = c()
  viols   = c()
  likes   = c()
  actions = c()
  varsets = c()
  
  
  
  # eta.old = X %*% b.lasso
  # eta.old[eta.old < eta.min] = eta.min
  # eta.old[eta.old > eta.max] = eta.max
  # mu.old  = mu.from.eta(eta.old, family = family, mu.min, mu.max)
  # 
  # like    = like.calc(X, family, ob.wt, mu = mu.old, y, alpha, lambda, beta = b.lasso, penalty = TRUE)
  
  likes   = c(likes, like1)
  
  # if (length(lamb) == 1 & lamb[1] == 0)
  # {	
  #   diff   = 0
  #   num.it = max.it + 1
  # }
  
  while(diff > tol & num.it < max.it)
  {
    setchange = 0
    
    # b.est.old = optim(par = b.lasso, fn = LL_Lasso_Wt_Comb,
    #               Xenv = Xenv,
    #               Xbias = Xbias,
    #               y = y, ob.wt = ob.wt,
    #               env_ix = env_ix,
    #               bias_ix = bias_ix,
    #               site.area = site.area,
    #               lasso = lambda, wt.vec = wt.vec, dat.type = dat.type, is.in = is.in,
    #               method = method, control = list(ndeps = rep(b.min, num_p), trace = 6, REPORT = 1))
    
    #function(Xenv, Xbias = NULL, y, beta, ob.wt = NULL, link = "logit", site.area = 1, dat.type = "PO", lasso = 0)
    
    #
    
    b.est = optim(par = b.lasso, fn = LL_Lasso_Wt_Comb_new, gr = if (grad == TRUE) scoreeq_new else NULL, 
                  Xenv = Xenv, 
                  Xbias = Xbias, 
                  y = y, ob.wt = ob.wt, link = link, coord = coord,
                  env_ix = env_ix,
                  bias_ix = bias_ix,
                  site.area = site.area, 
                  lasso = lambda, wt.vec = wt.vec, dat.type = dat.type, is.in = is.in,
                  method = method, control = list(ndeps = rep(b.min, num_p)))
    
    # b.est = optim(par = b.lasso[is.in], fn = LL_Lasso_Wt_Comb, 
    #               Xenv = subsetlist(Xenv, env_ix, keep = which(is.in == TRUE)), 
    #               Xbias = subsetlist(Xbias, bias_ix, keep = which(is.in == TRUE)), 
    #               y = y, ob.wt = ob.wt,
    #               env_ix = subsetix(env_ix, keep = which(is.in == TRUE)),
    #               bias_ix = subsetix(bias_ix, keep = which(is.in == TRUE)),
    #               site.area = site.area, 
    #               lasso = lambda[is.in], wt.vec = wt.vec, dat.type = dat.type,
    #               method = method, control = list(ndeps = rep(b.min, sum(is.in))))
    
    
    #sign.change                   = signs[is.in]*sign(b.est$par)
    sign.change                   = signs*sign(b.est$par)
    #sign.change[cumsum(is.in)[noshrink]] = 1
    sign.change[sign.change == 0] = 1
    sign.change[abs(b.est$par) < b.min] = -1
    sign.change[is.in == FALSE]   = 1
    sign.change[noshrink] = 1
    score.beta                    = b.est$par
    
    #score.beta                    = rep(0, length(b.lasso))
    #score.beta[is.in]             = b.est$par
    
    if (any(sign.change != 1) == TRUE)
    {
      #delta                       = b.est$par - b.lasso[is.in]
      delta                       = b.est$par - b.lasso
      tozero                      = abs(b.lasso[sign.change != 1]) - b.min
      #tozero                      = abs(b.lasso[is.in][sign.change != 1]) - b.min
      prop                        = min(tozero/abs(delta[sign.change != 1]))
      
      #prop                        = min(abs(b.lasso[is.in]/delta)[sign.change != 1])
      b.lasso                     = b.lasso + prop*delta
      #b.lasso[is.in]              = b.lasso[is.in] + prop*delta
      
      is.in[abs(b.lasso) <= b.min + b.min/(1e10)] = FALSE
      b.lasso[is.in == FALSE]     = 0
      signs                       = sign(b.lasso)
      setchange                   = 2*tol
      score.beta                  = b.lasso
      
      # update$eta                        = as.matrix(X[,is.in]) %*% b.lasso[is.in]
      # update$eta[update$eta < eta.min]  = eta.min
      # update$eta[update$eta > eta.max]  = eta.max
      # update$mu                         = mu.from.eta(update$eta, family = family, mu.min, mu.max)
    }
    
    #function(Xenv, Xbias = NULL, y, par, ob.wt = NULL, link = "logit", site.area = 1, dat.type = "PO", lasso = 0, env_ix, bias_ix, coord = c("x", "y"), wt.vec, is.in = rep(TRUE, length(par)))
    score                   = scoreeq_new(Xenv, Xbias, y, b.lasso, ob.wt, link = link, site.area = site.area, dat.type = dat.type, lasso = 0, env_ix, bias_ix, coord, wt.vec, is.in = rep(TRUE, num_p))
    
    #score                   = ppmscoreeq(y, X, b.lasso, ob.wt)
    #score                   = t(as.vector(update$deriv) * t(update$s.wt)) %*% (y - update$mu)
    score.lamb              = alpha*lambda + (1 - alpha) * lambda * as.vector(score.beta)
    score.lamb[lambda == 0] = 100000
    viol                    = as.vector(abs(score))/as.vector(score.lamb)
    bigviol                 = max(viol)
    
    if (any(sign.change != 1) != TRUE)
    {
      b.lasso[is.in]        = b.est$par[is.in]
      #b.lasso[is.in] = b.est$par
      signs          = sign(b.lasso)
      if (bigviol > 1 + 1.e-6 & is.in[viol == bigviol] == FALSE)
      {
        is.in[viol == bigviol][1]        = TRUE
        is.different[viol == bigviol][1] = TRUE
        #signs[viol == bigviol][1]        = sign(score[viol == bigviol][1])
        signs[viol == bigviol][1]        = -1*sign(score[viol == bigviol][1])
        setchange                        = 2*tol
      }
    }
    
    like    = b.est$value
    
    #like    = like.calc(X, family, ob.wt, mu = update$mu, y, alpha, lambda, beta = b.lasso, penalty = TRUE)
    #mu.old  = update$mu
    #eta.old = update$eta
    
    for (act in 1:length(sign(b.lasso)))
    {
      if (sign(b.lasso)[act] == 0 & sign(as.matrix(betas)[act,dim(as.matrix(betas))[2]]) != 0)
      {
        actions = rbind(actions, paste("Step ", dim(as.matrix(betas))[2] + 1, ": Delete variable ", act, sep = ""))
      }
      if (sign(b.lasso)[act] != 0 & sign(as.matrix(betas)[act,dim(as.matrix(betas))[2]]) == 0)
      {
        actions = rbind(actions, paste("Step ", dim(as.matrix(betas))[2] + 1, ": Add variable ", act, sep = ""))
      }
    }
    
    betas        = cbind(betas, b.lasso)
    scores       = cbind(scores, score)
    viols        = cbind(viols, viol)
    likes        = cbind(likes, like)
    varsets      = cbind(varsets, as.numeric(is.in))
    diff         = setchange + abs(likes[length(likes)] - likes[length(likes) - 1])
    num.it       = num.it + 1
  }
  
  if (error.flag == TRUE)
  {
    return(list(b = NA, mu = NA, e.df = NA, deviance = NA, likelihood = NA, GCV = NA, AIC = NA, BIC = NA, HQC = NA, AICc = NA, ls = NA, pen.likelihood = NA, bs = NA, s = NA, v = NA, actions = NA, flag = "Singular matrix"))
    cat(paste("Singular matrix error. No model fit.", "\n"))
    flush.console()
    stop
  }
  
  if (error.flag != TRUE)
  {
    
    #mu = exp(X %*% b.lasso)
    mu = list()
    psi = list()
    p_detect = list()
    for (i in 1:n.components)
    {
      mu[[i]] = exp(Xenv[[i]] %*% b.lasso[env_ix[[1]]])
      psi[[i]] = 1 - exp(-mu[[i]]*site.area)
      if (dat.type[i] == "PO")
      {
        p_detect[[i]] = NULL
      }
      else
      {
        lin_det  = Xbias[[i]] %*% b.lasso[bias_ix[[i]]]
        if (link == "logit")
        {
          p_detect[[i]] = expit(lin_det)
        }
        if (link == "cloglog")
        {
          p_detect[[i]] = clogloginv(lin_det)
          #p_detect[[i]] = 1 - exp(-exp(lin_det))
        }
      }
    }
    # 
    # eta          = as.matrix(X[,is.in]) %*% as.matrix(b.lasso[is.in])
    # eta[eta < eta.min] = eta.min
    # eta[eta > eta.max] = eta.max
    # mu           = mu.from.eta(eta, family = family, mu.min, mu.max)
    # # if (area.int != FALSE) # no shrinking of AI?
    # # {
    # #   is.in[dim(X)[2]] = FALSE
    # # }
    # wt           = irls.update(y, X, ob.wt, is.in, signs, eta, mu, alpha, lambda, beta.old = b.lasso, penalty = TRUE, family = family, mu.min, mu.max, eta.min, eta.max, tol)$wt
    # xwx          = wt %*% X[,is.in]
    # bpinv        = rep(0, length(b.lasso))
    # bpinv[is.in] = 1/abs(b.lasso[is.in])
    # ginv         = diag(bpinv)	
    # n.param      = length(b.lasso[is.in]) - 1 + area.int
    # 
    # if (sum(is.in) > 1)
    # {
    #   eff.df = sum(diag(solve(xwx + diag(as.vector(lambda[is.in])) %*% ginv[is.in, is.in]) %*% xwx)) + area.int
    # }
    # if (sum(is.in) == 1)
    # {
    #   eff.df = sum(diag(solve(xwx + diag(as.matrix(lambda[is.in])) %*% ginv[is.in, is.in]) %*% xwx)) + area.int
    # }	
    # 
    # n.pres = sum(y > 0)
    # 
    # unp.likelihood = like.calc(X, family, ob.wt, mu, y, alpha, lambda, beta = b.lasso, penalty = FALSE)
    # pen.likelihood = like.calc(X, family, ob.wt, mu, y, alpha, lambda, beta = b.lasso, penalty = TRUE)
    # 
    # if (family$family == "poisson")
    # {
    #   dev            = -2*(sum(ob.wt*(y*log(mu) - mu)) + sum(y > 0)) 
    # }
    # 
    # if (family$family == "binomial")
    # {
    #   dev            = -2*unp.likelihood
    # }
    # 
    # if (family$family == "gaussian")
    # {
    #   dev            = -2*(unp.likelihood - like.calc(X, family, ob.wt[y > 0], y[y > 0], y, alpha, lambda, beta = b.lasso, penalty = FALSE))
    # }
    # 
    # if (eff.df < n.pres)
    # {
    #   gcv     = dev/(n.pres*(1 - (eff.df)/n.pres)^2)
    # }
    # if (eff.df >= n.pres)
    # {
    #   gcv = NA
    # }
    # 
    # aic = -2*unp.likelihood + 2*(n.param + 1)
    # bic = -2*unp.likelihood + log(n.pres)*(n.param + 1)
    # hqc = -2*unp.likelihood + 2*(n.param + 1)*log(log(n.pres))
    # aicc = aic + 2*(n.param + 1)*(n.param + 2)/(n.pres - n.param - 2)
    # 
    # b.out = rep(NA, length(killed))
    # b.out[which(killed == FALSE)] = b.lasso
    # b.out[which(killed == TRUE)] = 0
    
    return(list(b = b.lasso, mu = mu, psi = psi, p_detect = p_detect, ls = likes, pen.likelihood = likes[length(likes)], bs = betas, s = scores, v = viols, actions = actions))
    #return(list(b = b.out, mu = mu, e.df = eff.df, deviance = dev, likelihood = unp.likelihood, GCV = gcv, AIC = aic, BIC = bic, HQC = hqc, AICc = aicc, ls = likes, pen.likelihood = likes[length(likes)], bs = betas, s = scores, v = viols, actions = actions))
    
  }
}

logit = function(pp) { log(pp) - log(1-pp) }

expit = function(eta) {1/(1+exp(-eta))}

cloglog = function(p) {log(-log(1 - p))}

clogloginv = function(lin) {1 - exp(-exp(lin))}

newenv.var = function(sp.xy, env.grid, env.scale, coord = c("X","Y"), file.name = NA)
{
  convert = FALSE
  if (any(lapply(env.grid, class) == "factor"))
  {
    convert  = TRUE
    out.grid = CatConvert(env.grid)
    env.grid = out.grid$X
  }
  x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
  y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
  x.back  = env.grid[,which(names(env.grid) == coord[1])]
  y.back  = env.grid[,which(names(env.grid) == coord[2])]
  x.col   = which(names(env.grid) == coord[1])
  y.col   = which(names(env.grid) == coord[2])
  var.col = setdiff(1:dim(env.grid)[2], c(x.col, y.col))
  s.res   = SpatRes(env.grid)
  
  sp.dat        = as.data.frame(matrix(NA, length(x.dat), length(var.col)))
  names(sp.dat) = names(env.grid[var.col])
  
  for (var in 1:length(var.col))
  {
    loop.scale = min(c(s.res$x.step, s.res$y.step))
    loc        = which(is.na(sp.dat[,var]))
    while(sum(is.na(sp.dat[,var])) > 0)
    {
      loc = which(is.na(sp.dat[,var]))
      sp.dat[loc, var] = newinterp(sp.xy[loc,], loop.scale, env.grid[,var.col[var]], env.grid, coord = c("X","Y"))
      loop.scale = loop.scale*2
    }
    cat(paste("Calculating species environmental data for variable:", names(sp.dat)[var], "\n"))
    flush.console()
  }
  
  sp.dat = data.frame(x.dat, y.dat, sp.dat)
  names(sp.dat)[1:2] = c("X", "Y")
  if (is.na(file.name) == FALSE)
  {
    save.name = paste(file.name, ".RData", sep = "")
    save(sp.dat, file = save.name)
    print(paste("Output saved in the file", save.name))
  }
  if (convert == TRUE)
  {
    sp.dat = list(X = sp.dat, cat.names = out.grid$cat.names)
  }
  sp.dat
}

sample.quad = function(env.grid, sp.scale, coord = c("X", "Y"), file = "Quad")
{
  convert = FALSE
  x.col   = which(names(env.grid) == coord[1])
  y.col   = which(names(env.grid) == coord[2])
  res     = SpatRes(env.grid, coord = coord)
  x.step  = res$x.step
  y.step  = res$y.step
  f.name = list()
  for(i in 1:length(sp.scale))
  {
    i.scale = sp.scale[i]
    x.o = min(env.grid[,x.col]) - floor(min(env.grid[,x.col])/x.step)*x.step
    y.o = min(env.grid[,y.col]) - floor(min(env.grid[,y.col])/y.step)*y.step
    if (x.o/x.step > 0.5)
    {
      x.o = x.o - x.step
    }	
    if (y.o/y.step > 0.5)
    {
      y.o = y.o - y.step
    }
    
    is.on.scale   = abs((env.grid[,x.col]/i.scale) - round(env.grid[,x.col]/i.scale) - x.o/i.scale) + abs((env.grid[,y.col]/i.scale) - round(env.grid[,y.col]/i.scale) - y.o/i.scale) < 1.e-8
    dat.quad      = env.grid[is.on.scale,]
    ### Get rid of machine error in coordinates
    dec.x = max(unlist(lapply(dat.quad[,x.col], DecimalCount)))
    dec.y = max(unlist(lapply(dat.quad[,y.col], DecimalCount)))
    dat.quad[,x.col] = unlist(lapply(dat.quad[, x.col], ZapCoord, dec.x))
    dat.quad[,y.col] = unlist(lapply(dat.quad[, y.col], ZapCoord, dec.y))
    #n.x           = floor((dat.quad[,x.col] - min(dat.quad[,x.col]))/x.step + 0.1)
    #n.y           = floor((dat.quad[,y.col] - min(dat.quad[,y.col]))/y.step + 0.1)
    #dat.quad[,x.col] = min(dat.quad[,x.col]) + n.x*x.step
    #dat.quad[,y.col] = min(dat.quad[,y.col]) + n.y*y.step
    
    if (is.na(file) == FALSE)
    {
      f.name[[i]] = paste(file, sp.scale[i], ".RData", sep = "")
      save(dat.quad, file = f.name[[i]])
      print(paste("Output saved in the file", f.name[[i]]))
    }
  }
  if (convert == TRUE)
  {
    dat.quad = list(X = dat.quad, cat.names = cat.names)
  }
  if (length(sp.scale) == 1)
    dat.quad
  else
    f.name
}

weights = function(sp.xy, quad.xy, coord = c("X", "Y"))
{
  sp.col   = c(which(names(sp.xy) == coord[1]), which(names(sp.xy) == coord[2]))
  quad.col = c(which(names(quad.xy) == coord[1]), which(names(quad.xy) == coord[2]))
  
  X.inc   = sort(unique(quad.xy[,quad.col[1]]))[2] - sort(unique(quad.xy[,quad.col[1]]))[1]
  Y.inc   = sort(unique(quad.xy[,quad.col[2]]))[2] - sort(unique(quad.xy[,quad.col[2]]))[1]
  quad.0X = min(quad.xy[,quad.col[1]]) - floor(min(quad.xy[,quad.col[1]])/X.inc)*X.inc
  quad.0Y = min(quad.xy[,quad.col[2]]) - floor(min(quad.xy[,quad.col[2]])/Y.inc)*Y.inc
  
  X = c(sp.xy[,quad.col[1]], quad.xy[,quad.col[1]])
  Y = c(sp.xy[,quad.col[2]], quad.xy[,quad.col[2]])
  
  round.X     = round((X - quad.0X)/X.inc)*X.inc
  round.Y     = round((Y - quad.0Y)/Y.inc)*Y.inc
  round.id    = paste(round.X, round.Y)
  round.table = table(round.id)
  wt          = X.inc*Y.inc/as.numeric(round.table[match(round.id, names(round.table))])
  
  wt
}

SpatRes = function(env.grid, coord = c("X", "Y"))
{
  x.col   = which(names(env.grid) == coord[1])
  y.col   = which(names(env.grid) == coord[2])
  x.uq    = sort(unique(env.grid[, x.col]))
  y.uq    = sort(unique(env.grid[, y.col]))
  n.dec   = max(unlist(lapply(x.uq, DecimalCount)))
  x.diff  = diff(x.uq)
  y.diff  = diff(y.uq)
  x.dec   = unlist(lapply(x.diff, DecimalCount))
  y.dec   = unlist(lapply(y.diff, DecimalCount))
  x.step  = min(floor(x.diff*10^max(x.dec) + 0.1))/(10^max(x.dec))
  y.step  = min(floor(y.diff*10^max(y.dec) + 0.1))/(10^max(y.dec))
  return(list(x.step = x.step, y.step = y.step))
}

ZapCoord = function(x, numdec = DecimalCount(x))
{
  x.out  = floor(x*10^numdec + 0.1)/(10^numdec)
  x.out
}

DecimalCount = function(x, max.dec = max(10, max(nchar(x))), tol = 1.e-1)
{
  digits  = 0:max.dec
  x.diff = (x - round(x, digits))/(10^(-1*(digits + 3)))
  num.dec = digits[min(which(abs(x.diff) < tol))]
  num.dec
}

newinterp = function(sp.xy, sp.scale, f, back.xy, coord = c("X","Y"))
{
  options(scipen = 999)
  x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
  y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
  x.back  = back.xy[,which(names(back.xy) == coord[1])]
  y.back  = back.xy[,which(names(back.xy) == coord[2])]
  
  grid    = data.table(x.back, y.back, f, key = c("x.back", "y.back"))
  
  ux    = sort(unique(x.back))
  uy    = sort(unique(y.back))
  
  x.col   = which(names(back.xy) == coord[1])
  y.col   = which(names(back.xy) == coord[2])
  
  x.step  = ux[2] - ux[1]
  y.step  = uy[2] - uy[1]
  
  x.o = min(back.xy[,x.col]) - floor(min(back.xy[,x.col])/x.step)*x.step
  y.o = min(back.xy[,y.col]) - floor(min(back.xy[,y.col])/y.step)*y.step
  
  x.1   = floor((x.dat - x.o)/sp.scale)*sp.scale + x.o
  y.1   = floor((y.dat - y.o)/sp.scale)*sp.scale + y.o
  x.2   = pmin(x.1 + sp.scale, max(ux))
  y.2   = pmin(y.1 + sp.scale, max(uy))
  
  w11   = (x.2 - x.dat)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
  w12   = (x.2 - x.dat)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))
  w21   = (x.dat - x.1)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
  w22   = (x.dat - x.1)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))
  
  f11   = grid[list(x.1, y.1)]$f
  f12   = grid[list(x.1, y.2)]$f
  f21   = grid[list(x.2, y.1)]$f
  f22   = grid[list(x.2, y.2)]$f
  
  c11 = 1 - is.na(f11)
  c12 = 1 - is.na(f12)
  c21 = 1 - is.na(f21)
  c22 = 1 - is.na(f22)
  
  env.wt.mat = cbind(f11*w11*c11, f12*w12*c12, f21*w21*c21, f22*w22*c22) 
  
  f.interp = apply(env.wt.mat, 1, sum, na.rm = TRUE)/(w11*c11 + w12*c12 + w21*c21 + w22*c22)
  f.interp
}

standardise.X = function(mat)
{
  X = scale(as.matrix(mat))
  dat.means = apply(as.matrix(mat), 2, mean, na.rm = TRUE)
  dat.sds   = apply(as.matrix(mat), 2, sd, na.rm = TRUE)
  return(list(X = X, dat.means = dat.means, dat.sds = dat.sds))
}

polynamesadd = function(X)
{
  Xnames = polynames(X)
  if (is.null(Xnames) == FALSE)
  {
    dimnames(X)[[2]][Xnames$ids] = Xnames$names
  }
  X
}

polynames = function(X)
{
  coef.names = dimnames(X)[[2]]
  which.poly = grep("poly\\(", coef.names)
  if (length(which.poly) > 0)
  {
    split1     = strsplit(coef.names[which.poly], "\\)")
    polyframe  =  plyr::ldply(split1, rbind)
    varframe   = plyr::ldply(strsplit(unlist(strsplit(as.character(polyframe[,1]), "poly\\("))[2*(1:(length(which.poly)))], ","), rbind)
    varframe   = varframe[,-dim(varframe)[2]]
    varframe   = data.frame(lapply(varframe, as.character), stringsAsFactors = FALSE)
    if (dim(varframe)[1] == 1)
    {
      varframe = t(varframe)
    }
    expframe   = plyr::ldply(strsplit(as.character(polyframe[,2]), "\\."), rbind)
    expframe   = data.frame(lapply(expframe, as.character), stringsAsFactors = FALSE)
    expframe[is.na(expframe)] = 0
    vframe = varframe[,1:dim(expframe)[2]]
    nameframe  = matrix(paste(as.matrix(vframe), "^", as.matrix(expframe), sep = ""), nrow(vframe), ncol(vframe))
    nameframe[expframe == "0"] = ""
    nameframe[expframe == "1"] = as.character(vframe[expframe == 1])
    nameframe = gsub(" ", "", nameframe)
    nameframe  = as.data.frame(nameframe)
    if (is.null(dim(nameframe)) == FALSE)
    {
      names.out = apply(nameframe, 1, function(row) paste(row[nzchar(row)], collapse = "*"))
    }
    if (is.null(dim(nameframe)) == TRUE)
    {
      names.out = nameframe
    }
    id.out    = which.poly
  }
  else
  {
    names.out = NULL
    id.out = NULL
  }
  return(list(names = names.out, ids = id.out))
}

point.interactions = function(env, pres, r, coord = c("X", "Y"), availability = NA)
{
  if (any(is.na(availability)))
  {
    availability = makeMask(env, coord)
  }	
  
  quad.x = env[,match(coord[1], names(env))]
  quad.y = env[,match(coord[2], names(env))]
  
  pres.x = pres[,match(coord[1], names(pres))]
  pres.y = pres[,match(coord[2], names(pres))]
  
  cat(paste("Calculating point interactions", "\n"))
  flush.console()
  occupied = matrix(0, dim(availability)[1], dim(availability)[2])
  rownames(occupied) = rownames(availability)
  colnames(occupied) = colnames(availability)
  
  x.mat = availability
  y.mat = availability
  
  for (i in 1:dim(x.mat)[1])
  {
    x.mat[i,] = as.numeric(colnames(availability))
  }
  for (i in 1:dim(y.mat)[2])
  {
    y.mat[,i] = as.numeric(rownames(availability))
  }
  
  grain = as.numeric(colnames(availability)[2]) - as.numeric(colnames(availability)[1])
  
  quad.int = rep(0, length(quad.x))
  
  for (i in 1:length(pres.x))
  {
    sub.col = which(as.numeric(colnames(occupied)) >= pres.x[i] - (r + grain) & as.numeric(colnames(occupied)) <= pres.x[i] + (r + grain))
    sub.row = which(as.numeric(rownames(occupied)) >= pres.y[i] - (r + grain) & as.numeric(rownames(occupied)) <= pres.y[i] + (r + grain))
    
    sub.occ = occupied[sub.row, sub.col]
    sub.x   = x.mat[sub.row, sub.col]
    sub.y   = y.mat[sub.row, sub.col]
    
    sub.occ[(sub.x - pres.x[i])^2 + (sub.y - pres.y[i])^2 < r^2] = sub.occ[(sub.x - pres.x[i])^2 + (sub.y - pres.y[i])^2 < r^2] + 1
    occupied[sub.row, sub.col] = sub.occ
    quad.cells           = which((quad.x - pres.x[i])^2 + (quad.y - pres.y[i])^2 <= (2*r)^2)
    quad.int[quad.cells] = quad.int[quad.cells] + 1
  }
  
  int.q = rep(0, length(quad.int))
  
  for (quad.i in which(quad.int > 0))
  {
    sub.col = which(as.numeric(colnames(occupied)) >= quad.x[quad.i] - (r + grain) & as.numeric(colnames(occupied)) <= quad.x[quad.i] + (r + grain))
    sub.row = which(as.numeric(rownames(occupied)) >= quad.y[quad.i] - (r + grain) & as.numeric(rownames(occupied)) <= quad.y[quad.i] + (r + grain))
    
    sub.occ  = occupied[sub.row, sub.col]
    sub.availability = availability[sub.row, sub.col]
    sub.x    = x.mat[sub.row, sub.col]
    sub.y    = y.mat[sub.row, sub.col]
    
    sub.cell = (sub.x - quad.x[quad.i])^2 + (sub.y - quad.y[quad.i])^2 <= r^2 & sub.availability > 0
    
    #int.q[quad.i] = sum(sub.occ[sub.cell] > 0)/sum(sub.cell)
    int.q[quad.i] = sum(sub.occ[sub.cell] > 0, na.rm = TRUE)/sum(sub.cell, na.rm = TRUE)
  }
  
  int.p = rep(0, length(pres.x))
  
  for (pres.i in 1:length(pres.x))
  {
    sub.col = which(as.numeric(colnames(occupied)) >= pres.x[pres.i] - (r + grain) & as.numeric(colnames(occupied)) <= pres.x[pres.i] + (r + grain))
    sub.row = which(as.numeric(rownames(occupied)) >= pres.y[pres.i] - (r + grain) & as.numeric(rownames(occupied)) <= pres.y[pres.i] + (r + grain))
    
    sub.occ  = occupied[sub.row, sub.col]
    sub.availability = availability[sub.row, sub.col]
    sub.x    = x.mat[sub.row, sub.col]
    sub.y    = y.mat[sub.row, sub.col]
    
    sub.cell = (sub.x - pres.x[pres.i])^2 + (sub.y - pres.y[pres.i])^2 <= r^2 & sub.availability > 0
    
    #int.p[pres.i] = sum(sub.occ[sub.cell] > 1)/sum(sub.cell)
    int.p[pres.i] = if (sum(sub.cell, na.rm = TRUE) == 0) 0 else sum(sub.occ[sub.cell] > 1, na.rm = TRUE)/sum(sub.cell, na.rm = TRUE)
    
  }
  
  interactions = c(int.p, int.q)
  interactions
}

makeMask = function(env, coord = c("X", "Y"))
{
  x.col = match(coord[1], names(env))
  y.col = match(coord[2], names(env))
  ux = sort(unique(env[,x.col]))
  uy = sort(unique(env[,y.col]))
  nx = length(ux)
  ny = length(uy)
  
  col.ref = match(env[,x.col], ux)
  row.ref = match(env[,y.col], uy)
  
  all.vec          = rep(0, max(row.ref)*max(col.ref))
  vec.ref          = (col.ref - 1)*max(row.ref) + row.ref
  all.vec[vec.ref] = 1
  mask.out         = matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux))
  mask.out
}

LL_Lasso_Wt_Comb_new = function(par, Xenv, Xbias, y, ob.wt, link = "logit", env_ix, bias_ix, site.area = 1, coord = c("x", "y"), lasso = 0, wt.vec, dat.type = "PO", is.in = rep(TRUE, length(par)))
{
  Xenv.use = subsetlist(Xenv, env_ix, keep = which(is.in == TRUE))
  Xbias.use = subsetlist(Xbias, bias_ix, keep = which(is.in == TRUE))
  env_ix.use = subsetix(env_ix, keep = which(is.in == TRUE))
  bias_ix.use = subsetix(bias_ix, keep = which(is.in == TRUE))
  #ix = mapply(c, env_ix.use, bias_ix.use)
  ix = lapply(Map(list, env_ix.use, bias_ix.use), unlist)
  n.components = length(ix)
  ll.comp = rep(0, n.components)
  for (i in 1:n.components)
  {
    ll.comp[i] = if (dat.type[i] == "PO") 
      wt.vec[i]*wtlassonegLL.po(beta = par[ix[[i]]], X.des = cbind(Xenv.use[[i]], Xbias.use[[i]]), 
                                ob.wt = ob.wt[[i]], y = y[[i]], lasso = 0)
    else wt.vec[i]*wtlassonegLL.occ(par = par[ix[[i]]], W = Xbias.use[[i]], 
                                    X = Xenv.use[[i]], y = y[[i]], site.area = site.area, link = link, lasso = 0)
    
    # ll.comp[i] = if (dat.type[i] == "PO") 
    #                   wt.vec[i]*wtlassonegLL.po(beta = par[ix[[i]]], X.des = cbind(Xenv[[i]], Xbias[[i]]), 
    #                                             ob.wt = ob.wt[[i]], y = y[[i]], lasso = 0)
    #              else wt.vec[i]*wtlassonegLL.occ(par = par[ix[[i]]], W = Xbias[[i]], 
    #                                              X = Xenv[[i]], y = y[[i]], site.area = site.area, lasso = 0)
  }
  maxf = sum(ll.comp) + sum(as.vector(lasso)*abs(par))
  maxf
}

subsetlist = function(mylist, myix, keep)
{
  keep.ix  = lapply(myix, function(x) x %in% keep) # list of indices to keep
  keep.out = if (length(keep.ix) == 1) list(mylist[[1]][,keep.ix[[1]], drop = FALSE]) else lapply(mapply(subset, lapply(mylist, t), keep.ix), t) # list of subsetted matrices
  keep.dim = unlist(sapply(lapply(keep.out, dim), function(x) x[2])) # number of columns in subsetted matrices
  out.list = list()
  for (i in 1:length(keep.dim))
  {
    if (is.null(keep.dim[i]))
    {
      out.list[i] = list(NULL)
    }
    else
    {
      if (keep.dim[i] > 0)
      {
        out.list[[i]] = keep.out[[i]]
      }
      else
      {
        out.list[i] = list(NULL)
      }
    }
    #out.list[[i]] = if (keep.dim[i] > 0) keep.out[[i]] else list(NULL)
  }
  out.list
}

subsetix = function(myix, keep)
{
  keep.ix  = lapply(myix, function(x) intersect(x, keep))
  keep.ix
}

wtlassonegLL.po = function(beta, X.des, ob.wt, y, lasso = 0) {
  # coord.col = match(coord, names(X.back))
  # back.xy   = X.back[, coord.col]
  # po.xy     = X.po[, coord.col]
  # X.back    = X.back[, -coord.col]
  # X.po      = X.po[, -coord.col]
  ## ob.wt     = weights(po.xy, back.xy, coord = coord)
  ## pres      = c(rep(1, dim(X.po)[1]), rep(0, dim(X.back)[1]))
  ## y         = pres/ob.wt
  
  # if (is.null(W.back) == FALSE)
  # {
  #   beta    = par[1:dim(X.back)[2]]
  #   alpha   = par[(dim(X.back)[2] + 1):length(par)]
  #   #	 X.des   = cbind(as.matrix(rbind(X.po, X.back)), as.matrix(rbind(W.po, W.back)))
  #   all.par = c(beta, alpha)
  # }
  # else
  # {
  #   beta = par
  #   all.par = beta
  #   #	X.des   = as.matrix(rbind(X.po, X.back))
  # }
  
  mu = exp(as.matrix(X.des) %*% beta)
  
  #logL.po = sum(ob.wt*(y*log(mu) - mu)) - sum(log(1:sum(y > 0))) - sum(alpha * as.vector(lambda)*abs(beta))
  logL.po = sum(ob.wt*(y*log(mu) - mu)) - sum(log(1:sum(y > 0))) - sum(as.vector(lasso)*abs(beta))
  (-1)*sum(logL.po)
}

wtlassonegLL.occ = function(par, W, X, y, site.area, link = "logit", lasso = 0){
  beta = par[1:dim(as.matrix(X))[2]]
  alpha = par[(dim(as.matrix(X))[2]+1):length(par)]
  if (is.null(W)) # added March 2019
  {
    alpha = NULL
    p = rep(0.5, nrow(y))
  }
  if (link == "logit" & is.null(W) == FALSE)
  {
    p = expit(as.matrix(W) %*% alpha)
  }
  if (link == "cloglog" & is.null(W) == FALSE)
  {
    p = clogloginv(as.matrix(W) %*% alpha)
  }
  #p = 1/(1 + exp(-as.matrix(W) %*% alpha))
  M = nrow(y)
  J = ncol(y) - apply(is.na(y), 1, sum)
  y.bin = data.frame(matrix(as.numeric(as.matrix(y)), M, ncol(y)))
  y.bin = apply(y.bin, 1, sum, na.rm = TRUE)
  
  #psi<-1/(1+exp(-X %*% beta))
  # cloglog(psi) = log(lambda); 
  # replace the logit link by the cloglog link function, then use the exponential to get lambda
  
  lambda = exp(as.matrix(X) %*% beta)
  psi = 1 - exp(-lambda*site.area)
  logL.occ = sum(log(dbinom(y.bin, J, p)*psi + ifelse(y.bin == 0, 1, 0)*(1 - psi))) - sum(as.vector(lasso)*abs(c(beta, alpha)))
  (-1)*sum(logL.occ)
}

scoreeq_new = function(Xenv, Xbias = NULL, y, par, ob.wt = NULL, link = "logit", site.area = 1, dat.type = "PO", lasso = 0, env_ix, bias_ix, coord = c("x", "y"), wt.vec, is.in = rep(TRUE, length(par)))
{
  if (class(y) == "list")
  {
    n.components = length(Xenv)
    num_env  = dim(Xenv[[1]])[2]
    num_bias = if (length(unlist(Xbias)) == 0) 0 else (unlist(lapply(lapply(Xbias, as.matrix), ncol)))
    num_p    = num_env + sum(num_bias)
    ix       = lapply(Map(list, env_ix, bias_ix), unlist)
    score_mat = matrix(0, num_p, n.components)
  }
  else
  {
    num_env  = dim(Xenv)[2]
    num_bias = if (is.null(Xbias) == FALSE) ncol(as.matrix(Xbias))
    else 0
    num_p    = num_env + num_bias
    n.components = 1
    ix       = list(1:num_p)
    Xenv     = list(Xenv)
    Xbias    = list(Xbias)
    y        = list(y)
    ob.wt    = list(ob.wt)
    score_mat = matrix(0, num_p, 1)
  }
  for (i in 1:n.components)
  {
    if (dat.type[i] == "PO")
    {
      X.i = cbind(Xenv[[i]], Xbias[[i]])
      score_mat[ix[[i]], i] = ppmscoreeq_new(y = y[[i]], X = X.i, par = par[ix[[i]]], ob.wt = ob.wt[[i]])
    }
    if (dat.type[i] == "Occ")
    {
      score_mat[ix[[i]], i] = occscoreeq_new(y = y[[i]], Xenv = Xenv[[i]], Xbias = Xbias[[i]], par = par[ix[[i]]], 
                                             ob.wt = ob.wt[[i]], link = link, site.area = site.area)
    }
  }
  #score_out = apply(score_mat, 1, sum)
  score_out = -1*(apply(score_mat, 1, sum) - lasso*sign(par))
  score_out[is.in == FALSE] = 0
  score_out
}

occscoreeq_new = function(Xenv, Xbias = NULL, y, par, ob.wt = NULL, link = "logit", site.area = 1)
{
  env_cov  = 1:dim(Xenv)[2]
  mu       = exp(Xenv %*% par[env_cov])
  psi      = 1 - exp(-mu*site.area)
  if (is.null(Xbias) == FALSE)
  {
    bias_cov = (max(env_cov) + 1):(max(env_cov) + dim(Xbias)[2])
    lin_det  = Xbias %*% par[bias_cov]
    if (link == "logit")
    {
      p_detect = expit(lin_det)
    }
    if (link == "cloglog")
    {
      p_detect = clogloginv(lin_det)
      #p_detect = 1 - exp(-exp(lin_det))
    }
  }
  
  M = nrow(y)
  J = ncol(y) - apply(is.na(y), 1, sum)
  y.bin = data.frame(matrix(as.numeric(as.matrix(y)), M, ncol(y)))
  y.bin = apply(y.bin, 1, sum, na.rm = TRUE)
  
  detected = which(y.bin > 0)
  
  env_pres = (site.area*mu[detected]*Xenv[detected,])/(exp(site.area*mu[detected]) - 1)
  env_abs  = (site.area*mu[-detected]*Xenv[-detected,]*((1 - p_detect[-detected])^J[-detected] - 1))/((1 - p_detect[-detected])^J[-detected]*(exp(site.area*mu[-detected]) - 1) + 1)
  
  score_env = apply(env_pres, 2, sum) + apply(env_abs, 2, sum)
  score_out = score_env
  # score_env = sum((site.area*mu[detected]*Xenv[detected,])/(exp(site.area*mu[detected]) - 1)) 
  #           + sum((site.area*mu[-detected]*Xenv[-detected,]*((1 - p_detect[-detected])^J[-detected] - 1))/((1 - p_detect[-detected])^J[-detected]*(exp(site.area*mu[-detected]) - 1) + 1))
  
  if (is.null(Xbias) == FALSE)
  {
    if (link == "logit")
    {
      bias_pres = Xbias[detected,]*((y.bin[detected] - J[detected])*exp(lin_det[detected]) + y.bin[detected])/(exp(lin_det[detected]) + 1)
      #bias_abs  = -1*(J[-detected]*psi[-detected]*Xbias[-detected,]*(1 - (1/exp(-lin_det[-detected]) + 1))^J[-detected]*exp(lin_det[-detected]))/((psi[-detected]*(1 - (1/exp(-lin_det[-detected]) + 1))^J[-detected] - psi[-detected] + 1)*(exp(lin_det[-detected]) + 1))
      #bias_abs  = -1*(J[-detected]*psi[-detected]*Xbias[-detected,]*(1 - (1/(exp(-lin_det[-detected]) + 1)))^J[-detected]*exp(lin_det[-detected]))/((psi[-detected]*(1 - (1/(exp(-lin_det[-detected]) + 1)))^J[-detected] - psi[-detected] + 1)*(exp(lin_det[-detected]) + 1))
      
      
      bias_abs  = -1*(J[-detected]*psi[-detected]*Xbias[-detected,]*(1 - p_detect[-detected])^J[-detected]*exp(lin_det[-detected]))/((psi[-detected]*(1 - p_detect[-detected])^J[-detected] - psi[-detected] + 1)*(exp(lin_det[-detected]) + 1))
      
      score_bias = if (is.null(dim(bias_pres))) sum(bias_pres) + sum(bias_abs) else apply(bias_pres, 2, sum) + apply(bias_abs, 2, sum)
      # score_bias = sum(Xbias[detected,]*((y.bin[detected] - J[detected])*exp(lin_det[detected]) + y.bin[detected])/(exp(lin_det[detected] + 1))) 
      #            + sum(-1*(J[-detected]*psi[-detected]*Xbias[-detected]*(1 - (1/exp(-lin_det[-detected]) + 1))^J[-detected]*exp(lin_det[-detected]))/((psi[-detected]*(1 - (1/exp(-lin_det[-detected]) + 1))^J[-detected] - psi[-detected] + 1)*(exp(lin_det[-detected]) + 1)))
    }
    
    if (link == "cloglog")
    {
      bias_pres = Xbias[detected,]*exp(lin_det[detected])*((y.bin[detected] - J[detected])*(exp(exp(lin_det[detected]))) + J[detected])/(exp(exp(lin_det[detected])) - 1)
      bias_abs  = J[-detected]*psi[-detected]*Xbias[-detected,] * exp(lin_det[-detected])/((psi[-detected] - 1)*exp(J[-detected]*exp(lin_det[-detected])) - psi[-detected])
      
      score_bias = if (is.null(dim(bias_pres))) sum(bias_pres) + sum(bias_abs) else apply(bias_pres, 2, sum) + apply(bias_abs, 2, sum)
      # edited March 2019
      
      # score_bias = sum(Xbias[detected,]*exp(lin_det[detected])*((y.bin[detected] - J[detected])*(exp(exp(lin_det[detected]))) + J[detected])/(exp(exp(lin_det[detected])) - 1))
      #            + sum(J[-detected]*psi[-detected]*Xbias[-detected,] * exp(lin_det[-detected])/((psi[-detected] - 1)*exp(J[-detected]*exp(lin_det[-detected])) - psi[-detected]))
    }
    score_out = c(score_env, score_bias)
  }
  score_out
}

ppmscoreeq_new = function(y, X, par, ob.wt, family = poisson())
{
  mu    = exp(as.matrix(X) %*% par)
  vari  = family$variance(mu)
  deriv = 1/vari
  weii  = ob.wt*1/(deriv*vari*deriv)
  Xw.s  = t(as.vector(weii) * t(t(X)))
  score = t(as.vector(deriv) * t(Xw.s)) %*% (y - mu)
  score
}

plotfit = function(fit, pred_data = NULL, xy = NULL, z = "intensity", model = NULL, source = NULL, link = "logit",
                   coord = c("X", "Y"), asp = "iso", ylab = "", xlab = "", col.regions = heat.colors(1024)[900:1], 
                   cuts = length(col.regions), cex = 1.4, main.text = NULL, cex.color = 1.4)
{
  # To do:
  # Allow users to input env data and coordinates for plotting
  
  beta = if (is.null(model)) fit$beta else fit$betas[,model]
  plot_pen = if (is.null(model)) fit$penalty else fit$penalty_vec[model]
  if (is.null(source))
  {
    if (z == "intensity" | z == "bias" | z == "occupancy")
    {
      source = min(which(fit$dat.type == "PO"))
    }
    if (z == "p_detect")
    {
      source = min(which(fit$dat.type == "Occ"))
    }
  }
  
  num_env  = dim(fit$X_env[[1]])[2]
  num_bias = unlist(lapply(lapply(fit$X_bias, as.matrix), ncol))
  num_p    = num_env + sum(num_bias)
  n.components = length(num_bias)
  bias_ind1 = num_env + 1 + cumsum(num_bias) - num_bias
  bias_ind2 = num_env + cumsum(num_bias)
  env_ix = rep(list(1:num_env), n.components)
  bias_ix = list()
  for (i in 1:n.components)
  {
    bias_ix[[i]] = bias_ind1[i]:bias_ind2[i]
  }
  
  if (fit$dat.type[source] == "PO")
  {
    Xdes = if (z == "intensity" | z == "occupancy") fit$X_env[[source]][fit$y[[source]] == 0,] else fit$X_bias[[source]][fit$y[[source]] == 0,]
    beta_sub = if (z == "intensity" | z == "occupancy") beta[env_ix[[source]]] else beta[bias_ix[[source]]]
    XY = fit$quad_xy
  }
  if (fit$dat.type[source] == "Occ")
  {
    Xdes = if (z == "intensity" | z == "occupancy") fit$X_env[[source]] else fit$X_bias[[source]]
    beta_sub = if (z == "intensity" | z == "occupancy") beta[env_ix[[source]]] else beta[bias_ix[[source]]]
    XY = fit$sp_xy[[source]]
  }
  
  if (z == "intensity" | z == "bias")
  {
    plot_z = exp(as.matrix(Xdes) %*% beta_sub)
  }
  if (z == "occupancy")
  {
    plot_z = 1 - exp(-exp(as.matrix(Xdes) %*% beta_sub)*fit$site.area) # need site.area from fit!
  }
  if (z == "p_detect")
  {
    lin_det  = as.matrix(Xdes) %*% beta_sub
    if (link == "logit")
    {
      plot_z = expit(lin_det)
    }
    if (link == "cloglog")
    {
      plot_z = clogloginv(lin_det)
      #plot_z = 1 - exp(-exp(lin_det))
    }
  }
  
  plot_x = XY[,match(coord[1], names(XY))]
  plot_y = XY[,match(coord[2], names(XY))]
  
  if (z == "p_detect")
  {
    z = "Detection Probability"
  }
  
  z_name = strsplit(z, " ")[[1]]
  z_name = paste(toupper(substring(z_name, 1, 1)), substring(z_name, 2),
                 sep = "", collapse = " ")
  
  if (is.null(main.text))
  {
    main.text = paste(z_name, " from Source ", source, "\n", "Model with penalty ", round(plot_pen, 3), sep = "")
  }
  
  levelplot(plot_z ~ plot_x + plot_y, asp = asp, 
            ylab = ylab, xlab = xlab, col.regions = col.regions, cuts = cuts, 
            main = list(main.text, cex = cex), 
            scales = list(y = list(draw = FALSE), x = list(draw = FALSE)), 
            cex = cex, colorkey = list(labels = list(cex = cex.color)))
}

predfit = function(fit, pred_data = NULL, xy = NULL, z = "intensity", model = NULL, source = NULL, link = "logit",
                   use = NULL, coord = c("X", "Y"), asp = "iso", ylab = "", xlab = "", col.regions = heat.colors(1024)[900:1], 
                   cuts = length(col.regions), cex = 1.4, main.text = NULL, cex.color = 1.4)
{
  # To do:
  # Allow users to input env data and coordinates for plotting
  
  beta = if (is.null(model)) fit$beta else fit$betas[,model]
  plot_pen = if (is.null(model)) fit$penalty else fit$penalty[model]
  if (is.null(source))
  {
    if (z == "intensity" | z == "bias" | z == "occupancy")
    {
      source = min(which(fit$dat.type == "PO"))
    }
    if (z == "p_detect")
    {
      source = min(which(fit$dat.type == "Occ"))
    }
  }
  
  num_env  = dim(fit$X_env[[1]])[2]
  num_bias = unlist(lapply(lapply(fit$X_bias, as.matrix), ncol))
  num_p    = num_env + sum(num_bias)
  n.components = length(num_bias)
  bias_ind1 = num_env + 1 + cumsum(num_bias) - num_bias
  bias_ind2 = num_env + cumsum(num_bias)
  env_ix = rep(list(1:num_env), n.components)
  bias_ix = list()
  for (i in 1:n.components)
  {
    bias_ix[[i]] = bias_ind1[i]:bias_ind2[i]
  }
  
  beta_sub = if (z == "intensity" | z == "occupancy") beta[env_ix[[source]]] else beta[bias_ix[[source]]]
  if (is.null(pred_data))
  {
    if (fit$dat.type[source] == "PO")
    {
      Xdes = if (z == "intensity" | z == "occupancy") fit$X_env[[source]][fit$y[[source]] == 0,] else fit$X_bias[[source]][fit$y[[source]] == 0,]
      #beta_sub = if (z == "intensity" | z == "occupancy") beta[env_ix[[source]]] else beta[bias_ix[[source]]]
      XY = fit$quad_xy
    }
    if (fit$dat.type[source] == "Occ")
    {
      Xdes = if (z == "intensity" | z == "occupancy") fit$X_env[[source]] else fit$X_bias[[source]]
      #beta_sub = if (z == "intensity" | z == "occupancy") beta[env_ix[[source]]] else beta[bias_ix[[source]]]
      XY = fit$sp_xy[[source]]
    }
    if (is.null(use) == FALSE)
    {
      Xcols = list()
      Xcols[[1]] = which(env_ix[[1]] %in% use)
      for (j in 1:length(bias_ix))
      {
        Xcols[[j + 1]] = which(bias_ix[[j]] %in% use)
      }
      
      Xdes = fit$X_env[[1]][fit$y[[1]] == 0, Xcols[[1]]]
      set_use = unlist(lapply(Xcols, length))
      for (i in 2:length(Xcols))
      {
        if (set_use[[i]] > 0)
        {
          Xdes = cbind(Xdes, fit$X_bias[[i - 1]][fit$y[[i - 1]] == 0, Xcols[[i]]])
        }
      }
      XY = fit$quad_xy
      beta_sub = beta[use]
    }
  }
  else
  {
    Xdes = pred_data
    XY = xy
  }
  
  if (z == "intensity" | z == "bias")
  {
    pred_z = exp(as.matrix(Xdes) %*% beta_sub)
  }
  if (z == "occupancy")
  {
    pred_z = 1 - exp(-exp(as.matrix(Xdes) %*% beta_sub)*fit$site.area) # need site.area from fit!
  }
  if (z == "p_detect")
  {
    lin_det  = as.matrix(Xdes) %*% beta_sub
    if (link == "logit")
    {
      pred_z = expit(lin_det)
    }
    if (link == "cloglog")
    {
      pred_z = clogloginv(lin_det)
      #pred_z = 1 - exp(-exp(lin_det))
    }
  }
  
  pred_z
}

plotpath = function(fit, colors = c("gold", "green3", "blue", "pink"), v.cut = NULL)
{
  best.models = apply(fit$criterion.matrix, 2, which.min) #See which fit optimises other criteria
  unique.beta = unique(best.models)
  names = rep(NA, length(unique.beta))
  for (i in 1:length(unique.beta))
  {
    names[i] = paste(names(best.models[best.models == unique.beta[i]]), collapse = "/")
  }
  
  v = if (is.null(v.cut)) 1:dim(fit$betas)[1] else setdiff(1:dim(fit$betas)[1], v.cut)
  
  betas = fit$betas[v,]
  
  min.y = min(apply(betas, 1, min))
  max.y = max(apply(betas, 1, max))
  
  lambdas = fit$penalty_vec
  
  if (lambdas[1] == 0)
  {
    lambda_d = diff(log(lambdas[2:3]))
    lambdas[1] = exp(log(lambdas[2]) - lambda_d)
  }
  
  
  plot(lambdas, betas[1,], log = "x", type = "l", ylim = c(min.y, max.y), xlab = "LASSO penalty", ylab = "Coefficients")
  for (i in 2:dim(betas)[1])
  {
    points(lambdas, betas[i,], type = "l")
  }
  
  for (i in 1:length(unique.beta))
  {
    abline(v = lambdas[unique.beta[i]], lwd = 3, col = colors[i])
  }	
  
  legend("topright", names, lwd = rep(3, length(unique.beta)), col = colors[1:length(unique.beta)])
}

criterion_curve = function(fit, criterion = "BIC")
{
  plot_y = fit$criterion.matrix[,match(criterion, names(fit$criterion.matrix))]
  lambdas = fit$penalty_vec
  
  if (lambdas[1] == 0)
  {
    lambda_d = diff(log(lambdas[2:3]))
    lambdas[1] = exp(log(lambdas[2]) - lambda_d)
  }
  min_pen = fit$penalty_vec[which.min(plot_y)]
  plot(lambdas, plot_y, xlab = "LASSO Penalty", ylab = criterion, type = "o", log = "x")
  points(lambdas[which.min(plot_y)], min(plot_y), col = "green3", pch = 19)
  legend("bottomright", legend = bquote("Min" ~ .(criterion): ~ lambda == .(min_pen) ~ phantom(x)  ), col = "green3", pch = 19)
}
