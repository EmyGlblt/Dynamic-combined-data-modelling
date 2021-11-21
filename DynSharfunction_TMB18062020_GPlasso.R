comb_lasso = function(env_formula, bias_formula = NULL, intercept_env = NULL, intercept_bias = NULL, 
                      quad_data = NULL, sp_data, sp_y, dat.type = "PO", coord = c("X", "Y"), sp_res = 1, 
                      penalty_vec = NULL, alpha = 1, gamma = 0, init.coef = NA, standardise = TRUE, criterion = "BIC",
                      family = "poisson", tol = 1.e-7, b.min = 1.e-6, 
                      max.it = 25, n.fits = 50, noshrink = rep(list(NULL), n.time), method = "BFGS", 
                      link = "logit", site.area = 1, obs.scores = NULL,
                      extinct_comp=c(NULL), coloni_comp=c(NULL),
                      area.int = FALSE, r = NULL, wt.vec = NULL, pen.min = 1.e-2, verbose = FALSE)
{
  
  # To do:
  # Automatic choice of penalty vector
  # output standardisation means and standard deviations
  # block CV
  
  n.time = length(sp_data)
  
  if( length(sp_data[[1]])!=length(sp_data[[2]]) ) stop('all data types are not dynamic')
  
  start.time = Sys.time()
  formula.out = list(env_formula, bias_formula)
  dat.use = dataprep(env_formula = env_formula, bias_formula = bias_formula, 
                     intercept_env = intercept_env, intercept_bias = intercept_bias, 
                     quad_data = quad_data, sp_data = sp_data, 
                     sp_y = sp_y, coord = coord, sp.res = sp_res, dat.type = dat.type, 
                     standardise = standardise, area.int = area.int, r = r, 
                     obs.scores = obs.scores, n.time = n.time)
  
  n_components = dat.use$n_components
  
  if(is.null(n.time)==FALSE){
    colname.bias=no.shrink = all_names = rep(list(c()), n.time)
    for (k in 1:n.time) {
      for (j in 1:n_components[[k]]) {
        colname.bias[[k]][[j]] = colnames(dat.use$X_bias[[k]][[j]])
      }
      no.shrink[[k]] =  c("Intercept", noshrink[[k]])
      all_names[[k]]  = c(colnames(dat.use$X_env[[k]][[1]]), colname.bias[[k]])
      
      no.shrink[[k]]   = which(all_names[[k]] %in% no.shrink[[k]])
    }
    
    noshrink = no.shrink
  }else{
    all_names  = c(colnames(dat.use$X_env[[1]]), unlist(lapply(dat.use$X_bias, colnames)))
    
    noshrink   = c("Intercept", noshrink)
    noshrink   = which(all_names %in% noshrink)
  }
  
  
  if (is.null(wt.vec)){
    wt.vec=NULL
    for (k in 1:n.time) {
      wt.vec[[k]] = rep(1, n_components[[k]])
    }
  }
  
  if(is.null(n.time)==FALSE){         # for a certain time period we will have the same number of of covariates for PO and occ data
    bias_p = bias_ind1 = bias_ind2 = rep(list(c()), n.time)
    env_p = num_p = c()
    env_ix = bias_ix = rep(list(list()), n.time)
    for (k in 1:n.time) {
      bias_p[[k]]    = unlist(lapply(dat.use$X_bias[[k]], function(x) dim(x)[[2]]))
      env_p[k]     = dim(dat.use$X_env[[k]][[1]])[2]
      
      num_p[k]     = env_p[k] + sum(bias_p[[k]])
      bias_ind1[[k]] = env_p[k] + 1 + cumsum(bias_p[[k]]) - bias_p[[k]]
      bias_ind2[[k]] = env_p[k] + cumsum(bias_p[[k]])
      
      
      for (i in 1:n_components[[k]])
      {
        env_ix[[k]][[i]]    = 1:env_p[k]
        bias_ix[[k]][[i]] = bias_ind1[[k]][i]:bias_ind2[[k]][i]
      }
    }
  }else{
    bias_p    = unlist(lapply(dat.use$X_bias, function(x) dim(x)[[2]]))
    env_p     = dim(dat.use$X_env[[1]])[2]
    
    num_p     = env_p + sum(bias_p)
    bias_ind1 = env_p + 1 + cumsum(bias_p) - bias_p
    bias_ind2 = env_p + cumsum(bias_p)
    env_ix    = rep(list(1:env_p), n_components)
    
    bias_ix   = list()
    for (i in 1:n_components)
    {
      bias_ix[[i]] = bias_ind1[i]:bias_ind2[i]
    }
  }
  
  
  if (is.na(init.coef)[1]){
    gamma = 0
  }
  
  if(is.null(n.time)==FALSE){
  # Dynamic path 
    adapt.weights = rep(list(c()), n.time)
    for (i in 1:n.time) {
      if (is.na(init.coef)[1])
      {
        gamma = 0
      }
      adapt.weights[[i]] = if (gamma == 0) rep(1, num_p[[i]]) else 1/abs(init.coef)^gamma
    }
    
    if (is.null(penalty_vec))
    {
      intmodel  = score.int_new(dat.use$X_env, dat.use$X_bias, dat.use$y_list, dat.use$ob_wt, n.time=n.time,
                                dat.type = dat.type, site.area = site.area, lasso = 0, env_ix = env_ix, 
                                bias_ix = bias_ix, coord = coord, wt.vec = wt.vec, noshrink = noshrink,
                                extinct_comp = extinct_comp, coloni_comp = coloni_comp, n_components=n_components)
      score0    = intmodel$score_int
      b.init    = intmodel$newpar
      
      num_startk = c(0, cumsum(num_p))
      num_endk   = cumsum(num_p)

      score0list = new.score = sub.score = pen.max = penalty_vec = list()
      for (i in 1:n.time) {
        score0list[[i]] = score0[(num_startk[[i]]+1):num_endk[[i]]]
        new.score[[i]]  = abs(score0list[[i]]/adapt.weights[[i]])
        sub.score[[i]]  = new.score[[i]][-noshrink[[i]]]
        pen.max[[i]]    = max(sub.score[[i]][is.infinite(sub.score[[i]]) == FALSE])
        
        if (is.na(pen.max[[i]]) == FALSE)
        {
          penalty_vec[[i]] = c(0, exp(seq(log(pen.min), log(pen.max[[i]] + 1.e-5), length.out = (n.fits - 1))))
        }
        
        if (is.na(init.coef[[i]])[1] == FALSE){
          init.coef[[i]][init.coef[[i]] == 0] = b.min
        }
        
      }

    }
    
  }else{
    
    adapt.weights = if (gamma == 0) rep(1, num_p[[1]]) else 1/abs(init.coef)^gamma ## see how it would work when we have different number of variable for the different years
    
    if (is.null(penalty_vec))
    {
      # prepare the data
      
      intmodel  = score.int_new(dat.use$X_env, dat.use$X_bias, dat.use$y_list, dat.use$ob_wt, dat.type = dat.type, site.area = site.area, lasso = 0, env_ix = env_ix, bias_ix = bias_ix, coord = coord, wt.vec = wt.vec, noshrink = noshrink, extinct_comp = extinct_comp, coloni_comp = coloni_comp)
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
    
    }
  
  # common for everything
  if (is.null(n.time)==TRUE) {
    penalty_vec0 = penalty_vec[2]
  }else{
    penalty_vec0 = map(penalty_vec, 2)
  }
    
  fit.0 = combsingle.lasso_new(y=dat.use$y_list, Xenv=dat.use$X_env, Xbias=dat.use$X_bias, 
                               lamb = penalty_vec0, ob.wt = dat.use$ob_wt, 
                               alpha = alpha, b.init = b.init, intercept = NA, family = family, 
                               tol = tol, gamma = gamma, init.coef = init.coef, mu.min = 1.e-16, 
                               mu.max = 1/mu.min, area.int = area.int, pen.max=pen.max,
                               interactions, max.it = 100, standardise = standardise, 
                               noshrink = noshrink, method = method, b.min = b.min, 
                               extinct_comp=extinct_comp, coloni_comp=coloni_comp,
                               link = link, site.area = site.area, dat.type = dat.type, 
                               wt.vec = wt.vec, grad = FALSE, n.time = n.time,
                               n_components=n_components)
  
  
    b.init = fit.0$b
    
    if (is.null(n.time)==TRUE){
      fitpath  = matrix(NA, num_p, length(penalty_vec))
      likepath = rep(NA, length(penalty_vec))
      
      for (i in 1:length(penalty_vec))
      {
        it.max = if (i == 1) 100 else max.it
        
        # function(y, Xenv, Xbias = NULL, lamb, ob.wt = rep(1, length(y)), alpha = 1, b.init = NA, intercept = NA, family = "gaussian", tol = 1.e-9, gamma = 0, init.coef = NA, mu.min = 1.e-16, mu.max = 1/mu.min, area.int = FALSE, interactions, max.it = 25, standardise = TRUE, noshrink = c(1), method = "BFGS", b.min = 1.e-4, link = "logit", site.area = 1, dat.type = "PO", wt.vec = NULL)
        
        fit.i = combsingle.lasso_new(dat.use$y_list, dat.use$X_env, dat.use$X_bias, lamb = penalty_vec[i], ob.wt = dat.use$ob_wt, 
                                     alpha = alpha, b.init = b.init, intercept = NA, family = family, 
                                     tol = tol, gamma = gamma, init.coef = init.coef, mu.min = 1.e-16, mu.max = 1/mu.min, area.int = area.int, 
                                     interactions, max.it = it.max, standardise = standardise, 
                                     noshrink = noshrink, method = method, b.min = b.min, n_components=n_components,
                                     extinct_comp=extinct_comp, coloni_comp=coloni_comp, n.time = n.time,
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
      
    }else{
      
      
        fitpath  = matrix(NA, sum(num_p), length(penalty_vec[[1]]))
        likepath = rep(NA, length(penalty_vec[[1]]))
        
        for (i in 1:length(penalty_vec[[1]]))
        {
          it.max = if (i == 1) 100 else max.it
          
          # function(y, Xenv, Xbias = NULL, lamb, ob.wt = rep(1, length(y)), alpha = 1, b.init = NA, intercept = NA, family = "gaussian", tol = 1.e-9, gamma = 0, init.coef = NA, mu.min = 1.e-16, mu.max = 1/mu.min, area.int = FALSE, interactions, max.it = 25, standardise = TRUE, noshrink = c(1), method = "BFGS", b.min = 1.e-4, link = "logit", site.area = 1, dat.type = "PO", wt.vec = NULL)
          
          fit.i = combsingle.lasso_new(dat.use$y_list, dat.use$X_env, dat.use$X_bias, lamb = map(penalty_vec, i), ob.wt = dat.use$ob_wt, 
                                       alpha = alpha, b.init = b.init, intercept = NA, family = family, 
                                       tol = tol, gamma = gamma, init.coef = init.coef, mu.min = 1.e-16, mu.max = 1/mu.min, area.int = area.int, 
                                       interactions, max.it = it.max, standardise = standardise, n_components=n_components,
                                       noshrink = noshrink, method = method, b.min = b.min, pen.max=pen.max,
                                       extinct_comp=extinct_comp, coloni_comp=coloni_comp, n.time = n.time,
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
        rownames(fitpath) = unlist(all_names)
        end.time = Sys.time()
        runtime = end.time - start.time
        
      }
    
  
    variable_set = abs(fitpath) >= b.min
    n_variables  = apply(variable_set, 2, sum)
    
    AICs = 2*likepath + 2*n_variables
    if(is.null(n.time)==TRUE){
      BICs = 2*likepath + log(sum(dat.use$y_n))*n_variables
      HQCs = 2*likepath + 2*(n_variables + 1)*log(log(sum(dat.use$y_n)))
      AICcs = AICs + 2*(n_variables + 1)*(n_variables + 2)/(sum(dat.use$y_n) - n_variables - 2)
    }else{
      BICs = 2*likepath + log(sum(sapply(dat.use$y_n, sum)))*n_variables
      HQCs = 2*likepath + 2*(n_variables + 1)*log(log(sum(sapply(dat.use$y_n, sum))))
      AICcs = AICs + 2*(n_variables + 1)*(n_variables + 2)/(sum(sapply(dat.use$y_n, sum)) - n_variables - 2)
    }
    
    criterion.matrix        = data.frame(AICs, AICcs, BICs, HQCs)
    names(criterion.matrix) = c("AIC", "AICc", "BIC", "HQC")
    
    meth.id   = paste(criterion, "s", sep = "")
    choice.id = max(which.min(get(meth.id)))
    
    beta.hat   = fitpath[,choice.id]
    like.hat   = likepath[choice.id]
    if(is.null(n.time)==TRUE){
      penalty    = penalty_vec[choice.id]
    }else{
      penalty = map(penalty_vec, choice.id)
    }
    
    
    if (is.null(n.time)==TRUE){
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
    }else{  ## Dynamic
      mu = rep(list(list()), n.time)
      bias = rep(list(list()), n.time)
      psi = rep(list(list()), n.time)
      p_detect = rep(list(list()), n.time)
      
      
      
      # to pass the rbind vector for b.lasso into a list
      num_envk = num_pk= list()
      num_biask= b.lassok = beta.hatk = rep(list(c()), n.time)
      for (i in 1:n.time) {
        num_envk[[i]]  = dim(dat.use$X_env[[i]][[1]])[2]
        num_biask[[i]] = unlist(lapply(lapply(dat.use$X_bias[[i]], as.matrix), ncol))
        num_pk[[i]]    = num_envk[[i]] + sum(num_biask[[i]])
      }
      num_startk = c(0, cumsum(num_pk))
      num_endk   = cumsum(num_pk)
      
      for (i in 1:n.time) {
        beta.hatk[[i]] = beta.hat[(num_startk[[i]]+1):num_endk[[i]]]
      }
      
      lin_det = list()
      for (k in 1:n.time) {
        for (i in 1:n_components[[k]])
        {
          mu[[k]][[i]] = exp(dat.use$X_env[[k]][[i]] %*% beta.hatk[[k]][env_ix[[1]][[1]]])
          psi[[k]][[i]] = 1 - exp(-mu[[k]][[i]]*site.area)
          if (dat.type[i] == "PO")
          {
            bias[[k]][[i]] = exp(dat.use$X_bias[[k]][[i]] %*% beta.hatk[[k]][bias_ix[[k]][[i]]])
            p_detect[[k]][[i]] = list(NULL)
          }
          else
          {
            lin_det[[k]]  = dat.use$X_bias[[k]][[i]] %*% beta.hatk[[k]][bias_ix[[k]][[i]]]
            bias[[i]] = list(NULL)
            if (link == "logit")
            {
              p_detect[[k]][[i]] = expit(lin_det[[k]])
            }
            if (link == "cloglog")
            {
              p_detect[[k]][[i]] = clogloginv(lin_det[[k]])
              #p_detect[[i]] = 1 - exp(-exp(lin_det))
            }
          }
        }
      }
    }
    
  end.time = Sys.time()
  runtime = end.time - start.time
  return(list(formula = formula.out, betas = fitpath, pen_likelihoods = likepath, penalty_vec = penalty_vec, 
              criterion.matrix = criterion.matrix, beta = beta.hat, pen_likelihood = like.hat, 
              penalty = penalty, mu = mu, bias = bias, psi = psi, p_detect = p_detect, pen.max = pen.max,
              X_env = dat.use$X_env, X_bias = dat.use$X_bias, ob_wt = dat.use$ob_wt, y = dat.use$y_list, 
              quad_data = quad_data, sp_data = sp_data, dat.type = dat.type, n.time=n.time,
              quad_xy = dat.use$quad_xy, sp_xy = dat.use$sp_xy, site.area = site.area, runtime = runtime))
}

dataprep = function(env_formula, bias_formula, intercept_env, intercept_bias, quad_data, sp_data, sp_y, coord, sp.res, dat.type, standardise, area.int, r, obs.scores = NULL, n.time=NULL)
{
  ob_wt = rep(list(list()), n.time)
  y_list = rep(list(list()), n.time)
  
  if (class(sp_data) != "list"){
    sp_data = list(sp_data)
  }
  
  #sp_data_xycol = lapply(lapply(sp_data, names), function(x) match(coord, x))
  
  
  if (is.null(n.time) == FALSE){
    ### Dynamic set up: multiple times
    n_components = quad.rep = list()
    y_n =rep(list(c()), n.time)
    sp_data_col = rep(list(list()), n.time)
    
    for (i in 1:n.time) {
      n_components[[i]]  = length(sp_data[[i]])  # should count within a year but set up is different so need to think about it more
      for (l in 1:n_components[[i]]) {
        sp_data_col[[i]][[l]] =  dim(sp_data[[i]][[l]])[2]

        y_n[[i]][l] = dim(sp_data[[i]][[l]])[1]
      }
      
      quad.rep[[i]] = replicate(n_components[[i]], quad_data[[i]], FALSE)
      
    }
    
    
    quad_vars =  rep(list(list()), n.time)
    need.interp = rep(list(list()), n.time)
    
    row_startk = row_endk = sp_rows = rep(list(c()), n.time)
    sp_rows = quad_rowsk = rep(list(list()), n.time)
    
    env_all = sp_data_rb = pres.i = list()
    
    for (i in 1:n.time){
      for (k in 1:n_components[[i]]) {
        if (is.null(quad_data) == FALSE){
          #colnames(quad.rep[[i]][[k]]) = colnames(quad_data[[i]][[k]])
          quad_coord = match(coord, names(quad.rep[[i]][[k]]))
          quad_vars[[i]][[k]]  = names(quad.rep[[i]][[k]][-quad_coord])
          
          need.interp[[i]][[k]] = setdiff(quad_vars[[i]][[k]], names(sp_data[[i]][[k]]))
          
          if (length(need.interp[[i]][[k]]) > 0) # change this to interpolate only needed variables?
          {
            sp_data[[i]][[k]] = newenv.var(sp.xy = sp_data[[i]][[k]], env.scale = sp.res, env.grid = quad.rep[[i]][[k]], coord = coord, file.name = "SpEnvData")
          }
        }
        
        quad.rep[[i]][[k]] = sample.quad(quad.rep[[i]][[k]], sp.res)
        
        if (dat.type[k] == "PO"){
          
          if (is.null(obs.scores[[i]][[k]]) == FALSE){
            ob_wt[[i]][[k]] = scoreweights(sp_data[[i]][[k]], quad.rep[[i]][[k]], obs.scores = obs.scores[[i]][[k]])
          }else(
            ob_wt[[i]][[k]] = weights(sp_data[[i]][[k]], quad.rep[[i]][[k]])
          )
          
          pres.i[[i]]     = rep(0, length(ob_wt[[i]][[k]]))
          pres.i[[i]][1:length(sp_y[[i]][[k]])] = sp_y[[i]][[k]]
          y_list[[i]][[k]] = pres.i[[i]]/ob_wt[[i]][[k]]
          
          # remove the null weights created from the list of components (occ data) to keep only the weights for PO
          #ob_wt[c(dat.type == "PO")]
        }
        else{
          ob_wt[[i]][[k]] = list(NULL)
          y_list[[i]][[k]] = sp_y[[i]][[k]]
        }
        
        #X_env = list()
        #X_bias = list()
        
      }
      
      
      # Combine PO and OCC for a certain year
      sp_data_rb[[i]] = do.call("rbind", sp_data[[i]])
     
      # Combine a certain year data (PO + OCC) with sampled quad
      env_all = Map(rbind, sp_data_rb, quad.rep[[1]])
      
      #env_all = if (any(dat.type[i] == "PO") == TRUE) rbind(do.call("rbind", sp_data[[i]]), quad.rep[[i]]) else do.call("rbind", sp_data[[i]])
      
      
      row_startk[[i]] = c(1, cumsum(y_n[[i]]) + 1)[1:n_components[[i]]]
      row_endk[[i]]   = cumsum(y_n[[i]])
      
      if (n_components[[i]] == 1){
        sp_rows[[i]]   = list(row_startk[[i]]:row_endk[[i]])
      }else{
        sp_rows[[i]] = mapply(function(x, y) list(x:y), row_startk[[i]], row_endk[[i]])
      }
       
      #quad_rows = (max(row_end) + 1):(dim(env_all)[1])  
      
      quad_rowsk[[i]] = (max(row_endk[[i]]) + 1):dim(env_all[[i]])[1]
      
      for (j in 1:n_components[[i]])
      {
        if (dat.type[j] == "PO")
        {
          sp_rows[[i]][[j]] = c(sp_rows[[i]][[j]], quad_rowsk[[i]])
        }
      }
      
    }
    
    mf_env = mt_env = list()
    X_env_all = list()
    
    # same variables for different datasets of the same time period
    for (k in 1:n.time) {
      mf_env[[k]]     = model.frame(env_formula, env_all[[k]])
      mt_env[[k]]     = attr(mf_env[[k]], "terms")
      
      attr(mt_env[[k]], "intercept") = 0
      X_env_all[[k]] = if (!is.empty.model(mt_env[[k]])) model.matrix(mt_env[[k]], mf_env[[k]], contrasts)else matrix(, NROW(X_env_all[[k]]), 0L)
    }
    
    data.i = rep(list(list()), n.time)
    bf.i = bt.i = list()
    X_env = X_bias= rep(list(list()), n.time)
    
    if (is.null(bias_formula) == FALSE)
    {
      for (i in 1:n.time) {
        for (j in 1:n_components[[i]]) {
          if (dat.type[j] == "PO")
          {
            data.i[[i]][[j]]      = rbind(sp_data[[i]][[j]], quad.rep[[i]][[j]])
          
            bf.i[[i]] = model.frame(bias_formula[[i]][[j]], data.i[[i]][[j]])
            bt.i[[i]] = attr(bf.i[[i]], "terms")
            
            attr(model.frame(bias_formula[[1]][[1]], data.i[[1]][[1]]), "terms")
            
            attr(bt.i[[i]], "intercept") = 0
            X_bias[[i]][[j]] = if (!is.empty.model(bt.i[[i]])) model.matrix(bt.i[[i]], bf.i[[i]], contrasts)
            else matrix(, NROW(sp_data[[i]][[j]]), 0L)
              
            
          }
          else
          {
            bf.i[[i]] = model.frame(bias_formula[[i]][[j]], sp_data[[i]][[j]])
            bt.i[[i]] = attr(bf.i[[i]], "terms")
            
            attr(bt.i[[i]], "intercept") = 0
            X_bias[[i]][[j]] = if (!is.empty.model(bt.i[[i]])) model.matrix(bt.i[[i]], bf.i[[i]], contrasts)
            else matrix(, NROW(sp_data[[i]][[j]]), 0L)
            
            
          }
        }
      }
      
    }
    
    
    if (standardise == TRUE)
    {
      X_env_all_stand = s_means_env = s_sds_env = X_env_s = X_env = rep(list(list()), n.time)
      X_env_all_stand = Map(standardise.X, X_env_all)
      
      for (k in 1:n.time) {
        s_means_env[[k]] = X_env_all_stand[[k]]$dat.means
        s_sds_env[[k]]   = X_env_all_stand[[k]]$dat.sds
        X_env_s[[k]]     = X_env_all_stand[[k]]$X
        X_env[[k]] = lapply(sp_rows[[k]], function(x) X_env_s[[k]][x,])
      }
      
      X_bias_stand_i = rep(list(list()), n.time)
      s_means_bias = s_sds_bias = rep(list(list()), n.time)
      s_means.bias = s_sds.bias = s_means = s_sds = list()
      
      if (is.null(bias_formula) == FALSE)
      {
        for (i in 1:n.time)
        {
          for (k in 1:n_components[[i]]) {
            X_bias_stand_i[[i]][[k]] = standardise.X(X_bias[[i]][[k]])
            s_means_bias[[i]][[k]]   = X_bias_stand_i[[i]][[k]]$dat.means
            s_sds_bias[[i]][[k]]  =  X_bias_stand_i[[i]][[k]]$dat.sds
            X_bias[[i]][[k]]    = X_bias_stand_i[[i]][[k]]$X 
          }
          
          s_means.bias[[k]] = as.data.frame(do.call(cbind, s_means_bias))[k,]
          s_sds.bias[[k]] = as.data.frame(do.call(cbind, s_sds_bias))[k,]
          
          s_means[[k]] = c(s_means_env[[k]], s_means.bias[[k]], recursive=TRUE)
          s_sds[[k]]   = c(s_sds_env[[k]], s_sds_bias[[k]], recursive=TRUE)
        }
        
      }
      
    }
    
    
    # add intercepts
    
    #intercept_env = 1
    #intercept_bias = list(NA, NA, 1)
    
    intercept_add = list()
    
    for (k in 1:n.time){
      
        if (is.null(intercept_env[[k]]) == FALSE)
        {
          X_env[[k]] = lapply(X_env[[k]], function(x) cbind(Intercept = 1, x))
          #X_env[[k]] = Map(cbind, Intercept = 1, X_env[[k]])
        }
      for (i in 1:n_components[[k]]) { 
        if (is.null(intercept_bias[[i]]) == FALSE)
        {
          intercept_add[[i]] = which(is.na(intercept_bias[[i]]) == FALSE)
          for (add.i in intercept_add[[i]])
          {
            X_bias[[i]][[add.i]] = cbind(Intercept = 1, X_bias[[i]][[add.i[[i]]]])
          }
        }
      }
    }
    
    
    X_env = lapply(X_env, polynamesadd)
    X_bias = lapply(X_bias, polynamesadd)
    
    # Try to get the correct row numbers
    # for (i in 1:n.time) {
    #     for (j in 1:n_components[[i]]) {
    #       if (dat.type[j] == "PO"){
    #         row.names(X_env[[i]][[j]]) = c(row.names(sp_data[[j]][[i]]),
    #                                      row.names(quad_data[[1]]))
    #         row.names(X_bias[[j]][[i]]) = c(row.names(sp_data[[j]][[i]]),
    #                                      row.names(quad_data[[1]]))
    #       }
    #     }
    #   
    # }
    
    int_i = rep(list(list()), n.time)
    if (area.int == TRUE)
    {
      for (i in 1:n.time) {
        for (k in 1:n_components[[k]])
        {
          if (dat.type[i] == "PO" & is.na(r[i]) == FALSE)
          {
            int_i[[i]][[k]] = point.interactions(env = quad.rep, pres = sp_data[[i]][[k]], r = r[i], coord = coord)
            if (standardise == TRUE)
            {
              int_i[[i]][[k]] = scale(int_i[[i]][[k]])
            }
            X_bias[[i]][[k]] = cbind(X_bias[[i]][[k]], Interaction = int_i[[i]][[k]])
            dimnames(X_bias[[i]][[k]])[[2]][dim(X_bias[[i]][[k]])[2]] = "Interaction"
          }
        }
      }
      
    }
    
    sp_xy = quad_xy = rep(list(list()), n.time)
    
    for (i in 1:n.time) {
      for (k in 1:n_components[[k]]) {
        sp_xy[[i]] = lapply(sp_data[[i]], function(x) {x[,match(coord, names(x))]})
        
        quad_xy[[i]][[k]] = quad.rep[[i]][[k]][,match(coord, names(quad.rep[[i]][[k]]))]  
      }
    }
    
    
    quad_data = quad.rep
    
    # X_env need to be re arrange to be in the form list [[1]] [[1]]po1, [[2]] po2,.. list [[2]] [[1]] Occ1, [[2]] Occ2 ..:
    # X_env.d = X_env.rea = list() 
    # 
    # for (i in 1:n.time) {
    #   for (k in 1:n_components[[k]]) {
    #     X_env.d[[i]] = X_env[[i]][[k]]
    #   }
    #   X_env.rea[[k]] = X_env.d
    #   
    # }
    # 
    # X_env = X_env.rea
    
    return(list(quad_data = quad_data, sp_data = sp_data, y_list = y_list, ob_wt = ob_wt, y_n = y_n, n_components = n_components, X_env = X_env, X_bias = X_bias, quad_xy = quad_xy, sp_xy = sp_xy))
    
    
    
  }else{
    
    ### Not dynamic set up
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
        
        if (is.null(obs.scores[[i]]) == FALSE){
          ob_wt[[i]] = scoreweights(sp_data[[i]], quad_data, obs.scores = obs.scores[[i]])
        }
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
    X_env_all = if (!is.empty.model(mt_env)) model.matrix(mt_env, mf_env, contrasts)
    else matrix(, NROW(X_env_all), 0L)
    
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
    return(list(quad_data = quad_data, sp_data = sp_data, y_list = y_list, ob_wt = ob_wt, y_n = y_n, n_components = n_components, X_env = X_env, X_bias = X_bias, quad_xy = quad_xy, sp_xy = sp_xy, n.time = n.time))
    
  }
  
  
}


scoreweights = function(sp.xy, quad.xy, coord = c("X", "Y"), obs.scores = NULL)
{
  if (is.null(obs.scores)){
    score.all = rep(1, (dim(sp.xy)[1]) + dim(quad.xy)[1])
  }else{
    score.all = c(obs.scores, rep(1, dim(quad.xy)[1]))
  }
  
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
  round.tab   = aggregate(data.frame(score.all), list(ID = round.id), sum)
  #round.table = table(round.id)
  #wt          = X.inc*Y.inc/as.numeric(round.table[match(round.id, names(round.table))])
  scorewt     = X.inc*Y.inc*score.all/round.tab$score.all[match(round.id, round.tab$ID)]
  scorewt
}

combsingle.lasso_new = function(y, Xenv, Xbias = NULL, lamb, ob.wt = ob.wt, alpha = 1, b.init = NA, intercept = NA, family = "gaussian", tol = 1.e-9, gamma = 0, init.coef = NA, mu.min = 1.e-16, mu.max = 1/mu.min, area.int = FALSE, interactions, max.it = 25, standardise = TRUE, noshrink = c(1), method = "BFGS", b.min = 1.e-4, link = "logit", site.area = 1, dat.type = "PO", wt.vec = NULL, grad = TRUE, extinct_comp=c(NULL), coloni_comp=c(NULL), n.time=NULL, pen.max=pen.max, n_components=NULL)
{
  # To do:
  # Fix behvaior with AI models
  # Other estimation procedures (IRLS, etc.?)
  
  if (class(y) == "list") # determine individual data sources and corresponding indices
  {
    
    if (is.null(n.time)==FALSE){
      
      num_env = num_p=n.components= bias_ind1= bias_ind2=env_ix=X = list()
      num_bias= rep(list(c()), n.time)
      ix = bias_ix = rep(list(list()), n.time)
      for (i in 1:n.time) {
        num_env[[i]]  = dim(Xenv[[i]][[1]])[2]
        num_bias[[i]] = unlist(lapply(lapply(Xbias[[i]], as.matrix), ncol))
        num_p[[i]]    = num_env[[i]] + sum(num_bias[[i]])
        n.components[[i]] = length(num_bias[[i]])
        bias_ind1[[i]] = num_env[[i]] + 1 + cumsum(num_bias[[i]]) - num_bias[[i]]
        bias_ind2[[i]] = num_env[[i]] + cumsum(num_bias[[i]])
        env_ix[[i]] = rep(list(1:num_env[[i]]), n.components[[i]])
        
        
        for (j in 1:n.components[[i]])
        {
          ix[[i]][[j]] = c(1:num_env[[i]], bias_ind1[[i]][j]:bias_ind2[[i]][j])
          bias_ix[[i]][[j]] = bias_ind1[[i]][j]:bias_ind2[[i]][j]
        }
        X[[i]] = Map(cbind, Xenv[[i]], Xbias[[i]]) 
      }
      
      if (is.null(wt.vec))
      {
        wt.vec = rep(list(c()), n.time)
        for (i in 1:n.time) {
          wt.vec[[i]] = rep(1, n.components[[i]])
        }
      }
      
    }else{ ##
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
    
    if (is.null(wt.vec))
    {
      wt.vec = rep(1, n.components)
    }
    
  }else{
    env_cov  = 1:dim(Xenv)[2]
    bias_cov = if (is.null(Xbias) == FALSE) (max(env_cov) + 1):(max(env_cov) + dim(Xbias)[2])
    else c()
    
    #bias_cov = (max(env_cov) + 1):(max(env_cov) + dim(Xbias)[2])
    
    X = cbind(Xenv, Xbias)
    num_p = length(c(env_cov, bias_cov))
  }
  
  
  error.flag = FALSE
  
  if (is.null(n.time)==FALSE){
    b.glm   = rep(list(c()), n.time)
    for (i in 1:n.time) {
      b.glm[[i]]   = rep(NA, num_p[[i]])
    }
    
    b.lasso = b.glm
    
    adapt.weights = lambda = killed = keep = ix.keep = list()
    lambda.start = rep(list(c()), n.time)
    lasso = rep(list(array()), n.time)
    lass.pen = lasso.store = list()
    
    for (i in 1:n.time) {
      adapt.weights[[i]] = (if(gamma == 0) rep(1, num_p[[i]]) # check area int
                       else 1/abs(init.coef)^gamma)
      
      lambda.start[[i]] = (if(length(lamb[[i]]) == 1) rep(lamb[[i]], num_p[[i]])
                           else rep(0, num_p[[i]]))
        
      lambda.start[[i]][noshrink[[i]]] = 0
      lambda[[i]] = as.array(lambda.start[[i]])
      
      lambda[[i]] = lambda[[i]] * abs(adapt.weights[[i]][1:length(lambda[[i]])])
      killed[[i]]  = is.infinite(lambda[[i]])
      keep[[i]] = which(killed[[i]] == FALSE)
      ix.keep[[i]] = lapply(ix[[i]], function(x) killed[[i]][x] == FALSE)
      
      #Xenv[[i]]  = subsetlist(X_env.rea[[i]], env_ix[[i]], keep[[i]])
      #Xbias[[i]] = subsetlist(X_bias.rea[[i]], bias_ix[[i]], keep[[i]])
      #b.lasso[[i]] = b.lasso[[i]][killed[[i]] == FALSE]
      #lambda[[i]]  = lambda[[i]][killed[[i]] == FALSE]
      
      # No subset method - have infinite for coeff lasso instead
      #Xenv[[i]]  = X_env[[i]]
      #Xbias[[i]] = X_bias[[i]]
      
      lass.pen[[i]] = pen.max[[i]]*1000000
      
      lasso.store[[i]] = lambda[[i]][killed[[i]] == TRUE]
      lambda[[i]][killed[[i]] == TRUE]  = lass.pen[[i]]
      
    }
 
  }else{
    adapt.weights = (if(gamma == 0) rep(1, num_p) # check area int
                     else 1/abs(init.coef)^gamma)
    
    b.glm   = rep(NA, num_p)
    b.lasso = b.glm
    
    lambda.start = (if(length(lamb) == 1) rep(lamb, num_p)
                    else rep(0, num_p))
    lambda.start[noshrink] = 0
    lambda = as.array(lambda.start)
    lambda = lambda * abs(adapt.weights[1:length(lambda)])
    
    lambda = lambda * abs(adapt.weights[1:length(lambda)])
    
    killed  = is.infinite(lambda)
    keep = which(killed == FALSE)
    
    ix.keep = lapply(ix, function(x) killed[x] == FALSE)
    
    Xenv  = subsetlist(Xenv, env_ix, keep)
    Xbias = subsetlist(Xbias, bias_ix, keep)
    
    b.lasso = b.lasso[killed == FALSE]
    lambda  = lambda[killed == FALSE]
  }

  
  #,,,,,,,,,,,,,,,,,,,,,, start from here

  
  if (any(is.na(b.init)) & is.na(intercept)) ## dynamic too
  {
    if (is.null(n.time)==TRUE){
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
    }else{
      
      num_env = num_p=n.components= bias_ind1= bias_ind2=env_ix=X = list()
      num_bias= rep(list(c()), n.time)
      ix = bias_ix = rep(list(list()), n.time)
      for (i in 1:n.time) {
        num_env[[i]]  = dim(Xenv[[i]][[1]])[2]
        num_bias[[i]] = unlist(lapply(lapply(Xbias[[i]], as.matrix), ncol))
        num_p[[i]]    = num_env[[i]] + sum(num_bias[[i]])
        n.components[[i]] = length(num_bias[[i]])
        bias_ind1[[i]] = num_env[[i]] + 1 + cumsum(num_bias[[i]]) - num_bias[[i]]
        bias_ind2[[i]] = num_env[[i]] + cumsum(num_bias[[i]])
        env_ix[[i]] = rep(list(1:num_env[[i]]), n.components[[i]])
        
        
        for (j in 1:n.components[[i]])
        {
          ix[[i]][[j]] = c(1:num_env[[i]], bias_ind1[[i]][j]:bias_ind2[[i]][j])
          bias_ix[[i]][[j]] = bias_ind1[[i]][j]:bias_ind2[[i]][j]
        }
      }
      
      
      # Keep only certain variables
      X_env.use=X_bias.use=list()
      for (i in 1:n.time) {
        is.in[[i]] = which(killed[[i]]==FALSE)
        
        # X_env.use[[i]]  = subsetlist(Xenv[[i]], env_ix[[i]], keep = is.in[[i]])
        # X_bias.use[[i]] = subsetlist(Xbias[[i]], bias_ix[[i]], keep = is.in[[i]])
        X_env.use[[i]]  = Xenv[[i]]
        X_bias.use[[i]] = Xbias[[i]]
        
      }
      
      y_po= X_env_po= X_bias_po= ob_wt_po=y_occ2=X_env_occ2=X_bias_occ2=rep(list(list()), n.time)
      y_occ=X_env_occ=X_bias_occ=list()
      
      for (t in 1:n.time) {
        for (i in 1:n_components[[t]]) {
          if (dat.type[i] == "PO") {
            y_po[[t]] = y[[t]][[i]]
            X_env_po[[t]] = Xenv[[t]][[i]]
            X_bias_po[[t]] = Xbias[[t]][[i]]
            ob_wt_po[[t]] = ob.wt[[t]][[i]]
          }else{
            y_occ2[[t]] = y[[t]][[i]] 
            X_env_occ2[[t]] = Xenv[[t]][[i]]
            X_bias_occ2[[t]] = Xbias[[t]][[i]]  # numeric for each occ year so need as.matrix
            #ob_wt.occ is NULL for occ
          }
          
        }
      }
      
      y_occ = abind(y_occ2, along=3)
      X_env_occ = abind(X_env_occ2, along=3)
      X_bias_occ = abind(X_bias_occ2, along=3)
      
      
      for (i in 1:n.time) {
        lass.pen[[i]] = pen.max[[i]]*1000000
        
        lasso.store[[i]] = lambda[[i]][killed[[i]] == TRUE]
        lambda[[i]][killed[[i]] == TRUE]  = lass.pen[[i]]
      }
      
      # Need to prepare the data again
      par=rep(list(rep(0, sum(killed[[1]] == FALSE))), n.time) # Need to be list first # can we have an issue if we have different dimension here?
      Oprep_TMB = TMB.optimprep(par=par, yocc=y_occ, Wocc=X_bias_occ, Xocc=X_env_occ, ypo=y_po, 
                                extinct_comp = extinct_comp, coloni_comp = coloni_comp, wt_vec = wt.vec,
                                Wpo=X_bias_po, Xpo=X_env_po, ob_wt.po=ob_wt_po, lasso=lambda, n.time=n.time,
                                site_area = site.area, E=length(extinct_comp), C=length(coloni_comp), dat.type)
      
      BEstartk = c(0, dim(Oprep_TMB$Xpobmat)[2])
      BExstartk = c(0, length(extinct_comp))
      Bcostartk = c(0, length(coloni_comp))
      BbiPOstartk = c(0, Oprep_TMB$Wocc_NB[1]) # doesn't work if we have different bias covar throught the years
      BbiOccstartk = c(0, dim(Oprep_TMB$Wpomat)[2])
      
      f1comb <- MakeADFun(data=list(K=Oprep_TMB$K, Xocc=X_env_occ, Wocc=Oprep_TMB$Wocc, site_area=site.area, 
                                    yocc=Oprep_TMB$yocc, R=Oprep_TMB$R, Jm=Oprep_TMB$Jm, S_NB=Oprep_TMB$S_NB, 
                                    Vocc_NB=Oprep_TMB$Vocc_NB, lasso_occ=Oprep_TMB$lasso_occ, extinctcomp=extinct_comp, 
                                    colonicomp=coloni_comp, Xpobmat=Oprep_TMB$Xpobmat, Xpoemat=Oprep_TMB$Xpoemat, 
                                    Xponmat=Oprep_TMB$Xponmat, Wpomat=Oprep_TMB$Wpomat, ypomat=Oprep_TMB$ypomat, 
                                    Wpos=Oprep_TMB$Wpos, Vpo_NB=Oprep_TMB$Vpo_NB, L_NB=Oprep_TMB$L_NB,
                                    row_startk=Oprep_TMB$row_startk, lassomat=Oprep_TMB$lassomat, lassostartk=Oprep_TMB$lassostartk,
                                    Lsize=Oprep_TMB$Lsize, ob_wtmat=Oprep_TMB$ob_wtmat, Wocc_NB=Oprep_TMB$Wocc_NB,
                                    BbiPOstartk=BbiPOstartk, BbiOccstartk=BbiOccstartk, 
                                    wtvecocc = Oprep_TMB$wtvecocc, wtvecpo = Oprep_TMB$wtvecpo, 
                                    BExstartk=BExstartk, Bcostartk=Bcostartk, BEstartk = BEstartk), 
                          parameters= list(b_init=Oprep_TMB$b_init), DLL = "TMB_comb")
      
      ## combine at each iteration
      f1comb.res <- optim(par=f1comb$par, fn=f1comb$fn, gr=f1comb$gr, method='BFGS', 
                          control = list(ndeps = rep(b.min, num_p[[1]]*n.time))) # control to allow very small coef to be 0
      
      b.lasso = f1comb.res$par
      like1 = f1comb.res$value
      
    }
    
    
  }
  
  if (any(is.na(b.init)) == FALSE) ## dynamic too
  {

    if (is.null(n.time)==TRUE){
      b.lasso = b.init[killed == FALSE]
      
      like1   = LL_Lasso_Wt_Comb_new(par = b.lasso, Xenv = Xenv, Xbias = Xbias, y = y, ob.wt = ob.wt,
                                   env_ix = env_ix, bias_ix = bias_ix, site.area = site.area,
                                   lasso = lambda[killed == FALSE], wt.vec = wt.vec, dat.type = dat.type, is.in = killed == FALSE)
    }else{
      
      # num_env = num_p=n.components= bias_ind1= bias_ind2=env_ix=X = list()
      # num_bias= rep(list(c()), n.time)
      # ix = bias_ix = rep(list(list()), n.time)
      # for (i in 1:n.time) {
      #   num_env[[i]]  = dim(Xenv[[i]][[1]])[2]
      #   num_bias[[i]] = unlist(lapply(lapply(Xbias[[i]], as.matrix), ncol))
      #   num_p[[i]]    = num_env[[i]] + sum(num_bias[[i]])
      #   n.components[[i]] = length(num_bias[[i]])
      #   bias_ind1[[i]] = num_env[[i]] + 1 + cumsum(num_bias[[i]]) - num_bias[[i]]
      #   bias_ind2[[i]] = num_env[[i]] + cumsum(num_bias[[i]])
      #   env_ix[[i]] = rep(list(1:num_env[[i]]), n.components[[i]])
      #   
      #   
      #   for (j in 1:n.components[[i]])
      #   {
      #     ix[[i]][[j]] = c(1:num_env[[i]], bias_ind1[[i]][j]:bias_ind2[[i]][j])
      #     bias_ix[[i]][[j]] = bias_ind1[[i]][j]:bias_ind2[[i]][j]
      #   }
      # }
      
      #killed.v = do.call(c, killed)
      #b.lasso = b.init[killed.v == FALSE]
      # Keep only certain variables
      X_env.use=X_bias.use=list()
      is.in=list()
      for (i in 1:n.time) {
        is.in[[i]] = killed[[i]]==FALSE
        
        #X_env.use[[i]]  = subsetlist(Xenv[[i]], env_ix[[i]], keep = which(is.in[[i]] == TRUE))
        #X_bias.use[[i]] = subsetlist(Xbias[[i]], bias_ix[[i]], keep = which(is.in[[i]] == TRUE))
        X_env.use[[i]]  = Xenv[[i]]
        X_bias.use[[i]] = Xbias[[i]]
        
      }
      
      # to pass the rbind vector for b.init into a list
      num_env = num_p= list()
      num_bias= b.lassokeep = rep(list(c()), n.time)
      for (i in 1:n.time) {
        num_env[[i]]  = dim(X_env.use[[i]][[1]])[2]
        num_bias[[i]] = unlist(lapply(lapply(X_bias.use[[i]], as.matrix), ncol))
        num_p[[i]]    = num_env[[i]] + sum(num_bias[[i]])
      }
      num_startk = c(0, cumsum(num_p))
      num_endk   = cumsum(num_p)
      
      if (class(b.init)=="list"){
        for (i in 1:n.time){
          #b.initlist[[i]] = b.init[(num_startk[[i]]+1):num_endk[[i]]]
          # subset for the one we keep
          b.lassokeep[[i]] = b.init[[i]]
          
        }
      }else{
        b.initlist = rep(list(c()), n.time)
        for (i in 1:n.time) {
          b.initlist[[i]] = b.init[(num_startk[[i]]+1):num_endk[[i]]]
          b.lassokeep[[i]] = b.initlist[[i]]
        }
      }
      

      y_po= X_env_po= X_bias_po= ob_wt_po=y_occ2=X_env_occ2=X_bias_occ2=rep(list(list()), n.time)
      y_occ=X_env_occ=X_bias_occ=list()
      
      for (t in 1:n.time) {
        for (i in 1:n_components[[t]]) {
          if (dat.type[i] == "PO") {
            y_po[[t]] = y[[t]][[i]]
            X_env_po[[t]] = Xenv[[t]][[i]]
            X_bias_po[[t]] = Xbias[[t]][[i]]
            ob_wt_po[[t]] = ob.wt[[t]][[i]]
          }else{
            y_occ2[[t]] = y[[t]][[i]] 
            X_env_occ2[[t]] = Xenv[[t]][[i]]
            X_bias_occ2[[t]] = Xbias[[t]][[i]]  # numeric for each occ year so need as.matrix
            #ob_wt.occ is NULL for occ
          }
          
        }
      }
      
      y_occ = abind(y_occ2, along=3)
      X_env_occ = abind(X_env_occ2, along=3)
      X_bias_occ = abind(X_bias_occ2, along=3)
      
      
      for (i in 1:n.time) {
        lass.pen[[i]] = pen.max[[i]]*1000000
        
        lasso.store[[i]] = lambda[[i]][killed[[i]] == TRUE]
        lambda[[i]][killed[[i]] == TRUE]  = lass.pen[[i]]
      }
      # Need to prepare the data again
      Oprep_TMB1 = TMB.optimprep(par=b.lassokeep, yocc=y_occ, Wocc=X_bias_occ, Xocc=X_env_occ, ypo=y_po, 
                                 extinct_comp = extinct_comp, coloni_comp = coloni_comp, wt_vec = wt.vec,
                                 Wpo=X_bias_po, Xpo=X_env_po, ob_wt.po=ob_wt_po, lasso=lambda, n.time=n.time,
                                 site_area = site.area, E=length(extinct_comp), C=length(coloni_comp), dat.type)
      
      BEstartk = c(0, dim(Oprep_TMB1$Xpobmat)[2])
      BExstartk = c(0, length(extinct_comp))
      Bcostartk = c(0, length(coloni_comp))
      BbiPOstartk = c(0, Oprep_TMB1$Wocc_NB[1]) # doesn't work if we have different bias covar throught the years
      BbiOccstartk = c(0, dim(Oprep_TMB1$Wpomat)[2])
      
      f1comb <- MakeADFun(data=list(K=Oprep_TMB1$K, Xocc=X_env_occ, Wocc=Oprep_TMB1$Wocc, site_area=site.area, 
                                    yocc=Oprep_TMB1$yocc, R=Oprep_TMB1$R, Jm=Oprep_TMB1$Jm, S_NB=Oprep_TMB1$S_NB, 
                                    Vocc_NB=Oprep_TMB1$Vocc_NB, lasso_occ=Oprep_TMB1$lasso_occ, extinctcomp=extinct_comp, 
                                    colonicomp=coloni_comp, Xpobmat=Oprep_TMB1$Xpobmat, Xpoemat=Oprep_TMB1$Xpoemat, 
                                    Xponmat=Oprep_TMB1$Xponmat, Wpomat=Oprep_TMB1$Wpomat, ypomat=Oprep_TMB1$ypomat, 
                                    Wpos=Oprep_TMB1$Wpos, Vpo_NB=Oprep_TMB1$Vpo_NB, L_NB=Oprep_TMB1$L_NB,
                                    row_startk=Oprep_TMB1$row_startk, lassomat=Oprep_TMB1$lassomat, lassostartk=Oprep_TMB1$lassostartk,
                                    Lsize=Oprep_TMB1$Lsize, ob_wtmat=Oprep_TMB1$ob_wtmat, Wocc_NB=Oprep_TMB1$Wocc_NB,
                                    BbiPOstartk=BbiPOstartk, BbiOccstartk=BbiOccstartk, 
                                    wtvecocc = Oprep_TMB1$wtvecocc, wtvecpo = Oprep_TMB1$wtvecpo,
                                    BExstartk=BExstartk, Bcostartk=Bcostartk, BEstartk = BEstartk), 
                          parameters= list(b_init=Oprep_TMB1$b_init), DLL = "TMB_comb")
      
      #No optim needed here
      b.lasso = f1comb$par
      like1 = f1comb$report()$comblik
      
    }
    
  }
  
  if (is.na(intercept) == FALSE) ## not dynamic - ask Ian
  {
    if (family$family == "poisson")
    {
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
  
  is.in = abs(b.lasso) > b.min
  
  if (is.null(n.time)==TRUE){   
    is.in[noshrink] = TRUE
    b.lasso[is.in == FALSE] = 0
    sign.change  = rep(-1, num_p)
    
  }else{ # Dynamic

    is.in.list =  list()
    for (i in 1:n.time) {
      is.in.list[[i]] = is.in[(num_startk[[i]]+1):num_endk[[i]]]
      is.in.list[[i]][noshrink[[i]]] = TRUE
    }
    is.in=do.call(c, is.in.list)
    b.lasso[is.in == FALSE] = 0
    
    sign.change = list()
    for (i in 1:n.time) {
      sign.change[[i]] = rep(-1, num_p[[i]])
    }

  }
  
  is.different = is.in
  diff   = 2*tol
  num.it = 0
  
  #++++++++++++++++++++++++++++++++++++++++++++ Start working from here
  
  signs = sign(b.lasso) # for dynamic set up each year is in the same vector following each other
  
  betas   = c(b.lasso)
  scores  = c()
  viols   = c()
  likes   = c()
  actions = c()
  varsets = c()
  
  likes   = c(likes, like1)
  
  while(diff > tol & num.it < max.it)
  {
    setchange = 0
    
    if (is.null(n.time)==TRUE){
      b.est = optim(par = b.lasso, fn = LL_Lasso_Wt_Comb_new, gr = if (grad == TRUE) scoreeq_new else NULL, 
                    Xenv = Xenv, 
                    Xbias = Xbias, 
                    y = y, ob.wt = ob.wt, link = link, coord = coord,
                    env_ix = env_ix,
                    bias_ix = bias_ix,
                    site.area = site.area, 
                    lasso = lambda, wt.vec = wt.vec, dat.type = dat.type, is.in = is.in,
                    method = method, control = list(ndeps = rep(b.min, num_p)))
      
      
    }else{
      # set up the variable we keep if it is dynamic process because the function doesn't do it directly
      
      # num_env = num_p=n.components= bias_ind1= bias_ind2=env_ix=X = list()
      # num_bias= rep(list(c()), n.time)
      # ix = bias_ix = rep(list(list()), n.time)
      # for (i in 1:n.time) {
      #   num_env[[i]]  = dim(Xenv[[i]][[1]])[2]
      #   num_bias[[i]] = unlist(lapply(lapply(Xbias[[i]], as.matrix), ncol))
      #   num_p[[i]]    = num_env[[i]] + sum(num_bias[[i]])
      #   n.components[[i]] = length(num_bias[[i]])
      #   bias_ind1[[i]] = num_env[[i]] + 1 + cumsum(num_bias[[i]]) - num_bias[[i]]
      #   bias_ind2[[i]] = num_env[[i]] + cumsum(num_bias[[i]])
      #   env_ix[[i]] = rep(list(1:num_env[[i]]), n.components[[i]])
      #   
      #   
      #   for (j in 1:n.components[[i]])
      #   {
      #     ix[[i]][[j]] = c(1:num_env[[i]], bias_ind1[[i]][j]:bias_ind2[[i]][j])
      #     bias_ix[[i]][[j]] = bias_ind1[[i]][j]:bias_ind2[[i]][j]
      #   }
      # }
      ## Need to transform data and subset for what we keep
      X_env.usek=X_bias.usek=list()
      for (i in 1:n.time) {
        #X_env.usek[[i]]  = subsetlist(Xenv[[i]], env_ix[[i]], keep = which(is.in.list[[i]] == TRUE))
        #X_bias.usek[[i]] = subsetlist(Xbias[[i]], bias_ix[[i]], keep = which(is.in.list[[i]] == TRUE))
        X_env.usek[[i]]  = Xenv[[i]]
        X_bias.usek[[i]] = Xbias[[i]]
        
      }
      
      # to pass the rbind vector for b.lasso into a list
      num_envk = num_pk= list()
      num_biask= b.lassokeep1 = rep(list(c()), n.time)
      for (i in 1:n.time) {
        num_envk[[i]]  = dim(X_env.use[[i]][[1]])[2]
        num_biask[[i]] = unlist(lapply(lapply(X_bias.use[[i]], as.matrix), ncol))
        num_pk[[i]]    = num_env[[i]] + sum(num_bias[[i]])
      }
      num_startk = c(0, cumsum(num_pk))
      num_endk   = cumsum(num_pk)
      
      for (i in 1:n.time) {
        b.lassokeep1[[i]] = b.lasso[(num_startk[[i]]+1):num_endk[[i]]]
      }
      
      
      y_po= X_env_po= X_bias_po= ob_wt_po=y_occ2=X_env_occ2=X_bias_occ2=rep(list(list()), n.time)
      y_occ=X_env_occ=X_bias_occ=list()
      
      for (t in 1:n.time) {
        for (i in 1:n_components[[t]]) {
          if (dat.type[i] == "PO") {
            y_po[[t]] = y[[t]][[i]]
            X_env_po[[t]] = Xenv[[t]][[i]]
            X_bias_po[[t]] = Xbias[[t]][[i]]
            ob_wt_po[[t]] = ob.wt[[t]][[i]]
          }else{
            y_occ2[[t]] = y[[t]][[i]] 
            X_env_occ2[[t]] = Xenv[[t]][[i]]
            X_bias_occ2[[t]] = Xbias[[t]][[i]]  # numeric for each occ year so need as.matrix
            #ob_wt.occ is NULL for occ
          }
          
        }
      }
      
      y_occ = abind(y_occ2, along=3)
      X_env_occ = abind(X_env_occ2, along=3)
      X_bias_occ = abind(X_bias_occ2, along=3)
      
      
      
      # Need to prepare the data again
      Oprep_TMB = TMB.optimprep(par=b.lassokeep1, yocc=y_occ, Wocc=X_bias_occ, Xocc=X_env_occ, ypo=y_po, 
                                extinct_comp = extinct_comp, coloni_comp = coloni_comp, wt_vec = wt.vec,
                                Wpo=X_bias_po, Xpo=X_env_po, ob_wt.po=ob_wt_po, lasso=lambda, n.time=n.time,
                                site_area = site.area, E=length(extinct_comp), C=length(coloni_comp), dat.type)
      
      BEstartk = c(0, dim(Oprep_TMB$Xpobmat)[2])
      BExstartk = c(0, length(extinct_comp))
      Bcostartk = c(0, length(coloni_comp))
      BbiPOstartk = c(0, Oprep_TMB$Wocc_NB[1]) # doesn't work if we have different bias covar throught the years
      BbiOccstartk = c(0, dim(Oprep_TMB$Wpomat)[2])
      
      # here we use lambda directly and not the prepared lassomat
      f.est <- MakeADFun(data=list(K=Oprep_TMB$K, Xocc=Oprep_TMB$Xocc, Wocc=Oprep_TMB$Wocc, site_area=site.area, 
                                    yocc=Oprep_TMB$yocc, R=Oprep_TMB$R, Jm=Oprep_TMB$Jm, S_NB=Oprep_TMB$S_NB, 
                                    Vocc_NB=Oprep_TMB$Vocc_NB, lasso_occ=Oprep_TMB$lasso_occ, extinctcomp=extinct_comp, 
                                    colonicomp=coloni_comp, Xpobmat=Oprep_TMB$Xpobmat, Xpoemat=Oprep_TMB$Xpoemat, 
                                    Xponmat=Oprep_TMB$Xponmat, Wpomat=Oprep_TMB$Wpomat, ypomat=Oprep_TMB$ypomat, 
                                    Wpos=Oprep_TMB$Wpos, Vpo_NB=Oprep_TMB$Vpo_NB, L_NB=Oprep_TMB$L_NB,
                                    row_startk=Oprep_TMB$row_startk, lassomat=Oprep_TMB$lassomat, lassostartk=Oprep_TMB$lassostartk,
                                    Lsize=Oprep_TMB$Lsize, ob_wtmat=Oprep_TMB$ob_wtmat, Wocc_NB=Oprep_TMB$Wocc_NB,
                                   BbiPOstartk=BbiPOstartk, BbiOccstartk=BbiOccstartk, 
                                   wtvecocc = Oprep_TMB$wtvecocc, wtvecpo = Oprep_TMB$wtvecpo, 
                                   BExstartk=BExstartk, Bcostartk=Bcostartk, BEstartk = BEstartk), 
                          parameters= list(b_init=Oprep_TMB$b_init), DLL = "TMB_comb")
      
      #Optim
      b.est <- optim(par=f.est$par, fn=f.est$fn, gr=f.est$gr, method='BFGS', 
                          control = list(ndeps = rep(b.min, num_p[[1]]*n.time))) # control to allow very small coef to be 0
      
      #b.lasso = b.est$par
     # like = b.est$value
    }

    
    sign.change                   = signs*sign(b.est$par)
    sign.change[sign.change == 0] = 1
    sign.change[abs(b.est$par) < b.min] = -1
    sign.change[is.in == FALSE]   = 1
    
    if (is.null(n.time)==TRUE){
      sign.change[noshrink] = 1
    }else{
      sign.change.L = list()
      for (i in 1:n.time) {
        sign.change.L[[i]] = sign.change[(num_startk[[i]]+1):num_endk[[i]]]
        sign.change.L[[i]][noshrink[[i]]] = 1
      }
      sign.change = do.call(c, sign.change.L)
    }
    
    
    score.beta  = b.est$par
    
    if (any(sign.change != 1) == TRUE)
    {
      delta                       = b.est$par - b.lasso
      tozero                      = abs(b.lasso[sign.change != 1]) - b.min
      prop                        = min(tozero/abs(delta[sign.change != 1]))
      
      b.lasso                     = b.lasso + prop*delta

      is.in[abs(b.lasso) <= b.min + b.min/(1e10)] = FALSE
      b.lasso[is.in == FALSE]     = 0
      signs                       = sign(b.lasso)
      setchange                   = 2*tol
      score.beta                  = b.lasso
      
    }
    
    if (is.null(n.time)==TRUE){
      score                   = scoreeq_new(Xenv, Xbias, y, b.lasso, ob.wt, link = link, site.area = site.area, dat.type = dat.type, lasso = 0, env_ix, bias_ix, coord, wt.vec, is.in = rep(TRUE, num_p))
    }else{
      NB_covar = dim(Oprep_TMB$b_init)[1]*dim(Oprep_TMB$b_init)[3]
      Grad1 = matrix(NA, nrow = NB_covar, ncol = b.est$counts[2])
      likcalc1=c()
      
      for (i in 1:b.est$counts[2]) {
        b.est <- optim(par=f.est$par, fn=f.est$fn, gr=f.est$gr, method='BFGS', control=list(maxit=i, ndeps = rep(b.min, num_p[[1]]*n.time)))
        likcalc1[i] = b.est$value
        ## for gradient value
        Grad1[,i] = f.est$gr(b.est$par)
      }
      
      score = apply(Grad1, 1, max)
      
      # to be able to work out everything using vectors ## lambda already a vector?
      if (class(lambda)=="list"){
        lambda = do.call(c, lambda)
      }else{
        lambda = lambda
      }
      
      # set up for infinite values
      lass.pen2 = unlist(pen.max)*1000000
      killed2 = unlist(killed)
      lambda[killed2 == TRUE]  = lass.pen2
      
      
      
    }

    
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
        
        if (is.null(n.time)==FALSE){
          for (i in 1:n.time) {
            lambda[which(lambda==pen.max[[i]])] = lasso.store[[i]]
          }
        }
      }
    }
    
    like    = b.est$value
    
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
  
  
  
  #,,, prepare the output results
  if (error.flag != TRUE)
  {
    
    #mu = exp(X %*% b.lasso)
    
    if (error.flag != TRUE)
    {
      if (is.null(n.time)==TRUE){
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
        }else{  ## Dynamic
          mu = rep(list(list()), n.time)
          psi = rep(list(list()), n.time)
          p_detect = rep(list(list()), n.time)
          
          
          # to pass the rbind vector for b.lasso into a list
          num_envk = num_pk= list()
          num_biask= b.lassok = rep(list(c()), n.time)
          for (i in 1:n.time) {
            num_envk[[i]]  = dim(X_env.use[[i]][[1]])[2]
            num_biask[[i]] = unlist(lapply(lapply(X_bias.use[[i]], as.matrix), ncol))
            num_pk[[i]]    = num_env[[i]] + sum(num_bias[[i]])
          }
          num_startk = c(0, cumsum(num_pk))
          num_endk   = cumsum(num_pk)
          
          for (i in 1:n.time) {
            b.lassok[[i]] = b.lasso[(num_startk[[i]]+1):num_endk[[i]]]
          }
          
          lin_det = list()
          for (k in 1:n.time) {
            for (i in 1:n.components[[k]])
            {
              mu[[k]][[i]] = exp(X_env.use[[k]][[i]] %*% b.lassok[[k]][env_ix[[1]][[1]]])
              psi[[k]][[i]] = 1 - exp(-mu[[k]][[i]]*site.area)
              if (dat.type[i] == "PO")
              {
                p_detect[[k]][[i]] = NULL
              }
              else
              {
                lin_det[[k]]  = X_bias.use[[k]][[i]] %*% b.lassok[[k]][bias_ix[[k]][[i]]]
                if (link == "logit")
                {
                  p_detect[[k]][[i]] = expit(lin_det[[k]])
                }
                if (link == "cloglog")
                {
                  p_detect[[k]][[i]] = clogloginv(lin_det[[k]])
                  #p_detect[[i]] = 1 - exp(-exp(lin_det))
                }
              }
            }
          }
        }
    }
    
    
    return(list(b = b.lasso, mu = mu, psi = psi, p_detect = p_detect, ls = likes, pen.likelihood = likes[length(likes)], bs = betas, s = scores, v = viols, actions = actions, n.time=n.time))
    #return(list(b = b.out, mu = mu, e.df = eff.df, deviance = dev, likelihood = unp.likelihood, GCV = gcv, AIC = aic, BIC = bic, HQC = hqc, AICc = aicc, ls = likes, pen.likelihood = likes[length(likes)], bs = betas, s = scores, v = viols, actions = actions))
    
  }
}


score.int_new = function(Xenv, Xbias = NULL, y, ob.wt = NULL, link = "logit", sp_res = 1, site.area = 1, dat.type = "PO", lasso = 0, env_ix, bias_ix, coord = c("x", "y"), wt.vec, noshrink, method = "BFGS", b.min = 1.e-5, extinct_comp, coloni_comp, n_components=NULL, n.time=NULL)
{
  
  if (class(y) == "list"){
    if (is.null(n.time)==FALSE){
      # Prepare data
      
      y.po= X_env.po= X_bias.po= ob_wt.po=y.occ2=X_env.occ2=X_bias.occ2=rep(list(list()), n.time)
      y.occ=X_env.occ=X_bias.occ=list()
      
      
        for (t in 1:n.time) {
          for (i in 1:n_components[[t]]) {
            if (dat.type[i] == "PO") {
              y.po[[t]] = y[[t]][[i]]
              X_env.po[[t]] = Xenv[[t]][[i]]
              X_bias.po[[t]] = Xbias[[t]][[i]]
              ob_wt.po[[t]] = ob.wt[[t]][[i]]
            }else{
              y.occ2[[t]] = y[[t]][[i]] 
              X_env.occ2[[t]] = Xenv[[t]][[i]]
              X_bias.occ2[[t]] = Xbias[[t]][[i]]  # numeric for each occ year so need as.matrix
              #ob_wt.occ is NULL for occ
            }
            
          }
        }
     
        y.occ = abind(y.occ2, along=3)
        X_env.occ = abind(X_env.occ2, along=3)
        X_bias.occ = abind(X_bias.occ2, along=3)
     
        
      
      num_env = num_bias=num_p=list()
      for (i in 1:n.time) {
        num_env[[i]]  = dim(Xenv[[i]][[1]])[2]
        num_bias[[i]] = unlist(lapply(lapply(Xbias[[i]], as.matrix), ncol))
        num_p[[i]]    = num_env[[i]] + sum(num_bias[[i]])
      }
    }else{
      num_env  = dim(Xenv[[1]])[2]
      num_bias = unlist(lapply(lapply(Xbias, as.matrix), ncol))
      num_p    = num_env + sum(num_bias)
    }
    
  }else{ # only for non dynamic so far when it is not a list
    num_p = dim(Xenv)[2] + dim(Xbias)[2]
  }
  
  if (is.null(n.time)==FALSE){
    par=beta0=list()
    for (t in 1:n.time) {
      par[[t]]=rep(0, num_p[[t]])
      beta0[[t]] = rep(0, length(y[[t]]))
      
      for (i in 1:length(dat.type))
      {
        if (dat.type[i] == "PO")
        {
          beta0[[t]][i] = sum(y[[t]][[i]] > 0)/(sum(y[[t]][[i]] == 0)*sp_res^2)
        }
        if (dat.type[i] == "Occ")
        {
          beta0[[t]][i] = sum(apply(y[[t]][[i]], 1, sum) > 0)/(dim(y[[t]][[i]])[1]*site.area)
        }
      }
      
      par[[t]][noshrink[[t]][which(noshrink[[t]] %in% env_ix[[t]][[1]])]] = log(mean(beta0[[t]]))
    }
    
    ## TMB part
    # Prepare data to use in the TMB script with this function
    Oprep_TMB = TMB.optimprep(par=par, yocc=y.occ, Wocc=X_bias.occ, Xocc=X_env.occ, ypo=y.po, 
                              extinct_comp = extinct_comp, coloni_comp = coloni_comp, wt_vec = wt.vec,
                              Wpo=X_bias.po, Xpo=X_env.po, ob_wt.po=ob_wt.po, lasso=0, n.time=n.time,
                              site_area = site.area, E=length(extinct_comp), C=length(coloni_comp), dat.type)

    BEstartk = c(0, dim(Oprep_TMB$Xpobmat)[2])
    BExstartk = c(0, length(extinct_comp))
    Bcostartk = c(0, length(coloni_comp))
    BbiPOstartk = c(0, Oprep_TMB$Wocc_NB[1]) # doesn't work if we have different bias covar throught the years
    BbiOccstartk = c(0, dim(Oprep_TMB$Wpomat)[2])
    
        ## Dynamic 
    model <- paste("#include <TMB.hpp>
                   #include <cmath>
				           #include <cstdlib>

                   template<class Type>
                   Type objective_function<Type>::operator() ()
                   {
                   // Define Data OCC
                   DATA_ARRAY(Xocc); 	// environment variable values
                   DATA_ARRAY(Wocc); 	// bias covariates
                   DATA_ARRAY(yocc); 	// survey occupancy state - sum per visits
                   DATA_SCALAR(site_area);
                   DATA_VECTOR(lasso_occ);
                   DATA_INTEGER(K);  // nb of time period
                   DATA_INTEGER(S_NB); // NB of sites
                   DATA_INTEGER(Vocc_NB); // NB of visits
                   DATA_VECTOR(extinctcomp); // extinction covar NB
                   DATA_VECTOR(colonicomp); // colonization covar NB
                   DATA_MATRIX(wtvecocc);  // Vector of source weights for OCC
                   DATA_MATRIX(Jm);
                   DATA_MATRIX(R);
                   DATA_IVECTOR(BEstartk); //vector of start NB covar for occ
                   DATA_IVECTOR(BExstartk); //vector of start NB covar for occ
                   DATA_IVECTOR(Bcostartk); //vector of start NB covar for occ
                   DATA_IVECTOR(BbiOccstartk); //vector of start NB covar for occ
                   DATA_INTEGER(Wocc_NB); //NB of bias covar for occ
                   
                   // Define Data PO
                   DATA_MATRIX(Xpobmat); 	// environment variable values
                   DATA_MATRIX(Xpoemat); 	// environment variable values - colonization
                   DATA_MATRIX(Xponmat); 	// environment variable values - extinction
                   DATA_VECTOR(ob_wtmat); //observation weight rbind vector per year
                   DATA_MATRIX(Wpomat); 	// bias covariates
                   DATA_MATRIX(ypomat); 	// survey po data
                   DATA_VECTOR(lassomat); // vector of rbind lasso coef
                   DATA_MATRIX(wtvecpo);  // Vector of source weights for PO
                   
                   DATA_INTEGER(Vpo_NB); // covar NB # of elements is the number of time period
                   DATA_IVECTOR(L_NB); // length enviro covar # of elements is the number of time period
                   DATA_INTEGER(Wpos); // NB of bias covar coef
                   DATA_IVECTOR(row_startk); // vector of X starting number for each year
                   DATA_IVECTOR(lassostartk); // vector of b coef starting number for each year
                   DATA_IVECTOR(Lsize); // vector of bcoef size for each year
                   DATA_IVECTOR(BbiPOstartk); //vector of start NB covar for po
                   
                   
                   // Define parameters
                   PARAMETER_ARRAY(b_init);
                   
                   // Some usefull info
                   vector<int> Xoccsize = Xocc.dim;	
                   int Xoccs = Xoccsize(1); // nb of var
                   vector<int> Woccsize = Wocc.dim; // nb of bias covar for OCC data
                   int Woccs = Woccsize(1);
                   int Extsize = extinctcomp.size();
                   int colsize = colonicomp.size();
                   int Xpos = Vpo_NB; // nb of var
                   
                   vector<int> dimbe(K);
                   for (int i=0; i<K; i++){
                   dimbe(i) = Xoccs-Extsize-colsize;
                   }
                   
                   vector<int> dimet(K);
                   for (int i=0; i<K; i++){
                   dimet(i) = Extsize;
                   }
                   
                   vector<int> dimnu(K);
                   for (int i=0; i<K; i++){
                   dimnu(i) = colsize;
                   }
                   
                   vector<int> dimal(K);
                   for (int i=0; i<K; i++){
                   dimal(i) = Woccs;
                   }
                   //REPORT(dimal);
                   
                   // Redefine b_init for occ only
                   vector<int> b_dim = b_init.dim ; // dimension of a1
                   array<Type> B_occ(Vpo_NB+Wocc_NB, 1, K) ;
                   for (int t=0; t<K; t++){
                   for (int j=0; j<Vpo_NB; j++){
                   B_occ(j,0,t)=b_init(j,0,t);
                   }
                   for (int j=0; j<Wocc_NB; j++){
                   B_occ(Vpo_NB+j,0,t)=b_init(Vpo_NB+Wpos+j,0,t);
                   }
                   }
                   REPORT(B_occ);
                   //REPORT(b_dim);
                   
                   
                   // Redefine b_init for PO only
                   array<Type> B_po(Vpo_NB+Wpos, 1, K);
                   for (int t=0; t<K; t++){
                   for (int j=0; j<Vpo_NB+Wpos; j++){
                   B_po(j,0,t)=b_init(j,0,t);
                   }
                   }
                   REPORT(B_po);
                   
                   // Separate Coefficients for OCC
                   matrix<Type> BE(Xoccs-Extsize-colsize, K);
                   matrix<Type> AL(Woccs, K);
                   matrix<Type> ETA(Extsize, K);
                   matrix<Type> NU(colsize, K);
                   
                   BE.setZero();
                   AL.setZero();
                   ETA.setZero();
                   NU.setZero();
                   
                   for (int i = 0; i < K; i++){
                   for (int j = 0; j < Xoccs-Extsize-colsize; j++){
                   BE(j,i) = B_occ(j,0,i);
                   }
                   for (int j = 0; j < Extsize; j++){
                   ETA(j,i) = B_occ(j+Xoccs-Extsize-colsize,0,i);
                   }
                   for (int j = 0; j < colsize; j++){
                   NU(j,i) = B_occ(j+Xoccs-Extsize,0,i);
                   }
                   for (int j = 0; j < Woccs; j++){
                   AL(j,i) = B_occ(j+Xoccs,0,i); // Need to start at 0
                   }
                   
                   }
                   
                   REPORT(BE);
                   REPORT(AL);
                   REPORT(ETA);
                   REPORT(NU);
                   
                   
                   // ***Loglikelihood function calculation for OCC DATA
                   // Environmental covariates
                   
                   matrix<Type> Xocc_time_beta(S_NB, K);
                   Xocc_time_beta.setZero();
                   
                   matrix<Type> Xocctime(S_NB, Xoccs-Extsize-colsize);
                   vector<Type> BEtime(Xoccs-Extsize-colsize);
                   
                   for (int t=0; t<K; t++){                         // we loop around the years to only consider Matrix vector product
                   for (int i=0; i<S_NB; i++){
                   for (int j=0; j<Xoccs-Extsize-colsize; j++){
                   Xocctime(i,j) = Xocc(i,j,t);
                   BEtime(j) = BE(j,t);
                   }
                   }
                   for (int i = 0; i < Xoccsize(0); i++){
                   for (int k = 0; k < Xoccs-Extsize-colsize; k++){
                   Xocc_time_beta(i,t) += BEtime(k)*Xocctime(i,k);  // We obtain a vector for each year
                   }
                   } 
                   }
                   
                   //REPORT(Xocc_time_beta)
                   
                   matrix<Type> LAMBDA=exp(Xocc_time_beta.array());
                   matrix<Type> L_time_siteA = LAMBDA*site_area;
                   matrix<Type> PSI = Type(1) - exp(Type(-1)*L_time_siteA.array());
                   
                   // BIAS covariates
                   
                   vector<int> Woccndim = Wocc.dim;
                   matrix<Type> Wocc_time_alpha(Woccndim(0), K);
                   Wocc_time_alpha.setZero();
                   
                   matrix<Type> Wocctime(S_NB, Woccs);
                   vector<Type> ALtime(Woccs);
                   
                   for (int t=0; t<K; t++){                         // we loop around the years to only consider Matrix vector product
                   for (int i=0; i<S_NB; i++){
                   for (int j=0; j<Woccs; j++){
                   Wocctime(i,j) = Wocc(i,j,t);
                   ALtime(j) = AL(j,t);
                   }
                   }
                   for (int i = 0; i < Woccndim(0); i++){
                   for (int k = 0; k < Woccndim(1); k++){
                   Wocc_time_alpha(i,t) += ALtime(k)*Wocctime(i,k);  // We obtain a vector for each year
                   }
                   } 
                   }
                   
                   //matrix<Type> P = Type(1)/(Type(1) + exp(Type(-1)*Wocc_time_alpha.array()));
                   matrix<Type> P = Type(1)-exp(Type(-1)*exp(Wocc_time_alpha.array()));

                   //REPORT(P);
                   //REPORT(Woccs);
                   
                   // Extinction process
                   matrix<Type> Xocce_time_eta(S_NB, K);
                   Xocce_time_eta.setZero();
                   
                   matrix<Type> Xoccetime(S_NB, Extsize);
                   vector<Type> ETAtime(Extsize);
                   
                   for (int t=0; t<K; t++){                         // we loop around the years to only consider Matrix vector product
                   for (int i=0; i<S_NB; i++){
                   for (int j=0; j<Extsize; j++){
                   Xoccetime(i,j) = Xocc(i,j+Xoccs-Extsize-colsize,t);
                   ETAtime(j) = ETA(j,t);
                   }
                   }
                   for (int i = 0; i < S_NB; i++){
                   for (int k = 0; k < Extsize; k++){
                   Xocce_time_eta(i,t) += ETAtime(k)*Xoccetime(i,k);  // We obtain a vector for each year
                   }
                   } 
                   }
                   
                   matrix<Type> EPSILON = Type(1)/(Type(1) + exp(Type(-1)*Xocce_time_eta.array()));
                   
                   // Colonization process
                   matrix<Type> Xoccc_time_nu(S_NB, K);
                   Xoccc_time_nu.setZero();
                   
                   matrix<Type> Xoccctime(S_NB, colsize);
                   vector<Type> NUtime(colsize);
                   
                   for (int t=0; t<K; t++){                         // we loop around the years to only consider Matrix vector product
                   for (int i=0; i<S_NB; i++){
                   for (int j=0; j<colsize; j++){
                   Xoccctime(i,j) = Xocc(i,j+Xoccs-Extsize,t);
                   NUtime(j) = NU(j,t);
                   }
                   }
                   for (int i = 0; i < S_NB; i++){
                   for (int k = 0; k < colsize; k++){
                   Xoccc_time_nu(i,t) += NUtime(k)*Xoccctime(i,k);  // We obtain a vector for each year
                   }
                   } 
                   }
                   
                   matrix<Type> GAMMA = Type(1)/(Type(1) + exp(Type(-1)*Xoccc_time_nu.array()));
                   
                   vector<Type> Jtmp(S_NB);
                   vector<Type> Rtmp(S_NB);
                   vector<Type> PSItmp(S_NB);
                   Jtmp.setZero();
                   Rtmp.setZero();
                   PSItmp.setZero();
                   
                   // State process
                   matrix<Type> ZBvec(S_NB,K);
                   vector<Type> ZBtmp(S_NB);
                   ZBvec.setZero();
                   
                   for (int i=0; i<K; i++){
                   for (int j=0; j<S_NB; j++){
                   Jtmp(j) = Jm(j,i);
                   Rtmp(j) = R(j,i);
                   PSItmp(j) = PSI(j,i);
                   
                   ZBtmp(j) = dbinom(Rtmp(j), Jtmp(j), PSItmp(j), false);
                   }
                   ZBvec.col(i) = ZBtmp;
                   }
                   
                   matrix<Type> muZ(S_NB, K);
                   muZ.setZero();
                   
                   muZ.col(0) = ZBvec.col(0);    // First year

                   vector<Type> ytmp(S_NB);
                   vector<Type> J2tmp(S_NB);
                   vector<Type> R2tmp(S_NB);
                   vector<Type> muZtmp(S_NB);
                   ytmp.setZero();
                   J2tmp.setZero();
                   R2tmp.setZero();
                   muZtmp.setZero();
                   
                   matrix<Type> Zcond(S_NB, K);
                   vector<Type> Ztmp(S_NB);
                   Zcond.setZero();
                   
                   for (int i=0; i<K-1; i++){
                   for (int j=0; j<S_NB; j++){
                   muZ(j,i+1) = ZBvec(j,i)*(Type(1)-EPSILON(j,i+1)) + (Type(1)-ZBvec(j,i))*GAMMA(j,i+1);
                   
                   J2tmp(j) = Jm(j,i);
                   R2tmp(j) = R(j,i);
                   muZtmp = muZ(j,i+1);
                   
                   Ztmp(j) = dbinom(R2tmp(j), J2tmp(j), muZtmp(j), false);
                   }
                   ZBvec.col(i+1) = Ztmp;
                   }
                   //REPORT(muZ);
                   //REPORT(ZBvec);
                  
          
                   // Observation process
                   vector<Type> y2tmp(Vocc_NB);
                   vector<Type> J3tmp(S_NB);
                   vector<Type> PtimeZtmp(S_NB);
                   y2tmp.setZero();
                   J3tmp.setZero();
                   PtimeZtmp.setZero();
                   
                   matrix<Type> Ycond(S_NB, K);
                   vector<Type> Ytmp(Vocc_NB);
                   Ycond.setZero();
                   //vector<Type> YsumP(S_NB);
                   matrix<Type> Yvisit(S_NB, Vocc_NB);
                   
                   for (int i=0; i<K; i++){
                   for (int v=0; v<Vocc_NB; v++){
                   for (int j=0; j<S_NB; j++){
                   
                   J3tmp(j) = Jm(j,i);
                   PtimeZtmp(j) = ZBvec(j,i)*P(j,i);
                   
                   y2tmp(v) = yocc(j,v,i);
                   Ytmp(v) = dbinom(y2tmp(v), J3tmp(j), PtimeZtmp(j), false);
                   
                   Yvisit(j,v) = Ytmp(v);
                   
                   }
                   }
                   Ycond.col(i) = Yvisit.rowwise().sum();   // sum accross visits
                   }
                   //REPORT(Ycond);
                   //REPORT(Zcond);
                   
                   matrix<Type> Ylog = log(Ycond.array());  
                   
                   vector<Type> liksumOCC(K);
                   vector<Type> liksumOCCw(K);
                   Type nLLocc(0);
                   
                  // Type lassoCoefocc(0);
                   //Type lassCoefSumOcc(0);
                   
                   // Group lasso
                   vector<Type> lassoenv = lasso_occ.segment(BEstartk(0), dimbe(0));
                   vector<Type> lassoext = lasso_occ.segment(BEstartk(0)+BExstartk(0), dimet(0));
                   vector<Type> lassocol = lasso_occ.segment(BEstartk(0)+BExstartk(0)+Bcostartk(0), dimnu(0));
                   vector<Type> lassobiasO = lasso_occ.segment(BEstartk(0)+BExstartk(0)+Bcostartk(0)+BbiOccstartk(0), dimal(0));
                   
                   vector<vector<Type> > BEv(K);
                   vector<vector<Type> > ETAv(K);
                   vector<vector<Type> > NUv(K);
                   vector<vector<Type> > ALv(K);
                   
                   for (int i=0; i<K; i++){
                   vector<Type> betmp = BE.col(i);
                   BEv(i)= abs(betmp.array());
                   
                   vector<Type> etatmp = ETA.col(i);
                   ETAv(i)= abs(etatmp.array());
                   
                   vector<Type> nutmp = NU.col(i);
                   NUv(i)= abs(nutmp.array());
                   
                   vector<Type> altmp = AL.col(i);
                   ALv(i)= abs(altmp.array());
                   }
                   
                   //REPORT(BEv);
                   

                   liksumOCC = Ylog.colwise().sum(); // sum accross sites 

                   for (int i = 0; i < K; i++){
                      liksumOCCw(i) = wtvecocc(i)*liksumOCC(i);
                   }
                   nLLocc = - liksumOCCw.sum(); //sum accross time
                   
                   //REPORT(lassoCoefocc);
                   //REPORT(liksumOCC);
                   //REPORT(liksumOCCw);                   
                   //REPORT(lassCoefSumOcc);
                   //REPORT(nLLocc);
                   
                   
                   // ***Likelihood calculation for PO DATA 

                   // Separate Coefficients for PO
                   matrix<Type> BE2(Xpos-Extsize-colsize, K);
                   matrix<Type> AL2(Wpos, K);
                   matrix<Type> ETA2(Extsize, K);
                   matrix<Type> NU2(colsize, K);
                   
                   BE.setZero();
                   AL.setZero();
                   ETA.setZero();
                   NU.setZero();
                   
                   for (int i = 0; i < K; i++){
                   for (int j = 0; j < Xpos-Extsize-colsize; j++){
                   BE2(j,i) = B_po(j,0,i);
                   }
                   for (int j = 0; j < Extsize; j++){
                   ETA2(j,i) = B_po(j+Xpos-Extsize-colsize,0,i);
                   }
                   for (int j = 0; j < colsize; j++){
                   NU2(j,i) = B_po(j+Xpos-Extsize,0,i);
                   }
                   for (int j = 0; j < Wpos; j++){
                   AL2(j,i) = B_po(j+Xpos,0,i); // Need to start at 0
                   }
                   
                   }

                   REPORT(BE2);
                   REPORT(ETA2);
                   REPORT(NU2);
                   REPORT(AL2);

                   // Prepare envrionmental matrix
                   // Separate the sub matrix for extinction, colonization, bias and other covariates
                   
                   // Basic covariates and MU contribution
                   
                   vector<matrix<Type> > Xpobt(K);
                   vector<Type> BEpotime(Xpos-Extsize-colsize);
                   vector<vector<Type> > Xpotbe(K);
                   
                   for (int t=0; t<K; t++){
                   matrix<Type> Xpob_block = Xpobmat.block(row_startk(t),0,L_NB(t),Xpos-Extsize-colsize);
                   Xpobt(t) = Xpob_block;
                   
                   for (int j=0; j<Xpos-Extsize-colsize; j++){
                   BEpotime(j) = BE2(j,t);
                   }
                   
                   Xpotbe(t) = exp(Xpobt(t)*BEpotime); // MU contribution
                   }
                   
                   //REPORT(BEtime);
                   //REPORT(Xpotbe);
                   
                   // Extinction covariates and MU contribution
                   
                   vector<matrix<Type> > Xpoet(K);
                   vector<Type> ETApotime(Extsize);
                   vector<vector<Type> > Xpoteta(K);
                   
                   for (int t=0; t<K; t++){
                   matrix<Type> Xpoe_block = Xpoemat.block(row_startk(t),0,L_NB(t),Extsize);
                   Xpoet(t) = Xpoe_block;
                   
                   for (int j=0; j<Extsize; j++){
                   ETApotime(j) = ETA2(j,t);
                   }
                   
                   Xpoteta(t) = Type(1)/(Type(1)+exp(Type(-1)*(Xpoe_block*ETApotime).array())); // MU contribution
                   
                   }
                   
                   //REPORT(Xpoteta);
                   
                   // Colonization covariates and MU contribution
                   
                   vector<matrix<Type> > Xpont(K);
                   vector<Type> NUpotime(colsize);
                   vector<vector<Type> > Xpotnu(K);
                   
                   for (int t=0; t<K; t++){
                   matrix<Type> Xpon_block = Xponmat.block(row_startk(t),0,L_NB(t),colsize);
                   Xpont(t) = Xpon_block;
                   
                   for (int j=0; j<colsize; j++){
                   NUpotime(j) = NU2(j,t);
                   }
                   
                   Xpotnu(t) =   Type(1)/(Type(1)+exp(Type(-1)*(Xpont(t)*NUpotime).array()));  // MU contribution
                   }
                   
                   //REPORT(Xpont);
                   
                   // Bias covariates and MU contribution
                   
                   vector<matrix<Type> > Wpot(K);
                   vector<Type> ALpotime(Wpos);
                   vector<vector<Type> > Wpotal(K);
                   
                   for (int t=0; t<K; t++){
                   matrix<Type> Wpob_block = Wpomat.block(row_startk(t),0,L_NB(t),Wpos);
                   Wpot(t) = Wpob_block;
                   
                   for (int j=0; j<Wpos; j++){
                   ALpotime(j) = AL2(j,t);;
                   }
                   
                   Wpotal(t) = exp(Wpot(t)*ALpotime);  // MU contribution
                   }
                   
                   //REPORT(Wpotal);
                   
                   // MU calculation
                   
                   vector<vector<Type> > MU(K);
                   for (int t=0; t<K; t++){
                   MU(t)= Xpotbe(t)*Xpoteta(t)*Xpotnu(t)*Wpotal(t);
                   }
                   
                   //REPORT(MU);
                   
                   vector<vector<Type> > ypot(K);
                   for (int t=0; t<K; t++){
                   vector<Type> ypo_block = ypomat.block(row_startk(t),0,L_NB(t),1);
                   ypot(t) = ypo_block;
                   }
                   
                   
                   // Observation weights
                   vector<vector<Type> > OBW(K); // Prepare observation weights
                   
                   for (int t=0; t<K; t++){
                   vector<Type> OBW_block = ob_wtmat.segment(row_startk(t), L_NB(t));
                   OBW(t) = OBW_block;
                   }
                   //REPORT(OBW);
                   
                   vector<vector<Type> >predypo(K);
                   vector<Type> predypo_S(K);
                   
                   for (int t=0; t<K; t++){
                   predypo(t) = OBW(t)*(ypot(t)*log(MU(t)) - MU(t));
                   predypo_S(t) = predypo(t).sum();
                   }
                   
                   Type predypo_sum(0);
                   predypo_sum = predypo_S.sum();
                   //REPORT(predypo);
                   //REPORT(predypo_sum);
                   
                   
                   // Partial likelihood
                   
                   vector<vector<Type> > partL(K);
                   vector<vector<Type> > partL_S(K);
                   vector<Type> partL_s2(K);
                   
                   for (int t=0; t<K; t++){
                    partL(t) = predypo(t)/predypo_sum;

                    partL_S(t) = log(abs(partL(t).array()));  
                    partL_s2(t) = wtvecpo(t)*sum(partL_S(t));  //sum accross sites
                   }
                   
                   Type partL_sum(0);
                   partL_sum = partL_s2.sum(); // sum accross years
                   //REPORT(partL_sum);
                   //REPORT(partL_s2);
                   //REPORT(partL);
                   
                   // Lasso coefficient penalties
                   // Group lasso
                    
                   vector<vector<Type> > BE2v(K);
                   vector<vector<Type> > ETA2v(K);
                   vector<vector<Type> > NU2v(K);
                   vector<vector<Type> > AL2v(K);
                   
                   for (int i=0; i<K; i++){
                   vector<Type> be2tmp = BE2.col(i);
                   BE2v(i)= abs(be2tmp.array());
                   
                   vector<Type> eta2tmp = ETA2.col(i);
                   ETA2v(i)= abs(eta2tmp.array());
                   
                   vector<Type> nu2tmp = NU2.col(i);
                   NU2v(i)= abs(nu2tmp.array());
                   
                   vector<Type> al2tmp = AL2.col(i);
                   AL2v(i)= abs(al2tmp.array());
                   }
                   
                   //REPORT(BEv);
                   
                   //Type lassoCoefPO(0);
                   
                   vector<Type> lassoenvpo = lassomat.segment(BEstartk(0), dimbe(0));
                   vector<Type> lassoextpo = lassomat.segment(BEstartk(0)+BExstartk(0), dimet(0));
                   vector<Type> lassocolpo = lassomat.segment(BEstartk(0)+BExstartk(0)+Bcostartk(0), dimnu(0));
                   vector<Type> lassobiasP = lassomat.segment(BEstartk(0)+BExstartk(0)+Bcostartk(0)+BbiPOstartk(0), dimal(0));
                   
                   
                   //lassCoefSumPO = sum(lassoCoefPO);
                   
                   // Likelihood calculation
                   
                   Type liksumPO(0);
                   Type nLLpo(0);
                   
                   liksumPO = partL_sum; 
                   nLLpo = - liksumPO;
                   
                   
                   //REPORT(lassCoefPO);
                   //REPORT(nLLpo);
                   
                   
                   Type AllLasso(0);
                   AllLasso = sum((lassoenvpo*BE2v).sum())+sum((lassoextpo*ETA2v).sum())+sum((lassocolpo*NU2v).sum())+sum((lassobiasO*ALv).sum())+sum((lassobiasP*AL2v).sum());
                   
                   Type comblik(0);
                   
                   comblik = nLLocc + nLLpo + AllLasso; 
                   REPORT(comblik);
                   
                   return comblik;
                   }
                   ")
    writeLines(model,"TMB_comb.cpp")
    compile("TMB_comb.cpp") 

    dyn.load(dynlib("TMB_comb"))
    
    f0comb <- MakeADFun(data=list(K=Oprep_TMB$K, Xocc=Oprep_TMB$Xocc, Wocc=Oprep_TMB$Wocc, site_area=site.area, 
                                  yocc=Oprep_TMB$yocc, R=Oprep_TMB$R, Jm=Oprep_TMB$Jm, S_NB=Oprep_TMB$S_NB, 
                                  Vocc_NB=Oprep_TMB$Vocc_NB, lasso_occ=Oprep_TMB$lasso_occ, extinctcomp=extinct_comp, 
                                  colonicomp=coloni_comp, Xpobmat=Oprep_TMB$Xpobmat, Xpoemat=Oprep_TMB$Xpoemat, 
                                  Xponmat=Oprep_TMB$Xponmat, Wpomat=Oprep_TMB$Wpomat, ypomat=Oprep_TMB$ypomat, 
                                  Wpos=Oprep_TMB$Wpos, Vpo_NB=Oprep_TMB$Vpo_NB, L_NB=Oprep_TMB$L_NB,
                                  row_startk=Oprep_TMB$row_startk, lassomat=Oprep_TMB$lassomat, lassostartk=Oprep_TMB$lassostartk,
                                  Lsize=Oprep_TMB$Lsize, ob_wtmat=Oprep_TMB$ob_wtmat, Wocc_NB=Oprep_TMB$Wocc_NB,
                                  BEstartk = BEstartk, BbiPOstartk=BbiPOstartk, BbiOccstartk=BbiOccstartk,
                                  BExstartk=BExstartk, Bcostartk=Bcostartk, 
                                  wtvecocc = Oprep_TMB$wtvecocc, wtvecpo = Oprep_TMB$wtvecpo), 
                        parameters= list(b_init=Oprep_TMB$b_init), DLL = "TMB_comb")
    
    
    ## combine at each iteration
    f0comb.res <- optim(par=f0comb$par, fn=f0comb$fn, gr=f0comb$gr, method='BFGS', 
                        control = list(ndeps = rep(b.min, num_p[[1]]*n.time))) # control to allow very small coef to be 0
    
    NB_covar = dim(Oprep_TMB$b_init)[1]*dim(Oprep_TMB$b_init)[3]
    bparam0 = Grad0 = matrix(NA, nrow = NB_covar, ncol = f0comb.res$counts[2])
    likcalc0=c()
    
    for (i in 1:f0comb.res$counts[2]) {
      f0comb.res <- optim(par=f0comb$par, fn=f0comb$fn, gr=f0comb$gr, method='BFGS', control=list(maxit=i, ndeps = rep(b.min, num_p[[1]]*n.time)))
      bparam0[,i] = f0comb.res$par
      likcalc0[i] = f0comb.res$value
      ## for gradient value
      Grad0[,i] = f0comb$gr(f0comb.res$par)
    }
    
    score_int = apply(Grad0, 1, max)
    newpar = f0comb.res$par
    
  }else{
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
    
    newpar = intercept_model$par
    
    # July 17 2018: behaves strangely with 0 as par for occ lasso with large site area. 
    # Better starting point than 0 (such as log(31/(73*site.area))) seems to work, but how to generalize?
    # Maybe an average intensity from PO plus occ?
    
    score_int = scoreeq_new(Xenv, Xbias = Xbias, y = y, par = intercept_model$par, ob.wt = ob.wt, 
                            link = link, site.area = site.area, dat.type = dat.type, 
                            lasso = 0, env_ix = env_ix, bias_ix = bias_ix, coord = coord, 
                            wt.vec = wt.vec, is.in = rep(TRUE, num_p))
  }
  
  
  
  return(list(score_int = score_int, newpar = newpar))
}


TMB.optimprep=function(par, yocc, Wocc, Xocc, ypo, Wpo, Xpo, ob_wt.po, lasso=0, site_area = site.area, extinct_comp = extinct_comp, coloni_comp = coloni_comp, E = length(extinct_comp), C = length(coloni_comp), n.time=n.time, wt_vec = wt.vec, dat.type=dat.type){
  # Prepare data and element to load into TMB function OCC
  n.times=dim(yocc)[3]
  
  if(is.null(n.times)==FALSE){
    K = dim(yocc)[3] # NB of years/time periods
  }else{
    K = NULL
  }
  
  Wocc_NB = dim(Wocc)[2]
  S_NB = dim(yocc)[1] # number of sites
  Vocc_NB = dim(yocc)[2] # number of visits
  
  R=matrix(NA, S_NB, K)
  Jm=matrix(NA, S_NB, K)
  
  wtvecpo = matrix(NA, ncol=length(which(dat.type=="PO")), nrow=n.times)
  wtvecocc = matrix(NA, ncol=length(which(dat.type=="Occ")), nrow=n.times)
  
  for (i in 1:n.times) {
    R[,i] = rep(1, S_NB) # NB of sites could be make general for different years. NB sites per year
    Jm[,i] = ncol(yocc[,,i]) - apply(is.na(yocc[,,i]), 1, sum) # NB of sites with different visits per year
    
    if(is.null(wt_vec)==FALSE){
      for (j in 1:length(which(dat.type=="PO"))) {
        wtvecpo[i,j] = wt_vec[[i]][j]
      }
      for (l in 1:length(which(dat.type=="Occ"))) {
        wtvecocc[i,l] = wt_vec[[i]][length(which(dat.type=="PO"))+l]
      }
    }else{
      for (j in 1:length(which(dat.type=="PO"))) {
        wtvecpo[i,j] = 1
      }
      for (l in 1:length(which(dat.type=="Occ"))) {
        wtvecocc[i,l] = 1
      }
    }
    
  }
  
  
  
  # set up coefficient for variables not to shrink and lasso
  noshrink=NULL
  
  noshrink   = c("Intercept", noshrink)
  colname.bias=c()
  for (k in 1:n.times) {
    colname.bias[k] = colnames(Wocc)
  }
  all_names  = c(colnames(Xocc[,,1]), colname.bias)
  
  noshrink   = which(all_names %in% noshrink)
  num_env = num_bias = num_p = lasso_occ.start = rep(list(list()), n.times)
  
  if (dim(Wpo[[1]])[2]==1){
    Wpo.pos = dim(Xpo[[1]])[2]+dim(Wpo[[1]])[2]
  }else{
    for (w in 1:dim(Wpo[[1]])[2]) {
      Wpo.pos = c(seq(from=dim(Xpo[[1]])[2]+1, to=dim(Xpo[[1]])[2]+Wpo[[1]][2], by=1))
    }
  }
  
  if (dim(Wocc)[2]==1){
    Wocc.pos = dim(Xpo[[1]])[2]+dim(Wpo[[1]])[2]+dim(Wocc)[2]
  }else{
    for (w in 1:dim(Wocc)[2]) {
      Wocc.pos = c(seq(from=dim(Xpo[[1]])[2]+Wpo[[1]][2]+1, to=dim(Xpo[[1]])[2]+Wpo[[1]][2]+dim(Wocc)[2], by=1))
    }
  }

  
  lasso_occ.start =  rep(list(list()), n.time)
  lasso_occ.list = rep(list(c()), n.time)
  
  lamb=lasso
  if (length(lamb)>1){
    if (class(lamb)=="list"){
      for (i in 1:K) {
        lasso_occ.start[[i]] = lamb[[i]][-Wpo.pos]
      }
      lasso_occ = unlist(lasso_occ.start)
    }else{
      lasso_occ = lamb[-Wpo.pos]
    }
  }else{
    for (i in 1:n.time) {
      num_env[[i]]   = dim(Xocc)[2]
      num_bias[[i]]  = dim(Wocc)[2]
      num_p[[i]]    = num_env[[i]] + sum(num_bias[[i]])
    }
    
    lasso_occ = (if(length(lamb) == 1) rep(lamb, num_p[[1]])
                 else rep(0, num_p[[1]]))
  }
  
  
  # Prepare data and element to load into TMB function PO
  n.time=length(ypo)
  K = length(ypo) # NB of years/time periods
  ob_wt=ob_wt.po
  
  # Need to standardize data:
 # for (i in 1:K) {
#    Xpo[[i]][,2:dim(Xpo[[i]])[2]] = scale(Xpo[[i]][,2:dim(Xpo[[i]])[2]]) # to not get the intercept scaled
#  }
#  Wpo = lapply(Wpo, scale)
  
  
  # set up coefficient for variables not to shrink
  noshrink = all_names = rep(list(NULL), n.time)
  
  colname.bias=c()
  no.shrink = rep(list(c()), n.time)
  for (k in 1:n.time) {
    colname.bias[k] = colnames(Wpo[[k]])
    no.shrink[[k]] =  c("Intercept", noshrink[[k]])
    
    all_names[[k]]  = c(colnames(Xpo[[k]]), colname.bias[[k]])
    no.shrink[[k]]   = which(all_names[[k]] %in% no.shrink[[k]])
  }
  
  noshrink = no.shrink
  lasso_po = rep(list(c()), n.time)
  
  num_env = num_bias = num_p = lassomat.start = rep(list(list()), n.time)
  if(length(lamb)>1){
    if (class(lamb)=="list"){
      for (i in 1:K) {
        lassomat.start[[i]] = lamb[[i]][-Wocc.pos]
      }
      lassomat = unlist(lassomat.start)
    }else{
      lassomat= lamb[-Wocc.pos]
    }
  }else{
    for (i in 1:n.time) {
      num_env[[i]]  = dim(Xpo[[i]])[2]
      num_bias[[i]] = unlist(dim(Wpo[[i]])[2])
      num_p[[i]]    = num_env[[i]] + sum(num_bias[[i]])
    }
    
    lassomat = (if(length(lamb) == 1) rep(lamb, num_p[[1]])
                else rep(0, num_p[[1]]))
    
  }
  
  
  # Prepare some usefool data before calculation in TMB
  # We cannot give TMB lists so I rbind every data
  b_init=array(NA, c(length(par[[1]]),1,n.time))
  for (k in 1:n.time) {
    b_init[,,k] = par[[k]]
  }
  
  b_n=lassostart= Lsize = c()
  for (i in 1:K) {
    b_n[i] = length(b_init[1:(dim(b_init)[1]-Wocc_NB),,i])
    Lsize[i] = b_n[i]
  }
  lassostartk = c(0, cumsum(b_n))[1:K] # need to be -1 to match TMB number system (starts at 0)
  lassoendk   = cumsum(b_n)
  
  
  Xpob=Xpoe=Xpon=list()
  
  for (i in 1:k) {
    Xpob[[i]] = Xpo[[i]][,-c(extinct_comp,coloni_comp)]
    Xpoe[[i]] = as.matrix(Xpo[[i]][,extinct_comp])
    Xpon[[i]] = as.matrix(Xpo[[i]][,coloni_comp])
  }
  
  # rbind the different list to read them in TMB
  
  Xpobmat = do.call(rbind, Xpob)
  Xponmat = do.call(rbind, Xpon)
  Xpoemat = do.call(rbind, Xpoe)
  
  Wpomat = do.call(rbind, Wpo)
  
  
  Xpo_n=Xpostart = row_startk= row_endk = c()
  for (i in 1:K) {
    Xpo_n[i] = dim(Xpo[[i]])[1]
  }
  row_startk = c(0, cumsum(Xpo_n))[1:K]
  row_endk   = cumsum(Xpo_n)
  
  
  for (i in 1:K) {
    ypo[[i]] = as.matrix(ypo[[i]]) 
  }
  
  ypomat = do.call(rbind, ypo)
  
  
  ob_wtmat = do.call(c, ob_wt)
  
 
  L_NB = c()
  for (i in 1:K) {
    #Vpo_NB[i] = dim(Xpo[[i]])[2]
    L_NB[i] = dim(Xpo[[i]])[1]
  }
  
  Vpo_NB= dim(Xpo[[1]])[2]
  Wpos = dim(Wpomat)[2]

return(list(K=K, Xocc=Xocc, Wocc=Wocc, site_area=site_area, yocc=yocc, R=R, Jm=Jm, 
            S_NB=S_NB, Vocc_NB=Vocc_NB, lasso_occ=lasso_occ, 
            extinct_comp=extinct_comp, coloni_comp=coloni_comp,
            Xpobmat=Xpobmat, Xpoemat=Xpoemat, Xponmat=Xponmat, Wpomat=Wpomat, 
            ypomat=ypomat, Wpos=Wpos, Vpo_NB=Vpo_NB, L_NB=L_NB, 
            wtvecocc = wtvecocc, wtvecpo = wtvecpo,
            row_startk=row_startk, lassomat=lassomat, lassostartk=lassostartk,
            Lsize=Lsize, Wocc_NB = Wocc_NB, ob_wtmat=ob_wtmat,b_init=b_init))
}


plotfit = function(fit, pred_data = NULL, xy = NULL, z = "intensity", model = NULL, source = NULL, link = "logit",
                   coord = c("X", "Y"), asp = "iso", ylab = "", xlab = "", col.regions = heat.colors(1024)[900:1], 
                   cuts = length(col.regions), cex = 1.4, main.text = NULL, cex.color = 1.4, Sscale=FALSE,
                   extinct_comp=NULL, coloni_comp=NULL)
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
  
  if (is.null(fit$n.time)==TRUE){
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
  }else{
    num_env = num_p=n.components= bias_ind1= bias_ind2=env_ix= list()
    num_bias= rep(list(c()), fit$n.time)
    bias_ix = rep(list(list()), fit$n.time)
    for (i in 1:fit$n.time) {
      num_env[[i]]  = dim(fit$X_env[[i]][[1]])[2]
      num_bias[[i]] = unlist(lapply(lapply(fit$X_bias[[i]], as.matrix), ncol))
      num_p[[i]]    = num_env[[i]] + sum(num_bias[[i]])
      n.components[[i]] = length(num_bias[[i]])
      bias_ind1[[i]] = num_env[[i]] + 1 + cumsum(num_bias[[i]]) - num_bias[[i]]
      bias_ind2[[i]] = num_env[[i]] + cumsum(num_bias[[i]])
      env_ix[[i]] = rep(list(1:num_env[[i]]), n.components[[i]])
      
      for (j in 1:n.components[[i]])
      {
        bias_ix[[i]][[j]] = bias_ind1[[i]][j]:bias_ind2[[i]][j]
      }
    }
  }
  
  num_start = c(1, cumsum(num_p)+1)[1:fit$n.time]
  num_end   = cumsum(num_p)
  
  if (is.null(fit$n.time)==TRUE){
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
  }else{
    if (fit$dat.type[source] == "PO"){
      Xdes = list()
      beta_list = beta_sub = rep(list(c()), fit$n.time)
      
      for (i in 1:fit$n.time) {
        if(is.null(extinct_comp)==TRUE & is.null(coloni_comp)==TRUE){
          Xdes[[i]] = if (z == "intensity" | z == "occupancy") fit$X_env[[i]][[source]][fit$y[[i]][[source]] == 0,] else fit$X_bias[[i]][[source]][fit$y[[i]][[source]] == 0,]
          beta_list[[i]] = beta[num_start[i]:num_end[i]]
          beta_sub[[i]] = if (z == "intensity" | z == "occupancy") beta_list[[i]][env_ix[[source]][[i]]] else beta_list[[i]][bias_ix[[i]][[source]]]
          
        }else{
          Xdes[[i]] = if (z == "intensity" | z == "occupancy") fit$X_env[[i]][[source]][fit$y[[i]][[source]] == 0,] else fit$X_bias[[i]][[source]][fit$y[[i]][[source]] == 0,]
          Xdes[[i]] = Xdes[[i]][,-c(extinct_comp, coloni_comp)]
          beta_list[[i]] = beta[num_start[i]:num_end[i]]
          beta_sub[[i]] = if (z == "intensity" | z == "occupancy") beta_list[[i]][env_ix[[i]][[source]]] else beta_list[[i]][bias_ix[[i]][[source]]]
          beta_sub[[i]] = beta_sub[[i]][-c(extinct_comp, coloni_comp)]
        }
        
      }
      
      XY = fit$quad_xy[[1]][[1]]
    }
    
    if (fit$dat.type[source] == "Occ")
    {
      Xdes = list()
      beta_list = beta_sub = rep(list(c()), fit$n.time)
      
      for (i in 1:fit$n.time){
        if(is.null(extinct_comp)==TRUE & is.null(coloni_comp)==TRUE){
          Xdes[[i]] = if (z == "intensity" | z == "occupancy") fit$X_env[[i]][[source]] else fit$X_bias[[i]][[source]]
          beta_list[[i]] = beta[num_start[i]:num_end[i]]
          beta_sub[[i]] = if (z == "intensity" | z == "occupancy") beta_list[[i]][env_ix[[i]][[source]]] else beta_list[[i]][bias_ix[[i]][[source]]]
          
        }else{
          Xdes[[i]] = if (z == "intensity" | z == "occupancy") fit$X_env[[i]][[source]] else fit$X_bias[[i]][[source]]
          Xdes[[i]] = Xdes[[i]][,-c(extinct_comp, coloni_comp)]
          beta_list[[i]] = beta[num_start[i]:num_end[i]]
          beta_sub[[i]] = if (z == "intensity" | z == "occupancy") beta_list[[i]][env_ix[[i]][[source]]] else beta_list[[i]][bias_ix[[i]][[source]]]
          
        }
        
      }
      
      XY = fit$sp_xy[[1]][[source]]
    }
  }
  
  
  if(is.null(fit$n.time)==TRUE){
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
  }else{
    plot_z = lin_det = list()
    for (i in 1:fit$n.time) {
      if (z == "intensity" | z == "bias")
      {
        plot_z[[i]] = exp(as.matrix(Xdes[[i]]) %*% beta_sub[[i]])
      }
      if (z == "occupancy")
      {
        plot_z[[i]] = 1 - exp(-exp(as.matrix(Xdes[[i]]) %*% beta_sub[[i]])*fit$site.area) # need site.area from fit!
      }
      
      if (z == "p_detect")
      {
        lin_det[[i]]  = as.matrix(Xdes[[i]]) %*% beta_sub[[i]]
        if (link == "logit")
        {
          plot_z[[i]] = expit(lin_det[[i]])
        }
        if (link == "cloglog")
        {
          plot_z[[i]] = clogloginv(lin_det[[i]])
          #plot_z = 1 - exp(-exp(lin_det))
        }
      }
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
  
  if(is.null(fit$n.time)==TRUE){
    if (is.null(main.text))
    {
      main.text = paste(z_name, " from Source ", source, "\n", "Model with penalty ", round(plot_pen, 3), sep = "")
    }
    
    levelplot(plot_z ~ plot_x + plot_y, asp = asp, 
              ylab = ylab, xlab = xlab, col.regions = col.regions, cuts = cuts, 
              main = list(main.text, cex = cex), 
              scales = list(y = list(draw = FALSE), x = list(draw = FALSE)), 
              cex = cex, colorkey = list(labels = list(cex = cex.color)))
    
  }else{
    #main.text = list()
    #if (is.null(main.text))
    #{
    #  for (i in 1:n.time) {
    #    main.text[[i]] = paste(z_name, " from Source ", source, "\n", "Model with penalty ", round(plot_pen[[i]], 3),
    #                           "\n", "Year", seq(1, n.time, 1)[i], sep = "")
    #  }
    #}
    
    if (Sscale==TRUE){
      minVal=min(unlist(plot_z))
      maxVal=max(unlist(plot_z))
      Med = median(unlist(plot_z))
      
      Lplot = list()
      Lplotname = rep(NA, fit$n.time)
      for (i in 1:fit$n.time){
        Lplot[[i]] = levelplot(plot_z[[i]] ~ plot_x + plot_y, asp = asp, 
                               ylab = ylab, xlab = xlab, col.regions = col.regions, #cuts = cuts, 
                               main = list(paste(z_name, " from Source ", source, "\n", "Model with penalty ", round(plot_pen[[i]], 3),
                                                 "\n", "years", i, sep = ""), cex = cex), 
                               at = unique(c(seq(minVal, Med, length=10),
                                             seq(Med, maxVal, length=30))),
                               scales = list(y = list(draw = FALSE), x = list(draw = FALSE)), 
                               cex = cex, colorkey = list(labels = list(cex = cex.color)))
        #Lplotname[i] = rep(paste("Lplot[[", i, "]]", sep = ""), n.time)
      }
    }else{
      Lplot = list()
      for (i in 1:fit$n.time){
        Lplot[[i]] = levelplot(plot_z[[i]] ~ plot_x + plot_y, asp = asp, 
                               ylab = ylab, xlab = xlab, col.regions = col.regions, cuts = cuts, 
                               main = list(paste(z_name, " from Source ", source, "\n", "Model with penalty ", round(plot_pen[[i]], 3),
                                                 "\n", "years", i, sep = ""), cex = cex), 
                               scales = list(y = list(draw = FALSE), x = list(draw = FALSE)), 
                               cex = cex, colorkey = list(labels = list(cex = cex.color)))
        #Lplotname[i] = rep(paste("Lplot[[", i, "]]", sep = ""), n.time)
      }
    }
    
    Lplot
  }
  
  
  
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
  
  if(is.null(fit$n.time)==TRUE){
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
    
  }else{
    num_env = num_p=n.components= bias_ind1= bias_ind2=env_ix= list()
    num_bias= rep(list(c()), fit$n.time)
    bias_ix = rep(list(list()), fit$n.time)
    for (i in 1:fit$n.time) {
      num_env[[i]]  = dim(fit$X_env[[i]][[1]])[2]
      num_bias[[i]] = unlist(lapply(lapply(fit$X_bias[[i]], as.matrix), ncol))
      num_p[[i]]    = num_env[[i]]-length(v.cut)/2 + sum(num_bias[[i]])
      n.components[[i]] = length(num_bias[[i]])
      bias_ind1[[i]] = num_env[[i]] + 1 + cumsum(num_bias[[i]]) - num_bias[[i]]
      bias_ind2[[i]] = num_env[[i]] + cumsum(num_bias[[i]])
      env_ix[[i]] = rep(list(1:num_env[[i]]), n.components[[i]])
      
      for (j in 1:n.components[[i]])
      {
        bias_ix[[i]][[j]] = bias_ind1[[i]][j]:bias_ind2[[i]][j]
      }
    }
    num_start = c(1, cumsum(num_p)+1)[1:fit$n.time]
    num_end   = cumsum(num_p)
    
    betas.l = rep(list(matrix()), fit$n.time)
    min.y = max.y = list()
    for (i in 1:fit$n.time) {
      betas.l[[i]] = betas[num_start[[i]]:num_end[[i]],]
      min.y[[i]] = min(apply(betas.l[[i]], 1, min))
      max.y[[i]] = max(apply(betas.l[[i]], 1, max))
    }
    
    lambdas = rep(list(c()), fit$n.time)
    lambda_d = list()
    for (i in 1:fit$n.time) {
      lambdas[[i]] = fit$penalty_vec[[i]]
      if (lambdas[[i]][1] == 0)
      {
        lambda_d[[i]] = diff(log(lambdas[[i]][2:3]))
        lambdas[[i]][1] = exp(log(lambdas[[i]][2]) - lambda_d[[i]])
      }
    }
    
    
    for (i in 1:fit$n.time) {
      plot(lambdas[[i]], betas.l[[i]][1,], log = "x", type = "l", ylim = c(min.y[[i]], max.y[[i]]), 
           xlab = "LASSO penalty", ylab = "Coefficients", xlim= c(min(lambdas[[i]]), max(lambdas[[i]])))
      for (j in 2:dim(betas.l[[i]])[1])
      {
        points(lambdas[[i]], betas.l[[i]][j,], type = "l")
      }
      
      for (j in 1:length(unique.beta))
      {
        abline(v = lambdas[[i]][unique.beta[j]], lwd = 3, col = colors[j])
      }	
      
      legend("topright", paste(names, "year", i, sep = ""), lwd = rep(3, length(unique.beta)), col = colors[1:length(unique.beta)])
    }
    
  }
  
  
}

criterion_curve = function(fit, criterion = "BIC")
{
  plot_y = fit$criterion.matrix[,match(criterion, names(fit$criterion.matrix))]
  
  if(is.null(fit$n.time)==TRUE){
    lambdas = fit$penalty_vec
    
    if (lambdas[1] == 0)
    {
      lambda_d = diff(log(lambdas[2:3]))
      lambdas[1] = exp(log(lambdas[2]) - lambda_d)
    }
    min_pen = fit$penalty_vec[which.min(plot_y)]
    
    plot(lambdas, plot_y, xlab = "LASSO Penalty", ylab = criterion, type = "o", log = "x")
    points(lambdas[which.min(plot_y)], min(plot_y), col = "green3", pch = 19)
    legend("topleft", legend = bquote("Min" ~ .(criterion): ~ lambda == .(min_pen) ~ phantom(x)  ), col = "green3", pch = 19)
    
  }else{
    lambdas=lambda_d=rep(list(c()),fit$n.time)
    min_pen=list()
    for (i in 1:fit$n.time) {
      lambdas[[i]] = fit$penalty_vec[[i]]
      
      if (lambdas[[i]][1] == 0)
      {
        lambda_d[[i]] = diff(log(lambdas[[i]][2:3]))
        lambdas[[i]][1] = exp(log(lambdas[[i]][2]) - lambda_d[[i]])
      }
      
      min_pen[[i]] = fit$penalty_vec[[i]][which.min(plot_y)]
      
      plot(lambdas[[i]], plot_y, xlab = "LASSO Penalty", ylab = criterion, type = "o", log = "x")
      points(lambdas[[i]][which.min(plot_y)], min(plot_y), col = "green3", pch = 19)
      legend("topleft", legend = bquote("Min" ~ .(criterion): ~ lambda == .(min_pen[[i]]) ~ phantom(x)  ), col = "green3", pch = 19)
      
    }
    
  }
  
}