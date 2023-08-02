# netQDA
An R package for Local Network Guided Quadratic Discriminant Analysis 

```{r setup}
library(devtools)
install_github("XuepingZhou/netQDA")
library(netQDA)
```

# Simulate Data
```{r echo= TRUE}
seed = 1111
N.grp = 200 # sample size 
p = 1000   # number of predictors, including noises 
p_info = 8  # number of informative fewatures 
p.noise = p - p_info # number of noises 
p.grp = 4   # each group/connected component contains 4 features 
rho = 0.7   # 
p1 = 4      # group 1, contains 4 informative features 
p2 = 2       # group 2, contains 2 informative features 
rho1 = 0.5   # group 1, informative features correlation
rho2 = 0.5    # group 2, informative features correlation
```

## precision matrix (omega) and covariance structure (cov) of connected component from star structure
```{r echo= TRUE}
omega1 <- matrix(0, p.grp, p.grp)
omega2 <- matrix(0, p.grp, p.grp)
omega1[1:p1, 1] <- rho1 ; omega1[1, 1:p1 ] <- rho1
omega2[ 1:p2, 1] <- rho2 ; omega2[1, 1:p2 ] <- rho2
diag(omega1) <- 1
diag(omega2) <- 1
cov1 <- solve(omega1)
cov2 <- solve(omega2)
```

## simulate informative features for the two groups
```{r echo= TRUE}
set.seed(seed)
X1.informative.1 = MASS::mvrnorm(n = N.grp,  mu = rep(0, p.grp) , Sigma =cov1 )
X1.informative.1 = scale( X1.informative.1 , scale = T,center = FALSE)
X1.informative.2 = MASS::mvrnorm(n = N.grp,  mu = rep(0, p.grp) , Sigma =cov1 )
X1.informative.2 = scale( X1.informative.2 , scale = T,center = FALSE)
X1.informative  <- cbind( X1.informative.1,    X1.informative.2)
X1.informative  <- scale(X1.informative ,
                           center = T, scale =  apply(X1.informative, 2, sd, na.rm = TRUE))
```

```{r echo= TRUE}
X2.informative.1 = MASS::mvrnorm(n = N.grp, mu = c(1.5,0,0,0), Sigma = cov2 )
X2.informative.1 = scale( X2.informative.1 , scale = T,center = FALSE)
X2.informative.2 = MASS::mvrnorm(n = N.grp, mu = c(1.5,0,0,0), Sigma = cov2 )
X2.informative.2 = scale( X2.informative.2, scale = T,center = FALSE)
X2.informative <- cbind( X2.informative.1,    X2.informative.2)
X =  rbind( X1.informative,  X2.informative)
```

## Simulate noise part  
```{r echo= TRUE}
X.noise.1 <-  matrix(rnorm(N.grp*(p.noise ), mean =0 , sd = 1), nrow = N.grp,
                     ncol= p.noise, byrow = T)
X.noise.2 <-  matrix(rnorm(N.grp*(p.noise ), mean =0 , sd = 1),
                     nrow = N.grp,
                     ncol= p.noise, byrow = T)
X.noise = rbind(X.noise.1 , X.noise.2 )
X  = cbind(X,  X.noise)
colnames(X) <- 1:dim(X)[2]

X.train = X[ c(1:100, 201:300), ]
X.new = X[ c(101:200, 301:400), ]
```

## Outcome class membership
```{r echo= TRUE}
Y.train = c(rep(1,100), rep(2,100))
Y.new = c(rep(1,100), rep(2,100))
```



# Train Model using netQDA
```{r echo= TRUE}
X = scale(X.train)
COV.COR <- individual_covcor(X.train, Y.train)
COV.list <- COV.COR[["COV.list"]]
COR.list <- COV.COR[["COR.list"]]
COV.scale.list <- COV.COR[["COV.scale.list"]]
try <-  netQDAtrain(X = X, y = Y.train,
                    COV.list=COV.list,
                    COV.scale.list=COV.scale.list,
                    COR.list=COR.list,
                    tau =2 , pair = c(1,2) , alpha = 0.5, nu = 4,  delta =8,
                    d = 1, nb = 5)
```

# Result of selected potential informative features 
```{r echo= TRUE}
try$screenset
```

# Classification of New Data


```{r echo= TRUE}
X.new <- scale( X.new ,  center = attr(X, 'scaled:center'),  scale = attr(X, 'scaled:scale'))
results <- netQDApredict(TrainModel= try, X.new)
table(results$predClass, Y.new)
pROC::roc(Y.new, results$predClass)
```




