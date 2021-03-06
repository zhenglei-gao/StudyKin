Trust Region Reflective?
========================================================
## To improve speed !

### Using Cpp

* https://github.com/hadley/devtools/wiki/Rcpp
* http://www.cplusplus.com/articles/1ywTURfi/
* http://en.wikipedia.org/wiki/Compatibility_of_C_and_C%2B%2B

### Better finite differece

<u>some text</u>

* Quoted from https://stat.ethz.ch/pipermail/r-help/2008-March/158228.html

If 'myfunc' is a vector function and can be vectorized in R, then it 
is even faster to use the following:

```r


grad.vec <- function(x, fn, ..., eps = sqrt(.Machine$double.neg.eps)) {
    x1 <- x + eps * pmax(abs(x), 1)
    x2 <- x - eps * pmax(abs(x), 1)
    (fn(x1, ...) - fn(x2, ...))/(x1 - x2)
}

grad.1 <- function(x, fn) {
    x <- sort(x)
    x.e <- head(embed(x, 2), -1)
    y.e <- embed(fn(x), 3)
    hh <- abs(diff(x.e[1, ]))
    apply(y.e, 1, function(z) (z[1] - z[3])/(2 * hh))
}

myfunc.1 <- function(x) {
    (exp(x) - x)/10
}

p0 <- rexp(1000)
system.time(for (i in 1:500) out1 <- grad.1(p0, myfunc.1))
```

```
##    user  system elapsed 
##    2.36    0.00    2.36
```

```r
system.time(for (i in 1:500) out2 <- grad.vec(p0, myfunc.1))
```

```
##    user  system elapsed 
##    0.12    0.00    0.13
```

```r

```


* http://www.mathworks.com/matlabcentral/answers/7898

## New Implememtaion and Test Results!


```r
library(KineticEval)
setwd("C:/Projects2013/KinEvalGit/")
data(schaefer07_complex_case)
model <- mkinmod.full(parent = list(type = "SFO", to = c("A1", "B1", "C1"), 
    sink = FALSE), A1 = list(type = "SFO", to = "A2"), B1 = list(type = "SFO"), 
    C1 = list(type = "SFO"), A2 = list(type = "SFO"), data = schaefer07_complex_case)

Fit2 <- NULL
alglist <- c("L-BFGS-B", "Marq", "Port", "spg", "solnp", "LM")
for (i in 1:6) {
    Fit2[[i]] <- mkinfit.full(model, plot = FALSE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}
## Save Fit2
names(Fit2) <- alglist
save(Fit2, file = "Fit_Schaefer.rda")


(lapply(Fit2, function(x) x$par))
unlist(lapply(Fit2, function(x) x$ssr))
##########################

##########################
load("ex1.rda")
load("ex1_a.rda")
# data(ex1_a)
alglist <- c("L-BFGS-B", "Marq", "Port", "spg", "solnp", "LM")
Fit <- NULL
for (i in 1:6) {
    Fit[[i]] <- mkinfit.full(ex1, plot = TRUE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}
names(Fit) <- alglist
(lapply(Fit, function(x) x$par))
kinplot(Fit[[6]])
unlist(lapply(Fit, function(x) x$ssr))
Fit_a <- NULL
for (i in 1:6) {
    Fit_a[[i]] <- mkinfit.full(ex1_a, plot = TRUE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}
unlist(lapply(Fit_a, function(x) x$ssr))

## Test on data andrew:
data(andrew)
mkinmodini <- mkinmod.full(Parent = list(type = "SFO", to = "Metab"), Metab = list(type = "SFO", 
    M0 = list(ini = 0, fixed = 0, lower = 0, upper = Inf)), data = andrew)
for (i in 1:6) {
    Fit[[i]] <- mkinfit.full(mkinmodini, plot = FALSE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}
names(Fit) <- alglist
(lapply(Fit, function(x) x$par))
unlist(lapply(Fit, function(x) x$ssr))

#####################
load("ex2.rda")
Fit2 <- NULL
alglist <- c("L-BFGS-B", "Marq", "Port", "spg", "solnp", "LM")
for (i in 1:5) {
    Fit2[[i]] <- mkinfit.full(ex2, plot = FALSE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}

```




## Trust Region Theory and Implementation

### Advantages:

The primary advantage of trust region methods is stability. 

*  requires gradient and Hessian of the objective function

* is guaranteed to converge to a point satisfying the first and second order necessary conditions for a local optimum (gradient zero and Hessian positive semidefinite for minimization or negative semidefinite for maximization). nlm and optim come with no guarantees and often perform badly on difficult problems. 

* handles restricted domains so long as the solution is in the interior of the domain. 

* does not do global optimization, but neither does anything else in R.

### General Ideas:

Suppose $f : R^n \to R$ is a function we wish to minimize. The basic idea of trust region is to approximate $f$ with a simpler function $q$ (a quadradtic model) at the neighborhood $N$(trust region) aound point $x$. 

$$q_k(p)=f(x_k)+\nabla{f(x_k)}^Tp+\frac{1}{2}p^T\nabla^2{f(x_k})p$$.

A trial step is caclculated by minimizing $q(s)$. This is the trust region supproblem,

$$\min \{q(s), s \in N\}$$

where $s\in N$ is defined by $\|D_is\|\leq\Delta$, $\Delta >0$ is the radius of the trusted region. So the subproblems is to solve

$$\min q_k(s) \text{ subject to } \|Ds\|\leq\Delta$$,

Suppose $p_k$ is a solution to the trust region subproblem. The adjustment of $\Delta_k$ is done as follows
$$
   \rho_k = \frac{f(x_k) - f(x_k + p_k)}{q_k(0) - q_k(p_k)}
$$
which is the actual decrease in the objective function $f$ in the step
compared to the predicted decrease using $q_k$.  If $\rho_k$ is small
or negative, then $q_k$ is a bad model at $x_k + p_k$ so the trial step should not be used and the trust region radius should be adjusted.

### Trust Region Subproblem(TRS)

A vector $p^*$ is a global solution to the trust region subproblem if and only if
$$

   \lVert p^* \rVert \le \Delta

$$
and there exists a scalar $\lambda \ge 0$ such that
$$
\begin{gather}
   (B + \lambda I) p^* = - g
   \label{eq:derivative}
   \\
   \lambda = 0 \quad \text{or} \quad \lVert p^* \rVert = \Delta
  
   \\
   B + \lambda I \ \text{is positive semidefinite}
   
\end{gather}
$$
where $B_k = \nabla^2 f(x_k)$. 

* Case $\lambda=0$: If $B$ is positive definite and
$p^* = - B^{-1} g$ satisfies $\|p^*\|\le \Delta$, then that is the solution.

* Case $\|p^*\| = \Delta$, Define
$$
\begin{equation} \label{eq:plambda}
   p(\lambda)
   =
   - (B + \lambda I)^{-1} g
   =
   - \sum_{j = 1}^{n} \frac{q_j^T g}{\lambda_j + \lambda} q_j
\end{equation}
$$
where $\lambda_j$ are the eigenvalues of $B$ and $q_j$ are the corresponding orthonormal eigenvectors (this is valid only when $\lambda \neq - \lambda_j$ for all $j$). Then
$$
   \lVert p(\lambda) \rVert^2
   =
   \sum_{j = 1}^{n} \left( \frac{q_j^T g}{\lambda_j + \lambda} \right)^2
$$
The analysis again splits into several cases.

  1. Let $\lambda_{\text{min}}$ denote the minimum eigenvalue of $B$. If
  $$ (q_j^T g) \neq 0, \qquad \text{for some $j$ such that
   $\lambda_j = \lambda_{\text{min}}$}
   $$
   then $p(\lambda)$ is a continuous, strictly decreasing function on the open interval $(- \lambda_{\text{min}}, \infty)$ and goes to $\infty$
as $\lambda \to - \lambda_{\text{min}}$ and to zero as $\lambda \to \infty$. Thus a $\lambda^*$ such that $\lVert p(\lambda^*) \rVert^2 = \Delta^2$ exists and is unique, and $p(\lambda^*)$ is the solution to the trust region subproblem.
  2. In the other case, $p(\lambda)$ is continuous at $\lambda = \lambda_{\text{min}}$ and now defines a continuous strictly decreasing function on
the closed interval $[- \lambda_{\text{min}}, \infty)$, and $\lambda$
must be in this interval in order for $B+\lambda I$ being positive definite.

If $\lVert p(- \lambda_{\text{min}}) \rVert^2 > \Delta^2$, then there is still a unique $\lambda^*$ in $(- \lambda_{\text{min}}, \infty)$
such that $\lVert p(\lambda^*) \rVert^2 = \Delta^2$, this $p(\lambda^*)$ is the solution to the TRS. Otherwise, we must have $\lambda = - \lambda_{\text{min}}$.Then $B + \lambda I$ is singular, and the $p(\lambda)$ defined above can no longer be used. Now solutions of $(B+\lambda_i)p=-g$ 
are non-unique, and we have
$$
\begin{equation} 
   p_{\text{hard}}(\tau)
   =
   - \sum_{\substack{j = 1 \\ \lambda_j \neq \lambda_{\text{min}}}}^n
   \frac{q_j^T g}{\lambda_j - \lambda_{\text{min}}} q_j
   +
   \sum_{\substack{j = 1 \\ \lambda_j = \lambda_{\text{min}}}}^n
   \tau_j q_j
\end{equation}
$$
is a solution of $(B+\lambda_i)p=-g$ for every vector $\tau$.
Since the first sum on the right hand side has
norm less than $\Delta$ by the falsity of $\lVert p(- \lambda_{\text{min}}) \rVert^2 > \Delta^2$, we can choose $\tau$ to make $\lVert p_{\text{hard}}(\tau) \rVert = \Delta$. The nonuniqueness is not an issue, because (at least as far as the subproblem is concerned) one solution is just as good as another.

* Steihuang(above procedure)
* Dogleg, STIR

### Practical Newton Methods - Trust-Region Newton Methods

* Newton-dogleg
* Newton-2D subspace 
* Newton- Nearly Exact Solution
* Trust Region Newton-(P)CG

### Finite Difference ?

* numerical gradients calculated by finite-difference approximation.

* http://terminus.sdsu.edu/SDSU/Math693a_f2012/?r=sched   (Best reference lecture notes for me on numerical optimization)

### Tricks in Matlab `lsqnonlin` implementation

The TRR in matlab involves the following steps:

1. Consider the trust region subproblem of finding trial step $s_i$ between the $i$-th point $x_i$ and the next point $x_{i+1}$ that minimizes $q_i(s)$. Formulate the 2-dimentional subspace $V$ with 2 spanning vectors $v_1$ and $v_2$. One vector is in the direction of the gradient $x_i$ and the other either in the approximate Gauss-Newton direction, i.e, a solution to $H\cdot v_2=-g$ or a direction of negative curvature($v_2^T\cdot H \cdot s_2 < 0$).

  For nonlinear least squares, suppose $f=\sum_i f_i(x)^2$ , let $F=(f_1(x), ..., f_n(x))$, then $f=\|F\|$. The structure of nonlinear least square problems can be exploited to enhance efficiency. In particular, the approximate Gauss-Newton direction, i.e. a solution $s$ to $\min\|Js+F\|^2_2$ is used to help define the 2-dimentional subspace $V$. In each iteration, the method of preconditioned conjugate gradients(PCG) is used to approximately solve the normal equation, i.e. $J^TJs=-J^TF$, although the normal equations are not explicitly fomed

2. Solve for the trial step $s_i$ of the 2-dimentional subproblem, i.e, solve the following equation: 
$$\min\{ \frac{1}{2}s^THs+s^Tg, \text{ subject to } \|D_is\|\leq\Delta_i\} $$
where $g$ is the gradient and $H$ is the Hessian matrix at the current point.

3. If for a predefined constant $0 < \tau < 1$, $f(x_i+s) < f(x_i)$, then $x_{i+1}$ becomes the current point, otherwise $x_{i+1}=x_i$.

4. Update $\Delta_i$. In particular, it is decreased if the trial step is not accepted.

These 4 steps are repeated until convergence, i.e., if $\nabla f(s_i)$ is below a chose tolerance level, the algorithm ends.

PCG uses a preconditioner matrix $P$ that makes the equation $Ax=b$ easier to solve by turning it into $P^{-1}(Ax-b)=0$

    lsqnonlin -> snls -> sfdnls, trdog
    
### Algorithm STIR(subspace trust region interior reflective)

Let $0<\mu<\eta<1, 0<\Lambda_l<\Lambda_u$, and $\gamma_1<1<\gamma_2$ be given. Let $x_0\in\text{int}(F)$, $\Delta_0<\Lambda_u$.

For $k=0,1,...$

1. Compute $f_k, g_k, D_k, H_k$, and $C_k$; define the quardratic model
$$
\psi_k(s)=g_k^Ts+\frac{1}{2}s^T(H_k+C_k)s.
$$

2. Compute a step $s_k$, with $s_k+x_k\in\text{int}(F)$, based on the subspace problem
$$
\min \psi_k(s) \text{ subject to } \|D_ks\|_2\leq\Delta_k, s\in S_k
$$
where the subspace $S_k$ is set up as below.

3. Compute 
$$
\rho_k=\frac{f(x_k+s_k)-f(x_k)+\frac{1}{2}s_k^TC_ks_k}{\psi_k(s_k)}
$$

4. If $\rho_k>u$ then set $x_{k+1}=x_k+s_k$. Otherwise set $x_{k+1}=x_k$.

5. Update $\Delta_k$ as specified below.

**Determin Subspace $S_k$**

**Updating Trust Region Size $\Delta_k$** 

1. If $\rho_k<\mu$ then set $\Delta_{k+1}\in (0, \gamma_1\Delta_k]$.
2. If $\rho_k\in(\mu,\eta)$ then set $\Delta_{k+1}\in [\gamma_1\Delta_k\, \Delta_k]$.
3. If $\rho_k\geq \eta$ then 
   * if $\Delta_k>\Lambda_l$ then set $\Delta_{k+1}\in$ either $(\gamma_1\Delta_k\, \Delta_k]$ or $[\Delta_k, \gamma_2\Delta_k]$.
   * otherwise, set $\Delta_{k+1}\in[\Delta_k, \min(\gamma_2\Delta_k,\Lambda_u)]$



### Some other's opinions

Lsqnonlin is a modern implementation of the Levenberg-Marquardt(LM) method, which is the most popular algorithm for nonlinear least squares.

Global convergence from arbitrary starting points is promoted via a trust region technique: at each iteration a quadratic surrogate model of the actual problem is constructed, and solved subject to a trust region - a ball around the current iterate, where the surrogate model is likely to be a good approximation of the true problem function.

There's also a more classical line-search implementation of LM in lsqnonlin - to run it you have to change an option (see the link below for details).

If the problem has bounds on the variables (say, all the parameters between zero and one), lsqnonlin uses an interior point method, which a modern approach to handling inequality constraints.

More general constraints other than bounds are handled by the solver fmincon.

If the problem is large scale (this doesn't happen very often in practice for nonlinear least squares), lsqnonlin has the option to use/estimate sparse Jacobians, or work with a Jacobian-times vector function instead of the Jacobian itself.

### Lectures and Notes

* http://www.mcs.anl.gov/~anitescu/CLASSES/2012/LECTURES/S310-2012-lect5.pdf
* http://www.cs.umd.edu/users/oleary/a607/2008/index.html

### Other's implementation

* https://bitbucket.org/carandraug/octave/src/13d1e9bfa362/scripts/optimization/lsqnonneg.m?at=default
* http://homepages.rpi.edu/~mitchj/pack.html
* https://github.com/dkogan/libdogleg



* http://mind.cog.jhu.edu/courses/680/octave/Installers/Octave/Octave.OSX10.6/Applications/MATLAB_R2009b.app/toolbox/optim/optim/

-------------------------

## Differential Equation and Jacobian!!

* http://leto.net/x/2008/11/5-minute-math-lesson-what-is-a.html

---------------------------

## Using Sparse Matrice in R


```r
library("Matrix")

m1 <- matrix(0, nrow = 1000, ncol = 1000)
m2 <- Matrix(0, nrow = 1000, ncol = 1000, sparse = TRUE)

object.size(m1)
# 8000200 bytes
object.size(m2)
# 5632 bytes

m1[500, 500] <- 1
m2[500, 500] <- 1

object.size(m1)
# 8000200 bytes
object.size(m2)
# 5648 bytes

```

The full matrix representation does not change in size because all of the zeros are being represented explicitly, while the sparse matrix is conserving that space by representing only the non-zero entries. With a simple subtraction, we find that adding one additional non-zero entry increases the storage requirement by 16 bytes. We can conclude that setting all of the entries to non-zero values would require 5632 + 16 * 1000 * 1000 bytes, which is 16,005,632 bytes or almost exactly twice the amount of storage required to use the full representation implemented by the matrix class.

Of course, the take away lesson is that sparse matrices are very efficient if your data is sparse and mildly wasteful if your data is not sparse.

Beyond simple initialization and assignment operations, we can perform quite a few other operations on objects of class Matrix, including vector multiplication, matrix addition and subtraction and transposition. In addition you can do binding operations using objects of class Matrix with the cBind and rBind functions.

An alternative to the Matrix package is the slam package by Kurt Hornik and others. The sparse matrices generated using this package can be noticeably smaller than those generated by the Matrix package in some cases.

* http://www.johnmyleswhite.com/notebook/2011/10/31/using-sparse-matrices-in-r/

---------------------------

## Least Square Special

Chain rule can be applied to reduce amount of computation.

The derivatives are calculated numerically using Richardson improvement
The the first order derivative with respect to $x_i$ is

$$f'_{i}(x) = (f(x_{1},.,x_{i}+d,.,x_{n}) - f(x_{1},.,x_{i}-d,.,x_{n}))/(2*d)$$

The second order derivative with respect to $x_i$ is

$$f''_{i}(x) = (f(x_{1},.,x_{i}+d,.,x_{n}) - 2 *f(x_{1},.,x_{n}) + f(x_{1},.,x_{i}-d,.,x_{n}))/(d^2)$$

The second order derivative with respect to $x_i, x_j$ is

$$f''_{i,j}(x) = (f(x_{1},.,x_{i}+d,.,x_{j}+d,.,x_{n}) - 2 *f(x_{1},.,x_{n}) +

f(x_{1},.,x_{i}-d,.,x_{j}-d,.,x_{n}))/(2*d^2) - (f''_{i}(x) + f''_{j}(x))/2$$

* http://neos-guide.org/content/nonlinear-least-squares
* http://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/viewer.htm#statug_nlin_sect020.htm


```r
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean, weights = varPower())
summary(fm1)
```



### using nls instead of optim directly.

The `nls` function works normally like the following:

     x <- 1:10
     y <- 2*x + 3                            # perfect fit
     yeps <- y + rnorm(length(y), sd = 0.01) # added noise
     nls(yeps ~ a + b*x, start = list(a = 0.12345, b = 0.54321))#

Because the model I use have a lot of parameters or I don't know beforehand what will be included in the parameter list, I want something like following

    tmp <- function(x,p) { p["a"]+p["b"]*x }
    p0 <- c(a = 0.12345, b = 0.54321)
    nls(yeps ~ tmp(x,p), start = list(p=p0))



```r
tmp <- function(x, coef) {
    a <- coef[1]
    b <- coef[2]
    a + b * x
}

x <- 1:10
y <- 2 * x + 3  # perfect fit
yeps <- y + rnorm(length(y), sd = 0.01)  # added noise
nls(yeps ~ a + b * x, start = list(a = 0.12345, b = 0.54321))  #
```

```
## Nonlinear regression model
##   model:  yeps ~ a + b * x 
##    data:  parent.frame() 
## a b 
## 3 2 
##  residual sum-of-squares: 0.00124
## 
## Number of iterations to convergence: 2 
## Achieved convergence tolerance: 5.06e-09
```

```r
nls(yeps ~ tmp(x, coef), start = list(coef = c(0.12345, 0.54321)))
```

```
## Nonlinear regression model
##   model:  yeps ~ tmp(x, coef) 
##    data:  parent.frame() 
## coef1 coef2 
##     3     2 
##  residual sum-of-squares: 0.00124
## 
## Number of iterations to convergence: 2 
## Achieved convergence tolerance: 5.06e-09
```

```r

tmp1 <- function(coef) {
    a <- coef[1]
    b <- coef[2]
    a + b * x
}
nls(yeps ~ tmp1(coef), start = list(coef = c(0.12345, 0.54321)))
```

```
## Nonlinear regression model
##   model:  yeps ~ tmp1(coef) 
##    data:  parent.frame() 
## coef1 coef2 
##     3     2 
##  residual sum-of-squares: 0.00124
## 
## Number of iterations to convergence: 2 
## Achieved convergence tolerance: 5.06e-09
```

Some R functions to do nonlinear least squares:


```r
## Example from R-help
y <- c(5.5199668, 1.5234525, 3.3557, 6.7211704, 7.4237955, 1.9703127, 4.3939336, 
    -1.4380091, 3.265018, 3.5760906, 0.2947972, 1.0569417)
x <- c(1, 0, 0, 4, 3, 5, 12, 10, 12, 100, 100, 100)
# Define target function as difference
f <- function(b) b[1] * (exp((b[2] - x)/b[3]) * (1/b[3]))/(1 + exp((b[2] - x)/b[3]))^2 - 
    y
x0 <- c(21.16322, 8.83669, 2.957765)
lsqnonlin(f, x0)  # ssq 50.50144 at c(36.133144, 2.572373, 1.079811)

# nls() will break down
nls(y ~ a * (exp((b - x)/c) * (1/c))/(1 + exp((b - x)/c))^2, start = list(a = 21.16322, 
    b = 8.83669, c = 2.957765), algorithm = "plinear")
# Error: step factor 0.000488281 reduced below 'minFactor' of 0.000976563
nls.lm(par = x0, fn = f)

```







```r
load("ex2.rda")
load("ex1.rda")
load("ex1_a.rda")
mkinmodini <- ex2
## Get the parametrization.
  inpartri <- mkinmodini$inpartri
  outpartri <- mkinmodini$outpartri
  ##

  ## mkinmodini is an object by mkinmod.full
  parms.ini <- mkinmodini$parms.ini
  state.ini <- mkinmodini$state.ini
  lower <- mkinmodini$lower
  upper <- mkinmodini$upper
  fixed_parms <- mkinmodini$fixed_parms
  fixed_initials <- mkinmodini$fixed_initials
  mod_vars <- names(mkinmodini$diffs)
  observed <-  mkin_wide_to_long(mkinmodini$residue,time='time')
  observed$err <-c(as.matrix(mkinmodini$weightmat))
  ## Subset dataframe with mapped (modelled) variables
  observed <- subset(observed, name %in% names(mkinmodini$map))
  ## Get names of observed variables
  ## NOTE HERE: the order may not be the same as the input mkinmod.full differential equations list. ## XXXXX TODO XXXX Reorder them maybe a good idea if the data is given from a data file while the mkinmod.full is defined not following the colnames order, although it is already taken care of in the cost(P) function to reorder the odeini using mod_vars
  obs_vars = unique(as.character(observed$name))


  ## Name the parameters if they are not named yet ## usually they are already names
  if(is.null(names(parms.ini))) names(parms.ini) <- mkinmodini$parms

  ## Name the inital parameter values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- mod_vars

  ## Parameters to be optimised
  parms.fixed <- parms.ini[fixed_parms]
  optim_parms <- setdiff(names(parms.ini), fixed_parms)
  parms.optim <- parms.ini[optim_parms]


  ## # ### ### ### ### ###
  state.ini.fixed <- state.ini[fixed_initials]
  optim_initials <- setdiff(names(state.ini), fixed_initials)
  state.ini.optim <- state.ini[optim_initials]
  state.ini.optim.boxnames <- names(state.ini.optim)
  state.ini.fixed.boxnames <- names(state.ini.fixed)
  if(length(state.ini.optim) > 0) {
    names(state.ini.optim) <- paste('M0',names(state.ini.optim),  sep="_")
  }
  if(length(state.ini.fixed) > 0) {
    names(state.ini.fixed) <- paste('M0',names(state.ini.fixed), sep="_")
  }
  
  eigen <- FALSE
  oldparms <- c(state.ini.optim,parms.optim)
  if (length(mkinmodini$map) == 1) {
    solution = "analytical"
  } else {
    if (is.matrix(mkinmodini$coefmat) & eigen) solution = "eigen"
    else solution = "deSolve"
  }
  ## always define mkindiff function since most of the time we will use it.
  mkindiff <- function(t, state, parms) {
    time <- t
    diffs <- vector()
    for (box in mod_vars)
    {
      diffname <- paste("d", box, sep="_")
      diffs[diffname] <- with(as.list(c(time,state, parms)),
                              eval(parse(text=mkinmodini$diffs[[box]])))
    }
    return(list(c(diffs)))
  }
#####################
ctr=kingui.control()
 #### Control parameters ####
  method <- ctr$method
  odesolver <- ctr$odesolver
  atol <- ctr$atol
  rtol <- ctr$rtol
  control <- ctr$control
  marqctr <- ctr$marqctr
  goMarq <- ctr$goMarq
  submethod <- ctr$submethod
  Hmethod1 <- ctr$Hmethod1
  Hmethod2 <- ctr$Hmethod2
####################
environment(kin_mod) <- environment()
oldparms <- c(state.ini.optim,parms.optim)
y <- observed$value
kin_mod(P=oldparms)-y

y <- observed$value
a <- nls(y ~ kin_mod(P,pnames=names(oldparms)),start=list(P=oldparms),lower=lower,upper=upper,algorithm="port")
summary(a)
f <- function(P){
    res <- observed$value-kin_mod(P)
    id <- which(is.na(res))
    return(res[-id])
  }
b <- nls.lm(par=oldparms,lower=lower,upper=upper,fn=f)
b1 <- nls.lm(par=oldparms,lower=lower,upper=upper,fn=f,control=nls.lm.control(maxiter=100))
compare_multi_kinmod(mkinmodini,t(coef(b)))
pnames <- names(oldparms)
obj1 <- function(P,lower,upper,...)
  {
    if(is.null(names(P))) names(P) <- pnames
    if((any(P<lower)) || (any(P>upper)) return(1e100) else return(0.5*sum((kin_mod(P,...)-y)^2,na.rm=TRUE))
  }
## utilities
gradfun_ls <- function(modelfun,par,obs,...)
{
  if(is.null(names(par))) names(par) <- pnames
  np <- length(par)
  resfun <- function(par,obs, ...){
    modelfun(par,...)-obs
  }
  res <- resfun(par,obs=obs)
  Jmat <- jacobian(modelfun,x=par)
  ## when there are NA's, not included in the Jaobian!!!
  id <- which(is.na(res))
  g <- t(Jmat[-id,])%*%res[-id]
  return(g)
}
gr1 <- function(P,...)
  {
    g<-gradfun_ls(modelfun=kin_mod,par=P,obs=observed$value)
    g
  }
c <- trust.optim(x=oldparms,fn=obj1,gr=gr1,method="SR1")
c.1 <- trust.optim(x=oldparms,fn=obj1,gr=gr1,method="BFGS")
## Next using trust.
obj <- function(par,Q=FALSE,...)
  {
    return(objfun_ls(modelfun=kin_mod,par,obs=observed$value,Q=Q,pnames=pnames...))
    
  }
d <- trust_ls(objfun=obj,parinit=oldparms,rinit=1,rmax=100)
###### Case 1: nls.lm better than nls
###### Case 2: trust.optim better than nls.lm
###### Case 3: STIR is the one for all

```



-------------------------


Stackoverflow Question
===========================================================
My original question was "Implementation of trust-region-reflective optimization algorithm in R". However, on the way of producing a reproducible example(thanks @Ben for his advice), I realize that my problem is that in Matlab, one function `lsqnonlin` is good(meaning no need to choose a good starting value, fast) enough for most cases I have, while in R, there is not such a one-for-all function. Different optmization algorithm works well in different cases. Different algorithms reach different solutions. The reason behind this may not be that the optimization algorithms in R is inferior to the trust-region-reflective algorithm in Matlab, it could also be related to how R handles Automatic Differentiation. This problem comes actually from interrupted work two years ago. At that time, Prof. John C Nash, one of the authors of the package **optimx** has suggested that there has been quite a lot of work for Matlab for Automatic Differentiation, which might be the reason that the Matlab lsqnonlin performs better than the optimization functions/algorithms in R. I am not able to figure it out with my knowledge. 

The example below shows some problems I have encountered(More reproducible examples are coming). To run the examples, first run `install_github("KineticEval","zhenglei-gao")`. You need to install package **mkin** and its dependencies and may also need to install a bunch of other packages for different optimization algorithms.

Basically I am trying to solve nonlinear least-squares curve fitting problems as described in the Matlab function `lsqnonlin` 's documentation (http://www.mathworks.de/de/help/optim/ug/lsqnonlin.html). The curves in my case are modeled by a set of differential equations. I will explain a bit more with the examples. Optimization algorithms I have tried including:

* Marq from `nls.lm`, the Levenburg-Marquardt
* Port from `nlm.inb`
* L-BGFS-B from `optim`
* spg from `optimx`
* `solnp` of package **Rsolnp**

I have also tried a few others but here are the ones I will often use. 

## Example 1: A simple case

I will give the R codes first and explain later.


```r
ex1 <- mkinmod.full(Parent = list(type = "SFO", to = "Metab", sink = TRUE, k = list(ini = 0.1, 
    fixed = 0, lower = 0, upper = Inf), M0 = list(ini = 195, fixed = 0, lower = 0, 
    upper = Inf), FF = list(ini = c(0.1), fixed = c(0), lower = c(0), upper = c(1)), 
    time = c(0, 2.8, 6.2, 12, 29.2, 66.8, 99.8, 127.5, 154.4, 229.9, 272.3, 
        288.1, 322.9), residue = c(157.3, 206.3, 181.4, 223, 163.2, 144.7, 85, 
        76.5, 76.4, 51.5, 45.5, 47.3, 42.7), weight = c(1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1)), Metab = list(type = "SFO", k = list(ini = 0.1, fixed = 0, 
    lower = 0, upper = Inf), M0 = list(ini = 0, fixed = 1, lower = 0, upper = Inf), 
    residue = c(0, 0, 0, 1.6, 4, 12.3, 13.5, 12.7, 11.4, 11.6, 10.9, 9.5, 7.6), 
    weight = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
ex1$diffs
Fit <- NULL
alglist <- c("L-BFGS-B", "Marq", "Port", "spg", "solnp")
for (i in 1:5) {
    Fit[[i]] <- mkinfit.full(ex1, plot = TRUE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}
names(Fit) <- alglist
(lapply(Fit, function(x) x$par))
kinplot(Fit[[2]])
unlist(lapply(Fit, function(x) x$ssr))
```

The output from the last line is:

    L-BFGS-B     Marq     Port      spg    solnp 
    5735.744 4714.500 5780.446 5728.361 4714.499 

Except for "Marq" and "solnp", the other algorithms did not reach the optimum.Besides, 'spg' method (also other methods like 'bobyqa') need too many function evaluations for such a simple case . Moreover, if I change the starting value and make `k_Parent=0.0058` (the optimum value for that parameter) instead of the random choosen `0.1`,  "Marq" cannot find the optimum any more! (Code provided below). I have also had datasets where "solnp" does not find the optimum. However, if I use `lsqnonlin` in Matlab, I haven't encountered any difficulties for such simple cases. 


```r
ex1_a <- mkinmod.full(Parent = list(type = "SFO", to = "Metab", sink = TRUE, 
    k = list(ini = 0.0058, fixed = 0, lower = 0, upper = Inf), M0 = list(ini = 195, 
        fixed = 0, lower = 0, upper = Inf), FF = list(ini = c(0.1), fixed = c(0), 
        lower = c(0), upper = c(1)), time = c(0, 2.8, 6.2, 12, 29.2, 66.8, 99.8, 
        127.5, 154.4, 229.9, 272.3, 288.1, 322.9), residue = c(157.3, 206.3, 
        181.4, 223, 163.2, 144.7, 85, 76.5, 76.4, 51.5, 45.5, 47.3, 42.7), weight = c(1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)), Metab = list(type = "SFO", k = list(ini = 0.1, 
    fixed = 0, lower = 0, upper = Inf), M0 = list(ini = 0, fixed = 1, lower = 0, 
    upper = Inf), residue = c(0, 0, 0, 1.6, 4, 12.3, 13.5, 12.7, 11.4, 11.6, 
    10.9, 9.5, 7.6), weight = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))

Fit_a <- NULL
alglist <- c("L-BFGS-B", "Marq", "Port", "spg", "solnp")
for (i in 1:5) {
    Fit_a[[i]] <- mkinfit.full(ex1_a, plot = TRUE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}
names(Fit_a) <- alglist
lapply(Fit_a, function(x) x$par)
unlist(lapply(Fit_a, function(x) x$ssr))
```

Now the output from last line is:

    L-BFGS-B     Marq     Port      spg    solnp 
    5653.132 4866.961 5653.070 5635.372 4714.499 

I will explain what I am optimising here. If you have run the above script and see the curves, we use a two-compartment model with first order reactions to describe the curves. The differential equations to express the model are:


```r
ex1$diffs
```

                                                                 Parent 
                                        "d_Parent = - k_Parent * Parent" 
                                                                   Metab 
    "d_Metab = - k_Metab * Metab + k_Parent * f_Parent_to_Metab * Parent" 

For this simple case, from the differential equations we can derive the equations to describe the two curves. The to be optimized parameters are $M_0,k_p, k_m, c=\mbox{FF_parent_to_Met} $ with the constraints $M_0>0,k_p>0, k_m>0, 1> c >0$.

$$
\begin{split}
            y_{1j}&= M_0e^{-k_pt_i}+\epsilon_{1j}\\
            y_{2j} &= cM_0k_p\frac{e^{-k_mt_i}-e^{-k_pt_i}}{k_p-k_m}+\epsilon_{2j}
            \end{split}
$$

Therefore we can fit the curve without solving differential equations. 


```r
BCS1.l <- mkin_wide_to_long(BCS1)
BCS1.l <- na.omit(BCS1.l)
indi <- c(rep(1, sum(BCS1.l$name == "Parent")), rep(0, sum(BCS1.l$name == "Metab")))
sysequ.indi <- function(t, indi, M0, kp, km, C) {
    y <- indi * M0 * exp(-kp * t) + (1 - indi) * C * M0 * kp/(kp - km) * (exp(-km * 
        t) - exp(-kp * t))
    y
}
M00 <- 100
kp0 <- 0.1
km0 <- 0.01
C0 <- 0.1
library(nlme)
result1 <- gnls(value ~ sysequ.indi(time, indi, M0, kp, km, C), data = BCS1.l, 
    start = list(M0 = M00, kp = kp0, km = km0, C = C0), control = gnlsControl())
# result3 <- gnls(value ~
# sysequ.indi(time,indi,M0,kp,km,C),data=BCS1.l,start=list(M0=M00,kp=kp0,km=km0,C=C0),weights
# = varIdent(form=~1|name)) Coefficients: M0 kp km C 1.946170e+02
# 5.800074e-03 8.404269e-03 2.208788e-01
```


Doing this way, the elapsed time is almost 0, and the optimum is reached. However, we do not always have this simple case. The model can be complex and solving the differential equations is needed. See example 2

## Example 2, a complex model

I worked on this dataset a long time ago and haven't got time to finish running the following script myself. (You might need hours to finish running.) 


```r
data(BCS2)
ex2 <- mkinmod.full(Parent = list(type = "SFO", to = c("Met1", "Met2", "Met4", 
    "Met5"), k = list(ini = 0.1, fixed = 0, lower = 0, upper = Inf), M0 = list(ini = 100, 
    fixed = 0, lower = 0, upper = Inf), FF = list(ini = c(0.1, 0.1, 0.1, 0.1), 
    fixed = c(0, 0, 0, 0), lower = c(0, 0, 0, 0), upper = c(1, 1, 1, 1))), Met1 = list(type = "SFO", 
    to = c("Met3", "Met4")), Met2 = list(type = "SFO", to = c("Met3")), Met3 = list(type = "SFO"), 
    Met4 = list(type = "SFO", to = c("Met5")), Met5 = list(type = "SFO"), data = BCS2)
ex2$diffs
Fit2 <- NULL
alglist <- c("L-BFGS-B", "Marq", "Port", "spg", "solnp")
for (i in 1:5) {
    Fit2[[i]] <- mkinfit.full(ex2, plot = TRUE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}
kinplot(Fit2[[5]])
names(Fit) <- alglist
(lapply(Fit, function(x) x$par))
unlist(lapply(Fit, function(x) x$ssr))

```


This is an example where you will see warning messages like:

    DLSODA-  At T (=R1) and step size H (=R2), the    
      corrector convergence failed repeatedly     
      or with ABS(H) = HMIN   
    In above message, R = 
    [1] 0.000000e+00 2.289412e-09

## Issues in R

* Different algorithms reach different solutions. 
* Missing global minimum
* Lower Efficiency

I have to adimit that this problem comes actually from interrupted work two years ago. At that time, Prof. John C Nash, one of the authors of the package `optimx` has suggested that there has been quite a lot of work for Matlab for Automatic Differentiation, which might be the reason that the Matlab `lsqnonlin` performs better than the optimization functions/algorithms in R. 

Thanks @Ben Bolke at Stackoverflow for his advice of giving a **concrete, reproducible** example of a problem where R's optimizers do significantly worse than `lsqnonlin` and reframe my original question of "can I implement algorithm xxyy in R?". This makes me sit down and digging out emails and examples and have a review of my previous work. 


```r
summary(cars)
```

```
##      speed           dist    
##  Min.   : 4.0   Min.   :  2  
##  1st Qu.:12.0   1st Qu.: 26  
##  Median :15.0   Median : 36  
##  Mean   :15.4   Mean   : 43  
##  3rd Qu.:19.0   3rd Qu.: 56  
##  Max.   :25.0   Max.   :120
```


You can also embed plots, for example:


```r
###
```


## Sources

* http://stackoverflow.com/questions/13069627/library-for-trust-region-reflective-algorithms-in-c
* Ipopt
* http://www.mathworks.de/de/help/optim/ug/constrained-nonlinear-optimization-algorithms.html
* http://stackoverflow.com/questions/5527145/convert-matlab-code-to-r
