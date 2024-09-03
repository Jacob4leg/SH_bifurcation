# Bifurcation and the curve of periodic orbits of Swift-Hohenberg equation
This is suplement code for the article "Continuation and bifurcations of periodic orbits and symbolic dynamics in the Swift-Hohenberg equation". The program uses CAPD library to prove the Theorem 2, i.e. the existence of curves $x_{\pm}  :  \left [0,\xi_* \right ]   \rightarrow   \mathbb{R}, \tilde{\xi} : X_* \rightarrow \mathbb{R}$, the existence of bifurcation point $\left (\xi^* , x^* \right )$, as well as the curves $x_{\pm}, \tilde{\xi}$ glue into one smooth curve. The program is divided into curve_x, curve_xi modules, which are responsible for determining the curves $x_{\pm}, \tilde{\xi}$ respectively. To compile the programs call

```make CAPDBINDIR=<path_to_CAPD_build_bin_directory>```

where the default path is

```CAPDBINDIR=~/CAPD/build/bin/```

The output of the program looks as follows:

**If the curves $x_{\pm}$ are beeing determined**

```xi=[], X=[]```

where in i-th iteration the program checks the existence of a unique smooth curve on the set xi $\times$ X. Program outputs every thousandth successful iteration.

**If the curve $\tilde{\xi}$ is beeing determined**

```XI=[], x=[], xi''(x)=[]```

where x is treated as the parameter. xi''(x) is the bound of the implicit derivative of the second order on the set XI $\times$ x. Here program outputs every hundredth successful iteration. In the halfway of determining $\tilde{\xi}$ curve, there appears the result of determining the existence of the maximum $\tilde{\xi}( x^* )$, which looks as follows

```maximum in XI=[], x=[]```

At the beginning and at the end of determining curve $\tilde{\xi}$ the program checks if the curve glues correctly with curves $x_{\pm}$ into one smooth curve. If that is the case, then the program outputs message

```threshold box in parameterization```
