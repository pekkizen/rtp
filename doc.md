Why sample mean of a product of n Uniform(0,1)
distributed independant random
numbers is not either $1/2^{n}$ or $1/e^{n}$.
------------------------------------



I got to this phenomenon when coding
[Fisher's method](https://en.wikipedia.org/wiki/Fisher%27s_method)
of combining p-values. The test statistics is logarithm  of product of the p-values = sum of logarithms of individual p-values.
Logarithm of product is much faster to calculate than sum of logarithms. In testing I noticed,
that product of more than 700 p-values start to go to zero with mean ~745 p-values.
Smallest non zero value in 64-bit arithmetic is 2^-1074.
From this it is simple to deduce, that going to zero would need a product of around 1074 p-values. Actually it is probabilitistic impossibility to get so far. For going to zero with 745  p-values,
we can solve rate parameter x from  $x^{-745} = 2^{-1074}$.

$$log(x) = -1074 * log(2) / -745 = 0.999 \cong 1 \Rightarrow x \cong e$$  

When two random variables are statistically independent, the expectation of their product is the product of their expectations, E.g.
[Wikipedia](https://en.wikipedia.org/wiki/Product_distribution).
This is also used as definition of independance. If $E[x\cdot y] = E[x]\cdot E[y]$,
then x and y are independant.

$$x_i \sim U(0, 1)$$
$$E[x_i] = 1/2$$  
$$x = \prod_{i=1}^{n}x_i$$  
From the independance assumption we get the first candidate for E[x]
$$E[x] = \prod_{i=1}^{n} E[x_i] = 1/2^{n}$$

[WolframAlpha](https://www.wolframalpha.com/input/?i=uniform+product+distribution) presents
Uniform Product Distribution with density function.
StackExchange has post around the subject.
[This is one](https://math.stackexchange.com/questions/659254/product-distribution-of-two-uniform-distribution-what-about-3-or-more).

$$f_u(x;n) = \left | log(x) \right |^{n-1} / (n-1)!$$
which gives also expected value $1/2^{n}$. I could not solve the indefinite integral and the result $1/2^{n}$ is based on numerical experiments on WolframAlpha site.
Also WolframAlpha refused to solve the integral because "Standard computation time exceeded".

$$E[x] = \int_{0}^{1} x\cdot f_u(x;n)  dx = 1/2^{n}$$

Using logaritms of x and $x_i$ leeds Gamma and Exponential distribution and $E[x] = 1/e^{n}$.
$$y = log(x) = \sum_{i=1}^{n} log(x_i)$$
$$log(x_i) \sim -Exp(1) \Rightarrow E[y] = \sum_{i=1}^{n} E[log(x_i)] = -n$$
$$y \sim -Gamma(n, 1) \Rightarrow E[y] = -n$$
$$E[x] = e^{E[y]} = 1/e^{n}$$

Expected number of variables $x_i$ needed in the product to reach value $e^{-\lambda}$ is  distributed $Poisson(\lambda)$.  
 
$$n \sim Poisson(\lambda)$$

$$f_p(n; \lambda) = \lambda^{n-1}\cdot e^{-\lambda} / (n-1)!$$

The Poisson density function can be expressed by Gamma
density, where the normal parameters have changed place and role. The dgamma density is regarded as a function of the shape parameter n.
$$R: dgamma(shape=n, x=\lambda) = dpois(n-1, \lambda)$$

Mean number of $x_i$ to reach any limit L is log(L), not log2(L). If the mean number of $x_i$ to reach a limit is k, it does't imply that the limit is the mean value of all products of lenth k.


For fixed number n we can integrate the mean by integrating
over all possible limit values $e^{-\lambda}$ over $\lambda$ range $1 - \infty$.
$$E[x] = \int_{1}^{\infty}e^{-\lambda}\cdot f_p(n,\lambda)  d\lambda = 1/2^{-n}\cdot \Gamma (n, 2) / (n-1)!$$

$$E[x] = 1/2^{n}\cdot GammaCDF(n, 2) \cong   1/2^{n}$$
Here we have kind of controversial results: general mean of a producs of length n is $1/2^{n}$ and product of any
lenth k reaching limit $L = e^{-\lambda}$ have mean $1/e^{\lambda}$.

[D. Knuth](https://en.wikipedia.org/wiki/Poisson_distribution#Random_drawing_from_the_Poisson_distribution) has even given a simple algorithm to generate  Poisson-distributed random numbers.  

```R
# Knuth's algorithm in R
rpoissonKnuth <- function(lambda) {
    L <- exp(-lambda)
    k <- 0
    p <- 1
    while (TRUE) {
        k <- k + 1
        p <- p * runif(1)
        if (p < L) {
            return(k - 1)
        }
    }
}
```