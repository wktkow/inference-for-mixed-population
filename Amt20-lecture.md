## Contents


## II Non-linear Models II.

- 0 References
- I Inference for Mixed Populations I.
- I.1 Introduction to Finite Mixtures I.
- I.2 Mean and Variance of Finite Mixtures .I.
- I.3 Inference For Finite Mixtures With Fixed Support Size .I.
- I.4 Fitting Finite Mixtures in R and SAS .I.
- I.5 Inference for Number of Support Points .I.
- I.6 Non-parametric Maximum Likelihood .I.
- I.7 Numerical Algorithms .I.
- I.8 Examples in CAMAN .I.
- I.9 Classification .I.
- I.10 Model Extensions .I.
- II Non-linear Models II.
- II.1 Non-Linear Mixed Models .II.
- II.2 Pharmacokinetic and Pharmacodynamic Models .II.


Chapter 0

References

- B ̈ohning D. (1999)Computer-assisted Analysis of Mixtures and Applications. Meta-analysis,
    Disease Mapping and Others. London: Chapman & Hall.
- Cressie, N.A.C. (1991)Statistics for Spatial Data. New York: John Wiley.
- Davidian, M. and Giltinan, D.M. (1995)Nonlinear Models for Repeated Measurement Data.
    London: Chapman & Hall.
- Diggle, P.J. (1983)Statistical Analysis of Spatial Point Patterns. Mathematics in Biology. London:
    Academic Press.
- Fahrmeir, L. and Tutz, G. (2001)Multivariate Statistical Modelling Based on Generalized Linear
    Models. Heidelberg: Springer-Verlag.

Advanced Modeling Techniques 1


- Fitzmaurice, G.M., Davidian, M., Verbeke, G., and Molenberghs, G.(2009). Longitudinal Data
    Analysis. Handbook. Hoboken, NJ: John Wiley & Sons.
- Lindsay B.G. (1995)Mixture Models: Theory, Geometry and Applications. NSF-CBMS regional
    conference series in probability and statistics, vol.5, Hayward: Institute of Mathematical Statistics.
- McLachlan G.J. and Basford K.E. (1988) Mixture Models. Inference and Applications to
    Clustering. New York: Marcel Dekker.
- McLachlan, G.J. and Peel, D. (2000)Finite Mixture Models. New York: John Wiley.
- Molenberghs, G. and Verbeke, G. (2005). Models for Discrete Longitudinal Data. New York:
    Springer-Verlag.
- Ripley, B.D. (1981)Spatial Statistics. New York: John Wiley.
- Titterington D.M., Smith A.F.M., & Makov U.E. (1985)Statistical Analysis of Finite Mixture
    Distributions. New York: John Wiley.

Advanced Modeling Techniques 2


- Verbeke, G. and Molenberghs, G. (2000)Linear mixed models for longitudinal data. New York:
    Springer-Verlag.

Advanced Modeling Techniques 3


Part I

Inference for Mixed Populations


Chapter I.

Introduction to Finite Mixtures

```
⊲ Snapper data
```
```
⊲ Unobserved heterogeneity
```
```
⊲ The cocktail example
```
```
⊲ Mixture of normals
```

I.1.1 Snapper Data

- Data set snapper.dat
- Length measurements (inches) of 256 snappers, with histogram:


- Histogram shows multi-modality which cannot easily be described by standard
    distributions.
- Biological interpretation:

```
⊲ Underlying categories correspond to age classes
⊲ Within each age class a rather ‘homogeneous’ distribution seems plausible
⊲ The relative heights of the modes give an indication of the proportion of the
population in that particular age class.
```

I.1.2 Unobserved Heterogeneity

- The multi-modality observed in the histogram suggests the presence of some
    underlying (latent) group structure
- In many cases, as in the snapper data, the group structure is not known or has not
    been recorded
- Let us assume that the population P of interest is composed of g sub-populations
    P 1 ,P 2 ,... ,Pg :

```
P = {P 1 ,P 2 ,... ,Pg }
```

- Each population Pj represents a proportion πj of the total population, ∑gj=1πj = 1
- Let X indicate from which population an observation has been sampled:

```
X =j ⇐⇒ Observation belongs to Pj
```
- The distribution of X is discrete with support { 1 , 2 ,... , g} and corresponding
    probabilities {π 1 , π 2 ,... , πg}:

### X ∼

```




```
```
1 2... g
π 1 π 2... πg
```
```




```
- X is latent, as it is not observed


- Let the density of the outcome Y in sub-population Pj be fj(y)
- The density of Y in the entire population P then equals:

```
f(y) = ∑
j
```
```
f(y|X = j)P(X = j) = ∑
j
```
```
πjfj(y)
```
- The distribution of Y is called a (finite) mixture with g components
- The densities f 1 (y),... , fg(y) often depend on (vectors of) (un-)known
    parameters θ 1 ,... , θg.
- The densities f 1 (y),... , fg(y) can be continuous, discrete, or a mixture of both
    types.


I.1.3 The Cocktail Example

- A mixture can be compared to a cocktail which is a stirred mixture of a number of
    ingredients, each representing a percentage of the cocktail:

```
Cocktail ←→ P ←→ Y ∼ f(y)
Ingredient ←→ Pj ←→ Y|X =j ∼ fj(y)
```
```
Relative proportion ←→ πj ←→ P(X =j)
```
- Research questions:
    ⊲ How many ingredients?
    ⊲ Which ingredients?
    ⊲ Relative proportions?


I.1.4 Snapper Data Revisited

- The four modes suggest a 4-component mixture
- A 4-component mixture of normals with equal variance has been fitted:

```
Y|X = j ∼ N(μj, σ^2 ) X ∼
```
```




```
### 1 2 3 4

```
π 1 π 2 π 3 π 4
```
```




```
- Equivalently, this can be written as:

```
Y|μ ∼ N(μ, σ^2 ) μ ∼
```
```




```
```
μ 1 μ 2 μ 3 μ 4
```
```
π 1 π 2 π 3 π 4
```
```




```

- Fitted model:

```
Y|μ ∼ N(μ, 0. 672 ) μ ∼
```
```




```
### 3 .43 5.32 7.60 10. 33

### 0 .12 0.53 0.27 0. 08

```




```

I.1.5 Mixture of Two Normals With Equal Variance

- Graphical representation of the mixture: π N(μ 1 , σ^2 ) + (1−π) N(μ 2 , σ^2 )


- Very flexible class of models:

```
⊲ Symmetric as well as skewed
⊲ Unimodal as well as multimodal
```
- If |μ 1 −μ 2 |/σ ≤ 2 then the mixture is unimodal for all π
- If |μ 1 −μ 2 |/σ > 2 then the modality of the mixture depends on π
- In general, the modes of the mixture are closer to each other than the modes of
    the components.
- Hence the number of components may not be graphically visible as was the case
    with the snapper data
- Also, the ‘appropriate’ g very much depends on the component densities fj(y)


Chapter I.

Mean and Variance of Finite Mixtures

```
⊲ General principle
```
```
⊲ Examples
```

I.2.1 General Principle

- Moments can easily be obtained using the latent variable representation:

### E(Y) = E[E(Y|X)]

```
Var(Y) = Var[E(Y|X)] + E[Var(Y|X)]
```
- Conditional moments E(Y|X) and Var(Y|X) directly follow from the component
    densities fj


I.2.2 Normals With Common Variance

```
Y|μ ∼ N(μ, σ^2 ) μ ∼
```
```




```
```
μ 1 μ 2 ··· μg
```
```
π 1 π 2 ··· πg
```
```




```
```
E(Y) = E(μ) = ∑
j
```
```
πjμj
```
```
Var(Y) = Var(μ) +E(σ^2 ) = Var(μ) +σ^2
```
### = ∑

```
j
```
```
πjμ^2 j −
```
```

∑
j
```
```
πjμj
```
```

^2 +σ 2
```

I.2.3 Normals With Common Mean

```
Y|σ^2 ∼ N(μ, σ^2 ) σ^2 ∼
```
```




```
```
σ 12 σ^22 ··· σg^2
```
```
π 1 π 2 ··· πg
```
```




```
```
E(Y) = E(μ) = μ
```
```
Var(Y) = Var(μ) +E(σ^2 ) = E(σ^2 )
```
### = ∑

```
j
```
```
πjσj^2
```

I.2.4 Normals With General Mean and Variance

```
Y|(μ, σ^2 ) ∼ N(μ, σ^2 ) (μ, σ^2 ) ∼
```
```




```
```
(μ 1 , σ^21 ) (μ 2 , σ^22 ) ··· (μg, σg^2 )
```
```
π 1 π 2 ··· πg
```
```




```
```
E(Y) = E(μ) = ∑
j
```
```
πjμj
```
```
Var(Y) = Var(μ) +E(σ^2 ) = Var(μ) +
```
```
∑
j
```
```
πjσ^2 j
```
### = ∑

```
j
```
```
πjμ^2 j −
```
```

∑
j
```
```
πjμj
```
```

^2 +∑
j
```
```
πjσj^2
```

I.2.5 Binomials

```
Y|p ∼ Bin(n, p) p ∼
```
```




```
```
p 1 p 2 ··· pg
π 1 π 2 ··· πg
```
```




```
```
E(Y) = E(np) = n∑
j
```
```
πjpj
```
```
Var(Y) = Var(np) +E[np(1−p)] = n^2 Var(p) +nE(p)−nE(p^2 )
```
```
= n(n−1)E(p^2 )−n^2 [E(p)]^2 +nE(p)
```
```
= n(n−1)∑
j
```
```
πjp^2 j −n^2
```
```

∑
j
```
```
πjpj
```
```

^2 +n∑
j
```
```
πjpj
```

I.2.6 Poissons

```
Y|λ ∼ Poisson(λ) λ ∼
```
```




```
```
λ 1 λ 2 ··· λg
```
```
π 1 π 2 ··· πg
```
```




```
```
E(Y) = E(λ) = ∑
j
```
```
πjλj
```
```
Var(Y) = Var(λ) +E(λ)
```
### = ∑

```
j
```
```
πjλ^2 j −
```
```

∑
j
```
```
πjλj
```
```

^2 +∑
j
```
```
πjλj
```

Chapter I.3

Inference For Finite Mixtures With Fixed Support Size

```
⊲ Introduction
```
```
⊲ EM algorithm
```
```
⊲ Example: Mixture of normals with general mean and variance
```
```
⊲ Properties and remarks
```

I.3.1 Introduction

- In this chapter, we will study how the parameters in a finite mixture distribution
    can be estimated using maximum likelihood estimation (MLE)
- We will consider the number g of components to be fixed (known)
- Let Y 1 ,... , YN be distributed as:

```
Yi ∼ π 1 fi 1 (yi) + π 2 fi 2 (yi) +... + πgfig(yi)
```
### =

```
∑g
j=1
```
```
πjfij(yi)
```
- fi 1 (yi),... , fig(yi) are the density functions of Yi in the g components of the
    mixture


- Often, we have that
    fij(yi) = fj(yi), for all j
assuming that all Yi follow the same distribution
- The index i now allows the Yi to have different distributions, e.g., to include
    covariates (see later)
- As before, we allow the densities fij(yi) to depend on unknown parameters, which
    are combined in the vector θ.
- This will often be explicitly denoted as fij(yi|θ)
- Further, let π be the vector of component probabilities: π′ = (π 1 ,... , πg)
- The vector ψ is the vector containing all unknown parameters in the model:
    ψ′ = (π′,θ′)


- The likelihood function equals:

```
L(ψ|y) =
```
```
N∏
i=1
```
```


```
```
∑g
j=1
```
```
πj fij(yi|θ)
```
```


```
```
where y′ = (y 1 ,... , yN) is the vector containing all observed response values.
```
- The corresponding log-likelihood equals:

```
ℓ(ψ|y) =
```
```
N∑
i=1
```
```
ln
```
```


```
```
∑g
j=1
```
```
πj fij(yi|θ)
```
```


```
- Maximizing ℓ(ψ|y) with respect to ψ in general requires numerical iterative
    procedures


- Also, the analytic expression of ℓ(ψ|y) suggests that numerical maximization will
    be far from straightforward
- For example, classical Newton-Raphson procedures would require calculation of
    first- and second-order derivatives of ℓ(ψ|y).
- An alternative procedure, especially convenient for mixture models, is the EM
    algorithm, the Expectation–Maximization algorithm
- EM is designed for MLE in situations with missing data
- Here, the underlying latent variable X, i.e., the component membership, will be
    considered missing.


I.3.2 EM Algorithm

3.2.1 Observed and Complete Data Likelihoods

- We define indicators Zij, i = 1,... , N, j = 1,... , g:

```
Zij =
```
```







```
```
1 if observation i belongs to component j
```
```
0 otherwise
```
- We then have that

```
P(Zij = 1) = πj
```

- The joint density of Yi and all associated Zij equals

```
fi(yi, Zi 1 = zi 1 ,... , Zig = zig)
```
```
= fi(yi | Zi 1 =zi 1 ,... , Zig = zig) × P(Zi 1 =zi 1 ,... , Zig = zig)
```
### =

```


```
```
∏g
j=1
```
```
[fij(yi|θ)]zij
```
```

×
```
```


```
```
∏g
j=1
```
```
πzjij
```
```

 =
```
```
∏g
j=1
```
```
[πjfij(yi|θ)]zij
```
- The joint likelihood function for the observed measurements y and for the vector
    z of all unobserved zij therefore equals:

```
L(ψ|y,z) =
```
```
N∏
i=1
```
```
∏g
j=1
```
```
[πjfij(yi|θ)]zij
```

- The corresponding log-likelihood function equals:

```
ℓ(ψ|y,z) =
∑N
i=1
```
```
∑g
j=1
```
```
zij{lnπj + lnfij(yi|θ)}
```
- Terminology:

```
L(ψ|y,z) : Complete data likelihood
ℓ(ψ|y,z) : Complete data log-likelihood
L(ψ|y) : Observed data likelihood
ℓ(ψ|y) : Observed data log-likelihood
```
- Note that maximizing ℓ(ψ|y,z) is much easier than maximizing the log-likelihood
    ℓ(ψ|y) of the observed data only.


- However, the obtained estimates would depend on the unobserved indicators zij.
- Compromise: Maximize the expected value of ℓ(ψ|y,Z), i.e., maximize

```
E[ℓ(ψ|y,Z) | y]
```
- An intuitive explanation is that the ‘missing’ observations zij are replaced by their
    expected values.


3.2.2 EM Algorithm

- The EM algorithm acts iteratively, in the sense that, starting from a ‘first guess
    estimate’ (starting value) ψ(1) for ψ, a series of estimates ψ(t) is constructed,
    which converges to the MLE ψ̂ of ψ:

```
ψ(1) →ψ(2) →.. .→ψ(t)→ψ(t+1)→.. .→ψ(∞)=ψ̂
```
- Given ψ(t), the updated estimate ψ(t+1) is obtained through one E step
    and one M step.
- E step: Calculation of

```
Q(ψ|ψ(t)) = E
```
```
[
ℓ(ψ|y,Z) | y,ψ(t)
```
```
]
```

- M step: Maximize Q(ψ|ψ(t)) with respect to ψ to obtain the updated estimate
    ψ(t+1).
- The procedure keeps iterating between the E step and the M step until
    convergence is attained, i.e., until

```
|ℓ(ψ(t+1)|y) − ℓ(ψ(t)|y)| < ε,
```
```
for some small, pre-specified, ε > 0.
```

3.2.3 The E Step

- Q(ψ|ψ(t)) is obtained from:
    Q(ψ|ψ(t)) = E

```
[
ℓ(ψ|y,Z) | y,ψ(t)
```
```
]
```
### = E

```


```
```


```
```
∑N
i=1
```
```
∑g
j=1
```
```
Zij[lnπj + lnfij(yi|θ)]
```
```


```
```
∣∣∣
∣∣∣
∣ y,ψ
```
```
(t)
```
```


```
### =

```
∑N
i=1
```
```
∑g
j=1
```
### E

```
[
Zij | y,ψ(t)
```
```
]
[lnπj + lnfij(yi|θ)]
```
- Hence, the E step only requires calculation of

### E

```
[
Zij | yi,ψ(t)
```
```
]
= P
```
```
(
Zij = 1 | yi,ψ(t)
```
```
)
=
```
```
fi(yi | Zij = 1) P(Zij = 1)
fi(yi|θ)
```
```
∣∣∣
∣∣∣
∣∣
ψ(t)
```
### =

```
πj fij(yi|θ)
∑
j
```
```
πjfij(yi|θ)
```
```
∣∣∣
∣∣∣
∣∣∣
∣∣∣
ψ(t)
```
```
= πij(ψ(t))
```

- πij(ψ(t)) is the posterior probability for observation i to belong to the jth
    component of the mixture
- From now on, πj will be called the prior probability for observation i to belong to
    the jth component of the mixture
- The E step reduces to calculating all posterior probabilities πij(ψ(t)),
    i = 1,... , N, j = 1,... , g.


3.2.4 The M Step

- The updated estimate ψ(t+1) is obtained from maximizing

```
Q(ψ|ψ(t)) =
∑N
i=1
```
```
∑g
j=1
```
```
πij(ψ(t)) [lnπj + lnfij(yi|θ)]
```
```
with respect to ψ′ = (π′,θ′).
```
- We first maximize with respect to π:
    ⊲ This requires maximization of

```
∑N
i=1
```
```
∑g
j=1
```
```
πij(ψ(t)) lnπj =
∑N
i=1
```
```
g∑− 1
j=1
```
```
πij(ψ(t)) lnπj +
∑N
i=1
```
```
πig(ψ(t)) ln
```
```

 1 −g∑−^1
j=1
```
```
πj
```
```


```
```
with respect to π 1 ,... , πg− 1
```

⊲ We set all first-order derivatives equal to zero:

### ∂

```
∂πj
```
### = 0 ⇔

```
∑N
i=1
```
```
πij(ψ(t))
πj(t+1)
```
### =

```
∑N
i=1
```
```
πig(ψ(t))
πg(t+1)
```
### ⇔

```
π(jt+1)
π(gt+1)
```
### =

```
∑N
i=1
```
```
πij(ψ(t))
∑N
i=1
```
```
πig(ψ(t))
```
⊲ This implies that

### 1 =

```
∑g
j=1
```
```
πj(t+1) =
```
```
∑g
j=1
```
```
π(gt+1)
∑N
i=1
```
```
πij(ψ(t))
∑N
i=1
```
```
πig(ψ(t))
```
### =

```
πg(t+1)
∑N
i=1
```
```
︷ ︸︸^1 ︷
∑g
j=1
```
```
πij(ψ(t))
∑N
i=1
```
```
πig(ψ(t))
```
### =

```
N πg(t+1)
∑N
i=1
```
```
πig(ψ(t))
```

```
⊲ Hence, πg(t+1) is given by
```
```
π(gt+1) =
```
```
∑N
i=1
```
```
πig(ψ(t))
N
⊲ It now also follows that all π(jt+1) are given by
```
```
πj(t+1) =
```
```
∑N
i=1
```
```
πij(ψ(t))
N
```
```
⊲ The updated mixture component probabilities are the average posterior
probabilities.
```
- Maximization with respect to θ requires maximization of

```
∑N
i=1
```
```
∑g
j=1
```
```
πij(ψ(t)) lnfij(yi|θ)
```

- In simple examples, this can be done analytically
- In general, however, this cannot be done analytically, and a classical maximization
    procedure, such as Newton-Raphson, is used.
- In such cases, the EM algorithm is double iterative, which can have serious
    consequences on the computation times.


I.3.3 Example: Normals With General Mean and Variance

```
Yi ∼
∑g
j=1
```
```
πjN(μj, σj^2 )
```
```
θ = (μ 1 ,... , μg, σ 12 ,... , σg^2 )
```
- The log-likelihood corresponding to the above model is

```
ℓ(ψ|y) =
∑N
i=1
```
```
ln
```
```



```
```
∑g
j=1
```
```
πj
```
### √^1

```
2 πσ^2 j
```
```
exp
```
```


−
```
### 1

```
2 σj^2
```
```
(yi−μj)^2
```
```



```
```



```

- This can be re-written as

```
ℓ(ψ|y) =
∑N
i=2
```
```
ln
```
```



```
```
∑g
j=2
```
```
πj √^1
2 πσj^2
```
```
exp
```
```

−^1
2 σ^2 j (yi−μj)
```
```
2
```
```

 + π 1 √^1
2 πσ 12 exp
```
```

−^1
2 σ 12 (yi−μ^1 )
```
```
2
```
```


```
```



```
```
+ ln
```
```



```
```
∑g
j=2
```
```
πj √^1
2 πσj^2
```
```
exp
```
```

−^1
2 σj^2 (y^1 −μj)
```
```
2
```
```

 + π 1 √^1
2 πσ 12 exp
```
```

−^1
2 σ 12 (y^1 −μ^1 )
```
```
2
```
```


```
```



```
- Taking μ 1 equal to y 1 , this becomes

```
ℓ(ψ|y) =
∑N
i=2
```
```
ln
```
```



```
```
∑g
j=2
```
```
πj √^1
2 πσj^2
```
```
exp
```
```

−^1
2 σj^2 (yi−μj)
```
```
2
```
```

 + π 1 √^1
2 πσ 12 exp
```
```

−^1
2 σ 12 (yi−y^1 )
```
```
2
```
```


```
```



```
```
+ ln
```
```



```
```
∑g
j=2
```
```
πj √^1
2 πσj^2
```
```
exp
```
```

−^1
2 σj^2 (y^1 −μj)
```
```
2
```
```

 + π 1 √^1
2 πσ 12
```
```



```

- The above expression converges to +∞ when σ 12 approaches zero
- Hence, this mixture model leads to infinite likelihoods.
- This can only be solved by keeping the component variances away from zero
- One way to do so is by assuming all the variances to be equal, i.e., σ^2 j =σ^2
- Indeed, σ^2 = 0 would then result in a discrete marginal mixture distribution with g
    components, which is not possible as soon as the number of distinct data points is
    larger than g.


- Sometimes a well fitting mixture of normals with general mean and variance can
    be obtained after convergence to a local maximum of the likelihood.
- However, since the solution is then NOT MLE, inference does not follow from
    standard likelihood theory.
- We will therefore only consider mixtures of normal distributions with common
    variance.


I.3.4 Properties and Remarks

3.4.1 Identifiability

- Consider the following mixture of 3 Poissons:

```
Y ∼ π 1 Poisson(λ 1 ) +π 2 Poisson(λ 2 ) +π 3 Poisson(λ 3 )
```
- The parameter vector ψ then equals: ψ′ = (π 1 , π 2 , π 3 , λ 1 , λ 2 , λ 3 )
- Note that likelihood value for
    ψ′ = (0. 1 , 0. 7 , 0. 2 , 1 , 5 ,7)
is exactly the same as the likelihood value for
ψ′ = (0. 1 , 0. 2 , 0. 7 , 1 , 7 ,5)


- In fact, any permutation of the elements in

```
{(λ 1 , π 1 ),(λ 2 , π 2 ),(λ 3 , π 3 )}
leads to the same likelihood value.
```
- In general, for g components, there are g! possible permutations of the mixture
    components, all yielding the same likelihood value, i.e., the likelihood has at least
    g! local maxima with the same likelihood value.
- This shows that, for finite mixtures of distributions of the same parametric family
    (e.g., mixture of Normals, Binomials, Poissons,... ), the vector ψ is not uniquely
    identified.
- One way to make ψ identifiable is by ordering the mixture components according
    to the corresponding component probabilities, e.g.,
       π 1 ≥ π 2 ≥... ≥ πg


3.4.2 Monotonicity Property of EM Algorithm

- It can be shown that an EM step cannot decrease the likelihood valueℓ(ψ|y), i.e.,

```
ℓ(ψ(t+1)|y) ≥ ℓ(ψ(t)|y), for all t
```
- This is called the monotonicity property of the EM algorithm
- It guarantees convergence of the iterative procedure
- Note that this does not guarantee convergence to a global maximum


3.4.3 Existence of Local Maxima

- Apart from the local maxima resulting from the non-identifiability problem, there
    may be local maxima yielding different likelihood values
- Example from B ̈ohning (p.66). Mixture of two normals with common variance:
- Obviously, the second and third set of estimates correspond to local maxima, as
    the first set of estimates yields a higher log-likelihood value
- This suggests that multiple sets of starting values should be used in practice


3.4.4 Convergence to a Ridge

- Consider fitting the mixture
    Y ∼ π N(μ 1 , σ^2 ) + (1−π) N(μ 2 , σ^2 )
of 2 normals with common variance, while the true distribution of Y is a single
normal, i.e.,
Y ∼ N(μ, σ^2 )
- We then have that the likelihood is maximized on a ridge of parameter values:
    μ 1 = μ 2 or π = 0 or π = 1
- The EM algorithm is capable of converging to some particular point on that ridge.
- This is not the case for many other, more classical, maximization algorithms.
- This is why the EM algorithm is especially convenient for mixture models


3.4.5 Convergence Rate

- Although the monotonicity property guarantees convergence, this convergence can
    be painfully slow
- Example from B ̈ohning (p.63). Mixture of three known distributions:
- With badly selected starting values, such slow convergence can lead to long
    computation times, especially when the M step in the algorithm requires iterative
    maximization (i.e., when the EM is double iterative).


Chapter I.4

Fitting Finite Mixtures in R and SAS

```
⊲ R-package CAMAN
```
```
⊲ SAS procedure FMM
```
```
⊲ Comparison of R with SAS
```
```
⊲ Example: Child data
```
```
⊲ Example: SIDS data
```

I.4.1 R-package CAMAN

- Developed by Peter Schlattmann, Johannes Hoehne, and Maryna Verba, based on
    C.A.MAN (D. B ̈ohning & P. Schlattmann)
- Contains several functions for mixture analyses
- Function ‘mixalg.EM’ is based on EM algorithm for fitting finite mixtures with
    fixed number of components
- Let’s fit a 4 component normal mixture to the snapper data
- Loading the data:

```
> load("c:/analysis/mixtureR/snapper.rdata")
```

- Data structure:

```
> snapper
length frequency
1 2.875 6
2 3.125 7
3 3.375 9
...................
40 12.625 1
```
- Fitting normal mixture with 4 components:

```
> em<-mixalg.EM(obs="length", weights="frequency", family="gaussian",
data=snapper, p=c(0.10,0.50,0.30,0.10), t=c(3,5,8,10))
```
- ‘obs=’ specifies the outcome
- ‘weights=’ specifies a replication factor (weight)


- ‘family=’ specifies the distribution in each component
- ‘data=’ specifies the data set
- ‘p=’ specifies starting values for the component probabilities,
    and indirectly the number of components
- The procedure automatically rescales the starting values in ‘p=’
- Hence the following specifications are equivalent:

```
p=c(0.10,0.50,0.30,0.10)
p=c(10,50,30,10)
```
- ‘t=’ specifies starting values for the component locations,
    and indirectly the number of components


- Generated output:

```
> em
Computer Assisted Mixture Analysis:
Data consists of 256 observations (rows).
The Mixture Analysis identified 4 components of a gaussian distribution:
DETAILS:
p mean
1 0.11755367 3.432325
2 0.53355806 5.319268
3 0.27207539 7.601072
4 0.07681288 10.334596
component variance: 0.447414325225651
Log-Likelihood: -505.7188 BIC: 1050.254
```
- Fitted model:

```
Y|μ ∼ N(μ, 0. 672 ) μ ∼
```
```




```
### 3 .43 5.32 7.60 10. 33

### 0 .12 0.53 0.27 0. 08

```




```

I.4.2 SAS Procedure FMM

- The same 4-component normal mixture can be fitted to the snapper data, using
    the following syntax:

```
proc fmm data=snapper;
model length = / dist=gaussian equate=scale k=4
parms(1 0.5,5 0.5, 9 0.5, 13 0.5);
probmodel / parms(0,1.6,1.1);
freq frequency;
run;
```
- ‘dist=’ specifies the distribution in each component
- ‘k=’ specifies the number of components in the mixture
- ‘parms(... )’ specifies starting values for all parameters in the component densities


- Here, this implies specification of the mean and variance in each of the normal
    components of the mixture
- ‘equate=’ specifies parameter constraints across the components. In our model,
    we restricted all component variances to be equal
- The location parameters of mixture components are specified using default link
    functions which can be changed using a ‘link=’ option:

```
Distribution Default link Parameterisation
NormalN(μ, σ^2 ) identity μ
BernoulliB(p) logit ln[p/(1−p)]
Binomial B(n, p) logit ln[p/(1−p)]
Exponential Exp(λ) log ln(λ)
PoissonP(λ) log ln(λ)
LognormalLN(μ, σ^2 ) identity μ
```

- The PROBMODEL statement is used to specify starting values for the component
    probabilities using the logit link function, which can be changed using a ‘link=’
    option
- For our 4-component model, this implies specification of the following starting
    values: 
       
       
       
       
       
       
       

```
π 1 = 0. 1
π 2 = 0. 5
```
```
π 3 = 0. 3
π 4 = 0. 1
```
```








```
### −→

```








```
```
ln[π 1 /π 4 ] = 0
ln[π 2 /π 4 ] = 1. 6
```
```
ln[π 3 /π 4 ] = 1. 1
ln[π 4 /π 4 ] = 0
```
```








```
- Table with fit statistics (selection):

```
-2 Log Likelihood 1011.4
Effective Parameters 8
Effective Components 4
```

- Estimates for mixture components:

```
Parameter Estimates for ’Normal’ Model
Standard
Component Parameter Estimate Error z Value Pr > |z|
1 Intercept 3.4323 0.1665 20.62 <.0001
2 Intercept 5.3193 0.07456 71.34 <.0001
3 Intercept 7.6011 0.1119 67.94 <.0001
4 Intercept 10.3346 0.1887 54.76 <.0001
1 Variance 0.4474 0.06084
2 Variance 0.4474 0.06084
3 Variance 0.4474 0.06084
4 Variance 0.4474 0.06084
```
- Due to the ‘equate=scale’ option, all component variances are equal


- Estimates for component probabilities:

```
Parameter Estimates for Mixing Probabilities
----------------Linked Scale---------------
Standard
Component Parameter Estimate Error z Value Pr > |z| Probability
1 Probability 0.4255 0.3286 1.30 0.1953 0.1176
2 Probability 1.9382 0.2640 7.34 <.0001 0.5336
3 Probability 1.2647 0.2820 4.48 <.0001 0.2721
```
- Hence the fitted model is given by:

```
Y|μ ∼ N(μ, 0. 672 ) μ ∼
```
```




```
### 3 .43 5.32 7.60 10. 33

### 0 .12 0.53 0.27 0. 08

```




```

- SAS easily allows plotting the fitted mixture density with individual component
    densities:
       ods graphics on;
       proc fmm data=snapper plots=density(bins=15);
       model length = / dist=gaussian equate=scale k=4
          parms(3 0.5,5 0.5, 8 0.5, 10 0.5) ;
       probmodel / parms(0,1.6,1.1);
       freq frequency;
       run;
       ods graphics off;


- Omission of the option ‘equate=scale’ requires fitting of a normal mixture with
    component-specific means as well as variances:

```
ods graphics on;
proc fmm data=snapper plots=density(bins=15);
model length = / dist=gaussian k=4 parms(3 0.5,5 0.5, 8 0.5, 10 0.5) ;
probmodel / parms(0,1.6,1.1);
freq frequency;
run;
ods graphics off;
```
- Table with fit statistics (selection):

```
Fit Statistics
-2 Log Likelihood 977.2
Effective Parameters 11
Effective Components 4
```
- As expected the likelihood is larger (ℓℓ = − 488. 6 instead of ℓℓ =− 505. 7 before)


- Note that, since the likelihood is unbounded for this model, the estimation
    procedure converged to a local maximum
- Hence the reported estimates are not MLE’s and the likelihood value cannot be
    used for formal model comparison based on LR’s
- Estimates for mixture components:

```
Parameter Estimates for ’Normal’ Model
Standard
Component Parameter Estimate Error z Value Pr > |z|
1 Intercept 3.2130 0.06178 52.00 <.0001
2 Intercept 5.2596 0.07330 71.76 <.0001
3 Intercept 7.4428 0.1110 67.05 <.0001
4 Intercept 8.6186 1.2434 6.93 <.0001
1 Variance 0.06314 0.02432
2 Variance 0.3937 0.08644
3 Variance 0.2060 0.1270
4 Variance 3.0934 1.7754
```

- Estimates for component probabilities:

```
Parameter Estimates for Mixing Probabilities
----------------Linked Scale---------------
Standard
Component Parameter Estimate Error z Value Pr > |z| Probability
1 Probability -0.7376 0.6574 -1.12 0.2619 0.0948
2 Probability 0.9956 0.7092 1.40 0.1604 0.5365
3 Probability -0.1513 1.0026 -0.15 0.8801 0.1704
```
- Fitted model:

```
Y ∼ 0. 09 N(3. 21 , 0. 252 ) + 0. 54 N(5. 26 , 0. 632 )
+ 0. 17 N(7. 44 , 0. 452 ) + 0. 20 N(8. 62 , 1. 762 )
while the previous model was equal to
Y ∼ 0. 12 N(3. 43 , 0. 672 ) + 0. 53 N(5. 32 , 0. 672 )
+ 0. 27 N(7. 60 , 0. 672 ) + 0. 08 N(10. 33 , 0. 672 )
```

- Comparison of both fitted densities:
- The model with component-specific variances results in a less smooth density but
    seems to better capture the long right tail


I.4.3 Comparison of R with SAS

- Advantages of SAS procedure FMM:
    ⊲ Many more distributions available for the mixture components
    ⊲ Component densities not restricted to one parametric family (see later)
    ⊲ More flexibility to add restrictions to parameters or to set parameters equal to
       pre-specified values (see later)
    ⊲ More flexibility to include covariates in components and/or component
       probabilities (see later)
- Advantages of R package CAMAN:
    ⊲ Based on EM which is more stable than optimization in SAS
    ⊲ Less sensitive to starting values
    ⊲ Starting values on original scale (no link functions)
    ⊲ Estimation of g is possible (see later)


I.4.4 Example: Child Data

- Data set child.dat, on 602 pre-school children in north-eastern Thailand
- The response of interest is the number of episodes of acute respiratory infection
    (fever, cough, running nose,... ), recorded within a 3-year period, with histogram:


- The Poisson distribution is often used in practice for describing count data
- Since the average number of episodes equals 4.45, we first try to approximate the
    above histogram by the Poisson(4.45):


- Obviously, a single Poisson distribution cannot account for the large percentage of
    children with no or almost no episodes
- This can also be observed from comparing the sample mean with the sample
    variance:

```
y = 4. 45 ≪ 20 .45 = s^2 y
```
- Hence, there is more variability in the data than what can be explained from a
    single Poisson distribution
- This phenomenon is called overdispersion
- One way to take into account the overdispersion is modelling underlying
    heterogeneity using a finite mixture


- A 4-component Poisson mixture will be used.
- R code:

```
em<-mixalg.EM(obs="counts", weights="frequency", family="poisson",
data=child,t=c(0.5,3,10,15), p=c(0.25,0.25,0.25,0.25))
```
- SAS code:

```
proc fmm data=child;
model counts = / dist=poisson k=4 parms(-0.7, 1.1, 2.3, 2.7);
probmodel / parms(0,0,0);
freq frequency;
run;
```
- Note that the parameters in the Poisson densities are now specified on a
    logarithmic scale:
       (0. 5 , 3 , 10 ,15) → (ln(0.5),ln(3),ln(10),ln(15))


- Fitted model (ll =− 1553. 81 ):

```
Y|λ ∼ Poisson(λ) λ ∼
```
```




```
### 0 .143 2.817 8.164 16. 156

### 0 .197 0.480 0.270 0. 053

```




```
- Graphical representation:


- As derived earlier, the mean and variance of the obtained mixture can be
    calculated as
       E(Y) =

```
∑
j
```
```
πjλj = 4. 45
```
```
Var(Y) = ∑
j
```
```
πjλ^2 j −
```
```

∑
j
```
```
πjλj
```
```


```
```
2
+∑
j
```
```
πjλj = 20. 44
```
```
which are very close to the observed average (y = 4. 45 ) and observed variance
(s^2 y = 20. 45 ), illustrating that the mixture has taken account of the overdispersion
in the data.
```
- Biological interpretation: Latent variable represents the health status:

```
Component λj πj Interpretation
1 0.143 0.197 almost always healthy
2 2.817 0.480 normal
3 8.164 0.270 above normal
4 16.156 0.053 high risk for infection
```

- Note that the first component is a Poisson distribution with mean λ= 0. 143 ,
    which assigns probability 0.867 to the value Y = 0.
- One may wonder how much worse the model would be if the first component
    would be fixed at λ= 0, representing subjects who never experience acute
    respiratory infections
- In SAS, this can easily be achieved by mixing a degenerate distribution at 0 with a
    finite mixture of 3 Poisson distributions.
- SAS code:

```
proc fmm data=child;
model counts = / dist=constant(0) k=1;
model + / dist=poisson k=3 parms(1.1, 2.3, 2.7);
probmodel / parms(0,0,0);
freq frequency;
run;
```

- Fitted model (ll =− 1554. 4 ):

```
Y ∼ 0. 161 10 + 0. 496 Poisson(2.572)
+0. 286 Poisson(7.905) + 0. 057 Poisson(15.960)
rather than the previous model
Y ∼ 0. 197 Poisson(0.143) + 0. 480 Poisson(2.817)
+0. 270 Poisson(8.164) + 0. 053 Poisson(16.156)
```
- The model is only slightly worse in terms of likelihood:
    ll = − 1554. 4 versus ll= − 1553. 8
- Note that a classical LR test does not apply due to a boundary null-hypothesis
    H 0 : λ 1 = 0


I.4.5 Example: SIDS Data

- SIDS: Sudden Infant Death Syndrome
- Numbers of reported SIDS cases in 100 North Carolina counties, during the period
    1974–1978:


- As the counties are not of the same size, we need to correct for the number of
    live-births in each county
- Data set sids.dat
- Data structure:

```
County Yi ni Ri =Yi/ni
1 13 4672 0.00278
2 0 487 0.00000
3 15 1570 0.00955
... ... ... ...
```
```
100 16 14484 0.00110
```

- Terminology:

```
⊲ Observed counts Yi, i = 1,... , 100
⊲ Number ni of ‘exposed’ children
⊲ Rates: Ri = Ynii
```
- Histogram of the 100 rates:


- The histogram suggests the presence of heterogeneity among the counties
- This suggests that the counties are clustered with respect to their SIDS risk
- A 3-component mixture of Binomial distributions will be used
- R code:

```
em<-mixalg.EM(obs="frequency", pop.at.risk="nrisk", family="binomial",
data=sids, t=c(0.001,0.002,0.005), p=c(1,1,1))
```
- SAS code:

```
proc fmm data=sids ;
model frequency/nrisk = / dist=binomial k=3 parms(-6.9,-6.2,-5.3);
probmodel / parms(0,0);
run;
```

- Note that the parameters in the Binomial densities are now specified on a logit
    scale:
    (0. 001 , 0. 002 , 0 .005) → (ln(0. 001 / 0 .999),ln(0. 002 / 0 .998),ln(0. 005 / 0 .995))
- Fitted model (ll =− 233. 70 ):

```
Yi|p ∼ Binomial(ni, p) p ∼
```
```




```
### 0 .0013 0.0021 0. 0042

### 0 .33 0.53 0. 14

```




```
- Note that the mixture distribution cannot be super-imposed on the histogram of
    the rates (as in previous examples) since the mixture is the distribution of the
    counts Yi (all having different distributions) rather than the distribution of the
    rates Ri
- Later, it will be shown how the above mixture can be used to create so-called
    disease maps where counties are grouped based on their associated SIDS risk.


- Because a Binomial(n, p) distribution with large n and small p can be well
    approximated by a Poisson(np) distribution, disease rates are often modelled using
    Poisson models
- In R, the model can be specified as:

```
em<-mixalg.EM(obs="frequency", pop.at.risk="nrisk", family="poisson",
data=sids, t=c(0.001, 0.002, 0.005), p=c(1, 1, 1))
```
- Since SAS models the mean of a Poisson distribution on a log scale, an offset
    needs to be used:

```
E(Y) = np = exp[ln(n) + ln(p)]
```
- ln(n) is the offset that needs to be added to the linear predictor that models ln(p)


- Furthermore, starting values for p in the various mixture components need to be
    specified on a logarithmic scale, rather than the logit scale as before:
       (0. 001 , 0. 002 , 0 .005) → (ln(0.001),ln(0.002),ln(0.005))
- However, due to the small values of p, ln[p/(1−p)] ≈ ln(p)
- SAS code:

```
data sids;
set sids;
offset=log(nrisk);
run;
proc fmm data=sids ;
model frequency = / dist=poisson offset=offset k=3 parms(-6.9,-6.2,-5.3);
probmodel / parms(0,0);
run;
```

- Fitted model (ll =− 234. 41 ):

```
Yi|p ∼ Poisson(nip) p ∼
```
```




```
### 0 .0013 0.0021 0. 0042

### 0 .33 0.53 0. 14

```




```
```
which yields the same model for the risk parameter p as under the original
Binomial model (ll =− 233. 70 ):
```
```
Yi|p ∼ Binomial(ni, p) p ∼
```
```




```
### 0 .0013 0.0021 0. 0042

### 0 .33 0.53 0. 14

```




```

Chapter I.5

Inference for Number of Support Points

```
⊲ Introduction
```
```
⊲ Mixture of 2 known distributions
```
```
⊲ Some results
```
```
⊲ Conclusions
```

I.5.1 Introduction

- All examples discussed so far acted conditional on the number g of mixture
    components
- In this chapter, we will take a closer look at the selection for g
- An obvious approach is to fit models with increasing g
- As example, we continue the analysis of the snapper data (data set snapper.dat),
    and we fit mixtures of normals with common variance, for a variety of values g:

```
Yi ∼
∑g
j=1
```
```
πjN(μj, σ^2 )
```

- Summary of results: g

```

μ^1 μ^2 ··· μg
π 1 π 2 ··· πg
```
```

```
 σ (^2) ℓ
1

^6.^10
1

 3.60 -527.2
2

^5 .59 9.^22
0 .86 0. 14

 2.00 -515.65
3

^5 .05 7.60 10.^49
0 .65 0.29 0. 06

 1.11 -512.00
4

^3 .43 5.32 7.60 10.^33
0 .12 0.53 0.27 0. 08

 0.45 -505.72
5

^3 .40 5.31 7.50 9.68 11.^99
0 .12 0.52 0.26 0.08 0. 02

 0.32 -493.50
6

^3 .39 5.30 8.37 7.34 9.85 12.^04
0 .12 0.51 0.05 0.23 0.06 0. 03

 0.29 -492.65


- Graphical representation:

```
g = 1 g = 2 g = 3
```
```
g = 4 g = 5 g = 6
```

- The one-component mixture equals the normal distribution with the sample mean
    and the sample variance as mean and variance:
       Yi ∼ N(y, s^2 y)
- The residual variance σ^2 decreases as more components are added to the mixture.
- This is to be expected from the previously derived result Var(Y) =Var(X) +σ^2
    and from the fact that adding support points for X increases the variability of X.
- Adding components to the mixture increases the maximized log-likelihood value
- Selecting g requires some measure to compare models with different g
- One possible measure is the difference in maximized log-likelihood

```
=⇒ LR test
```
- However, this is not a standard testing procedure, as will be shown in the next
    section.


I.5.2 Mixture of 2 Known Distributions

- Suppose that for a continuous response Y, it is of interest to test whether the
    density of Y equals f 1 , versus the alternative that the density is a mixture of f 1
    and a ‘contaminating density’ f 2 :
       
       
       
       
       
       
       

```
H 0 :f(y) = f 1 (y)
```
```
HA : f(y) = πf 1 (y) + (1−π)f 2 (y)
```
- Note that this is equivalent with
    
    
    
    
    
    
    

```
H 0 : g = 1
```
```
HA: g = 2
```

- Let ℓ 1 and ℓ 2 denote the maximized log-likelihood values under the
    one-component model and the two-component model, respectively
- A classical likelihood ratio (LR) test would be based on the test statistic

```
ξ = 2(ℓ 2 −ℓ 1 )
and the p-value would be calculated assuming that
```
```
ξ −→H^0 χ^21 , when N → +∞
```
- However, we have the following result:

```
P[ξ = 0] −→H^00. 5
```

- Proof:

```
⊲ The log-likelihood function under the mixture equals
```
### ℓ 2 (Y) =

```
∑
i
```
```
ln[πf 1 (Yi) + (1−π)f 2 (Yi)]
```
```
⊲ The first-order derivative evaluated at π = 1 equals
```
### ∂ℓ 2 (Y)

```
∂π
```
```
∣∣∣
∣∣∣
∣∣
π=1
```
### = ∑

```
i
```
```
[f 1 (Yi)−f 2 (Yi)]
[πf 1 (Yi) + (1−π)f 2 (Yi)]
```
```
∣∣∣
∣∣∣
∣∣
π=1
```
### = N − ∑

```
i
```
```
f 2 (Yi)/f 1 (Yi)
```
```
⊲ Note that
```
### E

```

f 2 (Yi)
f 1 (Yi)
```
```
∣∣∣
∣∣∣
∣∣ H^0
```
```

 = ∫[f
2 (y)/f 1 (y)]f 1 (y)dy = 1
```

⊲ It now immediately follows from the Central Limit Theorem that

### P

```



```
### ∂ℓ 2 (Y)

```
∂π
```
```
∣∣∣
∣∣∣
∣∣
π=1
```
### > 0

```
∣∣∣
∣∣∣
∣∣ H^0
```
```



```
### −→H^00. 5

⊲ The second-order derivative of ℓ 2 equals:

### ∂^2 ℓ 2 (Y)

```
∂π^2
```
### = −∑

```
i
```
```
[f 1 (Yi)−f 2 (Yi)]^2
[πf 1 (Yi) + (1−π)f 2 (Yi)]^2
```
### < 0

⊲ Hence, we have that ℓ 2 is concave, and that the first-order derivative at π = 1
is strictly positive with 50% chance


⊲ Graphically:

⊲ This implies that

```
P(π̂ = 1) −→H^00. 5
from which the result immediately follows
```

- The above property states that in half of the cases, the LR test statistic ξ will be
    exactly zero, for sufficiently large samples.
- Obviously the classical χ^2 -approximation is not valid
- The reason is that H 0 is on the boundary of the parameter space under HA
- This can be seen from rewriting the hypothesis as

```




```
```
H 0 : π = 1
HA : π < 1
```
- In general, the classical LR test cannot be used for testing for the number of
    components in a finite mixture
- The asymptotic distribution is to be derived for each testing problem separately


I.5.3 Some Results

- Testing H 0 : π = 1 in a mixture of two known densities:

```
ξ −→H^00. 5 χ^20 + 0. 5 χ^21
where χ^20 equals the discrete distribution with all probability mass at 0
```
- Testing H 0 : π = 1 in a mixture of two normals with common unknown variance:

```
ξ −→H^0 χ^22
```
- Testing H 0 : π = 1 in a mixture of two Binomials Bin(2, p 1 ) and Bin(2, p 2 ):

```
ξ −→H^00. 5 χ^20 + 0. 5 χ^21
```

- Testing H 0 : π = 1 in a mixture of two Poissons Poisson(λ 1 ) and Poisson(λ 2 ) with
    small λj (λj ∈[0, 0 .1]):

```
ξ −→H^00. 5 χ^20 + 0. 5 χ^21
```
- In general, simulation methods are used to study the asymptotic behavior of ξ
- If there is any convergence at all, it can be painfully slow
- Thode, Finch, & Mendell (Biometrics, 1988, 44:1195–1201):

```
“... one would need what is usually an infeasible large sample size
(N > 1000 ) for the use of large-sample approximations to be justified.”
```
- In practice, the null-distribution of ξ can be derived via bootstrap methods, which
    will not be discussed here any further.


I.5.4 Conclusions

- Inference on the number g of components in a finite mixture is far from
    straightforward
- One way to avoid the selection of g is to treat g as a parameter in the likelihood,
    and to estimate g from the available data

```
=⇒ Non-Parametric
Maximum
Likelihood
Estimation
(NPMLE)
```
- Note that this involves more than classical ML estimation theory as the number of
    parameters in the likelihood is not fixed: The number of support points as well as
    the number of associated component probabilities increases with g.


Chapter I.6

Non-parametric Maximum Likelihood

```
⊲ Introduction
```
```
⊲ Definition of NPMLE
```
```
⊲ Characterization of NPMLE
```
```
⊲ Examples
```

I.6.1 Introduction

- All examples discussed so far acted conditional on the number g of mixture
    components, i.e., on the number g of support points for the latent variable X
- In this chapter, it will be investigated how g can be estimated from the data
- As an example, we analyze the SIDS data set using mixtures of g Poisson
    distributions, for varying g:

```
Yi ∼
∑g
j=1
```
```
πjPoisson(pjni)
```
```
or equivalently
```
```
Yi|p ∼ Poisson(pni) p ∼
```
```




```
```
p 1 p 2 ··· pg
π 1 π 2 ··· πg
```
```




```

- Summary of results: g

```

p^1 p^2 ··· pg
π 1 π 2 ··· πg
```
```

 ℓ
```
```
1
```
```

^0.^0020
1
```
```

 -255.58
```
```
2
```
```

^0 .0016 0.^0035
0 .75 0. 25
```
```

 -237.28
```
```
3
```
```

^0 .0012 0.0021 0.^0042
0 .33 0.53 0. 14
```
```

 -234.41
```
```
4
```
```

^0 .0013 0.0021 0.0037 0.^0090
0 .32 0.52 0.15 0. 01
```
```

 -233.40
```
```
5
```
```

^0 .0013 0.0021 0.0037 0.0037 0.^0090
0 .32 0.52 0.11 0.04 0. 01
```
```

 -233.40
```
```
6
```
```

^0 .0013 0.0021 0.0037 0.0037 0.0037 0.^0090
0 .32 0.52 0.09 0.05 0.01 0. 01
```
```

 -233.40
```

- Once the latent variable X has 4 support points, the log-likelihood value cannot
    be increased anymore by including more support points.
- It can be shown that for any g > 4 , the maximized log-likelihood equals − 233. 40.
- This suggests using a 4-component mixture to describe the data.
- The resulting estimate for the distribution of the latent variable X is called a
    NPMLE: It maximizes the log-likelihood value over the class of all distributions
    for X.


I.6.2 Definition and Properties of NPMLE

- Let fi(yi|x) denote the density function (continuous or discrete) of Yi given the
    latent variable X
- Let G be the distribution function of X (continuous or discrete)
- G is called the mixing distribution.
- The marginal density of Yi equals
    fi(y|G) =

```
∫
fi(yi|x)dG(x)
where the integral becomes a sum in case X is discrete.
```
- The log-likelihood is then obtained as

```
ℓ(G) =
```
```
∑N
i=1
```
```
ln[fi(yi|G)]
```

- In many applications (with discrete responses) data incorporate replications such
    that ℓ(G) is of the form

```
ℓ(G) =
∑m
i=1
```
```
ωiln[fi(yi|G)]
```
```
where there are only m different values, each occurring ω 1 ,... , ωm times.
```
- Note that possible dependence on unknown parameters is suppressed in the above
    notation
- A NPMLE for G is any distribution function Ĝ for which ℓ(G) is maximized over
    the class of all distributions, i.e.,

```
ℓ(Ĝ) = maxG∈Γ ℓ(G)
```

- Note that the log-likelihood is maximized over the class Γ of all distributions, i.e.,
    discrete as well as continuous distributions.
- Property:

```
The log-likelihood ℓ(G) is concave in Γ
```
- This implies that ℓ(G) has a unique mode.
- Hence, we have that

```
A NPMLE Ĝ exists
```
- In many cases, Ĝ will be unique. However, it is not in general.


- Further, it can be shown that

```
Ĝ is discrete with at most m support points
```
- The discreteness of Ĝ allows to maximize ℓ(G) over the class of all discrete
    distributions only.
- Let Ω be the class of all discrete distributions. We then have that

```
ℓ(Ĝ) = maxG∈Γ ℓ(G) = maxG∈Ω ℓ(G)
```
- The upper bound for the number of support points is seldom sharp, i.e., there will
    often be (many) less support points than indicated by the bound.


I.6.3 Characterization of a NPMLE

6.3.1 Directional Derivative and Gradient Function

- For G and H in Ω, the directional derivative of ℓ(·) at G into the direction H is
    defined as

```
Φ(G, H) = lim
α→> 0
```
```
ℓ[(1−α)G+αH]−ℓ(G)
α
```

- Graphical interpretation:


- Note that

```
Φ(G, H) = lim
α→> 0
```
```
ℓ[(1−α)G+αH]−ℓ(G)
α
```
```
=
```
```
∂ ℓ[(1−α)G+αH]
∂α
```
```
∣∣∣
∣∣∣
∣∣
α=0
```
```
=
```
### ∂

```
∑
i
```
```
ln[(1−α)fi(yi|G) +αfi(yi|H)]
∂α
```
```
∣∣∣
∣∣∣
∣∣∣
∣α=0
```
```
= ∑
i
```
```
fi(yi|H) − fi(yi|G)
fi(yi|G)
```
```
= ∑
i
```
```
fi(yi|H)
fi(yi|G)
```
### − N


- For every discrete distribution H in Ω,

### H =

```




```
```
x 1 x 2 ··· xg
```
```
π 1 π 2 ··· πg
```
```




```
### ,

```
we have that
```
### Φ(G, H) = ∑

```
i
```
```
fi(yi|H)
fi(yi|G)
```
### − N = ∑

```
i
```
```
∑
j
```
```
πjfi(yi|xj)
```
```
fi(yi|G)
```
### − N

### = ∑

```
j
```
```
πj
```
```

∑
i
```
```
fi(yi|xj)
fi(yi|G)
```
```

 − N = N
```
```

∑
j
```
```
πjd(G, xj) − 1
```
```


```
```
with
```
```
d(G, x) =
```
### 1

### N

```
∑
i
```
```
fi(yi|x)
fi(yi|G)
```

- d(G, x) is called the gradient function of G, evaluated at x
- A NPMLE can now be defined as any Ĝ in Ω such that

```
Φ(G, Ĥ ) ≤ 0 , for all H inΩ
```
- Hence, the NPMLE will be characterized with the gradient function d(G, x)
- Three theorems are useful to check if a candidate estimate Ĝ is really NPML.


- Theorem 1:

```
Ĝ is NPMLE if and only if for all x, d(G, x̂ ) ≤ 1
```
- Theorem 2:

```
If Ĝ is a NPMLE, we have that d(G, x̂ ) = 1 for all support points x of Ĝ
```
```
d(G, x̂ ) is identically one if and only if Ĝ is not unique
```
- Theorem 3:

```
If all fi(yi|x), as functions of x, have unique modes in some interval
[a, b], then Ĝ can only have support points in the interval [a, b].
```

- Theorem 1 provides a tool to check whether Ĝ is a NPMLE
- Theorem 2 provides a tool to check uniqueness
- Theorem 3 allows to restrict attention to the interval [a, b]


I.6.4 Example: The SIDS Data

- Under the Poisson(pni) model, we have that

```
fi(yi|p) = exp(−pni)
```
```
(pni)yi
yi!
which, as a function of p, is uniquely maximized for p =yi/ni.
```
- Theorem 3 then implies that a NPMLE Ĝ will have support points between the
    smallest and the largest observed rate, i.e., in the interval [0, 0 .0096]
- Gradually increasing the number of components in a mixture suggested that the
    following 4-component mixture is NPML:

```
Yi|p ∼ Poisson(pni) p ∼ G =
```
```




```
### 0 .0013 0.0021 0.0037 0. 0090

### 0 .33 0.51 0.15 0. 01

```




```

- The gradient function d(G, p) for the above mixing distribution G can be
    obtained using the SAS code:
       data test; set sids;

```
data test; set test;
do x=0 to 0.01 by 0.0001; output; end;
```
```
data test; set test;
gradient=pdf(’POISSON’,y,x*n)/(0.3263*pdf(’POISSON’,y,0.001254*n)
+ 0.5124*pdf(’POISSON’,y,0.002081*n)
+ 0.1505*pdf(’POISSON’,y,0.003747*n)
+ 0.0108*pdf(’POISSON’,y,0.009013*n));
```
```
proc sort data=test; by x;
```
```
proc means data=test; var gradient; by x;
output out=out;
```
```
data out; set out; if _STAT_=’MEAN’;
```

title h=2.5 ’Gradient SIDS’;
proc gplot data=out;
plot gradient*x / nolegend haxis=axis1 vaxis=axis2 cvref=blue vref=1;
symbol c=red i=join w=5 l=1 mode=include;
axis1 label=(h=2 ’Risk p’) value=(h=1.5) order=(0 to 0.01 by 0.002)
minor=none w=5 ;
axis2 label=(h=2 angle=90 ’Gradient’) value=(h=1.5) order=(0.4 to 1.1 by 0.1)
minor=none w=5 ;
run;quit;


- Result:

# ✁

# ✄☎✆✂

# ✞✝

# ✟✠✡

# ✟✠☛

# ✟✠☞

# ✟✠✌

# ✟✠✍

# ✟✠✎

# ✏✠✟

# ✏✠✏

# ✑✒✓✔✕

# ✟✠✟✟✟✟✠✟✟✖✟✠✟✟✡✟✠✟✟☞✟✠✟✟✍✟✠✟✏✟

# ✗✘✙✚✛✜✢✣✤✥✦✤

- G is NPMLE because d(G, p) ≤ 1 in the interval [0, 0 .0096]
- Note also that d(G, p) = 1 for p in { 0. 0013 , 0. 0021 , 0. 0037 , 0. 0090 }
- The obtained estimate is unique as the gradient function is not identically one.


I.6.5 Example: Accident Data

- We now consider the number of accident claims during one year, out of 9461
    insurance policies issued by La Royal Belge Insurance Company:

```
Count yi: 0 1 2 3 4 5 6 7
Frequency ωi: 7840 1317 239 42 14 4 4 1
```
- These data have been analyzed frequently in the statistical literature:
    ⊲ Thyrion (Astin Bulletin, 1960, 1:142–162)
    ⊲ Simar (Annals of Statistics, 1976, 4:1200–1209)
    ⊲ Carlin & Louis (Chapman & Hall, 1996)
    ⊲ B ̈ohning (Chapman & Hall, 1999)


- Under the Poisson(λ) model, we have that

```
fi(yi|λ) =
```
```
e−λ λyi
yi!
which, as a function of λ, is uniquely maximized for λ = yi.
```
- Theorem 3 then implies that a NPMLE Ĝ will have support points between the
    smallest and the largest observation, i.e., in the interval [0,7]
- Simar and Carlin & Louis report that a NPMLE is given by

```
Y|λ ∼ Poisson(λ) λ ∼ G =
```
```




```
### 0 .089 0.580 3.176 3. 669

### 0 .7600 0.2362 0.0037 0. 0002

```




```
- The reported maximized log-likelihood value equals ℓ = − 5341. 5310.


- The gradient function d(G, p) for the above mixing distribution G equals:

# ✧ ★

# ✪✫ ✬✩

# ✮✭

# ✯ ✰ ✱ ✲

# ✯ ✰ ✱ ✳

# ✯ ✰ ✱ ✱

# ✴ ✰ ✯ ✯

# ✴ ✰ ✯ ✴

# ✴ ✰ ✯ ✵

# ✯ ✴ ✵ ✶ ✷

# ✸ ✹ ✺ ✻ ✼ ✽ ✾ ✿ ❀ ❁ ❁ ✼ ✻ ✽ ✾ ✿ ❂ ✺ ✿ ✺

- Hence, G is not NPMLE because d(G, p) > 1 in the neighborhood of 0, and does
    not reach 1 at all support points


- B ̈ohning reports that a NPMLE is given by

```
Y|λ ∼ Poisson(λ) λ ∼ G =
```
```




```
### 0 .000 0.336 2. 545

### 0 .4184 0.5730 0. 0087

```




```
- The corresponding maximized log-likelihood value now equals ℓ =− 5340. 7040 ,
    which is indeed larger than the value reported by Simar and Carlin & Louis, and
    the gradient function equals:

# ❃ ❄

# ❆❅

# ❇ ❈

# ❊❉

# ❋ ● ❍ ❍ ❍

# ■ ● ❋ ❋ ❋

# ■ ● ❋ ❋ ■

# ❋ ■ ❏ ❑ ▲

# ▼ ◆ ❖ P ◗ ❘ ❙ ❚ ❯ ❱ ❱ ◗ P ❘ ❙ ❚ ❲ ❖ ❚ ❖


Chapter I.7

Numerical Algorithms

```
⊲ CAMAN approach to NPMLE
```
```
⊲ Vertex exchange method (VEM)
```
```
⊲ Stopping rule
```

I.7.1 CAMAN Approach to NPMLE

7.1.1 Introduction

- The theorems discussed earlier imply that, under mild regularity conditions, a
    NPMLE Ĝ is discrete, with support ‘within the range of the data’.
- Moreover, the gradient d(G, x̂ ) of Ĝ in any direction x is never larger than 1, and
    equals 1 for all support points of Ĝ.
- In general, although Ĝ will be discrete, the support could be large, especially for
    continuous data, or when analyzing rates.
- One approach could be to start the EM algorithm discussed earlier, for a mixture
    model with a ‘sufficiently large’ number of components g.


- The idea is then that, if g is larger than the support size of Ĝ, some of the
    estimated support points will coincide, or some support points will get weight
    (probability) zero.
- Coinciding support points have already been obtained earlier in the analysis of the
    SIDS data:

```
g
```
```

p^1 p^2 ··· pg
π 1 π 2 ··· πg
```
```

 ℓ
```
```
4
```
```

^0 .0013 0.0021 0.0037 0.^0090
0 .32 0.52 0.15 0. 01
```
```

 -233.40
```
```
5
```
```

^0 .0013 0.0021 0.0037 0.0037 0.^0090
0 .32 0.52 0.11 0.04 0. 01
```
```

 -233.40
```
```
6
```
```

^0 .0013 0.0021 0.0037 0.0037 0.0037 0.^0090
0 .32 0.52 0.09 0.05 0.01 0. 01
```
```

 -233.40
```

- Indeed, as was seen later, the NPMLE Ĝ was unique, and contained only 4
    support points.
- However, if the EM algorithm is started from a mixing distribution with many
    support points (e.g., g = 25), g will often be much too large, leading to extremely
    slow convergence of the algorithm.
- In CAMAN, this is solved by splitting up the estimation procedure in two different
    phases
- The procedure will be illustrated in a specific example where the NPMLE consists
    of 3 support points


- Graphical representation of final solution:

# ❳❨❩

# ❬❭❬

# ❪❫❪❴❵

# ❛❜❝❝❞❡❢


7.1.2 Phase 1 of CAMAN

- Let [a, b] be the interval that will contain all support points (Theorem 3)
- A large grid Λ ={a= x 1 ,... , xL= b} is specified as ‘first guess’ for the support
    points of Ĝ:

# ❣ ❤

# ✐ ❥❦

# ❧ ♠❧

# ♥♦♥ ♣ q

# r s t t ✉ ✈ ✇

- The log-likelihood ℓ(G) is maximized over this grid, i.e., over all probability
    distributions with support Λ


- This only requires estimation of the corresponding weights {π 1 ,... , πL} and
    possibly also parameters in the models fi(yi|x) (e.g., the variance σ^2 in the
    normal model).
- This can be done using any of the numerical methods which will be discussed in
    the next sections.
- If L was chosen sufficiently large, Phase 1 results in many grid points with zero
    weight, while grid points in the region of the final solution receive positive weight:

# ①②

# ③④⑤

# ⑥⑦⑥

# ⑧⑨⑧⑩❶

# ❷❸❹❹❺❻❼


- From now on, it will be assumed that the number g of support points in Ĝ is at
    most the number of grid points with positive weight.
- This assumption is only justified if the number L of grid points was chosen
    sufficiently large, i.e., if the grid can be assumed to be a good approximation to
    the support of Ĝ.


7.1.3 Phase 2 of CAMAN

- All grid points with positive weight in Phase I, together with the corresponding
    weights, are used as starting values for the EM algorithm, discussed earlier
- So, in this phase, weights as well as support points are estimated, while in Phase
    1 this was only the case for the weights.
- Note also that, in this phase, the number of support points is also kept fixed, but
    this number is now (much) smaller than in Phase I.
- In case the number of grid points with positive weight in Phase 1 was still larger
    than the number g of support points in Ĝ, than some estimated support points
    will coincide, or will receive zero weight.
- Coinciding points are combined. Points with zero weight are left out.


- In practice, the effect of Phase 2 is that the grid points with positive support
    obtained in Phase 1 converge to the final solution with combined weights:

# ❽❾

# ❿➀➁

# ➂➃➂

# ➄➅➄➆➇

# ➈➉➊➊➋➌➍

- In order to assure that the so-obtained estimate Ĝ is really NPMLE, it should be
    checked that the gradient function d(G, x̂ ) is never larger than 1, and equals 1 for
    all support points of Ĝ.


I.7.2 Vertex Exchange Method (VEM)

- The first phase in the CAMAN approach to the calculation of a NPMLE Ĝ is the
    fitting of a finite mixture model, with large but fixed support for the mixing
    distribution G.
- Several algorithms have been proposed but the vertex exchange method (VEM),
    implemented in CAMAN, is the most efficient one so far.
- All algorithms maximize the mixture log-likelihood ℓ(G) over the class of all
    discrete mixing distributions with support equal to (a subset of)
    Λ = {x 1 ,... , xL}.
- Let G be a current guess for the final solution. Improving G implies reducing the
    weight of some support points, while increasing the weight of others.


- General idea:

```
G can be improved by moving weight from
a ‘bad’ support point to a ‘good’ one
```
- VEM will replace G by

```
G−απ-Gx- +απ-Gx+ = G+απ-(Gx+−Gx-),
```
```
for some α ∈[0,1] and support points x- and x+, and with Gx representing the
degenerate distribution with support x.
```
- All support points of G, except x- and x+, keep their original weights, while a
    proportion α of the weight π- of the ‘bad’ support point x- is moved to a ‘good’
    support point x+.
- α is called the step-length


- The ‘optimal’ choice for x-, x+ and α is the one which maximizes the gain in
    log-likelihood, i.e., which maximizes

```
ℓ[G+απ-(Gx+ −Gx−)] −ℓ[G]
```
- First-order Taylor approximation for small α yields

```
ℓ[G+απ-(Gx+ −Gx-)] −ℓ[G]
```
```
≈ α
```
```
∂ℓ[G+απ-(Gx+−Gx-)]
∂α
```
```
∣∣∣
∣∣∣
∣∣α=0
```
```
= α
```
```
∂ℓ{(1−α)G+α[G+π-(Gx+ −Gx-)]}
∂α
```
```
∣∣∣
∣∣∣
∣∣
α=0
= α Φ[G, G+π-(Gx+ −Gx-)]
```
```
= α
```
```

∑
i
```
```
fi(yi | G+π-(Gx+ −Gx-))
fi(yi|G)
```
### − N

```


```

```
= α
```
```



```
```
∑
i
```
```
fi(yi|G)−π-fi(yi|x-) +π-fi(yi|x+)
fi(yi|G)
```
### − N

```



```
```
= α
```
```

N − π-∑
i
```
```
fi(yi|x-)
fi(yi|G)
```
```
+ π-∑
i
```
```
fi(yi|x+)
fi(yi|G)
```
### − N

```


```
```
= α N π- [d(G, x+) − d(G, x-)]
```
- Obviously, the best choice for x+ is the support point in Λ for which the gradient
    function is maximal
- Further, the best choice for x- is the support point in Λ for which the gradient
    function is minimal


- Given G, an updated version of G is obtained from the following algorithm:

```
⊲ Find x+ ∈Λ which maximizes d(G, x)
⊲ Find x- ∈ Λ which minimizes d(G, x)
⊲ Find α which maximizes
ℓ[G+απ-(Gx+ −Gx-)] −ℓ[G]
over α
⊲ Replace G by G+απ-(Gx+ −Gx-)
```
- This algorithm is repeated until convergence.
- Note that maximization with respect to α usually requires iterative optimization
    procedures


- Graphical representation of selection of x- and x+:


- Let G(t) be any sequence created by the above algorithm, and let Ĝ be NPMLE,
    then one can show that

```
ℓ[G(t)] −→ ℓ[Ĝ], monotonously
```
```
provided the grid Λ = {x 1 ,... , xL} is sufficiently dense.
```

I.7.3 Stopping Rule

- The numerical algorithms discussed earlier all result in a sequence
    {
       G(1), G(2),... , G(t), G(t+1),...

```
}
```
```
which converges to a NPMLE Ĝ.
```
- In practice, one needs a stopping rule to decide when the iterative process is
    terminated, i.e., which G(t) will be considered to be sufficiently close to Ĝin order
    to be acceptable as NPMLE.
- We know from the first theorem of the characterization of NPMLE’s that
    d(G, x̂ ) ≤ 1 ,∀x.


- An obvious stopping rule is then to select a small value ε > 0 , and to stop the
    iterative procedure at the smallest t for which

```
d
```
```
(
G(t), x
```
```
)
≤ 1 + ε, for all x
```
- For VEM, it is sufficient to stop the estimation process as soon as

```
d
```
```
(
G(t), x+
```
```
)
≤ 1 + ε
which needs to be calculated anyway.
```
- In CAMAN, the ε is specified as the accuracy
- In order to guarantee that the iterative procedure will eventually stop, a maximal
    number of iteration steps is also required.


Chapter I.8

Examples in CAMAN

```
⊲ Child data
```
```
⊲ Snapper data
```

I.8.1 Child Data

- We will now illustrate the use of CAMAN for the estimation of the NPMLE for
    the Child data.
- Recall that the response of interest is the number of episodes of acute respiratory
    infection (fever, cough, running nose,... ), recorded within a 3-year period, with
    histogram:


- The model equals

```
Y|λ ∼ Poisson(λ) λ ∼
```
```




```
```
λ 1 λ 2 ··· λg
π 1 π 2 ··· πg
```
```




```
```
with unknown number g of components in the mixture.
```
- In CAMAN, three related procedures are available:

```
⊲ mixalg.VEM: Phase 1
⊲ mixalg.EM: Phase 2 (discussed before)
⊲ mixalg: Phase 1 & Phase 2
```
- Phase 1 can be performed using the following syntax:

```
vem<-mixalg.VEM(obs="counts", weights="frequency", family="poisson",
data=child, acc=10^(-8),numiter=5000, startk=50)
```

- The option ‘acc=’ specifies the accuracy ε used in the stopping rule
- The option ‘numiter=’ specifies the maximum number of iteration steps
- The option ‘startk=’ specifies the number of grid points in the initial grid
    Λ = {x 1 ,... , xL}
- Results:

```
Computer Assisted Mixture Analysis (VEM):
Data consists of 602 observations (rows).
The VEM-algorithm identified 8 grid points with positive support
p t
1 0.1151572058 0.0000000
2 0.0929999812 0.4897959
3 0.0627187858 2.4489796
4 0.4099903721 2.9387755
5 0.0608071940 7.8367347
6 0.2051281574 8.3265306
7 0.0522385437 16.1632653
8 0.0009597599 16.6530612
Log-Likelihood: -1553.883 BIC: 3741.392
```

- Out of the initial 50 grid points, only 8 received positive weight after the VEM
    run. All other grid points received weight 0.
- The maximized log-likelihood value after finalizing Phase 1 of the CAMAN
    procedure equals − 1553. 883
- Positive weight is often given to neighboring points:

```
Grid value Weight
··· ···
7. 84
8. 33
```
### 0. 0608

### 0. 2051

### ··· ···

### 16. 16

### 16. 65

### 0. 0522

### 0. 0010

### ··· ···


- Phase 2 can be performed using the ‘mixalg.EM’ procedure, with the results from
    the VEM run as input
- Alternatively, Phase 1 & Phase 2 can be performed jointly using:

```
npml<-mixalg(obs="counts", weights="frequency", family="poisson",
data=child, acc=10^(-8),numiter=50000, startk=50)
```
- Results:

```
Computer Assisted Mixture Analysis:
Data consists of 602 observations (rows).
The Mixture Analysis identified 5 components of a poisson distribution:
p lambda
1 7.740583e-06 0.0000000
2 1.969225e-01 0.1433966
3 4.799752e-01 2.8172852
4 2.692583e-01 8.1641705
5 5.383626e-02 16.1558261
Log-Likelihood: -1553.81 BIC: 3165.223
```

- Two of the original 8 support points have converged to the values 2. 8173 , 8. 1642 ,
    and 16. 1558.
- The final solution has 5 distinct support points only.
- CAMAN allows to combine identical support points, i.e., support points which
    differ less than a limit, which can be pre-specified using an additional option
    ‘limit=’:
       ⊲ Support points with weights less than ‘limit’ are deleted
       ⊲ Support points less than ‘limit’ different are combined
- The default limit as well as some additional information about input parameters
    and results can be obtained with:

```
summary(npml)
```

- Results:

```
number of VEM-iterations done: 7445
alteration within the final VEM-iteration step: 9.9405e-0 9
number of EM-iterations done: 34306
alteration within the final EM-iteration step: 9.997101e- 09
User-defined parameters:
max number of iterations: 50000
limit for combining components: 0.01
threshold for converging: 1e-08
number of grid points (startk): 50
```
- The maximized log-likelihood after Phase 2 equals − 1553. 81 which is only a minor
    improvement compared to the approximation obtained from the first phase
    (log-likelihood − 1553. 883 ).
- This minor further increase in log-likelihood required 34306 additional EM-steps,
    further illustrating the slow convergence of the EM-algorithm


- This suggests that the 8-point result from the first phase was already a (very)
    good approximation of the full NPMLE Ĝ.
- This also explains the results from Phase 1:

```
Grid value Weight
··· ···
7. 84
8. 33
```
```
0. 0608
0. 2051
```
```

= 0.^2659 ≈^0.^2693 at^8.^1642
··· ···
16. 16
16. 65
```
```
0. 0522
0. 0010
```
```

= 0.^0532 ≈^0.^0538 at^16.^1558
··· ···
```
- The neighboring points with positive weights suggest that the NPMLE has a
    support point somewhere between the two neighbors, with weight approximately
    equal to the sum of the weights of the neighbors.


- The fitted model can be summarized as:

```
Y|λ ∼ Poisson(λ) λ ∼
```
```




```
### 0 0.143 2.817 8.164 16. 156

### 0 .001 0.197 0.480 0.269 0. 054

```




```
- Graphical representation:


I.8.2 Snapper Data

- We will now illustrate the use of CAMAN for the estimation of the NPMLE for
    the Snapper data.
- Recall that the response of interest is the length, with histogram:


- The model equals

```
Y|μ ∼ N(μ, σ^2 ) μ ∼
```
```




```
```
μ 1 μ 2 ··· μg
π 1 π 2 ··· πg
```
```




```
```
with unknown number g of components in the mixture.
```
- A NPML estimate can be obtained in CAMAN using following syntax:

```
npml<-mixalg(obs="length", weights="frequency", family="gaussian",
data=snapper, acc=10^(-8),startk=50)
```
- The result equals:

```
Computer Assisted Mixture Analysis:
Data consists of 256 observations (rows).
The Mixture Analysis identified 1 components of a gaussian distribution:
```

```
DETAILS:
p mean
1 1 6.103516
component variance: 0
Log-Likelihood: -563.7445 BIC: 1133.034
```
- Only 1 component is identified, and the component variance is set equal to 0,
    indicating a problem.
- Indeed the specified model is unidentified, which follows from

```
Var(Y) = σ^2 + Var(μ)
```
- The only term identified from the data is Var(Y), estimated by the sample
    variance s^2 y = 3. 60
- How to split the total variance into within-component variability σ^2 and
    between-comonent variability Var(μ) is to be decided by the user.


- When NPMLE is the objective, the within-component variance σ^2 needs to be
    pre-specified using a ‘var=’ option.
- Obviously, σ^2 should be set equal to some value less than s^2 y = 3. 60
- NPMLE for σ^2 = 3:

```
p mean
1 0.9304275 5.826143
2 0.0695725 9.812964
Log-Likelihood: -519.7317 BIC: 1056.099
```
- NPMLE for σ^2 = 2:

```
p mean
1 0.81525492 5.504159
2 0.15510549 8.370224
3 0.02963959 10.727394
Log-Likelihood: -514.9862 BIC: 1057.698
```

- NPMLE for σ^2 = 1:

```
p mean
1 0.07929415 3.975117
2 0.57967082 5.207934
3 0.26631879 7.544679
4 0.06047305 9.793590
5 0.01424319 11.787159
Log-Likelihood: -510.9503 BIC: 1071.807
```
- NPMLE for σ^2 = 0. 2 :

```
p mean
1 0.114842014 3.336916
2 0.381770078 5.069379
3 0.164837218 5.981665
4 0.209089544 7.444038
5 0.050440157 8.489854
6 0.052514104 9.747340
7 0.009905078 10.497398
8 0.004406573 11.187444
9 0.012195233 12.226948
Log-Likelihood: -488.2221 BIC: 1070.712
```

- Summary:

```
σ^2 ℓ ̂g
3 − 519. 73 2
2 − 514. 99 3
1 − 510. 95 5
0.2 − 488. 22 9
```
- The number of components ranges from 2 to 9
- The maximized log-likelihood values show considerable variation
- This suggests quite different fits of the models to the observed data


- This is also seen in the resulting mixture densities:

```
σ^2 = 3 σ^2 = 2
```
```
σ^2 = 1 σ^2 = 0. 2
```

- Specifying a large value for σ^2 results in a mixture with a small number of
    components with much variability, and therefore in a (very) smooth mixture
    density.
- Specifying a small value for σ^2 results in a mixture with a large number of
    components with little variability, and therefore in a (very) non-smooth mixture
    density.
- From this perspective, the above procedure can be viewed as a method for density
    estimation, with smoothness parameter σ^2.
- Note that the fact that g and σ^2 are not simultaneously identified from the data is
    due to absence of a variance-mean link in the Gaussian family.


- For other models (Poisson, Binomial,... ), this problem does not occur since the
    variance is immediately tied to the mean:

```
Y ∼Poisson(λ) ⇒ Var(Y) = λ =E(Y)
Y ∼Binomial(n, p) ⇒ Var(Y) = np(1−p) = E(Y)[n−E(Y)]/n
which shows that there is no additional parameter which can be used to tune the
variability within the mixture components.
```

Chapter I.9

Classification

```
⊲ Introduction
```
```
⊲ Posterior probabilities
```
```
⊲ Examples
```
```
⊲ Cluster analysis versus discriminant analysis
```

I.9.1 Introduction

- We re-consider the Child data, where the response of interest is the number of
    episodes of acute respiratory infection (fever, cough, running nose,... ), recorded
    within a 3-year period
- The NPMLE for the mixing distribution in a Poisson model was earlier found to be
    (ll= − 1553. 81 ):

```
Y|λ ∼ Poisson(λ) λ ∼
```
```




```
### 0 0.143 2.817 8.164 16. 156

### 0 .001 0.197 0.480 0.269 0. 054

```




```
- The extremely small weight for the first component motivates a 4-component
    mixture.


- The fitted model, obtained before, then becomes (ll =− 1553. 81 ):

```
Y|λ ∼ Poisson(λ) λ ∼
```
```




```
### 0 .143 2.817 8.164 16. 156

### 0 .197 0.480 0.270 0. 053

```




```
- Note the minor change in log-likelihood (< 0. 01 ).
- Graphical representation:


- As discussed before, a biological interpretation can be given to the components of
    the mixture:

```
Component λj πj Interpretation
1 0.143 0.197 almost always healthy
2 2.817 0.480 normal
3 8.164 0.270 above normal
4 16.156 0.053 high risk for infection
```
- Once a mixture model has been fitted, one might be interested in classifying
    observations in the different mixture components, i.e., in deciding what
    component of the mixture a specific observation is most likely to belong to.
- In practice this is often done based on posterior probabilities.


I.9.2 Posterior Probabilities

- Consider the following finite mixture model for the response Y of interest:

```
Yi ∼ π 1 fi 1 (yi) + π 2 fi 2 (yi) +... + πgfig(yi) =
∑g
j=1
```
```
πjfij(yi)
```
- fi 1 (yi),... , fig(yi) are the density functions of Yi (possibly depending on unknown
    parameters θ) in the g components of the mixture
- As in the discussion of the EM algorithm, we define indicators Zij, i= 1,... , N,
    j = 1,... , g:

```
Zij =
```
```




```
```
1 if observation i belongs to component j
```
```
0 otherwise
```

- We then have that
    P(Zij = 1) = πj
- The component probabilities πj are therefore often called prior probabilities, in
    the sense that they express how likely the ith subject is to belong to component j,
    without taking into account the observed response value yi for that observation.
- The posterior probability for observation i to belong to the jth component then
    equals
       πij = P(Zij = 1 | yi)

### =

```
fi(yi | Zij = 1) P(Zij = 1)
fi(yi)
```
```
=
```
```
πj fij(yi)
∑
j
```
```
πjfij(yi)
```

- πij expresses how likely the ith subject is to belong to component j, taking into
    account the observed response value yi for that observation.
- In practice, the posterior probabilities depend on the unknown parameters
    π 1 ,... , πg and θ, but once the mixture model has been fitted, these parameters
    can be replaced by their estimates.
- A natural classification rule now immediately follows:

```
Classify observation i into component j
if and only if
πij = maxk {πik},
```
```
i.e., classify into the component to which observation i is most likely to belong
```

I.9.3 Example: Child Data

- We will now classify the child data using the posterior probabilities corresponding
    to the 4-component Poisson model obtained earlier
- In CAMAN, the procedures ‘mixalg.EM’ and ‘mixalg’ automatically calculate
    posterior probabilities and perform classifications
- The results are saved in the attributes ‘prob’ and ‘classification’:

```
em<-mixalg.EM(obs="counts",weights="frequency",family="poisson",data=child,
t=c(0.5,3,10,15), p=c(0.25,0.25,0.25,0.25), acc=10^(-20))
round(cbind(child,em@prob,em@classification),digits=4)
```

- Result:

```
counts frequency 1 2 3 4 em@classification
0 120 0.8557 0.1439 0.0004 0.0000 1
1 64 0.2310 0.7631 0.0059 0.0000 2
2 69 0.0148 0.9635 0.0216 0.0000 2
3 72 0.0007 0.9382 0.0610 0.0000 2
4 54 0.0000 0.8413 0.1585 0.0002 2
5 35 0.0000 0.6463 0.3529 0.0007 2
6 36 0.0000 0.3863 0.6112 0.0025 3
7 25 0.0000 0.1779 0.8156 0.0066 3
8 25 0.0000 0.0690 0.9165 0.0146 3
9 19 0.0000 0.0246 0.9457 0.0298 3
10 18 0.0000 0.0084 0.9335 0.0581 3
11 18 0.0000 0.0027 0.8878 0.1094 3
12 13 0.0000 0.0009 0.8032 0.1959 3
13 4 0.0000 0.0002 0.6743 0.3254 3
14 3 0.0000 0.0001 0.5115 0.4885 3
15 6 0.0000 0.0000 0.3460 0.6540 4
16 6 0.0000 0.0000 0.2110 0.7890 4
17 5 0.0000 0.0000 0.1190 0.8810 4
18 1 0.0000 0.0000 0.0639 0.9361 4
19 3 0.0000 0.0000 0.0334 0.9666 4
20 1 0.0000 0.0000 0.0171 0.9829 4
21 2 0.0000 0.0000 0.0087 0.9913 4
23 1 0.0000 0.0000 0.0022 0.9978 4
24 2 0.0000 0.0000 0.0011 0.9989 4
```

- In SAS, posterior probabilities and classification results can be saved
    in an output data set:

```
proc fmm data=child;
model y= / dist=poisson k=4 parms(-0.7, 1.1, 2.3, 2.7);
probmodel / parms(0,0,0);
output out=out posterior class;
freq w;
run;
```
- Result:

```
y w Post_1 Post_2 Post_3 Post_4 Class
0 120 0.8557 0.1439 0.0004 0.0000 1
1 64 0.2310 0.7631 0.0059 0.0000 2
```
. .. ...... ...... ...... ...... 2
5 35 0.0000 0.6463 0.3529 0.0007 2
6 36 0.0000 0.3863 0.6113 0.0025 3
. .. ...... ...... ...... ...... 3
14 3 0.0000 0.0001 0.5115 0.4885 3
15 6 0.0000 0.0000 0.3460 0.6539 4
. .. ...... ...... ...... ...... 4
24 2 0.0000 0.0000 0.0011 0.9989 4


- All zero counts are classified in the first component, i.e., in the component which
    can be interpreted as the group of children which are almost always healthy.
- Counts in the range [1,5] are classified in the second component, i.e., in the
    component which can be interpreted as the group of children with normal risk for
    acute respiratory infection.
- Counts in the range [6,14] are classified in the third component, i.e., in the
    component which can be interpreted as the group of children with increased risk
    for acute respiratory infection.
- Counts which are at least 15 are classified in the last component, i.e., in the
    component which can be interpreted as the group of children with high risk for
    acute respiratory infection.


- The proportion of children classified in each component equals:

```
Component # Children Proportion ̂πj
1 120 0.20 0.197
2 294 0.49 0.480
3 161 0.27 0.270
4 27 0.04 0.053
602 1 1
```
- Note that these proportions are very close to the estimated component
    probabilities π̂j.


I.9.4 Example: SIDS Data

- We re-consider the SIDS data, with a Poisson model for the observed rates of
    SIDS in 100 counties in North-Carolina
- The NPMLE obtained before is given by:

```
Yi|p ∼ Poisson(p ni) p ∼
```
```




```
### 0 .0013 0.0021 0.0037 0. 0090

### 0 .33 0.51 0.15 0. 01

```




```
- The fitted mixture can now be used to classify the 100 counties in either one of
    the 4 mixture components


- We obtain the following classification results:

```
Component # Counties Proportion ̂πj
1 24 0.24 0.32
2 64 0.64 0.52
3 11 0.11 0.15
4 1 0.01 0.01
100 1 1
```
- In practice, such classifications are often used to create so-called disease maps,
    i.e., geographical maps in which the different regions are represented according to
    their risk for certain ‘diseases’.


- For the SIDS example, this results in the following map:
- Classification legend for the disease map:

```
Component # Counties π̂j Color Component # Counties ̂πj Color
1 24 0.33 green 3 11 0.15 blue
2 64 0.51 orange 4 1 0.01 black
```

I.9.5 Example: Iris Data

- So far, classification was done based on a fitted mixture model, i.e., the model
    was used to classify observations in detected clusters.
- Mixture models can also be used to classify observations in known groups.
- As an example, we consider Fisher’s Iris data set, and restrict attention to the
    Versicolor and Virginica species:

```
Versicolor Virginica
```

- We focus on the outcome Sepal Length (=Y):
- We have data for 50 flowers of both types:

```
Type Number Mean Stand.Dev. Variance
```
```
Versicolor 50 59.36 5.16 26.64
Virginica 50 65.88 6.36 40.43
```

- Histogram in both samples:

# ➎➏➐➑➒➓➔→➔➐

# ↕ ➣↔↔↕↔↔➙↕➙↔➛↕➛↔➜↕➜↔

# ↕➝↕↔

# ↕➝➞↕

# ↕➝➞↔

# ↕➝➟↕

# ↕➝➟↔

# ↕➝➠↕

# ↕➝➠↔

# ↕➝➣↕

# ➡➢➤

# ➦➥➤➢

# ➧➤➨

# ➩➫➭➯➲➳➫➵➸➺➻

# ➼➽➾➚➽➪➽➶➹

# ➷ ➘➴➴➷➴➴➬➷➬➴➮➷➮➴➱➷➱➴

# ➷✃➷➴

# ➷✃❐➷

# ➷✃❐➴

# ➷✃❒➷

# ➷✃❒➴

# ➷✃❮➷

# ➷✃❮➴

# ➷✃➘➷

# ❰ÏÐ

# ÒÑÐÏ

# ÓÐÔ

# ÕÖ×ØÙÚÖÛÜÝÞ

- The distribution for the sepal length of a randomly selected flower from one of the
    two species is:
       Y ∼π N(59. 36 , 5. 162 ) + (1−π) N(65. 88 , 6. 362 )


- π is the probability that the flower is of the Versicolor type, i.e. the proportion of
    Versicolor flowers in the general population of Versicolor and Virginica flowers.
- The mixture model can be used to classify new flowers, provided an ‘estimate’ for
    π is available.
- Note that classification is then based on a 2-component mixture which is not
    fitted as such to the available data.
- Since it was decided by design to select 50 flowers of each type, π cannot be
    estimated from the data set at hand.
- If the sample would have been a random sample of 100 flowers, which happens to
    contain 50 flowers of each type, π could be set equal to 0.5
- As before, classification is based on posterior probabilities.


- In our example, the ith flower would be classified into the Versicolor group if

```
πi 1 ≥ πi 2 ⇔ πi 1 =
```
```
π fi 1 (yi)
πfi 1 (yi) + (1−π)fi 2 (yi)
```
### ≥ 0. 5

- The posterior probabilities for both groups, as functions of the sepal length,
    assuming equal prior probability for both groups (i.e., π = 0. 5 ) are:


- Flowers with sepal length not larger than 63 mm are classified into the ‘Versicolor’
    group, otherwise they are classified into the ‘Virginica’ group.
- If we use this cut-off value to classify the flowers in our data set, we obtain the
    following classification table:

```
Reality
```
```
Versicolor Virginica
```
```
Versicolor 39 19
```
```
Virginica 11 31
```
```
50 50
```
```
Classification
```

- The above table can be used to estimate the error rate for the classification of
    future flowers:
       P(Flower wrongly classified) = 0. 5 ×P(Flower wrongly classified|of Versicolor type)

```
+ 0. 5 ×P(Flower wrongly classified |of Virginica type)
```
```
= 0. 5 ×^1150 + 0. 5 ×^1950 = 0. 30
```
- Hence, it is to be expected that, using the derived cut-off value, 30% of all flowers
    would be wrongly classified.
- To illustrate that the cut-off value and hence also the error rate highly depends on
    the prior probabilities π and 1 −π, we repeat the calculations assuming that there
    are three times as many Versicolor flowers as Virginica flowers, i.e., for π = 0. 75.


- The posterior probabilities for both groups, as functions of the sepal length, now
    become:
- Flowers with sepal length smaller than 68 mm are classified into the ‘Versicolor’
    group, otherwise they are classified into the ‘Virginica’ group.


- As was to be expected, more flowers will be classified into the ‘Versicolor’ group
- If we use this cut-off value to classify the flowers in our data set, we obtain the
    following classification table:

```
Reality
```
```
Versicolor Virginica
```
```
Versicolor 47 33
```
```
Virginica 3 17
```
```
50 50
```
```
Classification
```
- In comparison to our first analysis, many more Virginica flowers are wrongly
    classified, while the Versicolor flowers are now much better classified.


- The above table can again be used to estimate the error rate for the classification
    of future flowers:
       P(Flower wrongly classified) = 0. 75 ×P(Flower wrongly classified|of Versicolor type)

```
+ 0. 25 ×P(Flower wrongly classified|of Virginica type)
```
```
= 0. 75 × 503 + 0. 25 ×^3350 = 0. 21
```
- Hence, only 21% of all flowers are now expected to be wrongly classified
- Note that the above estimates for the error rates are likely to be over-optimistic as
    they are obtained from testing the discriminant rule with the same observations as
    those used to construct the rule.
- More realistic estimates can be obtained using ‘training’ and ‘test’ data, or using
    cross-validation. This will not be discussed here any further.


I.9.6 Cluster Analysis versus Discriminant Analysis

- Using the Child data and the SIDS data, it has been illustrated how observations
    can be classified in clusters which were detected using NPMLE’s.
- In those analyses, the first step was to check whether there is underlying
    heterogeneity. Afterwards, the observations were classified into the different
    mixture components.
- This is called cluster analysis
- There also exist other approaches to cluster analysis, which are not based on finite
    mixtures


- Often, as was the case for the Iris data, one is interested in finding ‘optimal’ rules
    for classifying observations (possibly future observations) in known groups.
- This is called discriminant analysis
- There also exist other approaches to discriminant analysis, which are not based on
    finite mixtures


Chapter I.10

Model Extensions

```
⊲ Introduction
```
```
⊲ Case study: MMSE data
```

I.10.1 Introduction

- Mixture models can be used to describe latent heterogeneity
- In all examples so far, interest was in describing the distribution of a single
    outcome Y
- Mixture models can be incorporated in statistical models as well, to account for
    heterogeneity not explained by covariates included in the model.
- This will be illustrated in a Binomial regression model, but equally well is
    applicable in other contexts
- Due to the flexibility of the SAS procedure FMM, all analyses will be performed in
    SAS


I.10.2 Case Study: MMSE Data

- We consider data from 58 elderly hip fracture patients, treated at the University
    Hospital Gasthuisberg in Leuven, between September 16, 1996, and February 28,
    1997.
- Of interest is the Mini Mental State Examination (MMSE) score.
- The MMSE is the number of correctly answered questions, out of 30.
- High MMSE values indicate good cognitive functioning, while low MMSE values
    indicate bad cognitive functioning.
- Of interest is the relation between age and MMSE one day after hip surgery


- Descriptive statistics:

```
Outcome Mean Stand.Dev. Minimum Maximum
```
```
Age 78.71 8.20 65 95
MMSE 18.88 8.32 0 30
```
- A natural model is a Binomial logistic model (Model 1):

```
MMSEi ∼ Binomial(30, pi), ln
```
```



```
```
pi
1 −pi
```
```


= β 0 +β 1 Agei
```
- Results:

```
Standard Wald 95% Confidence Wald
Parameter DF Estimate Error Limits Chi-Square Pr > ChiSq
Intercept 1 6.8862 0.5435 5.8210 7.9514 160.54 <.0001
age 1 -0.0801 0.0068 -0.0933 -0.0668 140.29 <.0001
Scale 0 1.0000 0.0000 1.0000 1.0000
```

- Allowing the scale parameter to deviate from one (Pearson χ^2 method) yields:

```
Standard Wald 95% Confidence Wald
Parameter DF Estimate Error Limits Chi-Square Pr > ChiSq
Intercept 1 6.8862 1.5138 3.9192 9.8531 20.69 <.0001
age 1 -0.0801 0.0188 -0.1170 -0.0432 18.08 <.0001
Scale 0 2.7853 0.0000 2.7853 2.7853
NOTE: The scale parameter was estimated by the square root ofPearson’s Chi-Square/DOF.
```
- There is strong evidence for overdispersion.
- Clinicians hypothesize that part of the overdispersion can be explained from the
    fact that some patients are neuro-psychiatric at admission, while others are not.
- Correction for the neuro-psychiatric status is only possible if it was recorded
- Alternatively, mixture models can be used as an attempt to account for this
    heterogeneity in the population studied.


- A 2-component mixture with component-specific regression coefficients (Model 2):

```
MMSEi ∼ π Binomial(30, p 1 i) + (1−π) Binomial(30, p 2 i)
```
```
ln
```
```



```
```
p 1 i
1 −p 1 i
```
```


 = β 10 +β 11 Agei, ln
```
```



```
```
p 2 i
1 −p 2 i
```
```


 = β 20 +β 21 Agei
```
- The model assumes that, at each age, the population consists of two
    sub-populations:
       ⊲ The proportion of each sub-population does not depend on age
       ⊲ The probability to correctly answer an MMSE item is age-specific
       ⊲ The relation between age and correctly answering an MMSE item is different
          for both sub-populations


- SAS code:

```
proc fmm data=test ;
model mmse/n = age / dist=binomial k=2;
run;
```
- Relevant SAS output:

```
Parameter Estimates for ’Binomial’ Model
Standard
Component Effect Estimate Error z Value Pr > |z|
1 Intercept 8.6137 0.7686 11.21 <.0001
1 age -0.09505 0.009294 -10.23 <.0001
2 Intercept 11.4677 1.7088 6.71 <.0001
2 age -0.1624 0.02319 -7.00 <.0001
```
```
Parameter Estimates for Mixing Probabilities
----------------Linked Scale---------------
Standard
Effect Estimate Error z Value Pr > |z| Probability
Intercept 1.2197 0.3294 3.70 0.0002 0.7720
```

- A simplified model is obtained by assuming the age effects to be the same for
    both sub-populations (Model 3):

```
MMSEi ∼ π Binomial(30, p 1 i) + (1−π) Binomial(30, p 2 i)
```
```
ln
```
```

 p 1 i
1 −p 1 i
```
```

 = β
10 +β 1 Agei, ln
```
```

 p 2 i
1 −p 2 i
```
```

 = β
20 +β 1 Agei
```
- The model assumes that, at each age, the population consists of two
    sub-populations:
       ⊲ The proportion of each sub-population does not depend on age
       ⊲ The probability to correctly answer an MMSE item is age-specific
       ⊲ The relation between age and correctly answering an MMSE item is the same
          for both sub-populations


- SAS code:

```
proc fmm data=test ;
model mmse/n = age / dist=binomial k=2;
restrict age 1, age -1;
run;
```
- The RESTRICT statement allows specification of linear equality or inequality
    constraints:
       ⊲ Fixing a parameter at a particular value
       ⊲ Equating parameters in different components in a mixture
       ⊲ Imposing order conditions on parameters
       ⊲ Specifying contrasts among parameters
- The above RESTRICT statement is equivalent to:

```
restrict age 1, age -1 = 0;
```

- Restrictions for effects in specific mixture components are separated by commas.
- Many options possible, see SAS help function
- Relevant SAS output:

```
Parameter Estimates for ’Binomial’ Model
Standard
Component Effect Estimate Error z Value Pr > |z|
1 Intercept 9.1069 0.9879 9.22 <.0001
1 age -0.1007 0.01206 -8.35 <.0001
2 Intercept 6.9031 0.9396 7.35 <.0001
2 age -0.1007 0.01206 -8.35 <.0001
```
```
Parameter Estimates for Mixing Probabilities
----------------Linked Scale---------------
Standard
Effect Estimate Error z Value Pr > |z| Probability
Intercept 1.1527 0.3239 3.56 0.0004 0.7600
```

- Allowing the mixture weights π and 1 −π in Model 2 to depend on age can be
    obtained as follows (Model 4):

```
MMSEi ∼ πi Binomial(30, p 1 i) + (1−πi) Binomial(30, p 2 i)
```
```
ln
```
```



```
```
p 1 i
1 −p 1 i
```
```


 = β 10 +β 11 Agei, ln
```
```



```
```
p 2 i
1 −p 2 i
```
```


 = β 20 +β 21 Agei
```
```
ln
```
```

 πi
1 −πi
```
```

 = α 0 +α 1 Age
i
```
- The model assumes that, at each age, the population consists of two
    sub-populations:
       ⊲ The proportion of each sub-population depends on age
       ⊲ The probability to correctly answer an MMSE item is age-specific
       ⊲ The relation between age and correctly answering an MMSE item is different
          for both sub-populations


- SAS code:

```
data test;
set test;
agec=age-80;
run;
```
```
proc fmm data=test ;
model mmse/n = agec / dist=binomial k=2;
probmodel age / parms(1.2 0);
run;
```
- Good starting values are needed for the parameters in the model for the
    component weights
- In order to be able to use the results from Model 2 as starting values, the age
    covariate was centered at 80 years in the model for πi


- Relevant SAS output:

```
Parameter Estimates for ’Binomial’ Model
Standard
Component Effect Estimate Error z Value Pr > |z|
1 Intercept 8.7175 0.7851 11.10 <.0001
1 age -0.09629 0.009485 -10.15 <.0001
2 Intercept 11.7760 1.7591 6.69 <.0001
2 age -0.1663 0.02362 -7.04 <.0001
```
```
Parameter Estimates for Mixing Probabilities
Standard
Effect Estimate Error z Value Pr > |z|
Intercept 1.2644 0.3397 3.72 0.0002
agec 0.03735 0.04282 0.87 0.3831
```

- Summary of results:

```
Effect Parameter Model 1 Model 2 Model 3 Model 4
Component 1: Intercept β 10 6.886 (0.544) 8.614 (0.769) 9.107 (0.988) 8.718 (0.785)
Age β 11 -0.080 (0.007) -0.095 (0.009) -0.101 (0.012) -0.096 (0.009)
Component 2: Intercept β 20 11.468 (1.709) 6.903 (0.940) 11.776 (1.759)
Age β 21 -0.162 (0.023) -0.101 (0.012) -0.166 (0.024)
Weight 1: Probability π 0.772 0.760
Intercept α 0 1.220 (0.329) 1.153 (0.324) 1.264 (0.340)
Age α 1 0.037 (0.043)
Deviance: − 2 ℓℓ 662.0 412.9 421.3 412.1
```
- Of all models fitted, Model 2 is best supported by the data


- Interpretation:

```
⊲ Two sub-populations representing 77% and 23% of the population, respectively
⊲ The proportion of each sub-population does not change with age
⊲ The probability of correctly answering an MMSE item decreases with age
⊲ This decrease is significantly steeper in the second sub-population than in the
first
```
- Fitted probabilities to correctly answer an MMSE item:

```
ln
```
```



```
```
p 1 i
1 −p 1 i
```
```


 = 8.^614 −^0.^095 Agei (Population 1)
```
```
ln
```
```



```
```
p 2 i
1 −p 2 i
```
```


 = 11.^468 −^0.^162 Agei (Population 2)
```

- Graphical representation:

# àáß

# á

# âã

# äåæ

# çè

# êé

# ëìë

# ëìí

# ëìî

# ëìï

# ëìð

# ñìë

# òóôõöô÷øùú

# ïûüûðûýû

# þÿ✁✂✄☎✂✁✆✝✞✟✠ ñ ✡üü ☛☞☎✂✁✆✝✞✟✠ í ✡í ✌☛☞

- At any age, subjects in population 1 are more likely to correctly answer MMSE
    items than subjects in population 2


- An informal check whether the sub-populations detected truly correspond to
    patients who are (not) neuro-psychiatric, observed proportions are calculated for
    both groups, and for specific age intervals:

```
Not neuro-psychiatric Neuro-psychiatric
Age range Average MMSE P(MMSE item= 1) Average MMSE P(MMSE item= 1)
[65,75] 25.88 25.88/30=0.86 17.00 17.00/30=0.57
]65,85] 21.58 21.58/30=0.72 13.40 13.40/30=0.45
]85,95] 14.91 14.91/30=0.50 4.33 4.33/30=0.14
```
- These observed proportions can now be graphically compared to the fitted
    probabilities to correctly answer an MMSE item


- Result:

# ✎ ✏✍

# ✏

# ✑✒

# ✓ ✔ ✕

# ✖ ✗

# ✙✘

# ✚ ✛ ✚

# ✚ ✛ ✜

# ✚ ✛ ✢

# ✚ ✛ ✣

# ✚ ✛ ✤

# ✥ ✛ ✚

# ✦ ✧ ★ ✩ ✪ ★ ✫ ✬ ✭ ✮

# ✣ ✯ ✰ ✯ ✤ ✯ ✱ ✯

# ✲ ✳ ✴ ✵ ✶ ✷ ✸ ✴ ✶ ✵ ✹ ✺ ✻ ✼ ✴ ✽ ✥ ✾ ✰ ✰ ✿ ❀ ✸ ✴ ✶ ✵ ✹ ✺ ✻ ✼ ✴ ✽ ✜ ✾ ✜ ❁ ✿ ❀

# ❂ ✴ ✽ ❃ ✽ ❄ ✵ ✳ ✴ ❂ ❄ ✵ ✳ ✴

- Success rates for patients who are not neuro-psychiatric are well described by the
    first mixture component


- Success rates for neuro-psychiatric patients are less well described by the second
    mixture component
- This is also observed when cross-classifying the neuro-status with the classification
    based on posterior probabilities:

```
Yes No
```
```
Component 1 9 36 45
```
```
Component 2 9 4 13
```
```
18 40
```
```
Neuro-psychiatric?
```
```
Classification
```

- Also, when the model is corrected for neuro-status, there still is evidence for the
    presence of two mixture components:

```
proc fmm data=test; proc fmm data=test;
model mmse/n = age neuro neuro*age model mmse/n = age neuro neuro*age
/dist=binomial k=1; /dist=binomial k=2;
run; run;
```
- Relevant SAS output for 1-component mixture:

```
Fit Statistics
-2 Log Likelihood 542.4
Parameter Estimates for ’Binomial’ Model
Standard
Effect Estimate Error z Value Pr > |z|
Intercept 7.4909 0.6669 11.23 <.0001
age -0.08238 0.008214 -10.03 <.0001
NEURO -1.5639 1.3462 -1.16 0.2453
age*NEURO 0.004425 0.01682 0.26 0.7925
```

- Relevant SAS output for 2-component mixture:

```
Fit Statistics
-2 Log Likelihood 373.3
Parameter Estimates for ’Binomial’ Model
Standard
Component Effect Estimate Error z Value Pr > |z|
1 Intercept 6.1356 1.0411 5.89 <.0001
1 age -0.05973 0.01303 -4.58 <.0001
1 NEURO -0.3562 1.8429 -0.19 0.8467
1 age*NEURO -0.00833 0.02303 -0.36 0.7176
2 Intercept 5.9603 1.7478 3.41 0.0006
2 age -0.07600 0.01978 -3.84 0.0001
2 NEURO 4.5602 3.2644 1.40 0.1624
2 age*NEURO -0.08705 0.04215 -2.07 0.0389
```
```
Parameter Estimates for Mixing Probabilities
----------------Linked Scale---------------
Standard
Effect Estimate Error z Value Pr > |z| Probability
Intercept 0.9591 0.4082 2.35 0.0188 0.7229
```

- Conclusion:

```
The mixture components do not coincide with the neuro-psychiatric
and non-neuro-psychiatric subpopulations
```
- Alternative conclusion:

```
The neuro-psychiatric status does not entirely explain
the presence of mixture components
```

Part II

Non-linear Models

Advanced Modeling Techniques II.1


Chapter II.1

Non-Linear Mixed Models

```
⊲ From linear to non-linear models
```
```
⊲ Orange Tree Example
```
```
⊲ Song Bird Example
```
Advanced Modeling Techniques II.2


