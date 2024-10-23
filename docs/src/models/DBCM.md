# DCBM
## Model description
The Directed Binary Configuration Model (DBCM) is a maximum-entropy null model for directed networks. It is based on the idea of fixing the in- and outdegree sequence of the network, i.e., the number of incoming and outgoing edges incident to each node, and then randomly rewiring the edges while preserving the degree sequence. The model assumes that the edges are unweighted and that the network is simple, i.e., it has no self-loops or multiple edges between the same pair of nodes [[1](#1),[2](#2)]. 

We define the parameter vector as ``\theta = [\alpha ; \beta]``, where ``\alpha`` and ``\beta`` denote the parameters associated with the out- and indegree respectively.

| Description                   | Formula |
| --------------------------    | :-------------------------------------------------------------------------------- |
| Constraints                   | `` \forall i: \begin{cases} k_{i, out}(A^{*}) = \sum_{j=1}^{N} a^{*}_{ij} \\ k_{i, in}(A^{*}) = \sum_{j=1}^{N} a^{*}_{ji} \end{cases} ``|
| Hamiltonian                   | `` H(A, \Theta) = H(A, \alpha, \beta) = \sum_{i=1}^{N} \alpha_i k_{i,out}(A) +  \beta_i k_{i,in}(A)`` |
| Factorized graph probability  | `` P(A \| \Theta) = \prod_{i=1}^{N}\prod_{j=1, j \ne i}^{N} p_{ij}^{a_{ij}} (1 - p_{ij})^{1-a_{ij}}``  |
| $\langle a_{ij} \rangle$      | `` p_{ij} = \frac{e^{-\alpha_i - \beta_j}}{1+e^{-\alpha_i - \beta_j}}`` |
| Log-likelihood                | `` \mathcal{L}(\Theta) = -\sum_{i=1}^{N} \left[ \alpha_i k_{i,out}(A^{*}) +  \beta_i k_{i,in}(A^{*}) \right] - \sum_{i=1}^{N} \sum_{j=1, j\ne i}^{N} \ln \left( 1+e^{-\alpha_i - \beta_j} \right) ``|
| $\langle a_{ij}^{2} \rangle$  | `` \langle a_{ij} \rangle`` |
| $\langle a_{ij}a_{ts} \rangle$| `` \langle a_{ij} \rangle \langle a_{ts} \rangle`` |
| $\sigma^{*}(X)$               | ``\sqrt{\sum_{i,j} \left( \sigma^{*}[a_{ij}] \frac{\partial X}{\partial a_{ij}}  \right)^{2}_{A = \langle A^{*} \rangle} + \dots }`` |
| $\sigma^{*}[a_{ij}]$          | ``\frac{\sqrt{e^{-\alpha_i - \beta_j}}}{1+e^{-\alpha_i - \beta_j}} ``   |





## Creation
```julia
using Graphs
using MaxEntropyGraphs

# define the network
G = SimpleDiGraph(rhesus_macaques())

# instantiate a UBCM model
model = DBCM(G)
```

## Obtaining the parameters
```julia
# solve using default settings
solve_model!(model)
```

## Sampling the ensemble
```julia
# generate 10 random instance of the ensemble
rand(model, 10)
```

## Model comparison
```julia
# compute the AIC  
AIC(model)
```

## Counting network motifs
The count of a specific network motif can be computed by using `M{motif_number}`. The motif numbers match the patterns shown on the image below. So for example if you want to compute the number of reciprocated triangles, you would use `M13(model)`.

```@example
HTML("""<object type="image/png" data=$(joinpath(Main.buildpath, "../assets/directed_motifs.png"))  alt="Motif image not found" style="width: 100%;" ></object> """) # hide
```
[source](https://snap-stanford.github.io/cs224w-notes/preliminaries/motifs-and-structral-roles_lecture)


```julia
# Compute the number of occurences of M13
M13(model)
```

_References_

```@raw html
<ul>
<li>
<a id="1">[1]</a> 
Squartini, Tiziano and Garlaschelli, Diego. <!--  author(s) --> 
<em>"Maximum-Entropy Networks: Pattern Detection, Network Reconstruction and Graph Combinatorics"</em> <!--  title --> 
Springer-Verlag GmbH; 1st ed. 2017 edition (25 Dec. 2017). <!--  publisher(s) --> 
<a href="https://link.springer.com/book/10.1007/978-3-319-69438-2">https://link.springer.com/book/10.1007/978-3-319-69438-2</a>
</li>
<li>
<a id="2">[2]</a> 
Squartini, Tiziano and Garlaschelli, Diego. <!--  author(s) --> 
<em>"Analytical maximum-likelihood method to detect patterns in real networks"</em> <!--  title --> 
2011 New J. Phys. 13 083001. <!--  publisher(s) --> 
<a href="https://iopscience.iop.org/article/10.1088/1367-2630/13/8/083001">https://iopscience.iop.org/article/10.1088/1367-2630/13/8/083001</a>
</li>
</ul>
```

