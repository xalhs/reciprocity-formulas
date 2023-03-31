# reciprocity-formulas

This is an algorithm I made accompanying my Master's thesis. It is used to calculate topological invariants of manifolds given some parameters.

The main point of these algorithms is to calculate and verify a reciprocity formula, the formula looks like this:

$$
\frac{1}{|\det(H)|^{n-l/2}}\sum_{x \in (\mathbb{Z}^r/H(\mathbb{Z}^r))^n}\exp\left(\pi i x^T (L\otimes H^{-1}) x\right)
=\frac{1}{|\det(A)|^{m - r/2}}
\exp\left(\frac{i \pi}{4}\sigma(K)\sigma(L)\right)
\sum_{y \in (\mathbb{Z}^l/A(\mathbb{Z}^l))^m}
\exp\left(-\pi i y^T (K\otimes A^{-1}) y\right)  \qquad (1)
$$ 

Where $L$ and $K$ are bilinear forms, i.e. symmetric matrices (also $K$ is even which means its diagonal elements are even) and $A$ and $H$ are the "invertible parts"
of $L$ and $K$. Also $n$ is the dimension of $L$, $m$ the dimension of $K$, $l$ is the rank of $L$ and $r$ is the rank of $K$.

however, its use can be much simpler than that:

At the base level it can compute gauss sums like: $$\sum_{n=0}^{p-1} \mathrm{exp}(i \pi q n^2 /p)$$

It can also do the following operations:
- Isolate the invertible part of a symmetric matrix 
- Calculate the cokernel of a symmetric matrix on the field of integers (i.e. the quotient group $\mathbb{Z}^n/L\mathbb{Z}^n$)
- Split the above group into a direct sum of cyclic groups
- Calculate a linking form that corresponds to a given (linking) matrix
- Calculate the integer linking matrix of a lens space L(p,q) (using Saveliev's method)
- Perform the second Kirby move on symmetric (linking) matrices

# Theory


## Homology

very briefly, we can create a manifold by acting on a set of curves inside $S^3$. These curves are partially described by what is called the linking matrix, which is a symmetric matrix. Now there is a way of obtaining the homology if the manifold by looking at the linking matrix. Suppose that $L$ is an $n\times n$ linking matrix then the first homology group can be written as:

$$ 
H_1 = \langle G_1 ,..., G_n | L\cdot \mathbb{G} = 0 \rangle,
$$

where $\mathbb{G} = (G_1 , ... G_n)^T$ and $G_1,...,G_n$ are each a generator of $\mathbb{Z}$. In other words the first homology group will be given by the cokernel: $\mathbb{Z}^n/L\mathbb{Z}^n$

One thing this algorithm does is calculate the homology by calculating the above cokernel. It does that in two steps, first it isolates the invertible part from the non invertible part of the matrix using only the allowed moves (Kirby moves), bringing it into a form $A\oplus 0_{n-l}$ where $A$ is a symmetric $l \times l$ invertible matrix. Then it calculates the cokernel of $A$: $\mathbb{Z}^l/A\mathbb{Z}^l$. 

more will be added on the above explanation

## Reciprocity formulas

The following is a topological invariant of our manifold

$$
\sum_{y \in (\mathbb{Z}^l/A(\mathbb{Z}^l))^m}
\exp\left(-\pi i y^T (K\otimes A^{-1}) y\right)
$$

where $A$ as we've seen before is the invertible part of the linking matrix and $K$ is a symmetrized coupling constant matrix. The main result of my Master's thesis is that we can relate the above invariant, through formula (1) to another Reshetikhin-Turaev-like invariant shown on the left hand side of the equation. The main function of this algorithm takes as input the linking matrix $L$ and the (even) symmetrized coupling constant matrix $K$ and outputs the left and right hand side of equation (1) in order to confirm it. Of course the formula is proven rigorously by Deloup and Turaev but it is nice to have numerical evidence.

More information will be added soon.
References to all of these are in my Master's Thesis and possibly a future paper.
