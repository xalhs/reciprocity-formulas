# reciprocity-formulas

This is an algorithm I made accompanying my Master's thesis. It is used to calculate topological invariants of manifolds given some parameters.

The main point of these algorithms is to calculate and verify a reciprocity formula, the formula looks like this:

$$
\frac{1}{|\det(H)|^{n-l/2}}\sum_{x \in (\mathbb{Z}^r/H(\mathbb{Z}^r))^n}\exp\left(\pi i x^T (L\otimes H^{-1}) x\right)
=\frac{1}{|\det(A)|^{m - r/2}}
\exp\left(\frac{i \pi}{4}\sigma(K)\sigma(L)\right)
\sum_{y \in (\mathbb{Z}^l/A(\mathbb{Z}^l))^m}
\exp\left(-\pi i y^T (K\otimes A^{-1}) y\right)  
$$

Where $L$ and $K$ are bilinear forms, i.e. symmetric matrices (also $K$ is even which means its diagonal elements are even) and $A$ and $H$ are the "invertible parts"
of $L$ and $K$. Also $n$ is the dimension of $L$, $m$ the dimension of $K$, $l$ is the rank of $L$ and $r$ is the rank of $K$.

however, its use can be much simpler than that:

At the base level it can compute gauss sums like: $$\sum_{n=0}^{p-1} \mathrm{exp}(i \pi q n^2 /p)$$

It can also do the following operations:
- Isolate the invertible part of a symmetric matrix 
- Calculate the cokernel of a symmetric matrix (for invertible matrices it's the quotient group $\mathbb{Z}^n/L\mathbb{Z}^n$
- Split the above group into a direct sum of cyclic groups
- Calculate a linking form that corresponds to a given (linking) matrix

More information will be added soon.
