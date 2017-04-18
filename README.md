rook_decomp
===========

Software to find decompositions of 2 K_{n^2} into copies of the rook graph R(n,n)


The process is explained for the example where n=6.
K_36 has 36 vertices, and 630 edges, while R(6,6) has 36 vertices and 180 edges.
2 * 630 / 180 = 7, so we have to find 7 distinct copies of the rook graph in 2
K_36.

Label the vertices of K_36 infinity, (1,1), (1,2), ..., (1,7), (2,1), ..., (5,7).
Place infinity to one side, and group the remaining vertices into 5 groups such that group X contains vertices of the form (X,Y) for Y in {1,...,7}.

Under the action of Z_7 on each group (the same action applied simultaneously to all groups), we can place edges into edge-orbits. For example, the edge between (1,1) and (1,2) is in the same orbit as the edge between (1,2) and (1,3).

If we can choose two edges from each distinct orbit such that this collection of edges forms the rook graph R(n,n), then under the action of Z_7 we have all required copies of the rook graph, and thus the required decomposition.

This piece of software attempts to find these edges. It does so by trying to place the vertices of K_36 into R(n,n) such that only 2 edges from each orbit are used at any one time.
