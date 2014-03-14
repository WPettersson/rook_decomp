rook_decomp
===========

Software to find decompositions of 2 K_n,n into copies of the rook graph R_n


The process is explained for the example whre n=6.
K_n,n has 36 vertices.  Label these infinity, (1,1), (1,2), ..., (1,7), (2,1), ..., (5,7).

Place infinity to one side, and group the remaining vertices into 5 groups such that group X contains vertices of the form (X,Y) for Y in {1,...,7}.

Under the action of Z_7 on each group (the same action applied to all groups), we can place edges into edge-orbits. For example, the edge between (1,1) and (1,2) is in the same orbit as the edge between (1,2) and (1,3).

If we can choose two edges from each distinct orbit such that this collection of edges forms the rook graph R_n, then under the action of Z_7 we have all required copies of the rook graph, and thus the required decomposition.

This piece of software attempts to find these edges. It does so by trying to place the vertices of K_n,n into R_n such that only 2 edges from each orbit are used at any one time.
