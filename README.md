### Gaussian Elimination Parallel Program

Here's represented parallel program of Gaussian Elimination for solving systems of linear equations.
You can find more information about this method [here](https://en.wikipedia.org/wiki/Gaussian_elimination)

Program was tested both on PC and on cluster.

ImplementationThe algorithm has been implemented in C++, using the Message Passing Interface (MPI) API. Most of the code is a  straightforward  implementation  of  the  following  pseudocode 
~~~Cpp
All nodes:
  Read a subset of input data, [A|b], size N
Forward elimination: 
  For each row r[j], j in [1, N]:
  // Find pivot row
  rp = find_pivot(j)
  swap_rows(rp, r[j])
  
  // eliminate column j
  Node owning r[j]:
    Broadcast r[j]:
  All nodes:
  Receive r[j]
  For each local row r[i], i>j:
    Recompute row r[i]using r[j] (elem[i, j] will be 0 after this)
Back substitution:
  For each row r[j], j in [N, 1]:
    Node owning r[j]:
    Compute x[j] using partial solution p[j]
      Broadcast x[j]
    All nodes:
      Receive x[j]
      Compute partial solution p[j] using x[j]...x[N]
End:
  Root node:
    Receive x elements computed by each node, assemble solution vector x.
~~~

 mapping  communication  parts  to  the  MPI  Scatterv, Gatterv,  Broadcast  and  Recv  primitives.
