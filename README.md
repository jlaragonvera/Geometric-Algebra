# Geometric Algebra

CGAlgebra.m is a Mathematica package for the 5D Conformal Geometric Algebra.

This package contains declarations for calculations with Conformal
Geometric Algebra. Basis vectors {e_0, e_1, e_2, e_3, e_\infty} are 
denoted by e[0], e[1], e[2], e[3], e[\[Infinity]]. Geometric products
of basis elements are denoted as e[0,1,2] (=e_0 e_1 e_2), etc.
   
The results of any calculation is given in terms of the geometric product
of basis elements, that is, the outer (Grassman) product of basis
elements or multivectors is calculated by using OuterProduct[] and the
output is given in terms of geometric product of basis vectors.
        
Examples:     
                       
    The vector e_0 + 2 e_1 - a e_3 + e_\infty is written as
              
         e[0] + 2 e[1] - a e[3] + e[\[Infinity]]
              
    The multivector a + 5 e_1 + e_1e_2e_3 is written as
              
          a + 5 e[1] + e[1,2,3].

A tutorial can be downloaded from:
https://arxiv.org/abs/1711.02513
