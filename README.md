# Geometric Algebra

CGAlgebra.m 

is a Mathematica package for the 5D Conformal Geometric Algebra.

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


CliffordBasic.m

is a completely renewed but reduced version of Clifford.m package by 
G. Aragon-Camarasa, G. Aragon-Gonzalez, J.L. Aragon and 
M.A. Rodriguez-Andrade:
https://arxiv.org/abs/0810.2412

Using rule-base programming the algebra over R^{p,q} in arbitrary
dimensions is constructed as in A. Macdonald "An Elementary 
Construction of Geometric Algebra", Adv. Appl. Cliff. Alg. 12 (2002) 1-6.

In CliffordBasic, the j-th basis vector is denoted by e[j] and the 
geometric product of basis vectors, such as e1e3e4, as e[1,3,4].

Examples:     
                       
    The vector e_1 + 2 e_2 - a e_3 is written as
              
         A = e[2] + 2 e[2] - a e[3]
              
    The multivector a + 5 e_1 + e_1e_2e_3 is written as
              
         B = a + 5 e[1] + e[1,2,3].
     
    The geometric product AB is calculated as GeometricProduct[A,B]
    and yields:
    
         3 a e[2] - a^2 e[3] - 15 e[1,2] - a e[1,2] - 3 e[1,3] + 5 a e[1,3]
    
    which can be factored using GFactor[%]:
         
         3 a e[2] - a^2 e[3] + (-15-a) e[1,2] + (-3 + 5 a)e[1,3]
         
The signature of the R^{p,q} is set by $SetSignature={p,q}. If not specified,
the default value is:

         $SetSignature={20,0}        

