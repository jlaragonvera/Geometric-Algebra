# Geometric Algebra

*Mathematica* packages for Clifford (geometric) algebra calculations:

## CGAlgebra 

**CGAlgebra.m** is a *Mathematica* package for the 5D Conformal Geometric Algebra.

This package contains declarations for calculations with Conformal
Geometric Algebra. Basis vectors {e<sub>0</sub>, e<sub>1</sub>, e<sub>2</sub>, e<sub>3</sub>, e<sub>Infinity</sub>} are 
denoted by e[0], e[1], e[2], e[3], e[Infinity]. Geometric products
of basis elements are denoted as e\[0,1,2\] (= e<sub>0</sub> e<sub>1</sub> e<sub>2</sub>), etc.
   
The results of any calculation is given in terms of the geometric product
of basis elements, that is, the outer (Grassman) product of basis
elements or multivectors is calculated by using OuterProduct[] and the
output is given in terms of geometric product of basis vectors.
        
Examples:

The vector e<sub>0</sub> + 2 e<sub>1</sub> + e<sub>Infinity</sub> is written as:
              
         A = e[0] + 2 e[1] + e[\[Infinity]];
         
The multivector a + 5 e<sub>1</sub> +  e<sub>1</sub>e<sub>2</sub>e<sub>3</sub> is
              
         B = a + 5 e[1] + e[1,2,3];

The geometric product

        GeometricProduct[A,B]

yields:

      a e[0] + 2 a e[1] + a e[\[Infinity]] + 5 e[0,1] - 5 e[1,\[Infinity]] + e[0,1,2,3] - e[1,2,3,\[Infinity]]

The inner product (left contraction)

       InnerProduct[A,B]

yields

       10 + 2 e[2,3]

A tutorial can be downloaded from:
<https://arxiv.org/abs/1711.02513>


## Clifford

**clifford.m** is the most recent version of the package by G. Aragon-Camarasa, G. Aragon-Gonzalez, J.L. Aragon and M.A. Rodriguez-Andrade. A user guide (**CliffordUserGuide**) is available, as well as *Mathematica* palette (**CliffordPalette**). The fundamentals of the package are presented in:

<https://arxiv.org/abs/0810.2412>


## CliffordBasic

**CliffordBasic.m** is a completely renewed but reduced version of **clifford.m** package by G. Aragon-Camarasa, G. Aragon-Gonzalez, J.L. Aragon and M.A. Rodriguez-Andrade:
<https://arxiv.org/abs/0810.2412>

Using rule-base programming the algebra over R<sup>p,q</sup> in arbitrary dimensions is constructed as in A. Macdonald "An Elementary Construction of Geometric Algebra", Adv. Appl. Cliff. Alg. 12 (2002) 1-6.

In **CliffordBasic**, the j-th basis vector is denoted by e[j] and the geometric product of basis vectors, such as e<sub>1</sub>e<sub>3</sub>e<sub>4</sub>, as `e[1,3,4]`.

Examples:     

The vector e<sub>1</sub> + 2 e<sub>2</sub> - a e<sub>3</sub> is written as   
      
        A = e[2] + 2 e[2] - a e[3];
                  
The The multivector a + 5 e<sub>1</sub> + e<sub>1</sub> e<sub>2</sub> e<sub>3</sub> is written as           
              
         B = a + 5 e[1] + e[1,2,3];
         
The geometric product AB is calculated as 

        GeometricProduct[A,B]

yielding:

        3 a e[2] - a^2 e[3] - 15 e[1,2] - a e[1,2] - 3 e[1,3] + 5 a e[1,3]
    
which can be factored using

         GFactor[%]
         
         3 a e[2] - a^2 e[3] + (-15-a) e[1,2] + (-3 + 5 a)e[1,3]
         
The signature of the R<sup>p,q</sup> is set by `$SetSignature={p,q}`. If not specified,
the default value is:

         $SetSignature={20,0}        

