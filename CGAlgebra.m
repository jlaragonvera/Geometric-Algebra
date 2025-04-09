(* ::Package:: *)

(* Set up the Package Context. *)

(* :Title: CGAlgebra.m -- Conformal Geometric Algebra *)

(* :Authors: Alejandra Ortiz Duran and Jose L. Aragon *)

(* :Summary:
   This package contains declarations for calculations with Conformal
   Geometric Algebra. Basis vectors {e_0, e_1, e_2, e_3, e_\infty} are 
   denoted by e[0], e[1], e[2], e[3], e[\[Infinity]]. Geometric products
   of basis elements are denoted as e[0,1,2] (=e_0 e_1 e_2), etc.
   
   The results of any calculation is given in terms of the geometric product
   of basis elements, that is, the outer (Grassman) product of basis
   elements or multivectors is calculated by using OuterProduct[] and the
   output is given in terms of geometric product of basis vectors.
        
   Examples:     The vector e_0 + 2 e_1 - a e_3 + e_\infty is written as
                  e[0] + 2 e[1] - a e[3] + e[\[Infinity]].
                     
                 The multivector a + 5 e_1 + e_1e_2e_3 is written as
                  a + 5 e[1] + e[1,2,3].    
 *)
(* :Copyright: (c) 2017 by Jose L. Aragon *)

(* :Package Version: 2.0 *)

(* :Mathematica Version: 11.0 *)

(* :History:
   First version: 1.0 October 2017 by Jose L. Aragon
   
   Revision (October, 2021) J.L. Aragon and A. Ortiz Duran
        - GradeQ[] fixed.
        - Rotation[] default argument didn't work; fixed.
  Major revision (May, 2022) Alejandra Ortiz Duran and Jose L. Aragon
        - Several issues were found when working with the basis {e_0, e_1, e_2, e_3, e_\infty}
          by applying the general formulas for the outer and inner product:
          A_p^B_q = <Ap Bq>_{p+q}   and    A_p.B_q = <Ap Bq>_{q-p}.
          To fix this, the package was refactored and now works in R^{4,1}, with the basis
          {e1,e2,e3,e4,e5}, and using the changes of variables e4 = e0 - einf/2, e5 = e0 + einf/2
          and e0 = (e4 + e5)/2, einf = e5 - e4. 
*)

(* :Keywords: Clifford algebra, geometric algebra, conformal model *)

(* :Sources:
   
	1. K. Kanatani, 2015. Understanding Geometric Algebra. 
	   CRC Press
	2. L. Dorst, D. Fontijne, S. Mann, 2009. Geometric Algebra for Computer Science. 
	   Morgan Kaufmann
	3. A. Macdonald, 2011. Linear and Geometric Algebra.
	   CreateSpace Independent Publishing Platform
*)

(* :Examples:
   See CGAlgebraTutorial.pds for a general tutorial and examples.
*)

BeginPackage["CGAlgebra`"]

(* usage messages for the exported functions and the context itself *)
CGAlgebra::usage = "CGAlgebra.m is a package to resolve operations with Clifford Algebra in the Conformal Model."

e::usage = "e is used to denote the elements of the canonical basis of Euclidean vector
space where the Clifford Algebra is defined, thus e[i] is used as i-th basis
element (i=0,1,2,3,\[Infinity])"

Grade::usage = "Grade[m,r] returns the r-vector part of the multivector m."

GradeQ::usage = "GradeQ[m,r] returns True if m is a r-vector."

MultiplicationTable::usage = "MultiplicationTable[e] returns the multiplication table of symbol e."

GFactor::usage = "GFactor[m] groups terms with commom e[__]."

GeometricProduct::usage = "GeometricProduct[m1,m2,...] returns the Geometric Product of multivectors m1,m2,..."

OuterProduct::usage = "OuterProduct[m1,m2,...] returns the Outer (Grassmann) Product of multivectors m1,m2,..."

InnerProduct::usage = "InnerProduct[m1,m2] returns the Inner Product (Left Contraction) of multivectors m1 and m2"

Reversion::usage = "Reversion[m] returns the reversion of the multivector m"

Involution::usage = "Involution[m] returns the Grade Involution (\[Dagger]) of the multivector m"

Magnitude::usage = "Magnitude[m] returns the Magnitude (m^2) of the multivector m"

MultivectorInverse::usage = "MultivectorInverse[m] returns the Inverse (if it exists) of the multivector m"

I5::usage = "I5 is the pseudoscalar or G^{4,1}"

I5i::usage = "I5i is the inverse of the pseudoscalar I5"

Dual::usage = "Dual[m] returns the Dual of the multivector m"

ToVector::usage = "ToVector[v] transform elements of the conformal space of the kind v=xe[1]+ye[2]+ze[e] 
to the standard Mathematica notation {x,y,x}"

Rotation::usage = "Rotation[x,a,b,theta] Rotates the vector x by an angle theta (in radians), 
along the plane defined by a and b. The sense of the rotation is from a to b.
If no theta is given, the default value is the  angle between a and b."


Begin["`Private`"] (* Begin Private Context *) 

(* Returns True if either 0 or \[Infinity] is duplicated *)
DuplicateQ0Inf[i_] := Count[i, 0] >= 2 || Count[i, \[Infinity]] >= 2

(* Returns True if either 0 or \[Infinity] is duplicated and appears consecutively*)
Contiguous0Inf[i_] := SequenceCases[i, {0, 0}] != {} || SequenceCases[i, {\[Infinity], \[Infinity]}] != {}


(* Definition of the exported functions *)

(* GradeQ returns True if x contains only terms of grade r *)
GradeQ[a_, r_?NumberQ] := If[r === 0, True, False]                 /; FreeQ[a, e[__]]
GradeQ[e[i__], r_?NumberQ] := If[Length[{i}] === r, True , False]
GradeQ[(a_: 1) e[i__], r_?NumberQ] := GradeQ[e[i], r]               /; FreeQ[a, e[_]]
GradeQ[x_Plus, r_?NumberQ] := And @@ (GradeQ[#, r] & /@ Apply[List, x])

(* MultiplicationTable returns the multiplication table for Clifford products *)
MultiplicationTable[e_] := Grid[Transpose[Insert[Transpose[Prepend[Table[e[i, j], {i,{0,1,2,3,\[Infinity]}}, {j, {0,1,2,3,\[Infinity]}}], 
      Table[e[i], {i,{0,1,2,3,\[Infinity]}}]]], {"",e[0],e[1],e[2],e[3],e[\[Infinity]]}, 1]], Frame -> All]

(* GFactor groups terms with commom e[__] *)
GFactor[x_] := Collect[Expand[x], e[__]]

(* ToVector[v_] transform elements of the conformal space of the kind v=xe[1]+ye[2]+ze[e]
   to the standard Mathematica notation {x,y,x} *)
ToVector[v_]:= Table[Coefficient[v, e[i]], {i, 3}] /; GradeQ[v,1] && Coefficient[v,e[0]]==0 && Coefficient[v,e[\[Infinity]]]==0

(* Changes of variables *)

(* From (e[4], e[5]) to (e[0], e[Infinity]) *)
From45To0I[a_] := a                                    /; FreeQ[a, e2[__]] 
From45To0I[x_ + y_] := From45To0I[x] + From45To0I[y]
From45To0I[a_ e2[i__]] := a From45To0I[e2[i]]          /; FreeQ[a, e2[__]]
From45To0I[e2[i__]] := (Replace[e2[i], {e2[4] -> e[0] - e[\[Infinity]]/2, e2[5] -> e[0] + e[\[Infinity]]/2,e2[1] -> e[1],e2[2] -> e[2],e2[3] -> e[3] }] //Expand)/; Length[{i}] == 1
From45To0I[e2[i__]] := GeometricProduct @@ ((e2 /@ List @@ e2[i]) /. {e2[4] -> e[0] - e[\[Infinity]]/2, e2[5] -> e[0] + e[\[Infinity]]/2 ,e2[1] -> e[1],e2[2] -> e[2],e2[3] -> e[3]})  // Expand

(* From (e[0], e[Infinity]) to (e[4], e[5]) *)
FromOITo45[a_] := a                                    /; FreeQ[a, e[__]] 
FromOITo45[x_ + y_] := FromOITo45[x] + FromOITo45[y]
FromOITo45[a_ e[i__]] := a FromOITo45[e[i]]            /; FreeQ[a, e[__]]
FromOITo45[e[i__]] := (Replace[e[i], {e[0] -> (e2[4] + e2[5])/2, e[\[Infinity]] -> e2[5] - e2[4],e[1] -> e2[1],e[2] -> e2[2],e[3] -> e2[3]}]//Expand) /; Length[{i}] == 1
FromOITo45[e[i__]] := (GeometricProduct2@@ ((e /@ List @@ e[i]) /. {e[0] -> (e2[4] + e2[5])/2, e[\[Infinity]] -> e2[5] - e2[4] ,e[1] -> e2[1],e[2] -> e2[2],e[3] -> e2[3]})//Expand)  /; Length[{i}]>1

(*End of changes of variable*)

(*Begin of clifford's definitions*)

(*GFactor groups terms with commom e[__]*)
GFactor2[x_] := Collect[Expand[x],e2[__]]

(*Pseudoscalar function*)
Pseudoscalar[n_/;Element[n,Integers]&&n>0] := e2@@(Range[n])

(*HomogeneousQ function*)
HomogeneousQ[x_,r_/;Element[r,Integers]&&NonNegative[r]] := Simplify[Expand[x]===Grade2[x,r]]

(* Begin of the Clifford algebra of R^{p,q}; the basis vectors are denoted as e2[], and the functions are added the digit 2  *)

(* GFactor groups terms with commom e[__] *)
GFactor2[x_] := Collect[Expand[x],e2[__]]

(* Pseudoscalar function *)
Pseudoscalar[n_/;Element[n,Integers]&&n>0] := e2@@(Range[n])

(* HomogeneousQ returns True if the input is an homogeneous multivector *)
HomogeneousQ[x_,r_/;Element[r,Integers]&&NonNegative[r]] := Simplify[Expand[x]===Grade2[x,r]]

(* The RELATIONS of the Clifford algebra *)
e2[]:=1
e2[i_Integer,j__Integer] := e2[]                               /; i==j && i<=4 && EvenQ[Length[{i,j}]] && i>=0
e2[i_Integer,j__Integer] := e2[i]                              /; i==j && i<=4 && OddQ[Length[{i,j}]] && i>=0
e2[i_Integer,j__Integer] := (-e2[])^(Length[{i,j}]/2)          /; i==j && i>4 && i<=5&&EvenQ[Length[{i,j}]] && i>=0
e2[i_Integer,j__Integer] :=(-e2[])^((Length[{i,j}]-1)/2) e2[i] /; i==j && i>4 && i<=5&&OddQ[Length[{i,j}]] && i>=0
e2[i_Integer,j__Integer] :=0/;i==j&&Max[{i,j}]>5&&AllTrue[{i,j},NonNegative]
e2[i_Integer,j_Integer] :=-e2[j,i]                             /; i!=j && i>j && AllTrue[{i,j},NonNegative]
e2[i__Integer] := Signature[Ordering[{i}]] Apply[e2,Sort[{i}]] /; !OrderedQ[{i}] && AllTrue[{i},NonNegative]
e2[i__Integer] := (Times@@(Apply[e2,Gather[{i}],1]/. e2[_]->1))*e2@@(Cases[Apply[e2,Gather[{i}],1],Except[_Integer]]/. (e2[x_]|a_ e2[x_])->x)                                                                   /; OrderedQ[{i}] && !DuplicateFreeQ[{i}] && AllTrue[{i},Positive]

(* Begin Geometric Product Section *)
GeometricProduct2[_] := $Failed
GeometricProduct2[x_,y_] := GeometricProduct2[Expand[x],Expand[y]]          /; x=!=Expand[x] || y=!=Expand[y]
GeometricProduct2[x_,y_,z__] := Fold[GeometricProduct2,Expand[x],{Expand[y],z}]//Simplify
GeometricProduct2[x_,y_+z_] := GeometricProduct2[x,y]+GeometricProduct2[x,z]
GeometricProduct2[x_+y_,z_] := GeometricProduct2[x,z]+GeometricProduct2[y,z]
GeometricProduct2[a_,b_] := a b e2[]                                        /; FreeQ[a,e2[__]] && FreeQ[b,e2[__]]
GeometricProduct2[a_,(b_:1) e2[i__/;SubsetQ[Range[5],{i}]]] := a b e2[i]    /; FreeQ[a,e2[__]] && FreeQ[b,e2[__]]
GeometricProduct2[(a_:1) e2[i__/;SubsetQ[Range[5],{i}]],b_] := a b e2[i]    /; FreeQ[a,e2[__]] && FreeQ[b,e2[__]]
GeometricProduct2[(a_:1) e2[i__/;SubsetQ[Range[5],{i}]],(b_:1) e2[j__/;SubsetQ[Range[5],{j}]]] := a b e2[i,j]  /; FreeQ[a,e2[__]] && FreeQ[b,e2[__]]
(* End of Geometric Product Section *)

(* Begin Grade Section *)
Grade2[a_,r_/;Element[r,Integers]] := If[r===0,a,0]        /; FreeQ[a,e2[__]]
Grade2[x_,r_/;Element[r,Integers]] := Grade2[Expand[x],r]  /; x=!=Expand[x]
Grade2[x_,r_/;Element[r,Integers]] := 0                    /;r<0
Grade2[x_+y_,r_/;Element[r,Integers]] := Grade2[x,r]+Grade2[y,r]
Grade2[(a_:1) e2[i__/;SubsetQ[Range[5],{i}]],r_/;Element[r,Integers]]:=If[Length[{i}]===r,a e2[i],0]
(* End of Grade Section *)

(* Begin Inner Product Section *) 
InnerProduct2[_] := $Failed
InnerProduct2[x_,y_] := InnerProduct2[Expand[x], Expand[y]]         /; x=!=Expand[x] || y=!=Expand[y]
InnerProduct2[x_,y_+z_] := InnerProduct2[x,y]+InnerProduct2[x,z]
InnerProduct2[x_+y_,z_] := InnerProduct2[x,z]+InnerProduct2[y,z]
InnerProduct2[a_,b_] := 0                                          /; FreeQ[a,e2[__]] && FreeQ[b,e2[__]]
InnerProduct2[a_,(b_:1) e2[i__/;SubsetQ[Range[5],{i}]]] :=0        /; FreeQ[a,e2[__]]
InnerProduct2[(b_:1) e2[i__/;SubsetQ[Range[5],{i}]],a_] :=0        /; FreeQ[a,e2[__]]
(* A la Hestenes *)
InnerProduct2[(a_:1) e2[i__/;SubsetQ[Range[5],{i}]],(b_:1) e2[j__/;SubsetQ[Range[5],{j}]]] := Grade2[a b e2[i,j],Abs[Length[{i}]-Length[{j}]]] /; FreeQ[a,e2[__]] && FreeQ[b,e2[__]]
(* Dorst, Fontijne and Mann definition (with no absolute value as in Hestenes) *)
(* InnerProduct2[(a_:1) e2[i__/;SubsetQ[Range[5],{i}]], (b_:1) e2[j__/;SubsetQ[Range[5],{i}]]] := 0   /; FreeQ[a, e2[__]] && FreeQ[b, e2[__]] && Length[{i}] > Length[{j}] *)
(* InnerProduct2[(a_:1) e2[i__/;SubsetQ[Range[5],{i}]], (b_:1) e2[j__/;SubsetQ[Range[5],{j}]]] := Grade2[a b e2[i,j], Length[{j}] - Length[{i}]]  /; FreeQ[a, e2[__]] && FreeQ[b, e2[__]] && Length[{j}] >= Length[{i}] *)

(* End of Inner Product Section *)

(* Begin Outer Product Section *)
OuterProduct2[_] := $Failed
OuterProduct2[x_,y_] := OuterProduct2[Expand[x],Expand[y]]              /; x=!=Expand[x] || y=!=Expand[y]
OuterProduct2[x_,y_,z__] := Fold[OuterProduct2,x,{y,z}]//Simplify
OuterProduct2[x_,y_+z_] := OuterProduct2[x,y]+OuterProduct2[x,z]
OuterProduct2[x_+y_,z_] := OuterProduct2[x,z]+OuterProduct2[y,z]
OuterProduct2[a_,b_] := a b                                             /; FreeQ[a,e2[__]] && FreeQ[b,e2[__]]
OuterProduct2[a_,(b_:1) e2[i__/;SubsetQ[Range[5],{i}]]] := a b e2[i]    /; FreeQ[a,e2[__]]
OuterProduct2[(b_:1) e2[i__/;SubsetQ[Range[5],{i}]],a_] := a b e2[i]    /; FreeQ[a,e2[__]]
OuterProduct2[(a_:1) e2[i__/;SubsetQ[Range[5],{i}]],(b_:1) e2[j__/;SubsetQ[Range[5],{j}]]] := Grade2[a b e2[i,j],Length[{i}]+Length[{j}]]      /; FreeQ[a,e2[__]] && FreeQ[b,e2[__]]
(* End of Outer Product Section *)

(* End of the Clifford algebra of R^{p,q} section *)



(* The relations of the CONFORMAL GEOMETRIC ALGEBRA and its generaliations *)
e[] := 1
(* The next definitions takes care of the existence of a 0 after an \[Infinity], despite duplication *)
e[i___Integer, j___] := 0      /; Count[{i}, 0] >= 2 || (DuplicateQ0Inf[{j}] && (Max[Position[{j}, 0]] < FirstPosition[{j}, \[Infinity]][[1]] || Contiguous0Inf[{i}] ))
e[i_, j__] := e[]              /; i == j && EvenQ[Length[{i, j}]] && MemberQ[{1, 2, 3}, i] 
e[i_, j__] := e[i]             /; i == j && OddQ[Length[{i, j}]] && MemberQ[{1, 2, 3}, i]
e[\[Infinity], 0] := -2 - e[0,\[Infinity]]
e[\[Infinity], j_Integer] := -e[j, \[Infinity]]   /;  MemberQ[{1, 2, 3}, j]
e[i_, j_] := -e[j, i]          /; i != j && i > j && SubsetQ[{0, 1, 2, 3}, {i}]
e[i__] := Signature[Ordering[{i}]] e @@ Sort[{i}]   /; 
          !OrderedQ[{i}] && !DuplicateQ0Inf[{i}] && (! MemberQ[{i}, 0] || !MemberQ[{i}, \[Infinity]] ||FirstPosition[{i}, 0][[1]] < FirstPosition[{i}, \[Infinity]][[1]])
e[i__] := Module[{es = Cases[Apply[e, Gather[{i}], 1], Except[_Integer]]},
   Return[(Times @@ Join[Cases[Apply[e, Gather[{i}], 1], _Integer], Cases[{Times @@ es}, _Integer]]) e @@ Cases[es, e[x_] :> x]]
                ]           /; OrderedQ[{i}] && ! DuplicateFreeQ[{i}] 
e[i___Integer, \[Infinity], j__Integer, k___] := 
   Signature[Ordering[{j}]] (-2 e @@ Join[{i}, DeleteCases[Sort[{j}], 0], {k}] - e @@ Join[{i}, {0, \[Infinity]}, DeleteCases[Sort[{j}], 0], {k}])   /; MemberQ[{j}, 0]
e[i___Integer, \[Infinity], j__Integer, k___] := 0   /; Count[{i}, 0] >= 2  && ! MemberQ[{j}, 0] && FirstPosition[{k}, \[Infinity]][[1]] < FirstPosition[{k}, 0][[1]]


(* GeometricProduct returns the geometric product between two or more multivectors *)
GeometricProduct[x_, y_] := GeometricProduct[Expand[x], Expand[y]]   /;  x =!= Expand[x] || y =!= Expand[y]
GeometricProduct[x_, y_, z__] := Fold[GeometricProduct, Expand[x], {Expand[y], z}] // Simplify
GeometricProduct[x_, y_ + z_] := GeometricProduct[x, y] + GeometricProduct[x, z]
GeometricProduct[x_ + y_, z_] := GeometricProduct[x, z] + GeometricProduct[y, z]
GeometricProduct[a_, b_] := a b e[]    /; FreeQ[a, e[__]] && FreeQ[b, e[__]]
GeometricProduct[a_ , (b_: 1) e[i__ /; SubsetQ[{0, 1, 2, 3, \[Infinity]}, {i}]]] := a b e[i]         /; FreeQ[a, e[__]] && FreeQ[b, e[__]]
GeometricProduct[(a_: 1) e[i__ /; SubsetQ[{0, 1, 2, 3, \[Infinity]}, {i}]], b_] :=  a b e[i]         /; FreeQ[a, e[__]] && FreeQ[b, e[__]] 
GeometricProduct[(a_: 1) e[i__ /; SubsetQ[{0, 1, 2, 3, \[Infinity]}, {i}]], (b_: 1) e[j__ /; SubsetQ[{0, 1, 2, 3, \[Infinity]}, {j}]]] := 
              a b e[i, j]   /; FreeQ[a, e[__]] && FreeQ[b, e[__]]

(* Grade[m,r] extracts the term of grade r of the multivector m *)
Grade[a_, r_?NumberQ] := If[r === 0, a, 0]             /; FreeQ[a, e[__]]
Grade[x_, r_?NumberQ] := Grade[Expand[x], r]   /; x =!= Expand[x]
Grade[x_ + y_, r_?NumberQ] := Grade[x, r] + Grade[y, r]
Grade[(a_: 1) e[i__ /; SubsetQ[{0, 1, 2, 3, \[Infinity]}, {i}]], r_?NumberQ] := If[Length[{i}] === r, a e[i], 0]


(*OuterProduct returns the outer product between two or more multivectors*)
OuterProduct[_] := $Failed
OuterProduct[x_,y_] := OuterProduct[Expand[x],Expand[y]]        /; x=!=Expand[x] || y=!=Expand[y]
OuterProduct[x_,y_] := From45To0I[OuterProduct2[FromOITo45[x], FromOITo45[y]]//Expand] //Expand
OuterProduct[x_,y_,z__] := Fold[OuterProduct,x,{y,z}] //Simplify

(* InnerProduct returns the inner product between two or more multivectors *)
InnerProduct[_] := $Failed
InnerProduct[x_,y_] := InnerProduct[Expand[x],Expand[y]]        /; x=!=Expand[x] || y=!=Expand[y]
InnerProduct[x_,y_] := From45To0I[InnerProduct2[FromOITo45[x],FromOITo45[y]] //Expand]//Expand

(* Reversion returns the reversion of a multivector *)
Reversion[_] := $Failed
Reversion[x_] := Reversion[Expand[x]]                                                        /; x =!= Expand[x]
Reversion[a_] := a                                                                           /; FreeQ[a, e[__]]
Reversion[(a_: 1) e[i__ /; SubsetQ[{0, 1, 2, 3, \[Infinity]}, {i}]]] := a e @@ Reverse[{i}]  /; FreeQ[a, e[__]]
Reversion[x_ + y_] := Reversion[x] + Reversion[y]

(* Involution is the grade involution (\[Dagger]) of a multivector *)
Involution[_] := $Failed
Involution[x_] := Involution[Expand[x]]                                                           /; x =!= Expand[x]
Involution[a_] := a                                                                               /; FreeQ[a, e[__]]
Involution[(a_: 1) e[i__ /; SubsetQ[{0, 1, 2, 3, \[Infinity]}, {i}]]] := (-1)^Length[{i}] a e[i]  /; FreeQ[a, e[__]]
Involution[x_ + y_] := Involution[x] + Involution[y]

(* Magnitude returns the magnitude v^2 of a multivector v *)
Magnitude[x_] := Grade[GeometricProduct[Reversion[x], x], 0]

(* MultivectorInverse returns the inverse of a multivector, if it exists *)
MultivectorInverse[_] := $Failed
MultivectorInverse[x_] := MultivectorInverse[Expand[x]]                 /; x =!= Expand[x]
MultivectorInverse[x_] := Simplify[Reversion[x] / Magnitude[x]]      /; Magnitude[x] =!= 0

(* I5 = pseudoscalar; I5i = the inverse of I5. Dorst, Fontijne, Mann, 13.3.3 *)
 I5 := OuterProduct[e[0],OuterProduct[e[1],e[2],e[3]],e[\[Infinity]]]
 I5i := OuterProduct[e[0],MultivectorInverse[OuterProduct[e[1],e[2],e[3]]],e[\[Infinity]]]

(* Dual returns the dual of a multivector *)
Dual[x_] := -InnerProduct[x, I5]

(* Rotation[x,a,b,theta] Rotates the vector x by an angle theta (in radians), 
along the plane defined by a and b. The sense of the rotation is from a to b.
If no theta is given, the default value is the  angle between a and b.
*)
Rotation[x_, a_, b_, angle_:Automatic] := Module[{plane=OuterProduct[a, b]},
    If[angle === Automatic, theta=VectorAngle[ToVector[a],ToVector[b]], theta=angle];
	Return[GFactor[GeometricProduct[Cos[theta/2] - (plane/Magnitude[plane])*Sin[theta/2], x, Cos[theta/2] + (plane/Magnitude[plane])*Sin[theta/2]]]]
    ]  /;  GradeQ[x,1] && GradeQ[a,1] && GradeQ[b,1]

(* definitions for system functions *)
SetAttributes[e,NHoldAll]

(* Protect[ Evaluate[protected] ]   restore protection of system symbols *)


End[] (* End Private Context *)

(* Protecting exported symbols *)
Protect[e, GradeQ, MultiplicationTable, GFactor, GeometricProduct, Grade, OuterProduct, InnerProduct, Reversion, 
	     Involution, Magnitude, MultivectorInverse, I5, I5i, Dual]

EndPackage[]
