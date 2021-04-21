(* ::Package:: *)

(* ::Input:: *)
(*Quit*)


(* ::Title:: *)
(*Construction of Eq. Triangle correction functions*)


(* ::Chapter:: *)
(*Generating the 2D Dubiner polynomial basis*)


(* ::Input:: *)
(*a0[xi_,eta_]:=3*xi/(2-Sqrt[3]*eta)*)
(*b0[xi_,eta_]:=(2*Sqrt[3]*eta-1)/3*)
(**)
(*(*i[p_,v_,w_]:=w+(p+1)v-v(v-1)/2*)*)
(**)
(*Qv[v_,a_]:=JacobiP[v,0,0,a]/Sqrt[2/(2v+1)]*)
(*Qw[v_,w_,b_]:=JacobiP[w,2v+1,0,b]/Sqrt[2^(2v+2)/(2w+2v+2)]*)
(*phi[v_,w_,a_,b_]:=2*Qv[v,a]*Qw[v,w,b]*(1-b)^v/3^(1/4)*)
(**)
(*(* Geometry Parameters *)*)
(*pt1={-1,-1/Sqrt[3]};*)
(*pt2={1,-1/Sqrt[3]};*)
(*pt3={0,2/Sqrt[3]};*)
(*edg1=Line[{pt1,pt2}];*)
(*f1e[h_]:=h;*)
(*f1n[h_]:=-1/Sqrt[3];*)
(*edg2=Line[{pt2,pt3}];*)
(*f2e[h_]:=(1-h)/2;*)
(*f2n[h_]:=Sqrt[3]/2*h+Sqrt[3]/6;*)
(*edg3=Line[{pt3,pt1}];*)
(*f3e[h_]:=(-1-h)/2;*)
(*f3n[h_]:=-Sqrt[3]*h+Sqrt[3]/6;*)
(*tri0=Triangle[{pt1,pt2,pt3}];*)
(*cen0=RegionCentroid[tri0];*)
(**)
(*(* Select Polynomial Order *)*)
(*p=2;*)
(**)
(**)
(*vw=Select[Tuples[Range[0,p],2],Total[#]<=p&];*)
(**)
(*For[k=1,k<=(p+1)(p+2)/2,k++,{*)
(*Print[Extract[vw,k]],*)
(*v=Extract[vw,{k,1}],*)
(*w=Extract[vw,{k,2}],*)
(**)
(*Print[      Expand[   phi[v,w,a0[\[Xi],\[Eta]],b0[\[Xi],\[Eta]]]   ]      ],*)
(**)
(*GraphicsRow[{Plot3D[phi[v,w,a0[\[Xi],\[Eta]],b0[\[Xi],\[Eta]]], {\[Xi],-1,1}, {\[Eta],-1,1},PlotRange->Full],*)
(*                      Plot3D[phi[v,w,a0[\[Xi],\[Eta]],b0[\[Xi],\[Eta]]], {\[Xi],\[Eta]}\[Element]tri0                    ,PlotRange->Full]},ImageSize->Full]//Print,*)
(**)
(*testQv=Integrate[      Qv[v,a]^2                                                                            ,{a,-1,1}],*)
(*testQw=Integrate[      Qw[v,w,b]^2*(1-b)^(2v+1)                                  ,{b,-1,1}],*)
(*(*testSq=Integrate[      phi[v,w,a,b]^2                                      ,{a,-1,1},{b,-1,1}],*)*)
(*testTr=Integrate[      phi[v,w,a0[\[Xi],\[Eta]],b0[\[Xi],\[Eta]]]^2                          ,{\[Xi],\[Eta]}\[Element]tri0],*)
(**)
(*If[testQv!= 1,Print["Qv not normalized             ", testQv]],*)
(*If[testQw!= 1,Print["Qw not normalized             ", testQw]],*)
(*(*If[testSq\[NotEqual] 1,Print["Phi not normalized in square  ", testSq]],*)*)
(*If[testTr!= 1,Print["Phi not normalized?????       ", testTr]]*)
(*}]*)
(**)
(**)


(* ::Chapter:: *)
(*Solving the linear system for energy stability*)


(* ::Section:: *)
(*Construct the h*n polynomial over the edge*)


(* ::Input:: *)
(*sp=Solve[LegendreP[p+1,x]==0, x];*)
(*absc=x/.sp;*)
(*Print[StringForm["Flux Point locatons ``",absc]]*)
(*Print[]*)
(**)
(*(*For[f=1,f\[LessEqual]3,f++,{(* Loop through the faces of a tiangle *)*)*)
(*For[fp=1,fp<=p+1,fp++,{(* Loop through the FPs in a face *)*)
(*Print[StringForm["Flux Point ``:",fp]],*)
(*rhs={},*)
(**)
(*base= Transpose[{absc,Extract[IdentityMatrix[p+1],fp]}],(* Build the {h,y} reference points *)*)
(*hfkn[h]=InterpolatingPolynomial[base,h], (* Interpolate the h*n edge polynomial *)*)
(**)
(*For[i=1,i<=(p+1)(p+2)/2,i++,{(* Loop through the polynomial basis, each step is a eq in the system *)*)
(*vi=Extract[vw,{i,1}],*)
(*wi=Extract[vw,{i,2}],*)
(*AppendTo[rhs,Integrate[hfkn[h]*phi[vi,wi,a0[f2e[h],f2n[h]],b0[f2e[h],f2n[h]]],{h,-1,1}]], (* Integrate the RHS of the eq *)*)
(**)
(*(* If we chose cTriM \[Equal] 0 then the LHS is blank and we recover the DG scheme*)
(*lhs1={},*)
(*For[j=1,j\[LessEqual](p+1)(p+2)/2,j++,{(* Now we need to build the left side matrix, element by element *)*)
(*vj=Extract[vw,{j,1}],*)
(*wj=Extract[vw,{j,2}],*)
(*AppendTo[lhs1,Sum[D[phi[vi,wi,a0[\[Xi],\[Eta]],b0[\[Xi],\[Eta]]],{\[Xi],p-m+1},{\[Eta],m-1}]*D[phi[vj,wj,a0[\[Xi],\[Eta]],b0[\[Xi],\[Eta]]],{\[Xi],p-m+1},{\[Eta],m-1}],{m,1,p+1}]]*)
(*}],*)
(*Print[lhs1]*)*)
(*}],*)
(*Print[rhs],*)
(*Print[Plot3D[Sum[Extract[rhs,i]*phi[Extract[vw,{i,1}],Extract[vw,{i,2}],a0[\[Xi],\[Eta]],b0[\[Xi],\[Eta]]],{i,1,(p+1)(p+2)/2}], {\[Xi],\[Eta]}\[Element]tri0,PlotRange->Full]],*)
(*}]*)
(**)
(**)


(* ::Input:: *)
(**)
(**)


(* ::Input:: *)
(**)
(**)


(* ::Input:: *)
(* *)


(* ::Chapter:: *)
(*Plotting the Correction Functions*)
