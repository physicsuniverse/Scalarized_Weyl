Off[General::"spell", General::"spell1"];

TrigRules = {x_. Cos[\[Theta]_]^2 + x_. Sin[\[Theta]_]^2 -> x, 
   x_. Cot[\[Theta]_]^2 + y_. Csc[\[Theta]_]^2 :> y /; x + y === 0};

SetAttributes[FaSi, {Listable}]; FaSi[0] = 0; 
FaSi[x_SeriesData] := x /. x[[3]] -> FaSi[x[[3]]]; 
FaSi[x_] := Factor[x /. simpRules] /. simpRules;

SetAttributes[FacSimp, {Listable}]; FacSimp[0] = 0; 
FacSimp[x_Equal] := Equal[FacSimp[x[[1]] - x[[2]]], 0];
FacSimp[x_Rule] := Rule[x[[1]], FacSimp[x[[2]]]]; 
FacSimp[x_SeriesData] := x /. x[[3]] -> FacSimp[x[[3]]];
FacSimp[x_] := (Collect[coll[x] /. simpRules, FormVars, Factor] /. 
     simpRules) /; FormDegree[x] > 0;
FacSimp[x_] := Factor[x /. simpRules] /. simpRules;

zeroQ[x_] := Union[Flatten[{0, Normal[x]}]] === {0};

nonZeroL[x_] := Select[Union[Flatten[{Normal[x]}]], # =!= 0 &]

nonZeroN[x_] := Length[Select[Flatten[{Normal[x]}], # =!= 0 &]]

indepTerms[x_, opt_ : {}] := redPlus[redTimes[nonZeroL[x], opt]];

redTimes[x_, opt_ : {}] := 
 Block[{tmp$ = Select[x, Head[#1] === Times &], rst$}, 
  rst$ = Complement[x, tmp$]; Union[rst$, Map[numb2[#, opt] &, tmp$]]]

numb1[x_, opt_ : {}] := 1 /; NumericQ[x]; 
numb1[x_^y_., opt_ : {}] := 1 /; MemberQ[opt, x]; numb1[x_, opt_ : {}] := x; 
numb2[x_, opt_ : {}] := Times @@ Map[numb1[#, opt] &, Level[x, 1]]

redPlus[x_] := 
  Block[{tmp$ = Select[x, Head[#] === Plus &], rst$, i$, negi$}, 
   rst$ = Complement[x, tmp$]; i$ = 1; 
   While[i$ < Length[tmp$], negi$ = Select[tmp$, # === -tmp$[[i$]] &]; 
    If[negi$ =!= {}, tmp$ = Complement[tmp$, negi$]]; i$ = i$ + 1]; 
   Union[tmp$, rst$]];

RiemSym[h_] := (h[a_, a_, c_, d_] = 0; h[a_, b_, c_, c_] = 0;
   h[a_, b_, c_, d_] := -h[b, a, c, d] /; a > b; 
   h[a_, b_, c_, d_] := -h[a, b, d, c] /; c > d;
   h[a_, b_, c_, d_] := h[c, d, a, b] /; a > c || (a === c && b > d));

epsilon[a__] := Signature[{a}];

zeroTensor[rnk_] := Fold[Table, 0, List /@ Table[Dim, {rnk}]];

FuncRepRules[f_[y__], g_, ord_ : 2] := 
 FuncRepRules[f[y], g, ord] = 
  Block[{tmp$ = {f[y] -> g}, indL = indexList[Length[{y}], ord]},
          Do[AppendTo[tmp$, Derivative[Sequence @@ indL[[i]]][f][y] -> 
                  
      Expand[FaSi[D[g, Sequence @@ Transpose[{{y}, indL[[i]]}]]]]], {i, 
     Length[indL]}]; tmp$]

indexList[varNo_, ord_] := 
 Block[{tmp$, indL}, 
  tmp$ = NestList[RotateRight, Join[{1}, Table[0, {varNo - 1}]], varNo - 1]; 
  If[ord > 1, 
   Do[indL = tmp$; 
    Do[tmp$ = Join[tmp$, Map[Plus[indL[[i]], #] &, indL]], {i, 
      varNo}], {ord - 1}]]; Union[tmp$]]

MF[x_] := MatrixForm[x] /; Dim < 10; MF[x_] := x; DeclareForms[1, e[_]];

stPrint[x_, y_ : Null] := 
  If[y === Null, 
   Print[Style[x, FontSize -> 20, FontWeight -> "Bold", 
     FontColor -> RGBColor[1, 0, 1]]], 
   Print[Style[x, FontSize -> 16, FontSlant -> "Italic", 
     FontColor -> RGBColor[1, 0, 1]]]];

timePrint[x_, sT_] := Print[x <> " computed in ", TimeUsed[] - sT, " sec"];

SetAttributes[myCoef$, {Listable}]; 
myCoef$[x_SeriesData, y__] := 
 x /. (x[[3]] -> Map[Coefficient[#, y] &, x[[3]]]); 
myCoef$[x_, y__] := Coefficient[x, y];

RGtensors[gIN_, xIN_, opt___] := 
  Module[{frameOpt = 0, eIN, Ropt = True, Wopt = True, idMat, eFrame, dxRul, 
    Bmat, Amat, NP$ = False, de, a, b, c, id, k, 
    StrCon, \[Gamma]Udd, \[Gamma]ddd, gddd, Gamddd, rmn, tmp, RUd, Rg, gg, sT,
     Cflat, EinSp}, 
   Clear[GUdd, \[Omega]Ud, RUddd, Rdd, EUd, R, Wdddd, detg, rtAbsdetg, eta, 
    HStar, PS, PH];
   		
   		If[Transpose[gIN] =!= gIN, 
    Print[MatrixForm[gIN], " is not a symmetric matrix"]; Return[]];
    Dim = Length[xIN]; 
   If[Dim =!= Length[gIN], 
    stPrint["Metric and coordinate dimensions don't match", 1]; Return[]];
   		
   If[gIN === {{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, -1}, {0, 0, -1, 0}}, 
    NP$ = True]; eFrame = Array[e, {Dim}];
   (* frameOpt=0 => coordFrame, frameOpt=1 => explicitFrame, frameOpt=2 => 
   implicitFrame *)
   		
   	If[Length[{opt}] > 0, 
    If[Head[Flatten[{opt}][[1]]] === e, 
     If[Head[dxRuleList] === List && Head[deList] === List, frameOpt = 2; 
      eIN = eFrame, 
      stPrint["dxRuleList and/or deList have not been defined!", 1]; 
      Return[]]];
    	Which[! IntegerQ[{opt}[[1]][[1]]] && d[0] =!= 0, 
     stPrint["matrixEDCcode.nb must be evaluated before section 2 \
RG&TC-Code!", 1]; Return[],
     	Length[{opt}] === 2 && FormDegree[{opt}[[1]]] === 1, 
     If[frameOpt =!= 2, eIN = {opt}[[1]]; frameOpt = 1]; 
     If[{opt}[[2]][[1]] === 0, Ropt = False]; 
     If[{opt}[[2]][[2]] === 0, Wopt = False],
     	Length[{opt}] === 1 && IntegerQ[opt[[1]]], 
     If[opt[[1]] === 0, Ropt = False]; If[opt[[2]] === 0, Wopt = False],
     	Length[{opt}] === 1 && FormDegree[opt[[1]]] === 1, 
     If[frameOpt =!= 2, eIN = opt; frameOpt = 1],
     	True, stPrint["Wrong optional arguments!", 1]; Return[]]];
   		
   If[! FreeQ[{gIN, eIN, {deList}}, SeriesData],
    	tmp = 
     Select[Flatten[{gIN, eIN, {deList}}], Head[#] === SeriesData &][[1]][[
      1]]; If[MemberQ[xIN, tmp], de[1] = ToExpression["d" <> ToString[tmp]]; 
     de[2] = ToExpression["\[ScriptD]" <> ToString[tmp]]; 
     DeclareForms[1, \[ScriptD][_], de[1], de[2]]; 
     d[tmp] = de[1]; \[ScriptD][tmp] = de[2]; d[de[1]] = 0, d[tmp] = 0]];
   		
   If[frameOpt === 1, DeclareForms[1, \[ScriptD][_]]; 
    origFRin$ = eIN /. MapThread[Rule, {d /@ xIN, \[ScriptD] /@ xIN}]; 
    eTO$dx = MapThread[Rule, {e /@ Range[Dim], origFRin$}]];
   				
   If[frameOpt > 
      0 && (Union[FormDegree /@ eIN] =!= {1} || 
       Length[eIN] =!= Dim || ! 
        zeroQ[eIN /. {d[k_] -> 0, e[_] -> 0, 
           y_ | Bar[y_] :> 0 /; MemberQ[DifFormSymbols, y]}]), 
    stPrint["coframe has incorrect form", 1]; Return[]];	
   		
   If[Head[simpRules] === Rule || Head[simpRules] === RuleDelayed, 
    simpRules = {simpRules}];
   If[Head[simpRules] =!= List, simpRules = {}]; idMat = IdentityMatrix[Dim];
   
   Which[frameOpt === 0, Amat = idMat; Bmat = idMat; eIN = d /@ xIN; 
    dxRul = {},
    frameOpt === 1 && 
     zeroQ[eIN /. {d[k_] -> 0, 
        y_ | Bar[y_] :> 0 /; MemberQ[DifFormSymbols, y]}], 
    Bmat = Table[myCoef$[eIN[[a]], d[xIN[[b]]]], {a, Dim}, {b, Dim}];
    		If[! zeroQ[FaSi[Det[Bmat]]], sT = TimeUsed[]; 
     If[ByteCount[eIN] > 1500*Dim, stPrint["Computing d[coframe]...", 1]]; 
     Amat = FaSi[Inverse[Bmat]]; 
     dxRul = MapThread[Rule, {d[xIN], Amat . eFrame}]; 
     de = If[FreeQ[eIN, SeriesData], FacSimp[d[eIN] /. dxRul], 
       FacSimp[Wedge[FacSimp[FacSimp[d[Bmat]] /. dxRul], d[xIN] /. dxRul]]]; 
     timePrint["d[coframe]", sT]; tmp = d[xIN] /. dxRul; 
     MapThread[Set, {d[xIN], tmp}], 
     stPrint["coframe does not span space!", 1]; Return[]],
    frameOpt === 2, de = deList; dxRul = dxRuleList; tmp = d[xIN] /. dxRul; 
    MapThread[Set, {d[xIN], tmp}]; 
    Amat = Table[myCoef$[tmp[[a]], eIN[[b]]], {a, Dim}, {b, Dim}]];
   		
   k = {(-1 + x_Symbol)^b_. (1 + x_Symbol)^b_. -> (-1 + x^2)^
      b, (x_Symbol - y_Symbol)^b_. (x_Symbol + y_Symbol)^b_. -> (x^2 - y^2)^
      b}; gdd = gIN; coordList = xIN; Print["gdd = ", MF[gdd]]; 
   tmp = eIN . gIN . eIN; 
   Print["\!\(\*SuperscriptBox[\(ds\), \(2\)]\) = ", 
    If[Head[tmp] === SeriesData, FacSimp[tmp], 
      Collect[Expand[tmp], {d[a_]^2, d[a_] d[b_], e[a_]^2, e[a_] e[b_]}, 
       Factor]] //. k]; sT = TimeUsed[]; 
   gUU = FixedPoint[FaSi, Inverse[gIN]] //. k; 
   detg = FixedPoint[FaSi, Det[gIN]] //. k;
   Print["gUU = ", MF[gUU]]; Clear[k]; timePrint["gUU", sT];
   		
   rtAbsdetg = PowerExpand[Power[Factor[detg^2], 1/4]]; 
   SetAttributes[HStar, {Listable}]; If[NP$, rtAbsdetg = I]; 
   HStar[x_SeriesData] := x /. x[[3]] -> HStar /@ x[[3]]; 
   HStar[x_Plus] := FacSimp[Plus @@ HStar /@ x]; 
   HStar[x_*y_] := x*HStar[y] /; FormDegree[x] === 0; 
   HStar[x_] := rtAbsdetg Wedge @@ e /@ Range[Dim] /; FormDegree[x] === 0; 
   HStar[x_Wedge] := rtAbsdetg/detg /; Length[x] === Dim; 
   HStar[x_Wedge | x_e] := 
    HStar[x] = rtAbsdetg FacSimp[inn[x, Wedge @@ e /@ Range[Dim]]];		
   		
   eta[] := eta[] = rtAbsdetg Array[epsilon, Table[Dim, {Dim}]];
   eta[a__] := eta[a] = Raise[eta[], a]; sT = TimeUsed[];
   		
   If[frameOpt > 0, d[e[a_]] := de[[a]]; StrCon[c_, a_, a_] = 0; 
    StrCon[c_, a_, b_] := -StrCon[c, b, a] /; a > b; 
    StrCon[c_, a_, b_] := -myCoef$[de[[c]], e[a]\[Wedge]e[b]];
    \[Gamma]Udd = Array[StrCon, {Dim, Dim, Dim}]; \[Gamma]ddd = 
     gdd . \[Gamma]Udd, \[Gamma]Udd = 
     zeroTensor[3]; \[Gamma]ddd = \[Gamma]Udd];
   
   (** Tests **)
   If[frameOpt > 0, tmp = (d[gdd] /. dxRul) /. e[_] -> 0; 
    If[! zeroQ[tmp], stPrint["Warning: d[] on some symbols not defined!", 1]; 
     Print[indepTerms[FacSimp[tmp]]]]; 
    If[frameOpt === 2 && (ByteCount[deList] + ByteCount[dxRuleList])/Dim < 
       6000, tmp = FacSimp[d[deList]]; 
     If[! zeroQ[tmp], stPrint["Warning: deList not closed! d[deList]=", 1]; 
      Print[tmp]; If[! zeroQ[tmp /. d[_] -> 0], Return[]]]]];
   		
   Clear[pD]; 
   pD[x_, y_] := 
    pD[x, y] = 
     If[Amat === idMat, D[x, coordList[[y]]], 
      Sum[Amat[[k, y]]*D[x, coordList[[k]]], {k, Dim}]];	
   		
   gddd[c_, a_, b_] := gddd[c, b, a] /; a > b; 
   gddd[c_, a_, b_] := 
    gddd[c, a, b] = 
     FaSi[(pD[gdd[[a, c]], b] + pD[gdd[[b, c]], a] - pD[gdd[[a, b]], c])/2];
   If[frameOpt === 0, Gamddd = Array[gddd, {Dim, Dim, Dim}]; Clear[tmp]; 
    k = gUU . Gamddd; tmp[c_, a_, b_] := tmp[c, b, a] /; a > b; 
    tmp[c_, a_, b_] := tmp[c, a, b] = FaSi[k[[c, a, b]]]; 
    GUdd = Array[tmp, {Dim, Dim, Dim}],
    Gamddd = 
     FaSi[Array[
        gddd, {Dim, Dim, 
         Dim}] + (\[Gamma]ddd + Transpose[\[Gamma]ddd] - 
          Transpose[\[Gamma]ddd, {3, 2, 1}])/2]; 
    GUdd = If[
      Union[NumericQ /@ Flatten[gUU]] === {True} && nonZeroN[gUU] === Dim, 
      gUU . Gamddd, FaSi[gUU . Gamddd]]];
   		
   timePrint["Gamma", sT]; 
   If[frameOpt > 0, \[Omega]Ud = 
     Table[Sum[GUdd[[a, k, b]] e[k], {k, Dim}], {a, Dim}, {b, Dim}]];
   		
   sT = TimeUsed[]; RiemSym[rmn]; 
   rmn[a_, b_, c_, id_] := 
    rmn[a, b, c, id] = 
     FaSi[Sum[Amat[[k, c]] D[Gamddd[[a, id, b]], xIN[[k]]] - 
        Amat[[k, id]] D[Gamddd[[a, c, b]], xIN[[k]]] + 
        GUdd[[k, c, b]]*Gamddd[[k, id, a]] - 
        GUdd[[k, id, b]]*Gamddd[[k, c, a]] - 
        Gamddd[[a, k, b]]*\[Gamma]Udd[[k, c, id]], {k, Dim}]]; 
   Rdddd = Array[rmn, {Dim, Dim, Dim, Dim}]; timePrint["Riemann(dddd)", sT]; 
   sT = TimeUsed[];
   		If[zeroQ[Rdddd], stPrint["Flat Space!"]; Clear[Rdd, EUd, R, Wdddd]; 
    Return[]];
   		
   Rg = gUU . Rdddd;
   If[Ropt, Clear[rmn]; rmn[a_, b_, c_, c_] = 0; 
    rmn[a_, b_, c_, id_] := -rmn[a, b, id, c] /; c > id; 
    rmn[a_, b_, c_, id_] := rmn[a, b, c, id] = FaSi[Rg[[a, b, c, id]]]; 
    RUddd = Array[rmn, {Dim, Dim, Dim, Dim}]; timePrint["Riemann(Uddd)", sT]; 
    sT = TimeUsed[]; Rg = RUddd, stPrint["RUddd not computed", 1]];
   		
   Clear[tmp]; tmp[a_, b_] := tmp[b, a] /; a > b; 
   tmp[a_, b_] := tmp[a, b] = FaSi[Sum[Rg[[k, a, k, b]], {k, Dim}]]; 
   Rdd = Array[tmp, {Dim, Dim}];
   				
   RUd = gUU . Rdd; R = FaSi[Sum[RUd[[a, a]], {a, Dim}]]; 
   timePrint["Ricci", sT]; sT = TimeUsed[];
   		
   If[Wopt || Dim < 4, If[Dim > 3, If[zeroQ[Rdd], Wdddd = Rdddd,
      	gg = Transpose[Outer[Times, gdd, gdd], {1, 3, 2, 4}]; 
      gg = gg - Transpose[gg]; 
      Rg = Transpose[Outer[Times, Rdd, gdd], {1, 3, 2, 4}]; 
      Rg = Rg - Transpose[Rg]; Rg = Rg - Transpose[Rg, {1, 2, 4, 3}];
      Clear[rmn]; RiemSym[rmn]; 
      rmn[a_, b_, c_, id_] := 
       rmn[a, b, c, id] = 
        FaSi[Rdddd[[a, b, c, id]] - Rg[[a, b, c, id]]/(Dim - 2) + 
          R gg[[a, b, c, id]]/((Dim - 1) (Dim - 2))]; 
      Wdddd = Array[rmn, {Dim, Dim, Dim, Dim}]], Wdddd = zeroTensor[4]]; 
    timePrint["Weyl", sT]; 
    If[zeroQ[Wdddd], Cflat = True; 
     If[Dim =!= 3 || zeroQ[CF3d], stPrint["Conformally Flat"]]], 
    stPrint["Weyl tensor not computed", 1]];
   		
   If[zeroQ[Rdd], stPrint["Ricci Flat"]; EUd = Rdd,
    			sT = TimeUsed[]; EUd = FaSi[RUd - (1/2) R*idMat]; 
    timePrint["Einstein", sT];
    		If[zeroQ[EUd + (1/2 - 1/Dim) R*idMat], stPrint["Einstein Space"]; 
     EinSp = True]];
   
   If[(Dim > 2 && Cflat && EinSp) || (Dim === 2 && zeroQ[FaSi[covD[R]]]), 
    stPrint["Space of Constant Curvature!"]]; If[NP$, NPdefs];];

CF3d := (stPrint["Testing for 3-dim conformal flatness...", 1]; 
   Block[{tmp$ = covD[Rdd] - Outer[Times, gdd, covD[R]/2/(Dim - 1)]}, 
    FaSi[tmp$ - Transpose[tmp$, {1, 3, 2}]]]);

NPdefs := (PH[0, 0] = Rdd[[2, 2]]/2; PH[0, 1] = -Rdd[[2, 4]]/2; 
   PH[1, 0] = -Rdd[[2, 3]]/2; PH[2, 2] = Rdd[[1, 1]]/2; 
   PH[2, 1] = -Rdd[[3, 1]]/2; PH[1, 2] = -Rdd[[4, 1]]/2; 
   PH[0, 2] = Rdd[[4, 4]]/2; PH[2, 0] = Rdd[[3, 3]]/2; 
   PH[1, 1] = FaSi[(Rdd[[1, 2]] + Rdd[[3, 4]])/4]; \[CapitalLambda] = -R/2 4; 
   tmp$ = "NP variables PH[i,j], \[CapitalLambda] "; 
   If[Head[Wdddd] =!= Symbol, PS[0] = -Wdddd[[2, 4, 2, 4]]; 
    PS[1] = Wdddd[[2, 4, 2, 1]]; PS[2] = -Wdddd[[2, 4, 3, 1]]; 
    PS[3] = Wdddd[[2, 1, 3, 1]]; PS[4] = -Wdddd[[3, 1, 3, 1]]; 
    tmp$ = StringJoin[tmp$, "and PS[i] "]]; 
   stPrint[StringJoin[tmp$, "defined."], 1]);

inn[y_Wedge, x_] := inn[Rest[y], inn[First[y], x]];
inn[e[a_], x_ + y_] := inn[e[a], x] + inn[e[a], y];
inn[e[a_], x_*y_] := x*inn[e[a], y] /; FormDegree[x] === 0;
inn[e[a_], x_] := 0 /; FormDegree[x] === 0;
inn[e[a_], x_Wedge] := 
  inn[e[a], First[x]]*Rest[x] + (-1)^FormDegree[First[x]]*
    First[x]\[Wedge]inn[e[a], Rest[x]];
inn[e[a_], e[b_]] := gUU[[a, b]];

MoveIndex[x_, i_, g_] := Block[{rnk1 = TensorRank[x], indL},
       Which[i === 1, g . x, i === rnk1, x . g, True, 
    indL = Range[rnk1] /. {i -> 1, 1 -> i}; 
    Transpose[g . Transpose[x, indL], indL]]];

Lower[x_, i__] := 
  FacSimp[FoldList[MoveIndex[#1, #2, gdd] &, x, Flatten[{i}]][[-1]]];

Raise[x_, i__] := 
  FacSimp[FoldList[MoveIndex[#1, #2, gUU] &, x, Flatten[{i}]][[-1]]];

multiDot[x_, y_, indPr__] := 
  Block[{rnk1 = TensorRank[x], rnk2 = TensorRank[y], 
    indL = Sort /@ Transpose[{indPr}], ordPr = Sort[{indPr}]}, 
   If[Union[indL[[1]]] =!= indL[[1]] || indL[[1]][[-1]] > rnk1 || 
     Union[indL[[2]]] =!= indL[[2]] || indL[[2]][[-1]] > rnk2, errWarn; 
    Return[]]; 
   Which[Length[ordPr] === 1, FacSimp[dotOne[x, y, indPr]], 
    rnk1 === rnk2 && rnk1 === Length[ordPr], FacSimp[dotAll[x, y, ordPr]], 
    True, OKindCon[dotOne[x, y, ordPr[[1]]], 
     renumb[rnk1, ordPr[[1, 2]]] /@ Rest[ordPr]]]];

renumb[rnk_, j_][{a_, b_}] := 
  If[j > b, {a - 1, b + rnk - 1}, {a - 1, b + rnk - 2}];

dotOne[x_, y_, {i_, j_}] := 
 Block[{rnk1 = TensorRank[x], rnk2 = TensorRank[y], tmp$ = FormDegree[0]}, 
  If[tmp$ === 0 && FormDegree[x] > 0 && FormDegree[y] > 0, tmp$ = Wedge, 
   tmp$ = Dot]; 
  Transpose[
   tmp$[Transpose[x, Range[rnk1] /. {rnk1 -> i, i -> rnk1}], 
    Transpose[y, Range[rnk2] /. {j -> 1, 1 -> j}]], 
   Join[Range[i - 1], RotateRight[Range[i, rnk1 - 1], 1], 
    RotateLeft[Range[rnk1, rnk1 + j - 2], 1], 
    Range[rnk1 + j - 1, rnk1 + rnk2 - 2]]]]

dotAll[x_, y_, indPr_] := 
  Block[{s$, allInd, sumInd, tmp$ = FormDegree[0]}, 
   allInd = Map[s$, Transpose[indPr], {2}] /. 
     s$[i_] :> ToExpression["s$" <> ToString[i]]; 
   sumInd = Transpose[{allInd[[1]], Table[Dim, {Length[allInd[[1]]]}]}]; 
   If[tmp$ === 0 && FormDegree[x] > 0 && FormDegree[y] > 0, tmp$ = Wedge, 
    tmp$ = Times]; 
   Sum @@ {Unevaluated[
      tmp$[x[[Sequence @@ allInd[[1]]]], y[[Sequence @@ allInd[[2]]]]]], 
     Sequence @@ sumInd}];

Contract[x_, indPr__] := 
  Block[{rnk1 = TensorRank[x], indL = Flatten[{indPr}]}, 
   If[Length[Union[indL]] < Length[indL] || Length[indL] > rnk1 || 
     Sort[indL][[-1]] > rnk1, errWarn; Return[], 
    OKindCon[x, Sort[Sort /@ {indPr}]]]];

newPrList[x_] := Rest[x] - 1 /. k_ :> k - 1 /; k > x[[1, 2]] - 1;

OKindCon[x_, indPr_] := FacSimp[ConOne[x, indPr[[1]]]] /; Length[indPr] === 1;
 OKindCon[x_, indPr_] := 
 Block[{s$, allInd, sumInd}, 
   allInd = Range[TensorRank[x]] /. 
      Apply[Rule, 
       Apply[Sequence, Map[Transpose[{#, s$[#[[1]]] {1, 1}}] &, indPr], 1], 
       1] /. s$[i_] :> ToExpression["s$" <> ToString[i]]; 
   sumInd = Transpose[{Union[allInd], Table[Dim, {Length[allInd]/2}]}]; 
   FacSimp[Sum @@ {Unevaluated[x[[Sequence @@ allInd]]], 
      Sequence @@ sumInd}]] /; TensorRank[x] === 2 Length[indPr];
OKindCon[x_, indPr_] := OKindCon[ConOne[x, indPr[[1]]], newPrList[indPr]];

ConOne[x_, {i_, j_}] := 
  Block[{rnk1 = TensorRank[x], tr1List, tr2List}, 
   tr1List = (Range[rnk1] /. {j -> i}) /. k_ :> k - 1 /; k > j; 
   tr2List = Join[Range[i - 1], RotateRight[Range[i, rnk1 - 1], 1]]; 
   Apply[Plus, Transpose[Transpose[x, tr1List], tr2List], {rnk1 - 2}]];

errWarn := stPrint["Wrong tensor rank or number(s) of indices!", 1];

covDiv[x_, UpList_] := 
  covDiv[x, UpList] = 
   Module[{k = TensorRank[x], divInd, allInd, freeInd, hhD, s, lt, aInd}, 
    If[Length[Flatten[UpList]] > k || Sort[Flatten[UpList]][[-1]] > k || 
      Head[x] =!= List, errWarn; Return[]]; 
    If[Length[UpList] > 1, 
     divInd = Flatten[Select[UpList, Head[#] === List &]][[1]], 
     divInd = Flatten[{UpList}][[1]]];
    			allInd = (aInd /@ Range[k]) /. aInd[divInd] -> s; 
    freeInd = Select[allInd, # =!= s &]; Off[Part::"pspec"];
    lt = If[$VersionNumber < 
       3.5, {Sum[hhD[x[[Sequence @@ allInd]], s], {s, Dim}] + 
        Sum[If[MemberQ[Flatten[UpList], i], 
          GUdd[[allInd[[i]], s, j]] x[[
            Sequence @@ ReplacePart[allInd, j, i]]], -GUdd[[j, s, 
             allInd[[i]]]] x[[Sequence @@ ReplacePart[allInd, j, i]]]], {j, 
          Dim}, {i, k}, {s, Dim}], 
       Sequence @@ (Map[List[#, Dim] &, freeInd])}, {Unevaluated[
        Sum[hhD[x[[Sequence @@ allInd]], s] + 
          Sum[If[MemberQ[Flatten[UpList], i], 
            Sum[GUdd[[allInd[[i]], s, j]] x[[
               Sequence @@ ReplacePart[allInd, j, i]]], {j, Dim}], -Sum[
              GUdd[[j, s, allInd[[i]]]] x[[
                Sequence @@ ReplacePart[allInd, j, i]]], {j, Dim}]], {i, 
            k}], {s, Dim}]], Sequence @@ (Map[List[#, Dim] &, freeInd])}];
    If[k > 1, lt = Table @@ lt, 
     lt = Sum[pD[x[[s$]], s$], {s$, Dim}] + 
       Sum[GUdd[[s$, s$, j]] x[[j]], {j, Dim}, {s$, Dim}]]; On[Part::"pspec"];
     Return[FaSi[Evaluate[lt /. hhD -> pD]]]];

covD[x_List, y_List] := covD[x, y] = c$D[x, y]; covD[x_, y_List] := errWarn; 
covD[x_] := covD[x] = c$D[x];

c$D[x_, UpList_ : {}] := 
  Module[{k = TensorRank[x]}, 
   If[Length[UpList] > 
      0 && (Length[Flatten[UpList]] > k || Sort[Flatten[UpList]][[-1]] > k), 
    errWarn; Return[]]; 
   If[Head[x] === List, Array[cDcomp[x, k, UpList], Table[Dim, {k + 1}]], 
    FaSi[Table[pD[x, i], {i, Dim}]]]];

cDcomp[x_, rnk_, UpList_ : {}][y__] := 
  Block[{y$ = Drop[{y}, -1], z$ = {y}[[-1]]}, 
   FaSi[pD[x[[Sequence @@ y$]], z$] + 
     Sum[If[MemberQ[UpList, i], 
       GUdd[[y$[[i]], z$, j]] x[[
         Sequence @@ ReplacePart[y$, j, i]]], -GUdd[[j, z$, y$[[i]]]] x[[
         Sequence @@ ReplacePart[y$, j, i]]]], {j, Dim}, {i, rnk}]]];

Bianchi[0] := 
  Block[{tmp$ = covD[Rdddd]}, 
   FaSi[indepTerms[
     tmp$ + Transpose[tmp$, {1, 2, 4, 5, 3}] + 
      Transpose[tmp$, {1, 2, 5, 3, 4}]]]];

Bianchi[1] := 
  Block[{tmp$ = covD[Rdd]}, 
   If[Head[RUddd] === Symbol, RUddd = Raise[Rdddd, 1]]; 
   FaSi[indepTerms[covDiv[RUddd, {1}] + tmp$ - Transpose[tmp$, {1, 3, 2}]]]];

Bianchi[2] := indepTerms[covDiv[EUd, {1}]];

$PrePrint = If[ByteCount[#] > 200000, Short[#, 100], #] &;

Protect[cDcomp, CF3d, ConOne, c$D, dotAll, dotOne, epsilon, errWarn, FacSimp, 
  FaSi, indexList, inn, MF, MoveIndex, newPrList, NPdefs, numb1, numb2, 
  OKindCon, redPlus, redTimes, renumb, RiemSym, zeroTensor];

On[General::"spell", General::"spell1"];