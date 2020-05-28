(* ::Package:: *)

(* ::Section:: *)
(*Ensembling Random Features - Functions for Plots*)


(* ::Text:: *)
(*In this file, we give the main functions for plotting the generalisation error of a random features system as well as the terms appearing in the bias variance decomposition. There are also the ones needed for the various comparisons between ensembling and regularisation.*)
(*Throughout we define: \[Psi]1=P/D=(Number of features)/(Data dimension) and \[Psi]2=N/D=(Training set size)/(data dimension)*)


(*Please Add here the location of the file EnsembleRFErrorCode.m assumed that is in the same repository as tis file*)
Get[StringJoin[NotebookDirectory[],"Computations.m"]]


(* ::Subsubsection:: *)
(*Function to Generate Plot*)


ClearAll[genPlot];
Options[genPlot]:={"pts"->Null,"plotRange"->{{-1,2},{0,1.5}},"xaxisLabel"->None,"yaxisLabel"->None,"leg"->Null,"colors"->0,"legPos"->{.82,.7},"dashed"->0,"logy"->False};
genPlot[OptionsPattern[]]:= Module[{plot,pts=OptionValue["pts"],xaxisLabel=OptionValue["xaxisLabel"],
yaxisLabel=OptionValue["yaxisLabel"],colors=OptionValue["colors"],dashed=OptionValue["dashed"],leg=OptionValue["leg"],plotRange=OptionValue["plotRange"],
legPos=OptionValue["legPos"],logy=OptionValue["logy"],style,plotStyle,plotfun},
If[colors==0,
colors=ColorData[97,"ColorList"];
colors=colors[[1;;Length[pts]]];
];

If[dashed==0,
plotStyle=colors;
];

If[Length[dashed]>0,
plotStyle=Partition[Riffle[colors,dashed],2]
];
style={FontSize->20,FontWeight->Bold,FontColor->Black};
If[logy==True,plotfun=ListLogPlot,plotfun=ListLinePlot];
If[leg=={},
plot=Labeled[Show[ plotfun[pts,Joined->True, PlotRange->plotRange,PlotStyle->plotStyle,AxesLabel->{None,yaxisLabel},LabelStyle->style],ImageSize->Large],xaxisLabel,{{Bottom,Right}},LabelStyle->style],
plot=Labeled[Show[ plotfun[pts,Joined->True, PlotRange->plotRange,PlotStyle->plotStyle,AxesLabel->{None,yaxisLabel},LabelStyle->style],Graphics[Inset[LineLegend[colors,leg,LabelStyle->style],Scaled@legPos]],ImageSize->Large],xaxisLabel,{{Bottom,Right}},LabelStyle->style];
];
plot
]


(* ::Subsubsection:: *)
(*Functions to Plot Term by Term*)


ClearAll[plotTermbyTermEvol];
Options[plotTermbyTermEvol]:={"\[Psi]1"->1, "\[Psi]2"->1, "plot"->False, "NbPoints"->30, "minVal"-> -1, "maxVal"->2.0,"\[Lambda]"->0.0,"\[Psi]1"->0.0,
"\[Mu]1"-> 0.5, "\[Mu]Star"-> 0.3014051374945435`}
plotTermbyTermEvol::usage = "plotTermbyTermEvol returns the points for the different terms appering in the generalisation error.
The evolution is as a function of P/N=\!\(\*FormBox[\(\*FractionBox[\(Number\\\ of\\\ features\), \(Training\\\ set\\\ size\)]\\\ if\\\ no\\\ value\\\ is\\\ passed\\\ and\\\ \[Lambda]\\\ value\\\ is\\\ given . \\\\nTo\\\ plot\\\ as\\\ a\\\ function\\\ of\\\ \[Lambda], \\\ pass\\\ a\\\ value\\\ of\\\ P/N\\\ \(\(used\)\(.\)\(\\\ \\\ \)\)\),
TraditionalForm]\)
options are: \[Psi]2,k,lambda:0 (means that lambda evolves),\[Psi]1->0.0, Plot:True, NbPoints: 30,minVal: -1 , maxVal:2.0,
F1:1.0, FStar:0.0, \[Tau]:0.0, \[Mu]1: Relu, \[Mu]Star: ReLU";


plotTermbyTermEvol[OptionsPattern[]] := Module[{
\[Psi]1=OptionValue["\[Psi]1"],
\[Psi]2=OptionValue["\[Psi]2"],
\[Lambda]=OptionValue["\[Lambda]"],
plot,colors,xVals,plotBol=OptionValue["plot"],\[Mu]1=OptionValue["\[Mu]1"],\[Mu]Star=OptionValue["\[Mu]Star"],
minVal=OptionValue["minVal"],maxVal=OptionValue["maxVal"],NbPoints=OptionValue["NbPoints"]
,axesLabel,terms,
\[CapitalPsi]1,\[CapitalPsi]2,\[CapitalPsi]3,\[CapitalPsi]4,\[CapitalPsi]5,\[CapitalPsi]6,leg},

{\[CapitalPsi]1,\[CapitalPsi]2,\[CapitalPsi]3}={{},{},{}};
{\[CapitalPsi]4,\[CapitalPsi]5}={{},{}};
\[CapitalPsi]6={};
xVals=N@Subdivide[minVal,maxVal,NbPoints-1];

If[\[Psi]1==0.0,
{
Do[
{
terms=getVanillaTerms[\[Psi]2*10^x,\[Psi]2,\[Lambda],\[Mu]1,\[Mu]Star];
AppendTo[\[CapitalPsi]1,{x,terms[[1]]}];
AppendTo[\[CapitalPsi]2,{x,terms[[2]]}];
AppendTo[\[CapitalPsi]3,{x,terms[[3]]}];
terms=getMixedTerms[\[Psi]2*10^x,\[Psi]2,\[Lambda],\[Mu]1,\[Mu]Star];
AppendTo[\[CapitalPsi]4,{x,terms[[1]]}];
AppendTo[\[CapitalPsi]5,{x,terms[[2]]}];
terms=getPsi4[\[Psi]2*10^x,\[Psi]2,\[Lambda],\[Mu]1,\[Mu]Star];
AppendTo[\[CapitalPsi]6,{x,terms}];
},
{x,xVals}
];
axesLabel={"Log10[P/N]"};
}];

If[\[Lambda]==0.0,
{
Do[
{
terms=getVanillaTerms[\[Psi]1,\[Psi]2,10^x,\[Mu]1,\[Mu]Star];
AppendTo[\[CapitalPsi]1,{x,terms[[1]]}];
AppendTo[\[CapitalPsi]2,{x,terms[[2]]}];
AppendTo[\[CapitalPsi]3,{x,terms[[3]]}];
terms=getMixedTerms[\[Psi]1,\[Psi]2,10^x,\[Mu]1,\[Mu]Star];
AppendTo[\[CapitalPsi]4,{x,terms[[1]]}];
AppendTo[\[CapitalPsi]5,{x,terms[[2]]}];
terms=getPsi4[\[Psi]1,\[Psi]2,10^x,\[Mu]1,\[Mu]Star];
AppendTo[\[CapitalPsi]6,{x,terms}];
},
{x,xVals}
];
axesLabel={"Log10[\[Lambda]]"};
}];

If[\[Psi]2==0.0,
{
Do[
{
terms=getVanillaTerms[\[Psi]1,10^x,\[Lambda],\[Mu]1,\[Mu]Star];
AppendTo[\[CapitalPsi]1,{x,terms[[1]]}];
AppendTo[\[CapitalPsi]2,{x,terms[[2]]}];
AppendTo[\[CapitalPsi]3,{x,terms[[3]]}];
terms=getMixedTerms[\[Psi]1,10^x,\[Lambda],\[Mu]1,\[Mu]Star];
AppendTo[\[CapitalPsi]4,{x,terms[[1]]}];
AppendTo[\[CapitalPsi]5,{x,terms[[2]]}];
terms=getPsi4[\[Psi]1,10^x,\[Lambda],\[Mu]1,\[Mu]Star];
AppendTo[\[CapitalPsi]6,{x,terms}];
},
{x,xVals}
];
axesLabel={"Log10[N/D]"};
}];

{\[CapitalPsi]1,\[CapitalPsi]2,\[CapitalPsi]3,\[CapitalPsi]4,\[CapitalPsi]5,\[CapitalPsi]6}=Select[#,#[[2]]>0&]&/@{\[CapitalPsi]1,\[CapitalPsi]2,\[CapitalPsi]3,\[CapitalPsi]4,\[CapitalPsi]5,\[CapitalPsi]6};
colors=ColorData[97,"ColorList"];
colors=colors[[1;;6]];
leg={"\[CapitalPsi]1","\[CapitalPsi]2","\[CapitalPsi]3","\[CapitalPsi]4","\[CapitalPsi]5","\[CapitalPsi]6"};
If[plotBol,
plot=genPlot["pts"->pts,"xaxisLabel"->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)(\!\(\*FractionBox[\(P\), \(N\)]\))","yaxisLabel"->None,"leg"-> leg];
Return[{\[CapitalPsi]1,\[CapitalPsi]2,\[CapitalPsi]3,\[CapitalPsi]4,\[CapitalPsi]5,\[CapitalPsi]6,plot}];
];
{\[CapitalPsi]1,\[CapitalPsi]2,\[CapitalPsi]3,\[CapitalPsi]4,\[CapitalPsi]5,\[CapitalPsi]6}
 ];


(* ::Subsubsection:: *)
(*Functions to Plot Bias and Variance*)


ClearAll[plotBiasVarRatioEvol];
Options[plotBiasVarRatioEvol]:={"\[Psi]1"->1.0,"\[Psi]2"->1.0,"k"->2,"plot"->True,"logy"->False,
"NbPoints"->30, "minVal"-> -1, "maxVal"->2.0,"\[Lambda]"->0.0,"\[Psi]1"->0.0,"plotRange"->{Full,{0,1.5}},
"F1"-> 1.0, "FStar"-> 0.0, "\[Tau]"-> 0.0, "\[Mu]1"-> 0.5, "\[Mu]Star"-> 0.3014051374945435`,
"pts"->{0},"legPos"->{0.1,0.1},"leg"->{"Init","Noise","Sampling","Bias","Test Error"}}
plotBiasVarRatioEvol::usage = "plotTermbyTermEvol returns the points for the different terms appering in the generalisation error.
The evolution is as a function of P/N=\!\(\*FormBox[\(\*FractionBox[\(Number\\\ of\\\ features\), \(Training\\\ set\\\ size\)]\\\ if\\\ no\\\ value\\\ is\\\ passed\\\ and\\\ \[Lambda]\\\ value\\\ is\\\ given . \\\\nTo\\\ plot\\\ as\\\ a\\\ function\\\ of\\\ \[Lambda], \\\ pass\\\ a\\\ value\\\ of\\\ P/N\\\ \(\(used\)\(.\)\(\\\ \\\ \)\)\),
TraditionalForm]\)
options are: \[Psi]2,k,lambda:0 (means that lambda evolves),pts: if you have the various \[CapitalPsi] terms already computed or not,\[Psi]1->0.0, Plot:True, NbPoints: 30,minVal: -1 , maxVal:2.0,
F1:1.0, FStar:0.0, \[Tau]:0.0, \[Mu]1: Relu, \[Mu]Star: ReLU";  

plotBiasVarRatioEvol[OptionsPattern[]] := Module[{
\[Psi]1=OptionValue["\[Psi]1"],\[Psi]2=OptionValue["\[Psi]2"],k=OptionValue["k"],
\[Lambda]=OptionValue["\[Lambda]"], logy=OptionValue["logy"],
pts=OptionValue["pts"],plot,colors,xVals,plotBol=OptionValue["plot"],
F1=OptionValue["F1"],FStar=OptionValue["FStar"],
\[Tau]=OptionValue["\[Tau]"],\[Mu]1=OptionValue["\[Mu]1"],\[Mu]Star=OptionValue["\[Mu]Star"],
minVal=OptionValue["minVal"],maxVal=OptionValue["maxVal"],NbPoints=OptionValue["NbPoints"],
plotRange=OptionValue["plotRange"],legPos=OptionValue["legPos"],leg=OptionValue["leg"],
xaxesLabel,terms,
\[CapitalPsi]1,\[CapitalPsi]2,\[CapitalPsi]3,\[CapitalPsi]4,\[CapitalPsi]5,\[CapitalPsi]6,
resPsiV,resPsiM,varIn,varData,varNoise,bias,tot},
If[\[Psi]1==0.0&&\[Lambda]==0,Print["Please give a value of \[Lambda] or a value of \[Psi]1"];Return[$Failed];];
legPos={0.85,0.7};
If[\[Psi]1==0.0,xaxesLabel="\!\(\*SubscriptBox[\(Log\), \(10\)]\)(P/N)"];
If[\[Lambda]==0.0,xaxesLabel="Log10[\[Lambda]]"];
If[\[Psi]2==0.0,xaxesLabel="\!\(\*SubscriptBox[\(Log\), \(10\)]\)(N/D)"];

If[Length[pts]==0,
{\[CapitalPsi]1,\[CapitalPsi]2,\[CapitalPsi]3,\[CapitalPsi]4,\[CapitalPsi]5,\[CapitalPsi]6}=plotTermbyTermEvol["\[Psi]1"->\[Psi]1,"\[Psi]2"->\[Psi]2,"plot"->False, "NbPoints"->NbPoints, "minVal"-> minVal, "maxVal"->maxVal,"\[Lambda]"->\[Lambda],"\[Psi]1"->\[Psi]1, "\[Mu]1"-> \[Mu]1, "\[Mu]Star"-> \[Mu]Star],
{\[CapitalPsi]1,\[CapitalPsi]2,\[CapitalPsi]3,\[CapitalPsi]4,\[CapitalPsi]5,\[CapitalPsi]6}=pts
];
\[CapitalPsi]1=Select[\[CapitalPsi]1,#[[2]]>0.0&];
\[CapitalPsi]2=Select[\[CapitalPsi]2,#[[2]]>0.0&];
\[CapitalPsi]3=Select[\[CapitalPsi]3,#[[2]]>0.0&];
\[CapitalPsi]4=Select[\[CapitalPsi]4,#[[2]]>0.0&];
\[CapitalPsi]5=Select[\[CapitalPsi]5,#[[2]]>0.0&];
\[CapitalPsi]6=Select[\[CapitalPsi]6,#[[2]]>0.0&];
xVals=#[[1]]&/@\[CapitalPsi]1;
varIn=Table[{xVals[[i]],1/k F1^2 (\[CapitalPsi]2[[i]][[2]]-\[CapitalPsi]4[[i]][[2]])},{i,1,Length[\[CapitalPsi]4]}];
varNoise=Table[{xVals[[i]],(\[Tau]^2/k (\[CapitalPsi]3[[i]][[2]]-\[CapitalPsi]5[[i]][[2]])+\[Tau]^2 \[CapitalPsi]5[[i]][[2]])},{i,1,Length[\[CapitalPsi]5]}];
varData=Table[{xVals[[i]],F1^2 (\[CapitalPsi]4[[i]][[2]]-\[CapitalPsi]6[[i]][[2]])},{i,1,Length[\[CapitalPsi]6]}];
bias=Table[{xVals[[i]],F1^2 (1-2 \[CapitalPsi]1[[i]][[2]]+\[CapitalPsi]6[[i]][[2]])},{i,1,Length[\[CapitalPsi]6]}];
tot=Table[{xVals[[i]],bias[[i]][[2]]+varIn[[i]][[2]]+varData[[i]][[2]]+varNoise[[i]][[2]]},{i,1,Length[\[CapitalPsi]6]}];
If[plotBol,
colors=ColorData[97,"ColorList"];
colors=AppendTo[colors[[1;;4]],Black];
plot=genPlot["pts"->{varIn,varNoise,varData,bias,tot},"logy"->logy,"xaxisLabel"->xaxesLabel,"yaxisLabel"->None,"leg"-> leg,"colors"->colors,"plotRange"-> plotRange,"legPos"->legPos];
Return[{varIn,varNoise,varData,bias,tot,plot}];
];
{varIn,varNoise,varData,bias,tot}
   ];
