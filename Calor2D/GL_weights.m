function [xi, wi] = GL_weights(order)

switch (order)
  case 2
      xi = [-0.57735026918962576451 0.57735026918962576451 ];
      wi = [1.0000000000000000000 1.0000000000000000000 ];
  case 3
      xi = [0 -0.77459666924148337704 0.77459666924148337704 ];
      wi = [0.88888888888888888889 0.55555555555555555556 0.55555555555555555556 ];
  case 4
      xi = [-0.33998104358485626480 0.33998104358485626480 -0.86113631159405257522 0.86113631159405257522 ];
      wi = [0.65214515486254614263 0.65214515486254614263 0.34785484513745385737 0.34785484513745385737 ];
  case 5
      xi = [0 -0.53846931010568309104 0.53846931010568309104 -0.90617984593866399280 0.90617984593866399280 ];
      wi = [0.56888888888888888889 0.47862867049936646804 0.47862867049936646804 0.23692688505618908751 0.23692688505618908751 ];
  case 6
      xi = [0.66120938646626451366 -0.66120938646626451366 -0.23861918608319690863 0.23861918608319690863 -0.93246951420315202781 0.93246951420315202781 ];
      wi = [0.36076157304813860757 0.36076157304813860757 0.46791393457269104739 0.46791393457269104739 0.17132449237917034504 0.17132449237917034504 ];
  otherwise
      error ("invalid order value!"); 
end

return
%%
%% WOLFRAM ALPHA CODE TO GENERATE THIS FUNCTION
%% ADAPTED FROM
%% https://pomax.github.io/bezierinfo/legendre-gauss.html
%%
%$NumberMarks = False
%symboliclegendre[n_, x_] := Solve[LegendreP[n, x] == 0];
%legendreprime[n_, a_] := D[LegendreP[n, x], x] /. x -> a;
%weights[n_, x_] := 2/((1 - x^2) legendreprime[n, x]^2);
%
%(*how many terms should be generated*)
%h = 6;
%
%(*what numerical precision is desired?*)
%precision = 20;
%
%str = OpenWrite["lgvalues.txt", NumberMarks -> False];
%WriteString[str, "abcissae & weights\n\n"];
%WriteString[str, "switch (Order)\n"];
%Do[(*Write[str];Write[str,n];Write[str];*)
%  WriteString[str, "  case ", n, "\n  "];
%  nlist = symboliclegendre[n, x];
%  xnlist = Re[x /. nlist];
%  WriteString[str, "    xi = ["];
%  Do[
%   WriteString[str, N[Part[xnlist, i], precision], " "], {i, 
%    Length[xnlist]}];
%  WriteString[str, "];\n      wi = ["];
%  Do[WriteString[str, N[weights[n, Part[xnlist, i]], precision], 
%    " "], {i, Length[xnlist]}];
%  WriteString[str, "];\n"];,
%  {n, 2, h}];
%WriteString[str, 
%  "  otherwise\n      error (\"invalid order value!\"); \nendswitch\n\
%"];
%
%Close[str];