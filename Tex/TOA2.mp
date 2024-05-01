% TOA2.mp -- do not edit, generated automatically by "2024-04-04 Diagramm".tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 90*1.0pt#, 100*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vtop (__a);
vbottom (__b, __c);
vconnect ("plain", __b, __b1, __a);
vconnect ("plain", __a, __c1, __c);
vconnect ("boson,tension=0", __b1, __c);
vconnect ("boson,tension=0", __c1, __b);
vdot (__a, __b1, __b, __c1, __c);
vvertex ("lab=N,lab.dist=0.05w", __b);
vvertex ("lab=N,lab.dist=0.05w", __c);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
