% Chi1.mp -- do not edit, generated automatically by "2024-04-04 Diagramm".tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 220*1.0pt#, 40*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vtop (__a);
vbottom (__b);
vright (__c);
vconnect ("plain,left=0.1,label=$\mathcal {G}(\tau )$", __a, __c);
vconnect ("plain,left=0.1,label=$\mathcal {G}(\beta -\tau )$", __c, __b);
vdot (__a, __c);
vvertex ("lab=$\lambda _1$,lab.dist=0.05w", __a);
vvertex ("lab=$\lambda _1$,lab.dist=0.05w", __c);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
