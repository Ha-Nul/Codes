% propagator.mp -- do not edit, generated automatically by "2024-04-04 Diagramm".tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 50*1.0pt#, 120*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__a);
vright (__b);
vconnect ("double_arrow,label=$\mathcal {G}$", __a, __b);
vvertex ("lab=$\psi _a$,lab.dist=0.05w", __a);
vvertex ("lab=$\psi ^\dagger _b$,lab.dist=0.05w", __b);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
