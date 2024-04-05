% OCA.mp -- do not edit, generated automatically by "2024-04-04 Diagramm".tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 150*1.0pt#, 100*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__a);
vright (__b);
vconnect ("double_arrow,label=$\mathcal {G}$", __a, __a1);
vconnect ("double_arrow,label=$\mathcal {G}$", __a1, __b1);
vconnect ("double_arrow,label=$\mathcal {G}$", __b1, __b);
vconnect ("boson,left=0.7,tension=0.1", __a, __b1);
vconnect ("boson,left=0.7,tension=0.1", __a1, __b);
vvertex ("lab=$\psi _a$,lab.dist=0.05w", __a);
vvertex ("lab=$\psi ^\dagger _d$,lab.dist=0.05w", __b);
vvertex ("lab=$\psi _b$,lab.angle=120,lab.dist=0.05h", __a1);
vvertex ("lab=$\psi ^\dagger _c$,lab.angle=60,lab.dist=0.05h", __b1);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
