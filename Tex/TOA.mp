% TOA.mp -- do not edit, generated automatically by "2024-04-04 Diagramm".tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 110*1.0pt#, 30*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vtop (__q, __w, __d2, __e2, __r);
vbottom (__a, __b, __c, __d, __e);
vconnect ("plain,tension=1", __a, __b);
vconnect ("fermion,tension=1", __b, __c);
vconnect ("plain,tension=1", __c, __d);
vconnect ("plain,tension=1", __d, __e);
vconnect ("plain,tension=0", __d, __e2);
vconnect ("plain,tension=0", __e2, __e);
vconnect ("boson,tension=0,left=0.5", __a, __e2);
vvertex ("lab=$\psi _a$,lab.dist=0.05w", __a);
vvertex ("lab=$\psi _b$,lab.dist=0.05w", __c);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
