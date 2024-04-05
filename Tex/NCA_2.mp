% NCA_2.mp -- do not edit, generated automatically by "2024-04-04 Diagramm".tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 80*1.0pt#, 100*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__a);
vright (__b);
vconnect ("fermion", __a, __a1);
vconnect ("boson,left=0.7", __a, __b);
vconnect ("boson,left=0.7", __a1, __b1);
vconnect ("fermion", __b1, __b);
vconnect ("fermion", __a, __b);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
