% self_energy.mp -- do not edit, generated automatically by "2024-04-04 Diagramm".tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 80*1.0pt#, 120*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__a1);
vright (__b1);
vconnect ("fermion", __a1, __a2, __);
vconnect ("plain, left=1, tension=0.4", __a2, __b2);
vconnect ("plain, right=1, tension=0.4", __a2, __b2);
vconnect ("double_arrow", __b2, __b1);
vconnect ("phantom,label.dist=0,label=$\sum $", __a2, __b2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
