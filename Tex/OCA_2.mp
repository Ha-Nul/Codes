% OCA_2.mp -- do not edit, generated automatically by "2024-04-04 Diagramm".tex
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
vconnect ("fermion", __a, __a1');
vconnect ("fermion", __a1', __a2');
vconnect ("fermion", __a2', __b1');
vconnect ("fermion", __b1', __b2');
vconnect ("fermion", __b2', __a1);
vconnect ("fermion", __a1, __b1);
vconnect ("fermion", __b1, __a1'');
vconnect ("fermion", __a1'', __a2'');
vconnect ("fermion", __a2'', __b1'');
vconnect ("fermion", __b1'', __b2'');
vconnect ("fermion", __b2'', __b);
vconnect ("boson,left=0.7,tension=0.1", __a, __b1);
vconnect ("boson,left=0.9,tension=0.05", __a1', __b1');
vconnect ("boson,left=0.9,tension=0.05", __a2', __b2');
vconnect ("boson,left=0.7,tension=0.1", __a1, __b);
vconnect ("boson,left=0.9,tension=0.05", __a1'', __b1'');
vconnect ("boson,left=0.9,tension=0.05", __a2'', __b2'');
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
