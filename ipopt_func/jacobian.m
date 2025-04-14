function J = jacobian(x, d)
[hn, gn, dhn, dgn] = opf_consfcn(x, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il);
J = [dgn'; dhn'; d.A];
