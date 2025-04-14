function c = constraints(x, d)
[hn, gn] = opf_consfcn(x, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il);
if isempty(d.A)
    c = [gn; hn];
else
    c = [gn; hn; d.A*x];
end