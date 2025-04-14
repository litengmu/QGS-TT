function H = hessian(x, sigma, lambda, d)
lam.eqnonlin   = lambda(1:d.neqnln);
lam.ineqnonlin = lambda(d.neqnln+(1:d.niqnln));
H = tril(opf_hessfcn(x, lam, sigma, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il));
