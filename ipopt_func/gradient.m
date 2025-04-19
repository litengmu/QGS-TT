function df = gradient(x, d)

% optimality
%[f, df] = opf_costfcn(x, d.om);
% feasibility
f = 0;
df = zeros(length(x), 1);  % Zero vector, as the gradient of a constant is zero