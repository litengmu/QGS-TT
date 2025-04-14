function df = gradient(x, d)
    f = 0;
    df = zeros(length(x), 1);  % Zero vector, as the gradient of a constant is zero