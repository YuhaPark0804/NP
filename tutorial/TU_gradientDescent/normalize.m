function x_norm = normalize(x)
    maxs = max(x);
    mins = min(x);
    x_norm = (x-max(x)) ./ (max(x)-min(x));
end