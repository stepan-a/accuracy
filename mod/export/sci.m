function s = sci(x, digits)

% Formats a number in scientific notation for LaTeX math mode, e.g. 3.99\times10^{-3}.

if nargin<2
    digits = 2;
end

if x==0
    s = '0';
    return
end
if ~isfinite(x)
    s = '\mathrm{NaN}';
    return
end

e = floor(log10(abs(x)));
m = x/10^e;
if round(m*10^digits)==10^(digits+1)
    m = m/10;
    e = e+1;
end
if e==0
    s = sprintf('%.*f', digits, m);
else
    s = sprintf('%.*f\\times10^{%d}', digits, m, e);
end
