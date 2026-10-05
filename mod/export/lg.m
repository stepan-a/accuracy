function s = lg(x)

% Base ten logarithm of |x|, formatted with two decimals for a table.

if ~isfinite(x)
    s = '--';
elseif x==0
    s = '$-\infty$';
else
    s = sprintf('$%.2f$', log10(abs(x)));
end
