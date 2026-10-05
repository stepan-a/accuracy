function s = tex_name(s)

% Escapes a name (variable, equation) for use in LaTeX text.

s = strrep(s, '_', '\_');
