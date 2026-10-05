function write_macros(filename, macros)

% Writes LaTeX macros holding the numbers quoted in the slides.
%
% INPUTS
% - filename    [char]      output file
% - macros      [struct]    field names are the macro names (letters only),
%                           values are the LaTeX strings

fid = fopen(filename, 'w');
names = fieldnames(macros);
for i=1:length(names)
    fprintf(fid, '\\newcommand{\\%s}{%s}\n', names{i}, macros.(names{i}));
end
fclose(fid);
