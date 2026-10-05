function write_tabular(filename, colspec, header, rows)

% Writes a booktabs tabular.
%
% INPUTS
% - filename    [char]      output file
% - colspec     [char]      column specification, e.g. 'lrr'
% - header      [cell]      column titles (LaTeX)
% - rows        [cell]      one cell of strings (LaTeX) per row

fid = fopen(filename, 'w');
fprintf(fid, '\\begin{tabular}{%s}\n\\toprule\n', colspec);
fprintf(fid, '%s \\\\\n\\midrule\n', strjoin(header, ' & '));
for i=1:length(rows)
    fprintf(fid, '%s \\\\\n', strjoin(rows{i}, ' & '));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
