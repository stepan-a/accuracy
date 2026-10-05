function write_csv(filename, header, data)

% Writes a matrix as a CSV file with a header line, to be read by pgfplots.
%
% INPUTS
% - filename    [char]      output file
% - header      [cell]      column names
% - data        [double]    one column per name; non finite values are written as nan

fid = fopen(filename, 'w');
fprintf(fid, '%s\n', strjoin(header, ','));
fmt = [strjoin(repmat({'%.10g'}, 1, size(data,2)), ','), '\n'];
fprintf(fid, fmt, transpose(data));
fclose(fid);
