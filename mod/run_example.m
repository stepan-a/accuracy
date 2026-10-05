function run_example(name)

% Runs one example with the accuracy branch of Dynare, marks it as done, and
% records the commit of the branch that produced the results.
%
% INPUTS
% - name        [char]      name of the example: a .mod file of this directory,
%                           or 'irr' for the sequence of runs compared in
%                           compare_irr.mod
%
% REMARKS
% The location of Dynare is read from the environment variable DYNARE_ROOT,
% set by the Makefile. DYNARE_BUILD_DIR, when set, is honoured by Dynare itself.
% The commit of DYNARE_ROOT is written to results/dynare.hash, which is versioned
% with the results and shown in the footer of the slides.

root = getenv('DYNARE_ROOT');
if isempty(root)
    error('run_example: the environment variable DYNARE_ROOT is not set.')
end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'matlab'));
addpath(fullfile(here, 'export'));
cd(here);

if ~exist('results', 'dir')
    mkdir('results');
end

if strcmp(name, 'irr')
    % The contexts saved by the first three runs are loaded by the last one,
    % which also needs their +<fname> packages, hence a single session.
    runs = {'occbin_irr', 'ep_irr', 'sep_irr', 'sep_irr_h2', 'sep_irr_h4', 'compare_irr', 'sweep_irr'};
elseif strcmp(name, 'bs')
    runs = {'occbin_bs', 'ep_bs', 'sep_bs', 'sep_bs_h2', 'sep_bs_h4', 'sep_bs_o2', 'sep_bs_o3', 'compare_bs'};
else
    runs = {name};
end

for i=1:length(runs)
    dynare(runs{i});
end

[status, hash] = system(sprintf('git -C "%s" rev-parse --short HEAD', root));
if status || isempty(strtrim(hash))
    hash = 'unknown';
end
fid = fopen(fullfile('results', 'dynare.hash'), 'w');
fprintf(fid, '%s\n', strtrim(hash));
fclose(fid);

fclose(fopen(fullfile('results', [name '.stamp']), 'w'));
