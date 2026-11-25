% MAIN  Compute the four PUF metrics for a challenge-response dataset.
% Place puf_metrics.m on the MATLAB path or in the current folder.

% Point these paths to numeric matrices with observations in rows.
challengeFile = 'challenges.mat'; %or 'challenges.csv'
responseFile  = 'responses.mat'; %or 'responses.csv'

% 'rows'    : row i is one challenge and its complete response.
% 'entries' : C(i,j) and R(i,j) are paired scalar measurements on a grid.
pairing = 'rows';

C = readInput(challengeFile, 'C');
R = readInput(responseFile, 'R');

results = puf_metrics(C, R, pairing);
disp(results);

function X = readInput(filename, variableName)
    [~,~,ext] = fileparts(filename);
    if strcmpi(ext, '.mat')
        data = load(filename, '-mat');
        if ~isfield(data, variableName)
            error('Expected variable "%s" in %s.', variableName, filename);
        end
        X = data.(variableName);
    else
        X = readmatrix(filename);
    end
end
