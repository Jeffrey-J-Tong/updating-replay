pathOut = "D:\JeffreyTong\OneDrive\Study\BJMU_lab_2\data\rep_replay_2";
fileName = fullfile(pathOut, "_pbe_all.mat");
if ~exist(fileName, "file")
    pbe = DetectAllPbe(ops, raw, spikesHippo);
    save(fileName, 'pbe');
else
    load(fileName, 'pbe');
end
%%
fileName = fullfile(pathOut, "_decode.mat");
if ~exist(fileName, "file")
    [r, p, prob, maxPos, timePt] = DecodeAllPbe(ops, spikesHippo, pbe, ...
        rateMapLI, rateMapLO, rateMapRO, rateMapRI, ...
        tplIdxLI, tplIdxLO, tplIdxRO, tplIdxRI);
    save(fileName, 'r', 'p', 'prob', 'maxPos', 'timePt');
end
%%
function pbe = DetectAllPbe(ops, raw, spikes)
% detect pbe
mua = spikes.MUA;
dataRange = raw.DataRange([]);
pbe = PlaceCell.DetectPBE(mua, ops.ampSampleRate, dataRange, "threshold", [0.2 0.35]);
end
%%
function [r, p, prob, maxPos, timePt] = DecodeAllPbe(ops, spikes, pbe, rateMapLI, rateMapLO, rateMapRO, rateMapRI, tplIdxLI, tplIdxLO, tplIdxRO, tplIdxRI)
nPbe = size(pbe, 1);
% decode (parallel)
% params of decoder
bayesWin = 0.02;  % 20 ms
bayesStep = 0.01;  % 10 ms
posBinLI = ((1: 39) - 0.5)' / ops.binLength;
posBinLO = ((1: 55) - 0.5)' / ops.binLength;
posBinRO = ((1: 55) - 0.5)' / ops.binLength;
posBinRI = ((1: 39) - 0.5)' / ops.binLength;
% initialize outputs
r = zeros(nPbe, 4);
p = zeros(nPbe, 4);
maxPos = cell(nPbe, 1);
prob = cell(nPbe, 1);
timePt = cell(nPbe, 1);
% decode
disp('Performing Bayesian decoding.')
tic
parfor i = 1: nPbe
    % prepare loop variables
    probTmp = cell(1, 4);
    maxPosTmp = cell(1, 4);
    timePtTmp = cell(1, 4);
    rTmp = zeros(1, 4);
    pTmp = zeros(1, 4);
    pbeTimeRange = pbe(i, :)/ops.ampSampleRate;
    % decode based on LI templates
    [~, probTmp{1}, maxPosTmp{1}, timePtTmp{1}] = PlaceCell.BayesianDecode_v2(...
        spikes, posBinLI, rateMapLI, pbeTimeRange, ...
        bayesWin, bayesStep, "unitIdx", tplIdxLI);
    [~, rTmp(1), pTmp(1)] = PlaceCell.ReplayTest( ...
        timePtTmp{1}, maxPosTmp{1});
    % decode based on LO templates
    [~, probTmp{2}, maxPosTmp{2}, timePtTmp{2}] = PlaceCell.BayesianDecode_v2(...
        spikes, posBinLO, rateMapLO, pbeTimeRange, ...
        bayesWin, bayesStep, "unitIdx", tplIdxLO);
    [~, rTmp(2), pTmp(2)] = PlaceCell.ReplayTest( ...
        timePtTmp{2}, maxPosTmp{2});
    % decode based on RI templates
    [~, probTmp{3}, maxPosTmp{3}, timePtTmp{3}] = PlaceCell.BayesianDecode_v2(...
        spikes, posBinRO, rateMapRO, pbeTimeRange, ...
        bayesWin, bayesStep, "unitIdx", tplIdxRO);
    [~, rTmp(3), pTmp(3)] = PlaceCell.ReplayTest( ...
        timePtTmp{3}, maxPosTmp{3});
    % decode based on RO templates
    [~, probTmp{4}, maxPosTmp{4}, timePtTmp{4}] = PlaceCell.BayesianDecode_v2(...
        spikes, posBinRI, rateMapRI, pbeTimeRange, ...
        bayesWin, bayesStep, "unitIdx", tplIdxRI);
    [~, rTmp(4), pTmp(4)] = PlaceCell.ReplayTest( ...
        timePtTmp{4}, maxPosTmp{4});
    % output
    maxPos{i} = maxPosTmp;
    prob{i} = probTmp;
    timePt{i} = timePtTmp;
    r(i, :) = rTmp;
    p(i, :) = pTmp;
end
elapsedTime = toc;
disp(['Decoding takes ' num2str(elapsedTime) ' seconds.'])
end