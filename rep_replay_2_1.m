paths = IO.QuickLoadFolders("D:\JeffreyTong\OneDrive\Study\BJMU_lab_2\data\dirsToDeal.txt");
path = paths(1);
pathOut = "D:\JeffreyTong\OneDrive\Study\BJMU_lab_2\data\rep_replay_2";
%%
ops = CodeJeff.BuildOps(path);
raw = CodeJeff.BuildRawData(ops);
%%
path = char(path);
backslashIdx = strfind(paths, '\');
nameBase = path(backslashIdx(end)+1: end);
%%
fileName = fullfile(pathOut, "rateMap.mat");
digTrig = raw.DigitalInputs(ops.triggerChan);
spikesAll = SpikeData(ops.nameBase, ops.sorting, (1:16)', []);
spikesHippo = spikesAll.SelectTetrodes(ops.hippoTetrodes);
spikesHippo = spikesHippo.TypeIndex('pyr', 'good');
spikes = cell(size(raw.segments, 1), 1);
[spikes{:}] = spikesHippo.SplitSegments(raw.segments);
%%
tracksForc = CodeJeff.LoadTracks(ops, 'forc');
%%
[rateMapLI, rateMapLO, rateMapRO, rateMapRI, tplIdxLI, tplIdxLO, ...
    tplIdxRO, tplIdxRI] = TpltForc(ops, raw, digTrig, spikesHippo, ...
    spikes, tracksForc);
save(fileName, "rateMapLI", "rateMapLO", "rateMapRO", "rateMapRI", ...
    "tplIdxLI", "tplIdxLO", "tplIdxRO", "tplIdxRI");
%%
function [rateMapLI, rateMapLO, rateMapRO, rateMapRI, tplIdxLI, tplIdxLO, tplIdxRO, tplIdxRI] = TpltForc(ops, rawData, digTrig, spikesAll, spikes, tracksForc)
% initialize occupancy distribution & spike distribution
nUnits = size(spikesAll.units, 1);
occMapCell = cell(ops.numRegion, 1);
spkMapCell = cell(ops.numRegion, 1);
for i = 1: ops.numRegion
    numBinRegion = ops.armBinNum(i);
    spkMapCell{i} = zeros(numBinRegion, nUnits);
end
% accumulate occupancy distribution & spike distribution
numForced = size(ops.forcTrackFile, 1);
for i = 1: numForced
    % SpikeData segment
    spikeTmp = spikes{ops.forcSeg(i)};
    % correct by digitalin
    spikeTmp = spikeTmp.CorrectTrigger(digTrig(ops.forcTrigOrder(i)));
    % TrackJeff segment
    trackTmp = tracksForc{i};
    % trim whole track
    trackTmp = trackTmp.TrimTrack( ...
        rawData, ops.forcSeg(i), digTrig(ops.forcTrigOrder(i)));
    % compute length of maze arms (in pixel)
    armLengthsTmp = trackTmp.ComputeArmLengths;
    % spike histcounts
    spikeCountTmp = spikeTmp.SpikeHistcountFrame( ...
        trackTmp.frameSeries, ops.frameRate, ops.ampSampleRate);
    % speed threshold
    SpeedIdxTmp = trackTmp.SpeedIndex(ops.speedThreshold);
    % bin frame for each region
    % (double checked) it is ok to use trackProj here
    trackProjReorderTmp = [trackTmp.trackProj(:, [1 2 4 1 3 5])];
    for iRegion = 1: ops.numRegion
        numBinRegion = ops.armBinNum(iRegion);
        trackProjRegion = trackProjReorderTmp(:, iRegion);
        trackIdxRegion = trackTmp.regionIndex(:, iRegion) & SpeedIdxTmp;
        armLengthRegion = armLengthsTmp(iRegion);
        [occMapRegion, frameMapRegion] = trackTmp.BinFrame1D( ...
            trackProjRegion, trackIdxRegion, armLengthRegion, ...
            ops.armBinNum(iRegion));
        if isempty(occMapCell{iRegion})
            occMapCell{iRegion} = occMapRegion;
        else
            occMapCell{iRegion} = occMapCell{iRegion} + occMapRegion;
        end
        for iUnit = 1: nUnits
            spkCountUnit = spikeCountTmp(:, iUnit);
            for iBin = 1: numBinRegion
                % notice: spkCountUnit is not corrected
                numSpkBin = sum(spkCountUnit(frameMapRegion{iBin} + ...
                    ops.forcTrigFrame(i)));
                spkMapCell{iRegion}(iBin, iUnit) = ...
                    spkMapCell{iRegion}(iBin, iUnit) + numSpkBin;
            end
        end
    end
end
% concatenate occupancy distribution map
occMapLI = occMapCell{3};
occMapLO = [occMapCell{1}; occMapCell{2}];
occMapRO = [occMapCell{4}; occMapCell{5}];
occMapRI = occMapCell{6};
% concatenate spike distribution map
spkMapLI = spkMapCell{3};
spkMapLO = [spkMapCell{1}; spkMapCell{2}];
spkMapRO = [spkMapCell{4}; spkMapCell{5}];
spkMapRI = spkMapCell{6};
% calculate firing rate map - 1D
rateMapLI = ComputeRateMap(occMapLI, spkMapLI, ops.frameRate);
rateMapLO = ComputeRateMap(occMapLO, spkMapLO, ops.frameRate);
rateMapRO = ComputeRateMap(occMapRO, spkMapRO, ops.frameRate);
rateMapRI = ComputeRateMap(occMapRI, spkMapRI, ops.frameRate);
% calculate spatial info
infoLI = PlaceCell.CalculateSpatialInfo1D(rateMapLI, occMapLI);
infoLO = PlaceCell.CalculateSpatialInfo1D(rateMapLO, occMapLO);
infoRO = PlaceCell.CalculateSpatialInfo1D(rateMapRO, occMapRO);
infoRI = PlaceCell.CalculateSpatialInfo1D(rateMapRI, occMapRI);
infoIdxLI = infoLI.perSpike >= 0.25;
infoIdxLO = infoLO.perSpike >= 0.25;
infoIdxRO = infoRO.perSpike >= 0.25;
infoIdxRI = infoRI.perSpike >= 0.25;
% construct template for decoding
rateMapCat = [flip(rateMapLI, 1); flip(rateMapLO, 1); rateMapRO; rateMapRI];
[tplIdxLI, tplIdxLO, tplIdxRO, tplIdxRI] = ...
    CodeJeff.TemplateIndex(spikesAll, rateMapCat);
% only pyramidal cells -- [ADDITIONAL]
[~, pyIdx] = spikesAll.TypeIndex('pyr', 'good');
tplIdxLI = tplIdxLI & pyIdx & infoIdxLI;
tplIdxLO = tplIdxLO & pyIdx & infoIdxLO;
tplIdxRO = tplIdxRO & pyIdx & infoIdxRO;
tplIdxRI = tplIdxRI & pyIdx & infoIdxRI;
end
%%
function rateMap = ComputeRateMap(occMap, spkMap, frameRate)
rateMap = spkMap./occMap*frameRate;
% sometime there are missing data
rateMap = fillmissing(rateMap, 'spline', 1);
% smooth
% radius = 20;
sigma = 2;
radius = sigma * 6;
for i = 1: size(rateMap, 2)
    rateMap(:, i) = Plotting.GaussianFilter1D(rateMap(:, i), radius, sigma);
end
end