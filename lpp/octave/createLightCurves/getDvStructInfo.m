function [ output ] = getDvStructInfo( dvDir,kic,pn)
% Input: dvDir= directory location of the dv-outputs-0.mat file for your
% tce.
%  given a DVdirrectory get out the importantinformation and store it ina
%struct.
%   Susan Thompson December 2014



    dvfile=[dvDir,'/dv-outputs-0.mat'];
    load(dvfile);
    pstruct=outputsStruct.targetResultsStruct.planetResultsStruct;
     
    output=struct([]);
    output(1).period=pstruct(pn).allTransitsFit.modelParameters(11).value;
    output(1).kic=kic;
    output(1).pn=pn;
    output.phase=pstruct(pn).allTransitsFit.modelParameters(1).value;
    output.pulse=pstruct(pn).planetCandidate.trialTransitPulseDuration;
    output.dur=pstruct(pn).allTransitsFit.modelParameters(8).value;
    output.snr=pstruct(pn).allTransitsFit.modelFitSnr;
    output.maxmes=pstruct(pn).planetCandidate.maxMultipleEventSigma;
    output.maxses=pstruct(pn).planetCandidate.maxSingleEventSigma;
    output.wsmaxmes=pstruct(pn).planetCandidate.weakSecondaryStruct.maxMes;
    output.wsmaxphase=pstruct(pn).planetCandidate.weakSecondaryStruct.maxMesPhaseInDays;
    output.wsminmes=pstruct(pn).planetCandidate.weakSecondaryStruct.minMes;
    output.wsminphase=pstruct(pn).planetCandidate.weakSecondaryStruct.minMesPhaseInDays;
    
    output.whiteFlux=pstruct(pn).whitenedFluxTimeSeries.values;
    output.gapsind=pstruct(pn).whitenedFluxTimeSeries.gapIndicators;
    output.bjdTime=outputsStruct.targetResultsStruct.barycentricCorrectedTimestamps;
        
    

end

