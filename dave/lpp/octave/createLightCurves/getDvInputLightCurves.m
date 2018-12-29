function [ output ] = getDvInputLightCurves( dvDir )
%Return the PA and PDC light curves in an output struct.
%These are pulled from the DV inputs struct.
%

    clear inputsStruct;
    dvfile=[dvDir,'/dv-inputs-0.mat'];
    load(dvfile);
    pstruct=inputsStruct.targetStruct;
    
    output=struct([]);
    
    output(1).pdcTimeSeries=pstruct.correctedFluxTimeSeries.values;
    output.pdcGaps=pstruct.correctedFluxTimeSeries.gapIndicators;
    output.pdcFilled=pstruct.correctedFluxTimeSeries.filledIndices;
    
    output.paTimeSeries=pstruct.rawFluxTimeSeries.values;
    output.paGaps=pstruct.rawFluxTimeSeries.gapIndicators;
    %output.harmremTimeSeries=pstruct.harmonicFreeCorrectedFluxTimeSeries.values;
   % output.harmremGaps=pstruct.harmonicFreeCorrectedFluxTimeSeries.gapIndicators;
    %Previous commented out for transit inversion
    output.quarters=inputsStruct.dvCadenceTimes.quarters;
    
    output.kic=pstruct.keplerId;
    output.effTemp=pstruct.effectiveTemp;
    output.metal=pstruct.log10Metallicity;
    output.logg=pstruct.log10SurfaceGravity;

end

