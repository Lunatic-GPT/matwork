function reorder_TGSE_vfl_pss(  m_lEPIFactor,lSliceTurboFactor,lLinesToMeasureFullFT,lPartitionsToMeasure,TE,lEchoSpacing)

lDiscardedEchoes=0;

lMaxTrainLength_1Par = floor((lLinesToMeasureFullFT+m_lEPIFactor-1)/m_lEPIFactor);


lMaxTrainLengthBeforeKy0_1Par = floor((lLinesToMeasureFullFT-1)/2)-lMaxTrainLength_1Par*(m_lEPIFactor-1)/2; % need -1 because LinNo starts from lLinesToMeasureFullFT -1

lDeltaSkippedLines = (m_lEPIFactor+1)/2;  % the step number of lines to skip for the center echo segment

lNoValidTEs=floor(0.75*lMaxTrainLengthBeforeKy0_1Par/lDeltaSkippedLines);

m_lNoOfEchoSpacing4TE=floor((TE+lEchoSpacing/2)/lEchoSpacing)-lDiscardedEchoes;

% ok if lWantedNoOfEchoSpacingBefore Ky0 is not a multiple of lSliceTurboFactor; it is not final.

lApproxNoOfEchoSpacingBeforeKy0 = m_lNoOfEchoSpacing4TE - mod(floor(lPartitionsToMeasure/2),lSliceTurboFactor)-1;

%lSkippedLines is not the skipped line for partial Fourier but rather the number of leading echos skipped for shortenning TE

lSkippedLines = lMaxTrainLengthBeforeKy0_1Par - floor(lApproxNoOfEchoSpacingBeforeKy0/lSliceTurboFactor);

lSkippedLines = floor(lSkippedLines/lDeltaSkippedLines+0.5)*lDeltaSkippedLines;


if (lSkippedLines>(lNoValidTEs-1)*lDeltaSkippedLines)
    lSkippedLines = (lNoValidTEs-1)*lDeltaSkippedLines;
end
if (lSkippedLines<0)
    lSkippedLines = 0;
end
lTrainLengthBeforeKy0_1Par=lMaxTrainLengthBeforeKy0_1Par - lSkippedLines;

lSkippedLinesForPartialFT = lSkippedLines/lDeltaSkippedLines*m_lEPIFactor;

fprintf('TE = %f, lSkippedLines = %d; \nTrainLengthBeforeKy0_1Par = %d, SkippedLinesForPartialFT = %d\n',TE,lSkippedLines,lTrainLengthBeforeKy0_1Par,lSkippedLinesForPartialFT);
