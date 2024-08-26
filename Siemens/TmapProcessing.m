function [TMap,drifts,imagedata,NoiseMask]=TmapProcessing(imagedata,betaTgamma,TE,NoiseMask,drifts)
%imagedata is a X * Y * t array, for this example it is for a single, 2D
%slice only; angles in radians
%betaTgamma is the imaging frequency, or the field strength * gamma (gamma=42.58 MHz/T)
%TE is the echo time in seconds 
%noisemask, optional, is and X * Y array to show the areas in which you
%want to calculate the temperature
%drifts is a t * 1 array showing the drift in each successive image. You can use this to re-use previously calculated drifts

imagedata=-imagedata;

dodriftcorrect=1; %0 is no, 1 is yes, 2 is use loaded values
driftrejectthres=inf; %if drift is > than this in magnitude, reject slice. Set to inf to ignore this. Requires dodriftcorrect==1
domanualmask=1;

Tbase=37;
NumImg=size(imagedata,3);
alpha=-0.0101;
Kfact=1/(2*pi*alpha*betaTgamma*TE);
%Unwrapping limits
LowerUnwrappingLimit=-270/180*pi;
UpperUnwrappingLimit=270/180*pi;

TMap=zeros(size(imagedata));

if ~exist('NoiseMask','var')
    NoiseMask=MaskNoise(imagedata(:,:,1));
    NoiseMask=NoiseMask & MaskNoise(imagedata(:,:,2)); %mask any zeros on the first two slices.
end

if domanualmask
    manmask=manualmask(imagedata,NoiseMask);
else
    manmask=ones(size(imagedata(:,:,1)));
end

if ~exist('drifts','var')
    if dodriftcorrect==1
        driftinds=driftcorrect(imagedata,NoiseMask);
        drifts=zeros(NumImg,1);
    end
end

%Going through the regular images
for j=1:NumImg
    if(j==1)
        TMap(:,:,j)=Tbase;
    else
        PhPrev=imagedata(:,:,j-1);
        PhCurr=imagedata(:,:,j);
        PhDiff=(PhPrev-PhCurr).*NoiseMask;
        PhDiff=UnwrapPhase(PhDiff,LowerUnwrappingLimit,UpperUnwrappingLimit,2*pi);
        TMap(:,:,j)=TMap(:,:,j-1)+PhDiff*Kfact;
        
        
        if dodriftcorrect==1
            drift=mean(TMap(driftinds+(size(imagedata,1)*size(imagedata,2)*(j-1))))-Tbase;
            if abs(drift)<driftrejectthres
                TMap(:,:,j)=TMap(:,:,j)-(drift);
                drifts(j)=drift;
            else
                drifts(j)=NaN;
                TMap(:,:,j)=TMap(:,:,j-1);
            end
        elseif dodriftcorrect==2
            TMap(:,:,j)=TMap(:,:,j)-drifts(j);
        end
    end;

end;
if ~dodriftcorrect
    drifts=NaN;
end

for j=1:NumImg
    TMap(:,:,j)=TMap(:,:,j).*manmask;
end
TMap(TMap==0)=Tbase;

fprintf('Finished.\n');
return;

function PhDiffIm=UnwrapPhase(PhDiffIm,LowerUnwrappingLimit,UpperUnwrappingLimit,Period)

count=0;
while(~isempty(find(PhDiffIm<LowerUnwrappingLimit, 1)) && ~isempty(find(PhDiffIm>UpperUnwrappingLimit, 1)))
    count=count+1;
    PhDiffIm(PhDiffIm<LowerUnwrappingLimit)=PhDiffIm(PhDiffIm<LowerUnwrappingLimit)+Period;
    PhDiffIm(PhDiffIm>UpperUnwrappingLimit)=PhDiffIm(PhDiffIm>UpperUnwrappingLimit)-Period;
    if(count>1)
        fprintf('Existential crisis: why are you here? \n');
        keyboard;
    end;
end;

return

function mask=MaskNoise(refimage)
    mask=(refimage==0);
    mask=~mask;
return

function inds=driftcorrect(im,mask)
    figure(1);
    imagesc(squeeze(im(:,:,round(size(im,3)/2))).*mask);
    axis equal
    axis off
    colormap gray
    title('Select an unheated region for drift correction.');
    h = imrect;
    position = wait(h);
    position = floor(position); %perform appropriate pixel rounding
    delete(h);
    if position(3) < 3 || position(4) <3
        error('ROI must be at least 3x3.')
    end
    driftmask=zeros(size(mask));
    driftmask(position(2):position(2)+position(4)-1,position(1):position(1)+position(3)-1)=1;
    inds=find(driftmask==1);
    close all;
return

function manmask=manualmask(im,mask)
    figure(1);
    imagesc(squeeze(im(:,:,round(size(im,3)/2))).*mask);
    axis equal
    axis off
    colormap gray
    title('Select region to KEEP.');
    manmask = roipoly;
    close 1;
return