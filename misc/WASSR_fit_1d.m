function foff=WASSR_fit_1d(f,y)
%WASSR_fit(f,d)

            fitfunc=@(x) calcResidual(f,squeeze(y),x);
            
            fmin=min(y);
            
            if y(1)==fmin 
                foff=f(1);
                return;
            elseif y(end)==fmin
                foff=f(end); %|| ind==length(f)-1 
                return;
            end
            
            
            %%{
            [ymin,ind]=min(y);
            maxn=max(y(1:ind-1));
            maxp=max(y(ind+1:end));
            
            max2=min([maxn,maxp]);
            
            yt=y(1:ind);
            ft=f(1:ind);
            
            [yt,ind]=sort(yt);
            ft=ft(ind);
            [yt,ft]=ave_dup(yt,ft);
            
            dy=(max2+ymin)*0.5;
            f1=interp1(yt,ft,dy);
            
            
            yt=y(ind:end);
            ft=f(ind:end);
            
            [yt,ind]=sort(yt);
            ft=ft(ind);
            
            [yt,ft]=ave_dup(yt,ft);
            
            f2=interp1(yt,ft,dy);
            
            finit=(f1+f2)/2;
            %}
            
            %finit=mean(f(ind-1:ind+1));
            
            foff=fminsearch(fitfunc,finit);


function res=calcResidual(f2,y,x)


asym=z2asym_cutoff(f2,y,x);

res=sum(asym.^2)/length(asym);

function [yt2,ft2]=ave_dup(yt,ft)

            ind2=(diff(yt)==0);
            
            ft2=[];
            yt2=[];
            i=1;
            while i<=length(ind2)
                if ~ind2(i)
                  ft2(end+1)=ft(i);
                  yt2(end+1)=yt(i);
                  i=i+1;
                else
                    
                    
                    if ind2(i) && (i==length(ind2) || ~ind2(i+1))
                   ft2(end+1)=(ft(i)+ft(i+1))/2;
                   yt2(end+1)=yt(i);
                   i=i+2;
                 elseif ind2(i) && ind2(i+1) && ( i+1==length(ind2) ||~ind2(i+2))
                        ft2(end+1)=(ft(i)+ft(i+1)+ft(i+2))/3;
                   yt2(end+1)=yt(i);
                     i=i+3;
                 elseif ind2(i) && ind2(i+1) && ind2(i+2) && ( i+2==length(ind2) ||~ind2(i+3))
                        ft2(end+1)=(ft(i)+ft(i+1)+ft(i+2)+ft(i+3))/4;
                   yt2(end+1)=yt(i);
                     i=i+4;
                else
                    disp('more than 3 equal values found');
                end
                
            end
            
            if ~ind2(end)
                
                ft2(end+1)=ft(end);
                yt2(end+1)=yt(end);
            end
            
            
            
    
    



