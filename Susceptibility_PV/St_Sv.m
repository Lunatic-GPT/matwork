function [St0,Sv0]=St_Sv(TE,dchi)
% TE in ms;
% dchi in ppm
% 
        St0=100;
        T2s_v = dchi2T2(dchi*1e-6)*1000;
        T2s_t = 27; % for WM;
        Sv0=St0.*exp(-TE/T2s_v).*exp(TE/T2s_t)*1.05;% 1.05 is the blood-WM partition coefficient
       
        



