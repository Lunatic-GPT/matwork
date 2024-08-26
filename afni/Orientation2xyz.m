function [xyz,direction]=Orientation2xyz(orient)
% convert orient code in BrikInfo to xyz in dicom convention; 
%  x=R-L, y=A-P, and z=I-S (R,A,I are < 0; L,P,S are > 0).
%xyz: 1 is x; 2 is y; and 3 is z
% direction: +1; -1

for i=1:3
    
    switch orient(i,:)
        
        case 'LR'
            xyz(i)=1;
            direction(i)=-1;
        case 'RL'
             xyz(i)=1;
             direction(i)=1;
        case 'AP'
            xyz(i)=2;
            direction(i)=1;
        case 'PA'
            xyz(i)=2;
            direction(i)=-1;
        case 'IS'
            xyz(i)=3;
            direction(i)=1;
        case 'SI'
            xyz(i)=3;
            direction(i)=-1;
    
            
    end
end
