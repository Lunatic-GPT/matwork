function mnew=sudokuSolver(m)
% mnew=sudokuSolver(m)
% m is a 9 by 9 matrix with 0 in blank spots.
if isa(m,'char')
    m = load(m);
end
disp(m);
mnew = solver1(m);
disp('First pass');
fprintf('%d positions solved\n', length(find(mnew>0))-length(find(m>0)));

disp('Final pass...');
m_final = solver2(mnew);
disp('Solution:');
disp(m_final);

function m = solver2(m)
%% find positions with unique solution.
 ind = find(m==0);
 for i=1:length(ind)
     [x,y]=ind2sub([9,9],ind(i));
     s{i}= find_possible_values(m,x,y);
 end
 
 [m,status]=solve_loop(m,ind,s);
 
disp(status);
 
    function [m_new,status]=solve_loop(m,ind,s)
        if isempty(s) 
            m_new = m;
            status = 0;
            disp('dont reach here');
            return;
        end
        
        [x,y]=ind2sub([9,9],ind(1));    
           m_tmp = m;
           m_tmp(x,y) = s{1}(1);
            if check_solution(m_tmp)
                [m_new,status]=solve_loop(m_tmp,ind(2:end),s(2:end));
                if status <0 && length(s{1})>1
                  s{1}=s{1}(2:end);
                  [m_new,status] =  solve_loop(m,ind,s);
                   
                end
            elseif length(s{1})>1
               s{1}=s{1}(2:end);
               [m_new,status] =  solve_loop(m,ind,s);
            else
                status = -1;
                m_new = m;
                %error('should never reach here');
            end
    function [m_new,status]=solve_loop_old(m,ind,s)
        if isempty(s) 
            m_new = m;
            status = 0;
            disp('dont reach here');
            return;
        end
        
        [x,y]=ind2sub([9,9],ind(1));    
           m_tmp = m;
           m_tmp(x,y) = s{1}(1);
            if check_solution(m)
               [m_new,status]=solve_loop(m_tmp,ind(2:end),s(2:end));
               if status <0 && length(s{1})>1
                  s{1}=s{1}(2:end);
                  [m_new,status] =  solve_loop(m,ind,s);
                   
               end
            elseif length(s{1})>1
               s{1}=s{1}(2:end);
               [m_new,status] =  solve_loop(m,ind,s);
            else
                status = -1;
                m_new = m;
                %error('should never reach here');
            end
        
        function result = check_solution(m)
            
            result = true;
            for i=1:9
                a = m(:,i);
                a=a(a>0);
                if length(a)>length(unique(a))
                    result = false;
                end
                
                a = m(i,:);
                a=a(a>0);
                if length(a)>length(unique(a))
                    result = false;
                end
                
                
            end
             for i=1:3
                 for j=1:3
                   a = m(i*3-2:i*3,j*3-2:j*3);
                   a=a(:);
                   a=a(a>0);
                   if length(a)>length(unique(a))
                     result = false;
                    end
                 end
             end
             
            
        
function m = solver1(m)
%% find positions with unique solution.
 ind = find(m==0);
 for i=1:length(ind)
     [x,y]=ind2sub([9,9],ind(i));
     s{i}= find_possible_values(m,x,y);
 end
 
 for i=1:length(ind)
     
     [x,y]=ind2sub([9,9],ind(i));
     x2=floor((x-1)/3)*3+1;
     y2 = floor((y-1)/3)*3+1;
     
     [i_nb1,i_nb2,i_nb3] = find_neighbours(ind(i),ind);
  
     a1 = contain_unique(s{i},s(i_nb1)); 
     a2 = contain_unique(s{i},s(i_nb2));
     a3 = contain_unique(s{i},s(i_nb3));
     if a1>0
      m(x,y) = a1;
     end
     if a2>0
      m(x,y) = a2;
     end
     if a3>0
      m(x,y) = a3;
     end
 end
 
    function [i_nb1,i_nb2,i_nb3] = find_neighbours(ind0,ind)
     
        i_nb1 = [];
        i_nb2 = [];
        i_nb3 = [];
        [x0,y0]=ind2sub([9,9],ind0);
         [sx0,sy0]=square_index(x0,y0);
        for i=1:length(ind)
            if ind(i)==ind0
                continue;
            end
           [x,y]=ind2sub([9,9],ind(i));
           [sx,sy]=square_index(x,y);
           if (sx==sx0 && sy==sy0)
               i_nb1 = [i_nb1,i];
           end      
           if x==x0
               i_nb2 = [i_nb2,i];
           end
           if y==y0
               i_nb3 = [i_nb3,i];
           end
           
        end
     
        
    function a = contain_unique(s,s_nb)
            
        arr = [];
        for j=1:length(s_nb)
            arr = [arr,s_nb{j}];
        end
        arr = unique(arr);
        a = 0;
        for i=1:length(s)
            if ~any(s(i)==arr)
                a = s(i);
                return;
            end
        end
            
 
    function arr=find_possible_values(m,x,y)
        
        a1 = m(x,m(x,:)>0);
        a2 = m(m(:,y)>0,y);
        
        x2=floor((x-1)/3)*3+1;
        y2 = floor((y-1)/3)*3+1;
        
        a3 = m(x2:x2+2,y2:y2+2);
        
        a4 = unique([a1(:);a2(:);a3(:)]);
        
        arr = setdiff(1:9,a4);
        
        
        function [sx,sy] = square_index(x,y)
            
            
            sx=floor((x-1)/3);
            sy = floor((y-1)/3);
            