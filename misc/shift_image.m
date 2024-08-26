function shift_image(name_old,name_new,shift)
%shift_image(name_old,name_new,shift)
% shift: 1*3 vector.  positive to shift to larger indices and negative to
% smaller indices.  Same convention as in circshift.
    
 [a,info]=BrikLoad(name_old);
 if size(a,4)>1
     shift=[shift,0];
 end
 a=circshift(a,shift);
 Opt.OverWrite='y';
 Opt.prefix=name_new;
 history=sprintf('shift_image(%s,%s,%s)',name_old,name_new,num2str(shift));
 info.HISTORY_NOTE=[info.HISTORY_NOTE,history];
 WriteBrik(a,info,Opt);
