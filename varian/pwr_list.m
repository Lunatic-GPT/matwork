function  pwr_list()

pwr = 18:2:40;

pwr = repmat(pwr,[10,1]);

pwr=pwr(:);

pwr = [pwr,ones(size(pwr))];

save_mat_int(pwr,'pwr_list3.txt');
