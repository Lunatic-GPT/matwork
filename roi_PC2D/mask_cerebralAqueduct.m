function mask_cerebralAqueduct( m_large, dcm,prefix)
%  mask_cerebralAqueduct( m_large, dcm,prefix);
%m_large: a mask surround the aqueduct
% dcm: dcm dir name;
% prefix: output file name for aqueduct mask

mag=ri(dcm);

mag=mean(mag,4);
m_large=ri_d1(m_large);
y=mag(m_large>0);
ys=sort(y);
ys=ys(1:end-10);

roi2=mag>(mean(ys)+std(ys)*2);
roi=clusterize2(m_large&roi2);
roi=roi==1;

bound=bwmorph(roi,'remove');
[img,cm]=combine_over_under(mag,bound,[],[1,0,0],bound);
figure(10);
imshow(img,cm);
save(prefix, 'roi');