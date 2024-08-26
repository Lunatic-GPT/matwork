void crop(float *d, int l1, int l2,int l3, float *d2, int p1, int p2,int p3,int l1c,int l2c,int l3c)
//p1, p2, p3 is the starting position of the cropped region

if ((l1c+p1)>l1 || (l2c+p2)>l2 || (l3c+p3)>l3)
    fprintf("(l1c+p1)>l1 || (l2c+p2)>l2 || (l3c+p3)>l3 \n");
 
for (int i=0;i<l1c;i++)
    for (int j=0;j<l2c;j++)
        for (int k=0;k<l3c;k++)
        {
            ind1=s2i3(i+p1,j+p2,k+p3,l1,l2);
            ind2=s2i3(i,j,k,l1c,l2c);
            d2(ind2)=d(ind1);   
        }