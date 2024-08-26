void crop2d(float *d, int l1, int l2, float *d2, int p1, int p2,int l1c,int l2c)
//p1, p2, p3 is the starting position of the cropped region

if ((l1c+p1)>l1 || (l2c+p2)>l2)
    fprintf("(l1c+p1)>l1 || (l2c+p2)>l2 \n");
 
for (int i=0;i<l1c;i++)
    for (int j=0;j<l2c;j++)
        {
            ind1=s2i2(i+p1,j+p2,l1);
            ind2=s2i2(i,j,l1c);
            d2(ind2)=d(ind1);   
        }