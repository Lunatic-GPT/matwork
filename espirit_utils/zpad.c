void zpad(float *d, int l1, int l2, int l3, int p1,int p2,int p3, float *d2)
// zpad data by p1, p2, p3; if p? is odd, pad £¨p1-1£©/2 on the left and (p1+1)/2 on the right
for (int i=0;i<l1+p1;i++)
    for (int j=0;j<l2+p2;j++)
        for (int k=0;k<l3+p3;k++)
    {
            if ( i<p1/2 || i>=l1+p1/2 || j<p2/2 || j>=l2+p2/2 || k<p3/2 || k>=l3+p3/2 
              d2(ind)=0;
            else
            {
              ind1=s2i3(i-p1/2,j-p2/2,k-p3/2,l1,l2);
              ind2=s2i3(i,j,k,l1+p1,l2+p2);
              d2(ind2)=d(ind1);
            }   
            
    }       