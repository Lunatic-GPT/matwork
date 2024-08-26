void matMul(float *A, int rowsA, int colsA, float *B, int colsB,float *AB)

for (int i=0;i<rowsA;i++)
    for (int j=0;j<colsB;j++)
    {
        int iAB = s2i2(i,j,roswA);
        AB[iAB]=0;
        for (int k=0;k<colsA;k++)
        {
           int iA = s2i2(i,k,rowsA);
           int iB=s2i2(k,j,colsA);
           
           AB[iAB]+=A[iA]*B[iB];
         
            
        }
            
            
    }       