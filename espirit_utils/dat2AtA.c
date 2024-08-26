void dat2AtA(float *data, int sx,int sy, int nc,int *kSize,float *AtA)
// data: sx*sy*nc array
// AtA: (kSize[0]*kSize[1]*nc)*(kSize[0]*kSize[1]*nc);

//im2row(float *im, int sx,int sy,int nc, int* winSize, float *res)

int tsx=(sx-kSize[0]+1)*(sy-kSize[1]+1);
int tsy=kSize[0]*kSize[1];
int total =tsx*tsy*nc;

float *A=(float *)malloc(total*sizeof(float));
float *At=(float *)malloc(total*sizeof(float));

im2row(data,sx,sy,nc,kSize,A);

transpose(A,tsx,tsy*nc,At);
matMul(At,tsy*nc,tsx,At,tsy*nc,AtA);

free(A);
free(At);
        
        
//kernel = AtA;
//kernel = reshape(kernel,kSize(1),kSize(2),nc,size(kernel,2));
