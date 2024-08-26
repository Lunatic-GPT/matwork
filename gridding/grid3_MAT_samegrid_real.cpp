
#include "mex.h"

#include <math.h>
#include <stdlib.h>

#include "grid_utils_xp.c" 


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    unsigned long i;

      const char help_txt[] = "USAGE: (see grid3_MAT_2grid_real.c) \n\
 * grid_volume = grid3_MAT_2grid_real(0:data, 1:crds, 2:newcrds,effMtx) \n\
 * \n\
 * REQUIRED: \n\
 *  data: N-D double-precision matrix >=2D; Not complex\n\
 * \n\
 *  coords: coordinates for input data, fastest varying dimension is length 3 \n\
 *          trajectory coordinate points are typically scaled between -0.5 to 0.5 \n\
 * \n\
 *  effMtx: (integer) scalling factor for kernalRadius, kernelRadius = 4/effMtx.\n\
         ";

       

    /* Check for proper number of arguments */
    /* input: 
    *    REQUIRED: data, coords, weight, effMtx, numThreads
    */
    if (nrhs != 3) 
    {
        printf("%s",help_txt);
        mexErrMsgTxt("3 required inputs are: data, coords, effMtx.");
    }

    /* ouput: weights_out */
    if (nlhs != 1) 
    {
        printf("%s",help_txt);
        mexErrMsgTxt("Only 1 output arg is returned.");
    }

    /** check and retrieve input */
    /* PARAMS */
    int width       = (int) *mxGetPr(prhs[2]); /* grid size */

    /* DATA */
    printf("copy data\n");
        assert(prhs[0] != NULL);     /* check existence */
        assert(mxIsDouble(prhs[0])); /* check for type double */
    int nd = mxGetNumberOfDimensions(prhs[0]); /* get coordinate dimensions */
        assert( nd > 0 ); /* check for valid array size */
    const int *dims = mxGetDimensions(prhs[0]);
    unsigned long *dims_l = (unsigned long *) malloc(sizeof(unsigned long) * nd);
    for (i=0;i<nd;i++) dims_l[i] = (unsigned long) dims[i]; 
    dataArray_double *data = new_dataArray_double(nd,dims_l); /* alloc new dataArray_double */
        assert(data != NULL); /* make sure new mem is allocated */
    memcpy( data->data, mxGetPr(prhs[0]), sizeof(double)*(data->num_elem) ); /* copy data */
    free(dims_l);

 
    /* COORDS */
    printf("copy coords\n");
        assert(prhs[1] != NULL);     /* check existence */
        assert(mxIsDouble(prhs[1])); /* check for type double */
    nd = mxGetNumberOfDimensions(prhs[1]); /* get coordinate dimensions */
        assert( nd > 0 ); /* check for valid array size */
    dims = mxGetDimensions(prhs[1]);
        assert( dims[0] == 3 ); /* make sure the fastest varying dim holds an [x,y,z] triplet */
    dims_l = (unsigned long *) malloc(sizeof(unsigned long) * nd);
    for (i=0;i<nd;i++) dims_l[i] = (unsigned long) dims[i]; 
    dataArray_double *coords = new_dataArray_double(nd,dims_l); /* alloc new dataArray_double */
        assert(coords != NULL); /* make sure new mem is allocated */
    memcpy( coords->data, mxGetPr(prhs[1]), sizeof(double)*(coords->num_elem) ); /* copy coords */
    free(dims_l);
 

 
    /* check input data sizes */
    assert( coords->num_elem/3 == data->num_elem   );
  

    /* allocate output array */
    nd = mxGetNumberOfDimensions(prhs[0]); /* get coordinate dimensions */
    dims = mxGetDimensions(prhs[0]);
    dims_l = (unsigned long *) malloc(sizeof(unsigned long) * nd);
    for (i=0;i<nd;i++) dims_l[i] = (unsigned long) dims[i]; 
    
    dataArray_double *out = new_dataArray_double(nd,dims_l);
    free(dims_l);
    /* allocate kernel table */
    printf("allocate kernel\n");
    unsigned long dim_k[1];
    dim_k[0] = DEFAULT_KERNEL_TABLE_SIZE;
    dataArray_double *kern = (dataArray_double*) new_dataArray_double(1,dim_k);
    for(long i=0;i<kern->num_elem;i++) kern->data[i] = 0.;
    
    printf("load kernel\n");
    loadGrid3Kernel(kern);

    printf("End load kernel\n");
    /* radius FOV product */
    double rfp_d = DEFAULT_RADIUS_FOV_PRODUCT;
    double win_d = DEFAULT_WINDOW_LENGTH;
    double radiusFOVproduct=DEFAULT_RADIUS_FOV_PRODUCT;
    double windowLength = DEFAULT_WINDOW_LENGTH;
    for (long i=0;i<coords->num_elem;i++)
    coords->data[i]=coords->data[i]/windowLength;

   double kernelRadius           = radiusFOVproduct /width;  //in units of 1/width
    
    double   dist_multiplier = (kern->num_elem - 1)/kernelRadius/kernelRadius;
       
for (long i=0;i<data->num_elem;i++)
{    
       
    	double x = coords->data[i*3]; 
		double y = coords->data[i*3+1];
		double z = coords->data[i*3+2];
        double xb[2]={x-kernelRadius,x+kernelRadius};
        
        double yb[2]={y-kernelRadius,y+kernelRadius};
        
        double zb[2]={z-kernelRadius,z+kernelRadius};
        int kern_ind;
     if (i%2000==0) 
     {
        printf("%d/%d\n",i,data->num_elem);
        mexEvalString("drawnow;");
     }
	for (long j=i+1;j<data->num_elem;j++)
    { 
    	double x2 = coords->data[j*3]; 
		double y2 = coords->data[j*3+1];
		double z2 = coords->data[j*3+2];
        
    
        if(x2>=xb[1] || x2<=xb[0] || y2>=yb[1] || y2<=yb[0] || z2>=zb[1] || z2<=zb[0])
        continue;
        
    
       double dist_sqr=(x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2);
        
        kern_ind=int(dist_sqr*dist_multiplier+0.5);
        if(kern_ind>kern->num_elem-1)
            kern_ind=kern->num_elem-1;
        
     
        out->data[i]+=kern->data[kern_ind]*data->data[j];
        
        out->data[j]+=kern->data[kern_ind]*data->data[i];
        
    }	
}

       for (long i=0;i<data->num_elem;i++)
          out->data[i]+=data->data[i]*kern->data[0];

 
    /* Create an mxArray for the return argument */ 
    printf("copy output grid\n");
    
    nd = mxGetNumberOfDimensions(prhs[0]); /* get coordinate dimensions */
    int *odims = (int*) malloc( sizeof(int)*(nd) );

    dims = mxGetDimensions(prhs[0]);
    
    for(i=0;i<nd;i++)
        odims[i] = dims[i];
    /* matlab doesn't like odims to be static? */
    plhs[0] = mxCreateNumericArray(nd, odims, mxDOUBLE_CLASS, mxREAL);  
        assert(plhs[0] != NULL); /* check that mem was allocated */
    memcpy( mxGetPr(plhs[0]), out->data, (out->num_elem) * sizeof(double));
    free(odims);

    printf("free local memory\n");
    /* free temp data */
    free_dataArray_double(coords);
    free_dataArray_double(data);
    free_dataArray_double(out);

    /* free kernel table */
    free_dataArray_double(kern);
}


