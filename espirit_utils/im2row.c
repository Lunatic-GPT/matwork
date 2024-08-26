void im2row(float *im, int sx,int sy,int nc, int* winSize, float *res)
//res should have size ((sx-winSize[0]+1)*(sy-winSize[1]+1))*(winSize[0]*winSize[1])*nc;


int count=-1;
for (int y=0;y<winSize[1];y++)
{
    for (int x=0;x<winSize[0];x++)
    {
        count++;
        for (int i=0;i<sx-winSize[0]+1;i++)
            for (int j=0;j<sy-winSize[1]+1;j++)
                for (int k=0;k<nc;k++)
                {
                                        
                    int ind_i=s2i3(x+i,y+j,k,sx,sy);

                    int ind_tmp = s2i2(i,j,sx-winSize[0]+1);
                    int ind_o=s2i3(ind_tmp,count,k,(sx-winSize[0]+1)*(sy-winSize[1]+1),winSize[0]*winSize[1],nc);
                    
                    
                    res[ind_o]=im[ind_i];
                    
                /*    res(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:),...
                            (sx-winSize(1)+1)*(sy-winSize(2)+1),1,nc);*/
                    
                }
    }
    
}