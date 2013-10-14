# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        /* the algorithm is connect to matlab */
# include "math.h"
# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))
# define NEGMAX -1e10

/*
* Local Maximization.
* Input: Response map.
* Output: Local-maxed response map.
*/
int *int_vector(int n)
{
    int *v; 
    v = (int*) mxCalloc (n, sizeof(int));
    return v;
}

float *single_vector(int n)
{
    float *v; 
    v = (float*) mxCalloc (n, sizeof(float));    
    return v; 
}
/* Generating double matrix */
float **single_matrix(int m, int n)
{
    float **mat; 
    int i; 
    mat = (float**) mxCalloc(m, sizeof(float*)); 
    for (i=0; i<m; i++)
        mat[i] = single_vector(n); 
    return mat; 
}
/* Free matrix space */
int **int_matrix(int m, int n)
{
    int **mat;
    int i;
    mat = (int**) mxCalloc(m, sizeof(int*));
    for (i=0; i<m; i++)
        mat[i] = int_vector(n);
    return mat;
}

void free_matrix(void **a, int nrow, int ncol)
{
	int i;
	for(i = 0; i < nrow; i++)
		mxFree(a[i]);
	mxFree(a);
}



 /* variables */
 int nimage;                      /* number of images */
 int norient;                      /* number of orientations */
 int nAngle;
 int nLength;
 float **S2map;                /* filtered images */
 int h;                      /* halfsizes of filters */
 int sx, sy;                 /* current size of image */
 int ax, ay;                 /* half size of bounding box of learned model */ 
/* double *Mi, *Mx, *My, *Mm, *Mm1; */ /* storing selected bases */
 int L, ore;                 /* allowed ranges of shifting in location and orientation */
 int sub;                    /* subsampling */
 double *Sx1, *Sy1;           /* sizes of images */
 int **sinsh, **cossh;       /* store the shift to avoid repeated sin and cosin computation */
 double *maxi, *imaxi; 
 int nDims[4];

/* getting the index of matlab image, (x,y) location, (sx, sy) sizes */
__inline int px(int x, int y, int bx, int by)    
{            
	return (x + (y-1)*bx - 1); 
}
/*get the first index of the SUM2maps*/
__inline  int page(int iOrient, int iAngle, int iLength, int iImage){
	return iImage+iLength*nDims[1]+iAngle*nDims[2]+iOrient*nDims[3];
}


 /* store the shift values so that we do not need to repeat the sin and cos computation */ 
void storeshift()
{
    int ind, l;
    float theta; 
    
    sinsh = int_matrix(norient, L+L+1); 
    cossh = int_matrix(norient, L+L+1); 
    for (ind=0; ind<norient; ind++)        
    {
        theta = PI*ind/norient; 
        for (l=-L; l<=L; l++)
         {
            sinsh[ind][l+L] = floor(l*sub*sin(theta)+.5); 
            cossh[ind][l+L] = floor(l*sub*cos(theta)+.5); 
        }   
    }
}

/* for Gabor(x, y, orientation = ind), find local maximum in image i */
float shiftmax(int iOri,int iAng, int iLen,int iImg, int x, int y, int* rm)
{
	float m;
	int x1, y1, l, ml, mx, my, de, d, d1, here, mo, md, ths; 

	/*set to default location and orientation*/
	ml = 0; /*chnage in length*/
        mx = x; my= y;mo = iOri;  de =0;
	m = NEGMAX; 
	for (l=-L; l<=L; l++)   
	{
        x1 = x + cossh[iOri][l+L]; 
        y1 = y + sinsh[iOri][l+L];   /* shifting in normal direction */
        if ((x1>=1)&&(x1<=sx)&&(y1>=1)&&(y1<=sy))
        {
			here = px(x1, y1, sx, sy);
			for (de=-ore; de<=ore; de++)
            {
				d = de+iOri;
				d1 = d;
				if (d<0)
					d1 = d+norient;
				if (d>=norient)
					d1 = d-norient;   /* shifting in orientation */
				ths = page(d1,iAng,iLen,iImg); 
				if (S2map[ths][here]>m)
                {
					m = S2map[ths][here];   /* local maximization */
					ml = l; mx = x1; my = y1; mo = d1; md = de; 
                }
			}
        }
	}
	rm[0] = ml; rm[1] = mx; rm[2] = my; rm[3] = mo; rm[4] = md; 
	return(m); 
}


/* local maximization */
void Cmm(int i)
{
   int  ind, x, y, rm[5], j, o; 
   float mm; 
   int iOri,iAng,iLen;
   iOri = floor(maxi[0]+.5);
   iAng = floor(maxi[1]+.5);
   iLen = floor(maxi[2]+.5);
	
   x = floor(maxi[3]+.5); y = floor(maxi[4]+.5);      
   mm = shiftmax(iOri,iAng,iLen, i,x, y, rm);   /* local maximum of active Gabor */  
   imaxi[i*6+5] = mm; imaxi[i*6] = rm[3]; 
   imaxi[i*6+1] = iAng; imaxi[i*6+2] = iLen; 
   imaxi[i*6+3] = rm[1]; imaxi[i*6+4]=rm[2]; 
} 


/* mex function is used to pass on the pointers and scalars from matlab, 
   so that heavy computation can be done by C, which puts the results into 
   some of the pointers. After that, matlab can then use these results. 
   
   So matlab is very much like a managing platform for organizing the 
   experiments, and mex C is like a work enginee for fast computation. */

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
 int i, ind, j, c; 
 mxArray *f;  
	int iLen,iImg,iPage;
 
 c = 0; /* counter for input pointers and scalars */
 nimage = floor(mxGetScalar(prhs[c++])+.5); /*numImage*/
 norient = floor(mxGetScalar(prhs[c++])+.5);  /* number of orientations */
 nLength = floor(mxGetScalar(prhs[c++])+.5); /* line length*/
 nAngle = floor(mxGetScalar(prhs[c++])+0.5);

	/*parameters for page number*/
	nDims[0]= 1;
	nDims[1]= nDims[0]*nimage;
	nDims[2]= nDims[1]*nLength;
	nDims[3]= nDims[2]*nAngle;
	
	S2map = mxCalloc(nimage*norient*nLength*nAngle, sizeof(float*));   /* M1 maps */
	for (i=0; i<norient; i++){
		for (j=0; j<nAngle; j++){
			for(iLen = 0; iLen<nLength; ++iLen){
				for(iImg = 0 ; iImg<nimage; ++iImg){
					iPage = page(i,j,iLen,iImg);
					f = mxGetCell(prhs[c], iPage); 
					S2map[iPage] = mxGetPr(f);    /* get pointers to filtered images */       
				}
			}
		}
	}
	c++;

 h = floor(mxGetScalar(prhs[c++])+.5);    /* half size of filters */ 
 
 L = floor(mxGetScalar(prhs[c++])+.5);     /* range of shifting along normal direction of Gabor */
 ore = floor(mxGetScalar(prhs[c++])+.5);   /* range of shifting in orientation */
 sub = floor(mxGetScalar(prhs[c++])+.5);   /* sub-sampling of Gabor filters to improve speed */
 
 Sx1 = mxGetPr(prhs[c++]);
 Sy1 = mxGetPr(prhs[c++]);
 maxi = mxGetPr(prhs[c++]);
 imaxi = mxGetPr(prhs[c++]);
 
 storeshift();
 for (i=0; i<nimage; i++)
 {
	sx = floor(Sx1[i]+.5); 
	sy = floor(Sy1[i]+.5);
	  Cmm(i);      /* parallel local maximum pooling */
 }
}

