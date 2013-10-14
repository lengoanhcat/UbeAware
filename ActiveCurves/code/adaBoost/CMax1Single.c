/*perfoms operations in computing S1 maps*/
# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        
# include "math.h"
# include <mmintrin.h> /*header file of SIMD intrisics*/
# include <xmmintrin.h>/*header file of SIMD intrisics*/
 
# define PI 3.1415926
# define ROUND(x) (floor((x)+.5))
# define NEGMAX -1e10
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))

# define USE_SIMD /*enable SIMD */
/* Generating integer vector */
int *int_vector(int n)
{
    int *v; 
    v = (int*) mxCalloc (n, sizeof(int));
    return v; 
}
/* Generating integer matrix */
int **int_matrix(int m, int n)
{
    int **mat; 
    int i; 
    mat = (int**) mxCalloc(m, sizeof(int*)); 
    for (i=0; i<m; i++)
        mat[i] = int_vector(n); 
    return mat; 
}
/* Compute pixel index in the vector that stores image */
__inline int px(int x, int y, int lengthx, int lengthy)  /* the image is lengthx*lengthy */
{            
	return (x + (y-1)*lengthx - 1); 
}
/* variables */
int numOrient, locationShiftLimit, orientShiftLimit; /* key parameters */
int numImage; /* number of images */   
double *allSizex, *allSizey;  /* sizes of images */
float **S1maps, **M1maps;           
int  sx, sy; /* height, width of image */ 
int numShift, **xShift, **yShift, **orientShifted; /* stored shifts in x and y, and the shifted orientations */  
int h; /*half size of filter length*/
/* store all the shifts or perturbations in local maximum pooling */
void StoreShift()  
{
    int dLoc2,dOri2,currOri,currLoc;
    double alpha; 
   int orient,iShift;
    /* store all the possible shifts for all orientations */
    numShift = (locationShiftLimit*2+1)*(orientShiftLimit*2+1); 
    xShift = int_matrix(numOrient, numShift); 
    yShift = int_matrix(numOrient, numShift); 
    orientShifted = int_matrix(numOrient, numShift);
    for (orient=0; orient<numOrient; orient++)        
    {
	    iShift = 0;
	    for(dOri2 = -orientShiftLimit; dOri2<=orientShiftLimit; ++dOri2){
		currOri = orient + dOri2;
		if (currOri<0)
			currOri += numOrient; 
		else if (currOri>=numOrient)
			currOri -= numOrient; 

		alpha = PI*orient/numOrient;
		for(dLoc2 = -locationShiftLimit; dLoc2<=locationShiftLimit; ++dLoc2){
			xShift[orient][iShift] = ROUND(dLoc2*cos(alpha)); 
			yShift[orient][iShift] = ROUND(dLoc2*sin(alpha)); 
			orientShifted[orient][iShift] = currOri; 
			++iShift; 
		}	       
	    }
    
    }

}

/*set response to NEGMAX*/
void resetM1(int iImg)
{
	int  x, y; 
	int iOri,iAng,iLen,here;
	float *pS2,*pM2;
	for (iOri=0; iOri<numOrient; ++iOri)
		for (here =0; here<sx*sy;++here)
				M1maps[iOri*numImage+iImg][here]=NEGMAX;
}


/*local max operation, done in parallel*/
void Cmm2(int iImg)
{
	int  x, y; 
	int iOri,iAng,iLen;
	int iShift; /*index of deformations*/
	int here; /*1 dimension index of image*/
	int dx,dy,orient1,dHere; /*index for deformed entities*/
	float localMax,r;
	float *pS1,*pM1; /*pointer to specific S1 and M1 maps */
    #ifdef USE_SIMD
    int startx,starty,endx,endy;
    __m128 maxDest,maxSource;
    #endif
	for (iOri=0; iOri<numOrient; ++iOri){
		pM1 = M1maps[iOri*numImage+iImg];/*grab the current M1 map*/
		for (iShift = 0; iShift< numShift; ++iShift){
			dx = xShift[iOri][iShift];
			dy = yShift[iOri][iShift];
			orient1 = orientShifted[iOri][iShift]; 	
			pS1 = S1maps[orient1*numImage+iImg]; /*grab current S1 map*/
			dHere = dy*sx+dx;
			#ifndef USE_SIMD /*the two branches are doing the same operation*/
                    for (y = MAX(1+h,-dy+1); y<=MIN(sy-h,sy-dy); ++y){
                        for (x = MAX(1+h,-dx+1); x<=MIN(sx-h,sx-dx); ++x){
                            here = px(x,y,sx,sy);
                            localMax=pM1[here];
                            r = pS1[here+dHere];
                            pM1[here]=MAX(r,localMax);
                        }
                    }
                    
			#else /*this branch tries to compute them in parallel */
                    starty = MAX(1+h,-dy+1); startx = MAX(1+h,-dx+1);
                    endy = MIN(sy-h,sy-dy); endx = MIN(sx-h,sx-dx);
                    for(y = starty; y<=endy; ++y){
                        for (x =startx; x<=endx-3; x = x+4){
                            here = px(x,y,sx,sy);
                            maxDest = _mm_loadu_ps(&(pM1[here]));
                            maxSource=_mm_loadu_ps(&(pS1[here+dHere]));
                            maxDest=_mm_max_ps(maxDest,maxSource);
                            _mm_storeu_ps(&(pM1[here]),maxDest);                           
                        }
                        x = x-4;
                        for(;x<=endx;++x){
                            
                            here = px(x,y,sx,sy);
                            localMax=pM1[here];
                            r = pS1[here+dHere];
                            pM1[here]=MAX(r,localMax);
                        }
                    }
			#endif	

				}	
			}
}



/*grab memory pointers and get parameters*/
void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    int img, orient, c;
    mxArray *f;
    
    c = 0;
    numImage = ROUND(mxGetScalar(prhs[c++]));
    numOrient = ROUND(mxGetScalar(prhs[c++]));
    h = ROUND(mxGetScalar(prhs[c++]));
    
    S1maps = mxCalloc(numImage*numOrient, sizeof(float*));
    for (img=0; img<numImage; img++)
        for (orient=0; orient<numOrient; orient++)
    {
        f = mxGetCell(prhs[c], orient*numImage+img);
        S1maps[orient*numImage+img] = (float*)mxGetPr(f);
        }
    c++;
	
    M1maps = mxCalloc(numImage*numOrient, sizeof(float*));
    for (img=0; img<numImage; img++)
    {
        for (orient=0; orient<numOrient; orient++)
        {
            f = mxGetCell(prhs[c], orient*numImage+img);
            M1maps[orient*numImage+img] = (float*)mxGetPr(f);
        }
    }
    c++;
    locationShiftLimit = ROUND(mxGetScalar(prhs[c++])); /* parameter Lrange*/
    orientShiftLimit = ROUND(mxGetScalar(prhs[c++]));/*parameter Orange*/
    
    
    allSizex = mxGetPr(prhs[c++]);
    allSizey = mxGetPr(prhs[c++]);
    
    StoreShift();
    for (img=0; img<numImage; img++)
    {
        sx = ROUND(allSizex[img]);
        sy = ROUND(allSizey[img]);
        
        resetM1(img);
        Cmm2(img);
    }
}









