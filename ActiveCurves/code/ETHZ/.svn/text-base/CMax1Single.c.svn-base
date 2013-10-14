# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        
# include "math.h"

# include <mmintrin.h>
# include <xmmintrin.h>
 
# define PI 3.1415926
# define ROUND(x) (floor((x)+.5))
# define NEGMAX -1e10
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))

# define USE_SIMD
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
int numImage; /* number of resolutions */   
double *allSizex, *allSizey;  /* sizes of images at multiple resolutions */
float **S1map, **M1map;           
int sizex, sizey, sx, sy; /* MAX1 maps are smaller than SUM1 maps by subsample */ 
int numShift, **xShift, **yShift, **orientShifted; /* stored shifts in x and y, and the shifted orientations */      
/* store all the shifts or perturbations in local maximum pooling */
void StoreShift2()  
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
				M1map[iOri*numImage+iImg][here]=NEGMAX;
}


/*version 2 of local max*/
void Cmm2(int iImg)
{
	int  x, y; 
	int iOri,iAng,iLen;
	int iShift;
	int dx,dy,orient1,dHere,here;
	float localMax,r;
	float *pS1,*pM1;
    #ifdef USE_SIMD
    int startx,starty,endx,endy;
    __m128 maxDest,maxSource;
    #endif
	for (iOri=0; iOri<numOrient; ++iOri){
		pM1 = M1map[iOri*numImage+iImg];
		for (iShift = 0; iShift< numShift; ++iShift){
			dx = xShift[iOri][iShift];
			dy = yShift[iOri][iShift];
			orient1 = orientShifted[iOri][iShift]; 	
			pS1 = S1map[orient1*numImage+iImg];
			dHere = dy*sx+dx;
 #ifndef USE_SIMD
                    
                    for (y = MAX(1,-dy+1); y<=MIN(sy,sy-dy); ++y){
                        for (x = MAX(1,-dx+1); x<=MIN(sx,sx-dx); ++x){
                            here = px(x,y,sx,sy);
                            localMax=pM1[here];
                            r = pS1[here+dHere];
                            pM1[here]=MAX(r,localMax);
                            
                        }
                    }
                    
                    #else
                    starty = MAX(1,-dy+1); startx = MAX(1,-dx+1);
                    endy = MIN(sy,sy-dy); endx = MIN(sx,sx-dx);
                    for(y = starty; y<=endy; ++y){
                        for (x =startx; x<=endx-4; x = x+4){
                            here = px(x,y,sx,sy);
                            maxDest = _mm_loadu_ps(&(pM1[here]));
                            maxSource=_mm_loadu_ps(&(pS1[here+dHere]));
                            maxDest=_mm_max_ps(maxDest,maxSource);
                            _mm_storeu_ps(&(pM1[here]),maxDest);                           
                        }
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


void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    int img, orient, c;
    mxArray *f;
    
    c = 0;
    numImage = ROUND(mxGetScalar(prhs[c++]));
    numOrient = ROUND(mxGetScalar(prhs[c++]));
    
    S1map = mxCalloc(numImage*numOrient, sizeof(float*));
    for (img=0; img<numImage; img++)
        for (orient=0; orient<numOrient; orient++)
    {
        f = mxGetCell(prhs[c], orient*numImage+img);
        S1map[orient*numImage+img] = (float*)mxGetPr(f);
        }
    c++;
    M1map = mxCalloc(numImage*numOrient, sizeof(float*));
    for (img=0; img<numImage; img++)
    {
        for (orient=0; orient<numOrient; orient++)
        {
            f = mxGetCell(prhs[c], orient*numImage+img);
            M1map[orient*numImage+img] = (float*)mxGetPr(f);
        }
    }
    c++;
    locationShiftLimit = ROUND(mxGetScalar(prhs[c++]));
    orientShiftLimit = ROUND(mxGetScalar(prhs[c++]));
    
    
    allSizex = mxGetPr(prhs[c++]);
    allSizey = mxGetPr(prhs[c++]);
    
    
    
    StoreShift2();
    for (img=0; img<numImage; img++)
    {
        sizex = ROUND(allSizex[img]);
        sizey = ROUND(allSizey[img]);
        
        sx = floor((double)sizex);
        sy = floor((double)sizey);
        
        resetM1(img);
        
        
        Cmm2(img);
    }
}









