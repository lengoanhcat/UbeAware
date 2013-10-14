# include <stdio.h>
# include "mex.h"        /* the algorithm is connect to matlab */
# include "math.h"
# include "assert.h"
# include <mmintrin.h>
# include <xmmintrin.h>
# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))
# define NEGMAX -1e10
# define ROUND(x) (floor((x)+.5))


#define USE_SIMD

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


/* variables */
int nimage;                      /* number of images */
int norient;                      /* number of orientations */
int nLength;
int nAngle;
float **S2map;                /* filtered images */
float **M2map;                /* local maximum pooling for filtered images */
int h;                      /* halfsizes of filters */
int sx, sy;                 /* current size of image */
int ax, ay;                 /* half size of bounding box of learned model */
double *Mi, *Mx, *My, *Mm, *Mm1;  /* storing selected bases */
int L, ore;                 /* allowed ranges of shifting in location and orientation */

double *Sx1, *Sy1;           /* sizes of images */
int nDims[4];
float lowerBound;
int numShift, **xShift, **yShift, **orientShifted; /* stored shifts in x and y, and the shifted orientations */
int numOrient,locationShiftLimit, orientShiftLimit; /* key parameters */
/* getting the index of matlab image, (x,y) location, (sx, sy) sizes */
__inline int px(int x, int y, int bx, int by)
{
    return (x + (y-1)*bx - 1);
}

/*get the first index of the SUM2maps*/
__inline int page(int iOrient, int iAngle, int iLength, int iImage){
    return iImage+iLength*nDims[1]+iAngle*nDims[2]+iOrient*nDims[3];
}


/* store the shift values so that we do not need to repeat the sin and cos computation */
void StoreShift()
{
    int orient, i, j, shift, ci, cj, si, sj, os;
    double alpha;
    /* store all the possible shifts for all orientations */
    numShift = (locationShiftLimit*2+1)*(orientShiftLimit*2+1);
    xShift = int_matrix(numOrient, numShift);
    yShift = int_matrix(numOrient, numShift);
    orientShifted = int_matrix(numOrient, numShift);
    /* order the shifts from small to large */
    for (orient=0; orient<numOrient; orient++)
    {
        alpha = PI*orient/numOrient;
        ci = 0;
        
        for (i=0; i<=locationShiftLimit; i++)
            for (si=0; si<=1; si++)
        {
            cj = 0;
            for (j = 0; j<=orientShiftLimit; j++)
                for (sj=0; sj<=1; sj++)
            {
                shift = ci*(2*orientShiftLimit+1) + cj;
                xShift[orient][shift] = ROUND(i*(si*2-1)*cos(alpha));
                yShift[orient][shift] = ROUND(i*(si*2-1)*sin(alpha));
                os = orient + j*(sj*2-1);
                if (os<0)
                    os += numOrient;
                else if (os>=numOrient)
                    os -= numOrient;
                orientShifted[orient][shift] = os;
                
                if (j>0)
                    cj ++;
                else if (sj==1)
                    cj++;
                }
            if (i>0)
                ci ++;
            else if (si==1)
                ci++;
            }
    }
}


/* store the shift values so that we do not need to repeat the sin and cos computation */
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
            for(dLoc2 = locationShiftLimit; dLoc2>=-locationShiftLimit;--dLoc2){
                xShift[orient][iShift] = ROUND(dLoc2*cos(alpha));
                yShift[orient][iShift] = ROUND(dLoc2*sin(alpha));
                orientShifted[orient][iShift] = currOri;
                ++iShift;
            }
        }
        
    }
    
}





void resetM2(int iImg)
{
    int  x, y;
    int iOri,iAng,iLen,here;
    float *pS2,*pM2;
    for (iOri=0; iOri<norient; ++iOri){
        for (iAng = 0; iAng<nAngle;++iAng){
            for(iLen = 0; iLen<nLength; ++iLen){
                pM2=M2map[page(iOri,iAng,iLen,iImg)];
		pS2 = S2map[page(iOri,iAng,iLen,iImg)];	
                for (here =0; here<sx*sy;++here){
                    pM2[here]=NEGMAX;
		}
	
            }
        }
    }
    
}


void Cmm2(int iImg)
{
    int  x, y;
    int iOri,iAng,iLen;
    int iShift;
    int dx,dy,orient1,dHere,here;
    float localMax,r;
    float *pS2,*pM2;
    #ifdef USE_SIMD
    int startx,starty,endx,endy;
    __m128 maxDest,maxSource;
    #endif
    for (iOri=0; iOri<norient; ++iOri){
        for (iAng = 0; iAng<nAngle;++iAng){
            for(iLen = 0; iLen<nLength; ++iLen){
                pM2=M2map[page(iOri,iAng,iLen,iImg)];
                for (iShift = 0; iShift< numShift; ++iShift){
                    dx = xShift[iOri][iShift];
                    dy = yShift[iOri][iShift];
                    orient1 = orientShifted[iOri][iShift];
                    pS2 = S2map[page(orient1,iAng,iLen,iImg)];
                    dHere = dy*sx+dx;
                    
                    #ifndef USE_SIMD
                    
                    for (y = MAX(1+h,-dy+1); y<=MIN(sy-h,sy-dy); ++y){
                        for (x = MAX(1+h,-dx+1); x<=MIN(sx-h,sx-dx); ++x){
                            here = px(x,y,sx,sy);
                            localMax=pM2[here];
                            r = pS2[here+dHere];
                            pM2[here]=MAX(r,localMax);
                            
                        }
                    }
                    
                    #else
                    starty = MAX(1+h,-dy+1); startx = MAX(1+h,-dx+1);
                    endy = MIN(sy-h,sy-dy); endx = MIN(sx-h,sx-dx);
                    for(y = starty; y<=endy; ++y){
                        for (x =startx; x<=endx-3; x = x+4){
                            here = px(x,y,sx,sy);
                            maxDest = _mm_loadu_ps(&(pM2[here]));
                            maxSource=_mm_loadu_ps(&(pS2[here+dHere]));
                            maxDest=_mm_max_ps(maxDest,maxSource);
                            _mm_storeu_ps(&(pM2[here]),maxDest);                           
                        }
			x=x-4;
                        for(;x<=endx;++x){
                            here = px(x,y,sx,sy);
                            localMax=pM2[here];
                            r = pS2[here+dHere];
                            pM2[here]=MAX(r,localMax);
                        }
                    }
                    #endif
                }
            }
        }
    }
    
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
    int iLen,iImg,iPage;
    mxArray *f;
    
    c = 0; /* counter for input pointers and scalars */
    nimage = floor(mxGetScalar(prhs[c++])+.5); /*numImage*/
    norient = floor(mxGetScalar(prhs[c++])+.5);  /* number of orientations */
    numOrient=norient;
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
                    S2map[iPage] = (float*)mxGetPr(f);    /* get pointers to filtered images */
                }
            }
        }
    }
    c++;
    
    
    M2map = mxCalloc(nimage*norient*nLength*nAngle, sizeof(float*));   /* M1 maps */
    for (i=0; i<norient; i++){
        for (j=0; j<nAngle; j++){
            for(iLen = 0; iLen<nLength; ++iLen){
                for(iImg = 0 ; iImg<nimage; ++iImg){
                    iPage = page(i,j,iLen,iImg);
                    f = mxGetCell(prhs[c], iPage);
                    M2map[iPage] = (float*)mxGetPr(f);    /* get pointers to filtered images */
                }
            }
        }
    }
    c++;
    
    locationShiftLimit = floor(mxGetScalar(prhs[c++])+.5);     /* range of shifting along normal direction of Gabor */
    orientShiftLimit = floor(mxGetScalar(prhs[c++])+.5);   /* range of shifting in orientation */
    
    
    Sx1 = mxGetPr(prhs[c++]);
    Sy1 = mxGetPr(prhs[c++]);
    StoreShift2();
    /*mexPrintf("numShift: %d\n",numShift);*/
    for (i=0; i<nimage; i++)
    {
        sx = floor(Sx1[i]+.5);
        sy = floor(Sy1[i]+.5);
        /* parallel local maximum pooling */
        resetM2(i);
        Cmm2(i);
        
        /*  local maximum pooling */
        /*Cmm(i);     */
    }
}

