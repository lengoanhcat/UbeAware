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
int numOrient;                      /* number of orientations */
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
int numShift;
int **xShift, **yShift, **orientShifted; /* stored shifts in x and y, and the shifted orientations */
int  **curvatureShifted;
int numOrient;
int locationShiftLimit, orientShiftLimit; /* key activeness parameters */
int curvatureShiftLimit,lengthShiftLimit;
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
void StoreShift2()
{
	int dLoc2,dOri2,dCurve2,dLen2;
	int currOri,currLoc,currCurvature,currLen;
	double alpha;
	int orient,iShift,iCurve,iPage,iAngle;
	int numCurve;
	/* store all the possible shifts for all arcs */
	numShift = (locationShiftLimit*2+1)*(orientShiftLimit*2+1);
	numShift *= (curvatureShiftLimit*2+1);
	xShift = int_matrix(nAngle*numOrient, numShift);
	yShift = int_matrix(nAngle*numOrient, numShift);
	orientShifted = int_matrix(nAngle*numOrient, numShift);
	curvatureShifted = int_matrix(nAngle*numOrient,numShift);
	for (orient=0; orient<numOrient; orient++)
	{
		for(iAngle = 0; iAngle<nAngle; ++iAngle){
			iPage = orient*nAngle+iAngle;
			iShift = 0;
			for(dOri2 = -orientShiftLimit; dOri2<=orientShiftLimit; ++dOri2){
				currOri = orient + dOri2;
				if (currOri<0)
					currOri += numOrient;
				else if (currOri>=numOrient)
					currOri -= numOrient;
				for(dCurve2 = -curvatureShiftLimit; dCurve2<=curvatureShiftLimit; ++dCurve2){
					currCurvature = iAngle + dCurve2;
					currCurvature = MIN(currCurvature,nAngle-1);
					currCurvature = MAX(currCurvature,0);		
					alpha = PI*orient/numOrient;
					for(dLoc2 = locationShiftLimit; dLoc2>=-locationShiftLimit;--dLoc2){
						xShift[iPage][iShift] = ROUND(dLoc2*cos(alpha));
						yShift[iPage][iShift] = ROUND(dLoc2*sin(alpha));
						orientShifted[iPage][iShift] = currOri;
						curvatureShifted[iPage][iShift] = currCurvature;
						++iShift;
					}
				}
			}
		}
	}

}


void resetM2(int iImg)
{
	int  x, y;
	int iOri,iAng,iLen,here;
	float *pS2,*pM2;
	for (iOri=0; iOri<numOrient; ++iOri){
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
	int iShift,iPage;
	int dx,dy,dHere,here;
	int orient1,angle1;
	float localMax,r;
	float *pS2,*pM2;
	int dLen2,length1;
#ifdef USE_SIMD
	int startx,starty,endx,endy;
	__m128 maxDest,maxSource;
#endif
	h = 0;
	for (iOri=0; iOri<numOrient; ++iOri){
		for (iAng = 0; iAng<nAngle;++iAng){

			for(iLen = 0; iLen<nLength; ++iLen){
				pM2=M2map[page(iOri,iAng,iLen,iImg)];
				iPage = iOri*nAngle+iAng;
				for(dLen2= -lengthShiftLimit; dLen2<=lengthShiftLimit; ++dLen2){
					length1 = iLen+dLen2;
					if(length1<0) length1=0;
					if(length1>=nLength) length1 = nLength-1;

					for (iShift = 0; iShift< numShift; ++iShift){
						dx = xShift[iPage][iShift];
						dy = yShift[iPage][iShift];

						orient1 = orientShifted[iPage][iShift];
						angle1 = curvatureShifted[iPage][iShift];

						pS2 = S2map[page(orient1,angle1,length1,iImg)];
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
	numOrient = floor(mxGetScalar(prhs[c++])+.5);  /* number of orientations */
	h = ROUND(mxGetScalar(prhs[c++]));

	nLength = floor(mxGetScalar(prhs[c++])+.5); /* line length*/
	nAngle = floor(mxGetScalar(prhs[c++])+0.5);

	/*parameters for page number*/
	nDims[0]= 1;
	nDims[1]= nDims[0]*nimage;
	nDims[2]= nDims[1]*nLength;
	nDims[3]= nDims[2]*nAngle;

	S2map = mxCalloc(nimage*numOrient*nLength*nAngle, sizeof(float*));   /* M1 maps */
	for (i=0; i<numOrient; i++){
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


	M2map = mxCalloc(nimage*numOrient*nLength*nAngle, sizeof(float*));   /* M1 maps */
	for (i=0; i<numOrient; i++){
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
	curvatureShiftLimit = floor(mxGetScalar(prhs[c++])+.5); /*range of curvature shift*/
	lengthShiftLimit = floor(mxGetScalar(prhs[c++])+.5); /*ramge of length shift*/

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

