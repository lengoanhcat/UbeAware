
# include <mmintrin.h>
# include <xmmintrin.h>
# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        /* the algorithm is connect to matlab */
# include "math.h"
# include <assert.h>
# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))
# define NEGMAX -1e10
#define USE_SIMD

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
int nimage, norient, nbasis, sx, sy, h;                      
float **M1map, **S2map,**M2map;                              
int **lsinsh, **lcossh;
int locationShiftLimit,orientShiftLimit;
int numShift, **xShift, **yShift, **orientShifted; /* stored shifts in x and y, and the shifted orientations */
int cLengthLimit,cAngleLimit;
int **startDx,**startDy,**endDx,**endDy;
int **startO,**endO;
int nLength,nAngle;
float *pLambda;
int nCurves;
int nDeform;
double *selCurves;
double *maxLength;

/* getting the index of matlab image, (x,y) location, (sx, sy) sizes */
__inline int px(int x, int y, int bx, int by)    
{            
	return (x + (y-1)*bx - 1); 
}



/* store the shift values so that we do not need to repeat the sin and cos computation */
void StoreShift2()
{
    int dLoc2,dOri2,currOri,currLoc;
    double alpha;
    int orient,iShift;
    int numOrient;

    numOrient=norient;
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
                xShift[orient][iShift] = floor(dLoc2*cos(alpha)+0.5f);
                yShift[orient][iShift] = floor(dLoc2*sin(alpha)+0.5f);
                orientShifted[orient][iShift] = currOri;
                ++iShift;
            }
        }
        
    }
    
}



/*store the shift valuds for curve*/
void curve_shift()
{
    int iOri,iAng,iPage;
	int dAngle,l;
	float theta,dTheta;
	int nAngle = 2*cAngleLimit+1;
	startDx = int_matrix(norient*nAngle, nLength);
	startDy = int_matrix(norient*nAngle, nLength);
	
	endDx = int_matrix(norient*nAngle, nLength);
	endDy = int_matrix(norient*nAngle, nLength);
	
	startO = int_matrix(norient*nAngle, nLength);
	endO = int_matrix(norient*nAngle, nLength);
	
	for (iOri=0; iOri<norient; ++iOri){
		for (iAng = 0; iAng<nAngle;++iAng){
			dAngle = -cAngleLimit+iAng;
			iPage = iOri*nAngle+iAng;
			startDx[iPage][0] = 0;
			startDy[iPage][0] = 0;
			endDx[iPage][0]=0;
			endDy[iPage][0]=0;
			startO[iPage][0] = iOri;
			endO[iPage][0] = iOri;
			for (l=1; l<nLength; ++l){
				/*x-axis: from top-left to downside*/
				/*y-axis: from top-left to rightside*/
				startO[iPage][l] = (iOri + l*dAngle)%norient;
				if (startO[iPage][l]<0) startO[iPage][l]+=norient;
				theta = PI/norient*(iOri + l*dAngle);
				dTheta = dAngle*PI/norient;
				startDx[iPage][l] = startDx[iPage][l-1] - floor((h*1.8)*cos(dTheta/2.0f)*sin(theta-dTheta/2.0f)+.5);
				startDy[iPage][l] = startDy[iPage][l-1] + floor((h*1.8)*cos(dTheta/2.0f)*cos(theta-dTheta/2.0f)+.5);
				
				endO[iPage][l] = (iOri - l*dAngle)%norient;
				if (endO[iPage][l]<0)  endO[iPage][l]+= norient;
				theta = PI/norient*(iOri + norient- l*dAngle);
				endDx[iPage][l] = endDx[iPage][l-1] - floor((h*1.8)*cos(dTheta/2.0f)*sin(theta+dTheta/2.0f)+.5);
				endDy[iPage][l] = endDy[iPage][l-1] + floor((h*1.8)*cos(dTheta/2.0f)*cos(theta+dTheta/2.0f)+.5);
			}
		}
	}
	
	
	
}

int shiftRange(int x1,int y1, int x2,int y2,int *top_x, int *bot_x, int *left_y, int *right_y){
	int size_x,size_y;
	*top_x = MIN(x1,x2);
	*bot_x = MAX(x1,x2);
	*left_y= MIN(y1,y2);
	*right_y=MAX(y1,y2); 
	
	
	size_x = *bot_x-*top_x;
	size_y = *right_y-*left_y;
	if (size_x> sx ||size_y>sy) return -1;
	
	*top_x =MAX(1+h, 1 + h - *top_x);
	*bot_x =MIN(sx-h, sx - h - *bot_x);
	*left_y =MAX(1+h, 1 + h - *left_y);
	*right_y =MIN(sy-h,sy - h - *right_y);
	
	
	
	return 1;
}

/* the shared sketch algorithm */
void Csum2LC()
{
	int x, y, here;
	int dHere_s,dHere_e;
	int xStart,yStart,xEnd,yEnd;
        int iOri,iLen,numOrient;
	int iCache;
	float *pPage,*pPrevPage,*pStartM1,*pEndM1;
	int imgSize = sx*sy;
	int len0,ori0,ang0,ori1;
	int iCurve;
	int totalSketch = 0;
	float lambda_s,lambda_o;
	
        numOrient = norient;
	for (iCurve =0; iCurve<nCurves; ++iCurve) {
		len0=floor(selCurves[iCurve+nCurves*2]+0.5);
		ang0=floor(selCurves[iCurve+nCurves]+0.5);
		ori0=floor(selCurves[iCurve]+0.5);
		for(iOri=-orientShiftLimit; iOri<=orientShiftLimit;++iOri){
			/*adjust orientation*/
			ori1=ori0+iOri;
			if(ori1<0)
				ori1+=numOrient;
			else if(ori1>=numOrient)
				ori1-=numOrient;
			/*compute Sum1 for this orientation;*/
			iCache = ori1*nAngle+ang0;
			
			pPage = S2map[nCurves*(iOri+orientShiftLimit)+iCurve];
			/*For length = 0*/
			pPrevPage=M1map[ori1];
			for(here=0;here<imgSize;++here)
				pPage[here]=pLambda[totalSketch]*pPrevPage[here];
			
			/*For length = 1:len0*/
			for(iLen=1;iLen<=len0;++iLen){
				if(shiftRange(startDx[iCache][iLen],startDy[iCache][iLen],
							  endDx[iCache][iLen],endDy[iCache][iLen],
							  &xStart,&xEnd,&yStart,&yEnd)<0) 
					continue;
				dHere_s = startDy[iCache][iLen]*sx + startDx[iCache][iLen];
				dHere_e = endDy[iCache][iLen]*sx + endDx[iCache][iLen];
				pStartM1 = M1map[startO[iCache][iLen]];
				pEndM1= M1map[endO[iCache][iLen]];
				assert(xStart>0);
				assert(yStart>0);
				assert(xEnd<=sx);
				assert(yEnd<=sy);
				lambda_s= pLambda[totalSketch+iLen*2-1];
				lambda_o= pLambda[totalSketch+iLen*2];
				for (x=xStart; x<xEnd; ++x){
					for (y=yStart; y<yEnd; ++y){
						here=px(x,y,sx,sy);
						assert(here+dHere_s>0);
						assert(here+dHere_s<sx*sy);
						assert(here+dHere_e>0);
						assert(here+dHere_e<sx*sy);
						pPage[here]=pPage[here] 
						+ lambda_s*pStartM1[here+dHere_s]
						+ lambda_o*pEndM1[here+dHere_e]; 
					}
				}
			}
			
			
		}
		totalSketch+=2*len0+1;
	}
	
}


/*refresh the M2map to negmax*/
void resetM2()
{
	float *pM2;
	int imgSize = sx*sy;
	int here;
	int iDeform,iCurve;
	for(iDeform=0;iDeform<nDeform;++iDeform){
		for(iCurve=0;iCurve<nCurves;++iCurve){
			pM2 = M2map[iCurve+iDeform*nCurves];
			for(here=0;here<imgSize;++here)
				pM2[here]=NEGMAX;
		}
		
	}
	
    
}

/*shift: translation*/
/*deform: rotation*/
MAX2(){
	float *pMAX,*pCurr;
	int orient,dOri2,curOri;
	int dLoc2;
        int iCurve, iDeform,iShift;
        int dx,dy,here,dHere;
        int x,y;
	float r,localMax;
        #ifdef USE_SIMD
         int startx,starty,endx,endy;
         __m128 maxDest,maxSource;
        #endif
	
	for(iCurve=0;iCurve<nCurves;++iCurve){
		orient =selCurves[iCurve];/*orient is the fist column*/
		pMAX = M2map[iCurve];
		iShift = 0;iDeform=0;
		/*iDeform*/
		for(dOri2 = -orientShiftLimit; dOri2<=orientShiftLimit; ++dOri2){
			pCurr=S2map[iCurve+iDeform*nCurves];
			/*iShift*/
			for(dLoc2 = locationShiftLimit; 
				dLoc2>=-locationShiftLimit;--dLoc2){
				
				dx = xShift[orient][iShift];
				dy = yShift[orient][iShift];
				dHere = dy*sx+dx;
				
#ifndef USE_SIMD
				
				for (y = MAX(h+1,-dy+1); y<=MIN(sy-h,sy-dy); ++y){
					for (x = MAX(h+1,-dx+1); x<=MIN(sx-h,sx-dx); ++x){
						here = px(x,y,sx,sy);
						localMax=pMAX[here];
						r = pCurr[here+dHere];
						pMAX[here]=MAX(r,localMax);
						
					}
				}
				
#else
				starty = MAX(h+1,-dy+1); startx = MAX(h+1,-dx+1);
				endy = MIN(sy-h,sy-dy); endx = MIN(sx-h,sx-dx);
				for(y = starty; y<=endy; ++y){
					for (x =startx; x<=endx-3; x = x+4){
						here = px(x,y,sx,sy);
						maxDest = _mm_loadu_ps(&(pMAX[here]));
						maxSource=_mm_loadu_ps(&(pCurr[here+dHere]));
						maxDest=_mm_max_ps(maxDest,maxSource);
						_mm_storeu_ps(&(pMAX[here]),maxDest);                           
					}
					x = x-4;
					for(;x<=endx;++x){
						here = px(x,y,sx,sy);
						localMax=pMAX[here];
						r = pCurr[here+dHere];
						pMAX[here]=MAX(r,localMax);
					}
				}
#endif
				iShift++;
			}/*iShift*/
			iDeform++;
		}/*iDeform*/
		
	}/*iCurve*/
	
}



/* mex function is used to pass on the pointers and scalars from matlab, 
 so that heavy computation can be done by C, which puts the results into 
 some of the pointers. After that, matlab can then use these results. 
 
 So matlab is very much like a managing platform for organizing the 
 experiments, and mex C is like a work enginee for fast computation. */

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
	int o, i, j, c,iLen,iImg,iPage; 
	mxArray *f;  
	
	c = 0; /* counter for input pointers and scalars */
	
	norient = floor(mxGetScalar(prhs[c++])+.5);  /* number of orientations */
	nLength = floor(mxGetScalar(prhs[c++])+.5);  /* number of line segments */
	/*paratmeres for curve*/
	nAngle =floor(mxGetScalar(prhs[c++])+.5);
    cAngleLimit= (nAngle-1)/2;
	
	
	
	locationShiftLimit = floor(mxGetScalar(prhs[c++])+.5);  
	orientShiftLimit = floor(mxGetScalar(prhs[c++])+.5); 
	
	nDeform = 2*orientShiftLimit+1;

        nCurves = floor(mxGetScalar(prhs[c++])+.5);
        selCurves =mxGetPr(prhs[c++]);
	
	M1map = mxCalloc(norient, sizeof(float*));   /* M1 maps */
	for (i=0; i<norient; i++)
	{
		f = mxGetCell(prhs[c], i); 
		M1map[i] = (float*)mxGetPr(f);    /* get pointers to filtered images */       
	}
	c++;
	
	S2map = mxCalloc(nCurves*nDeform, sizeof(float*));   /* S2 maps */
	for (i=0; i<nDeform; i++){
		for (j=0; j<nCurves; j++){
			f = mxGetCell(prhs[c], i*nCurves+j); 
			S2map[i*nCurves+j] =(float*) mxGetPr(f);    /* get pointers to filtered images */       
		}
	}
	c++;
	
	
	M2map = mxCalloc(nCurves, sizeof(float*));   /* S2 maps */
	
	for (j=0; j<nCurves; j++){
		f = mxGetCell(prhs[c], j); 
		M2map[j] =(float*) mxGetPr(f);    /* get pointers to filtered images */       
	}
	c++;
	
        sx = floor(mxGetScalar(prhs[c++])+.5); 
	sy = floor(mxGetScalar(prhs[c++])+.5); 
	h =  floor(mxGetScalar(prhs[c++])+.5); 
	
	pLambda = mxGetPr(prhs[c++]);
	
	/*main functions*/
	curve_shift();
	StoreShift2();
	Csum2LC();     
	MAX2();
	
	free_matrix((void**)startDx,norient*nAngle, nLength);
	free_matrix((void**)startDy,norient*nAngle, nLength);
	
	free_matrix((void**)endDx,norient*nAngle, nLength);
	free_matrix((void**)endDy,norient*nAngle, nLength);
	
	free_matrix((void**)startO,norient*nAngle, nLength);
	free_matrix((void**)endO,norient*nAngle, nLength);
	
	
	free_matrix((void*)xShift,norient,numShift);
	free_matrix((void*)yShift,norient,numShift);
	free_matrix((void*)orientShifted,norient,numShift);
	
	
}
