/**************************************************
 Mex-C code performing S2 operation
 x-axis: from top-left to downside, start from 1
 y-axis: from top-left to rightside, start from 1
 **************************************************/
# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        /* the algorithm is connect to matlab */
# include "math.h"
# include <mmintrin.h> /*headers for SIMD intrinsics*/
# include <xmmintrin.h>

# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))
# define NEGMAX -1e10


#define USE_SIMD /*enable SIMD*/
/* Generating  int vector */
int *int_vector(int n)
{
    int *v; 
    v = (int*) mxCalloc (n, sizeof(int));
    return v;
}
/* Generating float vector */
float *single_vector(int n)
{
    float *v; 
    v = (float*) mxCalloc (n, sizeof(float));    
    return v; 
}
/* Generating float matrix */
float **single_matrix(int m, int n)
{
    float **mat; 
    int i; 
    mat = (float**) mxCalloc(m, sizeof(float*)); 
    for (i=0; i<m; i++)
        mat[i] = single_vector(n); 
    return mat; 
}
/* Generate int matrix */
int **int_matrix(int m, int n)
{
    int **mat;
    int i;
    mat = (int**) mxCalloc(m, sizeof(int*));
    for (i=0; i<m; i++)
        mat[i] = int_vector(n);
    return mat;
}

/*free matrix space*/
void free_matrix(void **a, int nrow, int ncol)
{
	int i;
	for(i = 0; i < nrow; i++)
		mxFree(a[i]);
	mxFree(a);
}



/* variables */
int nimage, norient, sx, sy, h;
float **M1maps, **S2maps;                              
int **startDx,**startDy,**endDx,**endDy;/*arc templates, for Gabor primitive position*/
int **startO,**endO;/*arc templates, for Gabor primitive orientation*/
int nLength,nAngle;/*number of arc lengthes and curvatures*/
int nDims[4];
double lambda,logZ; /*model parameter*/

/* getting the index of matlab image, (x,y) location, (sx, sy) sizes */
 __inline int px(int x, int y, int bx, int by)    
{            
	return (x + (y-1)*bx - 1); 
}

/*get the cell index of the SUM2maps*/
 __inline int page(int iOrient, int iAngle, int iLength, int iImage){
	return iImage+iLength*nDims[1]+iAngle*nDims[2]+iOrient*nDims[3];
}


/* store the constituent Gabor positon and orientaions of arcs centered at (0,0)
 for arcs centered at (x,y), just shift Gabor positions relatively*/
void arcShift() 
{
    int iOri,iAng,iPage;
	int dAngle,l;
	float theta,dTheta;
	int cAngleLimit = (nAngle-1)/2;

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

	
/*
 * decide pssible positions of an arc
 * x1,y1,x2,y2: bounding box of the current curve
 * top_x,bot_x,left_y,right_y: rectanguar area of possible positions
 */
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

/* perfoms Sum2 operations */
void Csum2LC(int iImg)
{
	int x, y, here;
	int dHere_s,dHere_e;
	int xStart,yStart,xEnd,yEnd;
	int iOri,iLen,iAng;
	int iPage,iCache;
	float *pPage,*pPrevPage,*pStartS2,*pEndS2;
	int imgSize = sx*sy;
	
#ifdef USE_SIMD
    __m128 sumDest,sumSource1,sumSource2,sumSource3;
#endif
	
	/*initialize SUM2maps for arcs of length 0*/
	for (iOri=0; iOri<norient;++iOri){
		for (iAng = 0; iAng<nAngle;++iAng){
			iPage = page(iOri, iAng, 0, iImg);
			pPage = S2maps[iPage];
			pPrevPage = M1maps[iOri*nimage+iImg];
			for(here = 0 ; here<sx*sy; ++here)
				pPage[here]=lambda*pPrevPage[here]-logZ;
	     }
	}
	/*compute scores for arc of length >0*/
	for (iOri=0; iOri<norient; ++iOri){
		for (iAng = 0; iAng<nAngle;++iAng){
			iCache = iOri*nAngle+iAng;
			for(iLen = 1; iLen<nLength; ++iLen){
				
				pPrevPage = S2maps[page(iOri,iAng,iLen-1,iImg)]; /*map in precedent layer*/
				pPage = S2maps[page(iOri,iAng,iLen,iImg)]; /*current map*/

				/*decide possible arc positions*/
				if(shiftRange(startDx[iCache][iLen],startDy[iCache][iLen],
							  endDx[iCache][iLen],endDy[iCache][iLen],
							  &xStart,&xEnd,&yStart,&yEnd)<0) 
						   continue; /*if not possible to happen on image of current size*/
				
				/*displacements of two Gabors at two ends*/
				dHere_s = startDy[iCache][iLen]*sx + startDx[iCache][iLen];
				dHere_e = endDy[iCache][iLen]*sx + endDx[iCache][iLen];
				pStartS2 = S2maps[page(startO[iCache][iLen],iAng,0,iImg)];
				pEndS2 = S2maps[page(endO[iCache][iLen],iAng,0,iImg)];
#ifndef USE_SIMD
				for (x=xStart; x<xEnd; ++x){
					for (y=yStart; y<yEnd; ++y){
						here=px(x,y,sx,sy);
						pPage[here]=pPrevPage[here] 
						+ pStartS2[here+dHere_s]
						+ pEndS2[here+dHere_e]; 
					}
				}
#else
				for (y=yStart; y<yEnd; ++y){
					for (x=xStart; x<xEnd-3; x+=4){
						here = px(x,y,sx,sy);
						sumDest = _mm_loadu_ps(&(pPage[here]));
						sumSource1 = _mm_loadu_ps(&(pPrevPage[here]));
						sumSource2 = _mm_loadu_ps(&(pStartS2[here+dHere_s]));
						sumSource3 = _mm_loadu_ps(&(pEndS2[here+dHere_e]));
						sumDest = _mm_add_ps(sumSource1,sumSource2);
						sumDest = _mm_add_ps(sumDest,sumSource3);
						_mm_storeu_ps(&(pPage[here]),sumDest);
					}
					x = x-4;
					for(; x<xEnd;++x){
						here = px(x,y,sx,sy);
						pPage[here]=pPrevPage[here] 
						+ pStartS2[here+dHere_s]
						+ pEndS2[here+dHere_e]; 
					}
				}
#endif
			}
		}
	}
	
}






/*Grab parameters and maps pointers*/
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
	int o, i, j, c,iLen,iImg,iPage; 
	mxArray *f;  
	
	c = 0; /* counter for input pointers and scalars */
	nimage = floor(mxGetScalar(prhs[c++])+.5);
	norient = floor(mxGetScalar(prhs[c++])+.5);  /* number of orientations */
	nLength = floor(mxGetScalar(prhs[c++])+.5);  /* number of arc lengthes */
	nAngle =floor(mxGetScalar(prhs[c++])+.5);   /* numer of arc curvatures */
	
	/*parameters for cell number*/
	nDims[0]= 1;
	nDims[1]= nDims[0]*nimage;
	nDims[2]= nDims[1]*nLength;
	nDims[3]= nDims[2]*nAngle;
	
	
	
	M1maps = mxCalloc(nimage*norient, sizeof(float*));   /* M1 maps */
	for (i=0; i<nimage; i++)
	{
		for (o=0; o<norient; o++)
		{  
			f = mxGetCell(prhs[c], o*nimage+i); 
			M1maps[o*nimage+i] = (float*)mxGetPr(f);    /* get pointers to filtered images */       
		}
	}
	c++;
	
	S2maps = mxCalloc(nimage*norient*nLength*nAngle, sizeof(float*));   /* S2 maps */
	for (i=0; i<norient; i++){
		for (j=0; j<nAngle; j++){
			for(iLen = 0; iLen<nLength; ++iLen){
				for(iImg = 0 ; iImg<nimage; ++iImg){
					iPage = page(i,j,iLen,iImg);
					f = mxGetCell(prhs[c], iPage); 
					S2maps[iPage] =(float*) mxGetPr(f);    /* get pointers to filtered images */       
				}
			}
		}
	}
	c++;
	
	sx = floor(mxGetScalar(prhs[c++])+.5); /*image height*/
	sy = floor(mxGetScalar(prhs[c++])+.5); /*image width*/
	h =  floor(mxGetScalar(prhs[c++])+.5); /*half length of Gabor Primitive*/
	
	lambda = mxGetScalar(prhs[c++]); /*model parameter*/
	logZ = mxGetScalar(prhs[c++]); /*model parameter*/
	
	/*main functions*/
	arcShift();
	for (i=0; i<nimage; i++)
		Csum2LC(i);     
	
	/*release memory*/
	free_matrix((void**)startDx,norient*nAngle, nLength);
	free_matrix((void**)startDy,norient*nAngle, nLength);
	
	free_matrix((void**)endDx,norient*nAngle, nLength);
	free_matrix((void**)endDy,norient*nAngle, nLength);
	
	free_matrix((void**)startO,norient*nAngle, nLength);
	free_matrix((void**)endO,norient*nAngle, nLength);
	
	
}