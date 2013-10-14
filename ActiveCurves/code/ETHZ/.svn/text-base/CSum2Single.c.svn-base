

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
float **M1map, **S2map;                              
int **lsinsh, **lcossh;


int cLengthLimit,cAngleLimit;
int **startDx,**startDy,**endDx,**endDy;
int **startO,**endO;
int nLength,nAngle;
int nDims[4];
double lambda,logZ;

/* getting the index of matlab image, (x,y) location, (sx, sy) sizes */
 __inline int px(int x, int y, int bx, int by)    
{            
	return (x + (y-1)*bx - 1); 
}

/*get the first index of the SUM2maps*/
 __inline int page(int iOrient, int iAngle, int iLength, int iImage){
	return iImage+iLength*nDims[1]+iAngle*nDims[2]+iOrient*nDims[3];
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
void Csum2LC(int iImg)
{
	int x, y, here;
	int dHere_s,dHere_e;
	int xStart,yStart,xEnd,yEnd;
	int iOri,iLen,iAng;
	int iPage,iCache;
	float *pPage,*pPrevPage,*pStartS2,*pEndS2;
	
	int nAngle = 2*cAngleLimit+1;
	int imgSize = sx*sy;
	/*initialize SUM2maps for curve of length 0*/
	for (iOri=0; iOri<norient;++iOri){
		for (iAng = 0; iAng<nAngle;++iAng){
			iPage = page(iOri, iAng, 0, iImg);
			pPage = S2map[iPage];
			pPrevPage = M1map[iOri*nimage+iImg];
			for(here = 0 ; here<sx*sy; ++here)
				pPage[here]=lambda*pPrevPage[here]-logZ;
	     }
	}

	for (iOri=0; iOri<norient; ++iOri){
		for (iAng = 0; iAng<nAngle;++iAng){
			iCache = iOri*nAngle+iAng;
			for(iLen = 1; iLen<nLength; ++iLen){
				pPrevPage = S2map[page(iOri,iAng,iLen-1,iImg)];
				pPage = S2map[page(iOri,iAng,iLen,iImg)];

				if(shiftRange(startDx[iCache][iLen],startDy[iCache][iLen],endDx[iCache][iLen],endDy[iCache][iLen],
						   &xStart,&xEnd,&yStart,&yEnd)<0) 
						   continue;
				dHere_s = startDy[iCache][iLen]*sx + startDx[iCache][iLen];
				dHere_e = endDy[iCache][iLen]*sx + endDx[iCache][iLen];
				pStartS2 = S2map[page(startO[iCache][iLen],iAng,0,iImg)];
				pEndS2 = S2map[page(endO[iCache][iLen],iAng,0,iImg)];
				assert(xStart>0);
				assert(yStart>0);
				assert(xEnd<=sx);
				assert(yEnd<=sy);
				for (x=xStart; x<xEnd; ++x){
					for (y=yStart; y<yEnd; ++y){
						here=px(x,y,sx,sy);
						assert(here+dHere_s>0);
						assert(here+dHere_s<sx*sy);
						assert(here+dHere_e>0);
						assert(here+dHere_e<sx*sy);
						pPage[here]=pPrevPage[here] 
						+ pStartS2[here+dHere_s]
						+ pEndS2[here+dHere_e]; 
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
	int o, i, j, c,iLen,iImg,iPage; 
	mxArray *f;  
	
	c = 0; /* counter for input pointers and scalars */
	nimage = floor(mxGetScalar(prhs[c++])+.5);
	norient = floor(mxGetScalar(prhs[c++])+.5);  /* number of orientations */
	nLength = floor(mxGetScalar(prhs[c++])+.5);  /* number of line segments */
	/*paratmeres for curve*/
	cAngleLimit =floor(mxGetScalar(prhs[c++])+.5);
	nAngle=2*cAngleLimit+1;
	/*parameters for page number*/
	nDims[0]= 1;
	nDims[1]= nDims[0]*nimage;
	nDims[2]= nDims[1]*nLength;
	nDims[3]= nDims[2]*(2*cAngleLimit+1);
	
	
	
	M1map = mxCalloc(nimage*norient, sizeof(float*));   /* M1 maps */
	for (i=0; i<nimage; i++)
	{
		for (o=0; o<norient; o++)
		{  
			f = mxGetCell(prhs[c], o*nimage+i); 
			M1map[o*nimage+i] = (float*)mxGetPr(f);    /* get pointers to filtered images */       
		}
	}
	c++;
	
	S2map = mxCalloc(nimage*norient*nLength*nAngle, sizeof(float*));   /* M1 maps */
	for (i=0; i<norient; i++){
		for (j=0; j<nAngle; j++){
			for(iLen = 0; iLen<nLength; ++iLen){
				for(iImg = 0 ; iImg<nimage; ++iImg){
					iPage = page(i,j,iLen,iImg);
					f = mxGetCell(prhs[c], iPage); 
					S2map[iPage] =(float*) mxGetPr(f);    /* get pointers to filtered images */       
				}
			}
		}
	}
	c++;
	
	sx = floor(mxGetScalar(prhs[c++])+.5); 
	sy = floor(mxGetScalar(prhs[c++])+.5); 
	h =  floor(mxGetScalar(prhs[c++])+.5); 
	
	lambda = mxGetScalar(prhs[c++]);
	logZ = mxGetScalar(prhs[c++]);
	
	
	/*main functions*/
	curve_shift();
	for (i=0; i<nimage; i++)
		Csum2LC(i);     
	
	
	free_matrix((void**)startDx,norient*nAngle, nLength);
	free_matrix((void**)startDy,norient*nAngle, nLength);
	
	free_matrix((void**)endDx,norient*nAngle, nLength);
	free_matrix((void**)endDy,norient*nAngle, nLength);
	
	free_matrix((void**)startO,norient*nAngle, nLength);
	free_matrix((void**)endO,norient*nAngle, nLength);
	
	
}







