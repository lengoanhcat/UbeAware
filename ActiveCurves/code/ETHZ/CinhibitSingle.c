/**************************************************
 Mex C Code for reproducing experiment 1
***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        /* the algorithm is connect to matlab */
# include "math.h"
# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define SIGN(x) ((x)>0? 1: -1)
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))
# define NEGMAX -1e10

double *double_vector(int n)
{
    double *v; 
    
    v = (double*) mxCalloc (n, sizeof(double));
    
    return v; 
}

int *int_vector(int n)
{
    int *v; 
    
    v = (int*) mxCalloc (n, sizeof(int));
    
    return v; 
}

double **double_matrix(int m, int n)
{
    double **mat; 
    int i; 
    
    mat = (double**) mxCalloc(m, sizeof(double*)); 
    for (i=0; i<m; i++)
        mat[i] = double_vector(n); 
    
    return mat; 
}

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


/* getting the index of matlab image, (x,y) location, (sx, sy) sizes */
__inline  int px(int x, int y, int bx, int by)    
{            
	return (x + (y-1)*bx - 1); 
}


 /* variables */
 int nimage, norient, nbasis, sx, sy, h, Lrange, Orange;                      
 float **S1map, **C, *ret, *sym, **Asym, **allsymbol;
 double   *maxi, *imaxi;                              
 int **lsinsh, **lcossh, **sinsh, **cossh;
 int nLength,nAngle;
int **startDx,**startDy,**endDx,**endDy;
int **startO,**endO;
/* store the shift values so that we do not need to repeat the sin and cos computation */ 
void curve_shift()
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


void storeshift()
{
    int o, l;
    float theta; 
    
    sinsh = int_matrix(norient, Lrange+Lrange+1); 
    cossh = int_matrix(norient, Lrange+Lrange+1); 
    for (o=0; o<norient; o++)        
    {
        theta = PI*o/norient; 
        for (l=-Lrange; l<=Lrange; l++)
         {
            sinsh[o][l+Lrange] = floor(l*sin(theta)+.5); 
            cossh[o][l+Lrange] = floor(l*cos(theta)+.5); 
        }   
    }
}
    
/* for Gabor(x, y, orientation = ind), find local maximum in image i */
float shiftmax(int i, int ind, int x, int y, int *rm) 
{
   float m;
   int x1, y1, l, ml, mx, my, de, d, d1, here, mo, md; 
 
   m = NEGMAX;  
   for (l=-Lrange; l<=Lrange; l++)   
     {
        x1 = x + cossh[ind][l+Lrange]; 
        y1 = y + sinsh[ind][l+Lrange];   /* shifting in normal direction */
        if ((x1>=1)&&(x1<=sx)&&(y1>=1)&&(y1<=sy))
        {
        here = px(x1, y1, sx, sy);
        for (de=-Orange; de<=Orange; de++)
            {
              d = de+ind;
              d1 = d;
              if (d<0)
                 d1 = d+norient;
              if (d>=norient)
                  d1 = d-norient;   /* shifting in orientation */
              if (S1map[d1*nimage+i][here]>m)
                {
                   m = S1map[d1*nimage+i][here];   /* local maximization */
                   ml = l; mx = x1; my = y1; mo = d1; md = de; 
                }
             }
        }
     }
   rm[0] = ml; rm[1] = mx; rm[2] = my; rm[3] = mo; rm[4] = md; 
   return(m); 
}


/* the Gabor(mx, my, orientation = mi) inhibits overlapping Gabors for image i */
void inhibit(int i, int mi, int mx, int my) 
{
   int x, y, o; 
   float *f, *fc; 
   
   for (o=0; o<norient; o++)   
     {
       f = S1map[o*nimage+i];   
       fc = C[mi+o*norient];   /* inhibition coefficients, zero = inhibit */
       for (x=MAX(1, mx-2*h); x<=MIN(sx, mx+2*h); x++)
         for (y=MAX(1, my-2*h); y<=MIN(sy, my+2*h); y++)
         {
          f[px(x, y, sx, sy)] *= fc[px(x-mx+2*h+1, y-my+2*h+1, 4*h+1, 4*h+1)];
         }
      }
}


/* draw the symbol of Gabor(x0, y0, orientation = ind) with intensity w */
void draw(float *symo, int x0, int y0, int ind, float w)
{
  int x, y; 
  float a; 
          
  for (x=x0-h; x<=x0+h; x++)
     for (y=y0-h; y<=y0+h; y++)
       {
        if ((x>=1)&&(x<=sx)&&(y>=1)&&(y<=sy))
        {
         a = allsymbol[ind][px(x-x0+h+1, y-y0+h+1, 2*h+1, 2*h+1)]*w; 
         if (symo[px(x, y, sx, sy)]<a)
             symo[px(x, y, sx, sy)] = a;
        }
       }
}

void Cshutdown()
{
    int i, ml, cx, cy, mi, l, al, sl, mx, my, rm[5]; 
    float mm; 
    int iOri,iAng,iLen,iCache;
	int len,ori,ang;
    ori = floor(maxi[0]+.5); 
    ang = floor(maxi[1]+.5);
	len = floor(maxi[2]+.5); 
    cx = floor(maxi[3]+.5); 
    cy = floor(maxi[4]+.5); 
    
    for(iLen = 0;iLen<=len;++iLen){
		iCache = ori*nAngle+ang;
        mx = cx+startDx[iCache][iLen];
		my = cy+startDy[iCache][iLen]; 
        mi = startO[iCache][iLen];
		draw(sym, mx, my, mi, 1);/*sqrt(maxi[5]));*/
		
        mx = cx+endDx[iCache][iLen];
		my = cy+endDy[iCache][iLen];
		mi = endO[iCache][iLen];
        draw(sym, mx, my, mi, 1); /*sqrt(maxi[5]));*/
    }
    
    for (i=0; i<nimage; i++)
    {
		ori = floor(imaxi[0+i*6]+.5); 
		ang = floor(imaxi[1+i*6]+.5);
		len = floor(imaxi[2+i*6]+.5); 
		cx = floor(imaxi[3+i*6]+.5); 
		cy = floor(imaxi[4+i*6]+.5);  
       for(iLen = 0;iLen<=len;++iLen){
		   
		iCache = ori*nAngle+ang;
		mx = cx+startDx[iCache][iLen];
		my = cy+startDy[iCache][iLen]; 
		mi = startO[iCache][iLen];
		mm = shiftmax(i, mi, mx, my, rm);
        /* draw(Asym[i], mx, my, mi, sqrt(imaxi[i*5])); */
        draw(Asym[i], rm[1], rm[2], rm[3], mm/6.0f); 
        inhibit(i, rm[3], rm[1], rm[2]);
		   
		mx = cx+endDx[iCache][iLen];
		my = cy+endDy[iCache][iLen]; 
		mi = endO[iCache][iLen];
		mm = shiftmax(i, mi, mx, my, rm);
		/* draw(Asym[i], mx, my, mi, sqrt(imaxi[i*5])); */
		draw(Asym[i], rm[1], rm[2], rm[3], mm/6.0f); 
		inhibit(i, rm[3], rm[1], rm[2]);
		   
        }
    }
}
        
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
 int o, j, i, c; 
 mxArray *f;  
 
 c = 0; /* counter for input pointers and scalars */
 nimage = floor(mxGetScalar(prhs[c++])+.5); 
 norient = floor(mxGetScalar(prhs[c++])+.5);  /* number of orientations */
 nLength = floor(mxGetScalar(prhs[c++])+.5); /* line length*/
 nAngle = floor(mxGetScalar(prhs[c++])+0.5);
 
 S1map = mxCalloc(nimage*norient, sizeof(float*));   /* M1 maps */
 for (i=0; i<nimage; i++)
 {
   for (o=0; o<norient; o++)
      {  
       f = mxGetCell(prhs[c], o*nimage+i); 
       S1map[o*nimage+i] = mxGetPr(f);    /* get pointers to filtered images */       
      }
 }
 c++;

 C = mxCalloc(norient*norient, sizeof(float*));    /* C: correlation/inhibition between filters */
 for (o=0; o<norient; o++)
     {  
       for (j=0; j<norient; j++)
        {
         f = mxGetCell(prhs[c], j*norient+o); 
         C[j*norient+o] = mxGetPr(f);         /* get correlation/inhibition */
        }   
     }
 c++; 
 
 sx = floor(mxGetScalar(prhs[c++])+.5); 
 sy = floor(mxGetScalar(prhs[c++])+.5); 
 h =  floor(mxGetScalar(prhs[c++])+.5); 
 Lrange = floor(mxGetScalar(prhs[c++])+.5);     /* range of shifting along normal direction of Gabor */
 Orange = floor(mxGetScalar(prhs[c++])+.5);   /* range of shifting in orientation */
 
 maxi= mxGetPr(prhs[c++]);   
 imaxi= mxGetPr(prhs[c++]);   
 
 sym= mxGetPr(prhs[c++]);   
 Asym = mxCalloc(nimage, sizeof(float*));   /* symbols of active Gabors for an individual image */
 for (i=0; i<nimage; i++)
    {
        f = mxGetCell(prhs[c], i);
        Asym[i] = mxGetPr(f);  
     }
 c++; 
 
 allsymbol = mxCalloc(norient, sizeof(float*));    
 for (o=0; o<norient; o++)
     {  
       f = mxGetCell(prhs[c], o); 
       allsymbol[o] = mxGetPr(f);  /* symbols of filters */         
     }
 c++; 
 
 storeshift(); 
 curve_shift(); 
 Cshutdown();     
}



 

                    