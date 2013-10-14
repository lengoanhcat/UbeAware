/**************************************************
 * Mex C Code for backtracking and inhibition of arc
 x-axis: from top-left to downside, start from 1
 y-axis: from top-left to rightside, start from 1
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
# define ROUND(x) (floor((x)+.5))

 

/*allocate int sized vector of size n */
int *int_vector(int n) {
    int *v;
    
    v = (int*) mxCalloc(n, sizeof(int));
    
    return v;
}


/*allocate int sized matrix of size m by n */
int **int_matrix(int m, int n) {
    int **mat;
    int i;
    
    mat = (int**) mxCalloc(m, sizeof(int*));
    for (i=0; i<m; i++)
        mat[i] = int_vector(n);
    
    return mat;
}

/*free memory for passed matrix a,
 the matrix is of size nrow by ncol*/
void free_matrix(void **a, int nrow, int ncol) {
    int i;
    for(i = 0; i < nrow; i++)
        mxFree(a[i]);
    mxFree(a);
}


/* getting the index of matlab image, (x,y) location, (bx, by) sizes */
__inline  int px(int x, int y, int bx, int by) {
    return (x + (y-1)*bx - 1);
}


/* variables */
int nimage, norient, nbasis, sx, sy, h, Lrange, Orange;
float **S1map, **C, *sym, **Asym, **allsymbol;
double   *maxi, *imaxi; /* prototype arc and deformed arc index */
int **lsinsh, **lcossh, **sinsh, **cossh; /*Gabor primitive deformation range*/
int nLength, nAngle; /*number of arc lengthes and curvatures*/
int **startDx, **startDy, **endDx, **endDy; /*arc templates, for Gabor primitive position*/
int **startO, **endO;/*arc templates, for Gabor primitive orientation*/
int is_inhibit;  /*if perform inhibition. if 0, only draw the arc*/
float *primResp; /*response of constinutent Gabor primitives */

/* store the constituent Gabor positon and orientaions of arcs centered at (0,0)
   for arcs centered at (x,y), just shift Gabor positions relatively*/
void arcShift() {
    int iOri, iAng, iPage;
    int dAngle, l;
    float theta, dTheta;
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

/*Store local deformation range of Gobor primitives, in order to save computation*/
void storeShift() {
    int o, l;
    double theta;
    
    sinsh = int_matrix(norient, Lrange+Lrange+1);
    cossh = int_matrix(norient, Lrange+Lrange+1);
    for (o=0; o<norient; o++) {
        theta = PI*o/norient;
        for (l=-Lrange; l<=Lrange; l++) {
            sinsh[o][l+Lrange] = ROUND(l*sin(theta));
            cossh[o][l+Lrange] = ROUND(l*cos(theta));
        }
    }
}

/* for Gabor(x, y, orientation = ind), find local maximum in image i */
float shiftmax(int i, int ind, int x, int y, int *rm) {
    float m;
    int x1, y1, l, ml, mx, my, de, d, d1, here, mo, md;
    mx = x; my=y;mo = ind;
     
    m = NEGMAX;
    for (l=-Lrange; l<=Lrange; l++) {
        x1 = x + cossh[ind][l+Lrange];
        y1 = y + sinsh[ind][l+Lrange];   /* shifting in normal direction */
        if ((x1>=1)&&(x1<=sx)&&(y1>=1)&&(y1<=sy)) {
            here = px(x1, y1, sx, sy);
            for (de=-Orange; de<=Orange; de++) {
                d = de+ind;
                d1 = d;
                if (d<0)
                    d1 = d+norient;
                if (d>=norient)
                    d1 = d-norient;   /* shifting in orientation */
                if (S1map[d1*nimage+i][here]>m) {
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
void inhibit(int i, int mi, int mx, int my) {
    int x, y, o;
    float *f, *fc;
    
    for (o=0; o<norient; o++) {
        f = S1map[o*nimage+i];
        fc = C[mi+o*norient];   /* inhibition coefficients, zero = inhibit */
        for (x=MAX(1, mx-2*h); x<=MIN(sx, mx+2*h); x++)
            for (y=MAX(1, my-2*h); y<=MIN(sy, my+2*h); y++) {
				if (fc[px(x-mx+2*h+1, y-my+2*h+1, 4*h+1, 4*h+1)]<1e-6)
            		f[px(x, y, sx, sy)] =0.0f ;
            }
    }
}


/* draw the symbol of Gabor(x0, y0, orientation = ind) with intensity w */
void draw(float *symo, int x0, int y0, int ind, float w) {
    int x, y;
    float a;
    
    for (x=x0-h; x<=x0+h; x++)
        for (y=y0-h; y<=y0+h; y++) {
        if ((x>=1)&&(x<=sx)&&(y>=1)&&(y<=sy)) {
            a = allsymbol[ind][px(x-x0+h+1, y-y0+h+1, 2*h+1, 2*h+1)];
			a = a * w;
            if (symo[px(x, y, sx, sy)]<a)         
		     symo[px(x, y, sx, sy)] = a;
        }
    }
}


/*find deformed positions of each Gabor Primive in an arc*/
void retrievAndShutdown() {
    int i, ml, cx, cy, mi, l, al, sl, mx, my, rm[5];
    float mm;
    int iOri, iAng, iLen, iCache;
    int len, ori, ang;    
    int gabor_index=0;
    int **posList; /*save index of each deformed Gabor*/
    
    /*get the curve index*/
    ori = floor(maxi[0]+.5);
    ang = floor(maxi[1]+.5);
    len = floor(maxi[2]+.5);
    cx = floor(maxi[3]+.5);
    cy = floor(maxi[4]+.5);
	
    posList =int_matrix(2*len+1,3);
    /*draw the shared template*/
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
    
    /*draw deformed template*/
    for (i=0; i<nimage; i++) {
        ori = floor(imaxi[0+i*6]+.5);
        ang = floor(imaxi[1+i*6]+.5);
        len = floor(imaxi[2+i*6]+.5);
        cx = floor(imaxi[3+i*6]+.5);
        cy = floor(imaxi[4+i*6]+.5);
        iCache = ori*nAngle+ang;
        
        gabor_index=0;
        mx = cx+startDx[iCache][0];
        my = cy+startDy[iCache][0];
        mi = startO[iCache][0];
        mm = shiftmax(i, mi, mx, my, rm);
        gabor_index = 0*2+i*(2*len+1);
        primResp[gabor_index]=mm;
        posList[gabor_index][0]=rm[1];
        posList[gabor_index][1]=rm[2];
        posList[gabor_index][2]=rm[3];
        
        draw(Asym[i], rm[1], rm[2], rm[3], mm/6.0f);
        for(iLen = 1;iLen<=len;++iLen){
			/*Gabor in one end*/
            iCache = ori*nAngle+ang;
            mx = cx+startDx[iCache][iLen];
            my = cy+startDy[iCache][iLen];
            mi = startO[iCache][iLen];
            mm = shiftmax(i, mi, mx, my, rm);
            gabor_index = iLen*2-1;
            primResp[gabor_index+i*(2*len+1)]=mm;
            posList[gabor_index][0]=rm[1];
            posList[gabor_index][1]=rm[2];
            posList[gabor_index][2]=rm[3];
            draw(Asym[i], rm[1], rm[2], rm[3], mm/6.0f);

            /*Gabor in another end*/
            mx = cx+endDx[iCache][iLen];
            my = cy+endDy[iCache][iLen];
            mi = endO[iCache][iLen];
            mm = shiftmax(i, mi, mx, my, rm);
            gabor_index = iLen*2;
            primResp[gabor_index+i*(2*len+1)]=mm;
            posList[gabor_index][0]=rm[1];
            posList[gabor_index][1]=rm[2];
            posList[gabor_index][2]=rm[3];            
            draw(Asym[i], rm[1], rm[2], rm[3], mm/6.0f);            
        }
        if( is_inhibit){
        for( gabor_index = 0; gabor_index <2*len+1; ++gabor_index){
            inhibit(i, posList[gabor_index][2],
                       posList[gabor_index][0],
                       posList[gabor_index][1]);
        }}
    }
    free_matrix((void**)posList,2*len+1,3);
}



/*Grab parameters and maps pointers*/
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    int o, j, i, c;
    mxArray *f;
    
    c = 0; /* counter for input pointers and scalars */
    nimage = floor(mxGetScalar(prhs[c++])+.5);  /*number of images*/
    norient = floor(mxGetScalar(prhs[c++])+.5);  /* number of orientations */
    nLength = floor(mxGetScalar(prhs[c++])+.5);  /* number of different arc length*/
    nAngle = floor(mxGetScalar(prhs[c++])+0.5); /*number of different arc curvature*/
    
    S1map = mxCalloc(nimage*norient, sizeof(float*));   /* S1 maps */
    for (i=0; i<nimage; i++) {
        for (o=0; o<norient; o++) {
            f = mxGetCell(prhs[c], o*nimage+i);
            S1map[o*nimage+i] = (float*)mxGetPr(f);    /* get pointers to filtered images */
        }
    }
    c++;
    
    C = mxCalloc(norient*norient, sizeof(float*));    /* C: correlation/inhibition between filters */
    for (o=0; o<norient; o++) {
        for (j=0; j<norient; j++) {
            f = mxGetCell(prhs[c], j*norient+o);
            C[j*norient+o] = (float*)mxGetPr(f);         /* get correlation/inhibition */
        }
    }
    c++;
    sx = floor(mxGetScalar(prhs[c++])+.5); /*height of image*/
    sy = floor(mxGetScalar(prhs[c++])+.5); /*width of image*/
    h =  floor(mxGetScalar(prhs[c++])+.5); /*half length of Gabor primitive*/
    Lrange = floor(mxGetScalar(prhs[c++])+.5);     /* range of shifting along normal direction of Gabor */
    Orange = floor(mxGetScalar(prhs[c++])+.5);   /* range of shifting in orientation */
    
    maxi= mxGetPr(prhs[c++]); /*index of prototype arc*/
    imaxi= mxGetPr(prhs[c++]); /*index of deformed arc, here it is the same as maxi*/
    
    primResp = (float*)mxGetPr(prhs[c++]); /*pointer to individual gabor primtive response*/
    
    sym= (float*)mxGetPr(prhs[c++]);
    Asym = (float**)mxCalloc(nimage, sizeof(float*));   /* symbols of active Gabors for an individual image */
    for (i=0; i<nimage; i++) {
        f = mxGetCell(prhs[c], i);
        Asym[i] = (float*)mxGetPr(f);
    }
    c++;
    
    allsymbol = mxCalloc(norient, sizeof(float*));
    for (o=0; o<norient; o++) {
        f = mxGetCell(prhs[c], o);
        allsymbol[o] =(float*) mxGetPr(f);  /* symbols of Gabor primitives */
    }
    c++;
	

    is_inhibit=1;/*inhibit by default, unless extra arguments are provided*/
    if(nrhs=c+1){
		is_inhibit=floor(mxGetScalar(prhs[c++])+.5);
    }
    
	/*pre-compute Gabor Primitive active range*/
    storeShift();
	/*pre-compute Gabor indexes of arcs centered at point (0,0)*/
    arcShift();

	/*main routine*/
    retrievAndShutdown();
	
	/*release pre-computed results*/
	free_matrix((void**)startDx,norient*nAngle, nLength);
	free_matrix((void**)startDy,norient*nAngle, nLength);
	
	free_matrix((void**)endDx,norient*nAngle, nLength);
	free_matrix((void**)endDy,norient*nAngle, nLength);
	
	free_matrix((void**)startO,norient*nAngle, nLength);
	free_matrix((void**)endO,norient*nAngle, nLength);
	free_matrix((void**)sinsh,norient,2*Lrange+1);
    free_matrix((void**)cossh,norient,2*Lrange+1);
}
