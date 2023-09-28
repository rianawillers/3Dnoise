/* 3Dnoise.c */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#define MAXPATHLEN    1024
#define MAX_FRAMES     100
#define MAXROW         100
#define MAXCOL         156

unsigned long ***ULImage;
float ***Image;
float ***tmpImage;
char  images_filename[MAXPATHLEN];
int   first_frame;
int   last_frame;
int   first_col;
int   last_col;
char  type[MAXPATHLEN];


int LoadImages();
int ProcessImages();
double GetS(int mincol, int maxcol, int maxrow, int maxframe);
double GetNT(int mincol, int maxcol, int maxrow, int maxframe);
double GetNVH(int mincol, int maxcol, int maxrow, int maxframe);
double GetNTV(int mincol, int maxcol, int maxrow, int maxframe);
double GetNTH(int mincol, int maxcol, int maxrow, int maxframe);
double GetNH(int mincol, int maxcol, int maxrow, int maxframe);
double GetNV(int mincol, int maxcol, int maxrow, int maxframe);
double GetNTVH(int mincol, int maxcol, int maxrow, int maxframe);
int oprDT(int mincol, int maxcol, int maxrow, int maxframe);
int oprDV(int mincol, int maxcol, int maxrow, int maxframe);
int oprDH(int mincol, int maxcol, int maxrow, int maxframe);
int oprSDT(int mincol, int maxcol, int maxrow, int maxframe);
int oprSDV(int mincol, int maxcol, int maxrow, int maxframe);
int oprSDH(int mincol, int maxcol, int maxrow, int maxframe);
float ***matrix3d(int nsl,int nsh,int nrl,int nrh,int ncl,int nch);
unsigned long ***matrix3dUL(int nsl,int nsh,int nrl,int nrh,int ncl,int nch);

/******************************************************************/
/*                                                                */
/* main                                                           */
/*                                                                */
/******************************************************************/

int main(int argc, char ** argv )
{

	if (argc < 6)
	{
		printf("\n\nUsage: 3Dnoise2008 filename first_frame last_frame first_column last_column type(float/uint)\n");
		printf("       Please note that counting starts from 1\n\n");
		exit(-1);
	}

	/* Get input filename */

	strcpy(images_filename,argv[1]);

	first_frame = atoi(argv[2]);
	last_frame  = atoi(argv[3]);
	first_col = atoi(argv[4]);
	last_col  = atoi(argv[5]);
	strcpy(type,argv[6]);

	if((last_frame-first_frame+1) > MAX_FRAMES){
		printf("\nAdjusting last_frame from %d ",last_frame);
		last_frame = first_frame + MAX_FRAMES;
		printf("to %d\n\n",last_frame);
	}
	if(first_col < 1) {
		printf("\nAdjusting first_col from %d ",first_col);
		first_col = 1;
		printf("to %d\n\n",first_col);
	}
	if(last_col > MAXCOL) {
		printf("\nAdjusting last_col from %d ",last_col);
		last_col = MAXCOL;
		printf("to %d\n\n",last_col);
	}

	printf("Input file name       : %s\n\n",images_filename);
	printf("Processing %3d frames : %d to %d (counting from 1!)\n",last_frame-first_frame+1,first_frame,last_frame);
	printf("Processing %3d rows   : %d to %d (counting from 1!)\n",MAXROW,1,MAXROW);
	printf("Processing %3d cols   : %d to %d (counting from 1!)\n",(last_col-first_col+1),first_col,last_col);
	if(strcmp(type,"float")==0) 
		printf("type = %s : Float images will be processed\n",type);
	else if(strcmp(type,"uint")==0) 
		printf("type = %s : Unsigned long int images will be processed\n",type);
	else
	{
		printf("Invalid image type - used f=float or i=unsigned long int\n");
		exit(-1);
	}

	/* Allocate memory for image arrays */

	ULImage = matrix3dUL(0,MAX_FRAMES-1,0,MAXROW-1,0,MAXCOL-1);
	Image = matrix3d(0,MAX_FRAMES-1,0,MAXROW-1,0,MAXCOL-1);
	tmpImage = matrix3d(0,MAX_FRAMES-1,0,MAXROW-1,0,MAXCOL-1);

	ProcessImages();

}

/******************************************************************/
/*                                                                */
/* LoadImages                                                     */
/*                                                                */
/* writes a sequence of images to the float image vector          */
/*                                                                */
/******************************************************************/

int LoadImages() 
{

	int  Cnt,i,j;
	int  cntI=0;
	FILE *fp;

	/* Clear the image arrays */

	for(Cnt=0;Cnt<MAX_FRAMES;Cnt++)
		for(i=0;i<MAXROW;i++)
			for(j=0;j<MAXCOL;j++)
				Image[Cnt][i][j] = 0.;

	/* Reload images from disk */

	if((fp = fopen(images_filename,"rb")) == NULL){
		printf("Unable to open input file\n");
		exit(-1);
	}

	for(Cnt=0;Cnt<last_frame;Cnt++)
	{
		for(i=0;i<MAXROW;i++)
		{
			for(j=0;j<MAXCOL;j++)
			{
				if(strcmp(type,"float")==0)
				{
					if (fread(&Image[cntI][i][j], sizeof(float),1, fp) ==0)
					{
						printf("\n\nCannot read image file - check the number of frames specified!\n\n");
						exit(-1);
					}
				}
				else
				{
					if (fread(&ULImage[cntI][i][j], sizeof(unsigned long),1, fp) ==0)
					{
						printf("\n\nCannot read image file - check the number of frames specified!\n\n");
						exit(-1);
					}
					Image[cntI][i][j] = (float)ULImage[cntI][i][j];

				}
			}
		}
		if(Cnt >= (first_frame-1)) cntI++;
	}

	fclose(fp);

	return cntI;
}



/******************************************************************/
/*                                                                */
/* ProcessImages                                                  */
/*                                                                */
/* the row axis is r ranging from 0 to MAXROW                     */
/* the col axis is c ranging from 0 to MAXCOL                     */
/* here we have to reload the images from the input file          */
/* many times, since each new noise component we calculate        */
/* destroys the original data set.                                */
/*                                                                */
/******************************************************************/

int ProcessImages()
{

	double sum;              
	double s,stvh,stv,sth,svh,st,sv,sh;
	int mincol,maxcol;

	mincol = first_col - 1;
	maxcol = last_col - 1;

	LoadImages();

	s=GetS(mincol,maxcol,MAXROW,last_frame-first_frame+1);

	sum=0.0;

	LoadImages();  
	st=GetNT(mincol,maxcol,MAXROW,last_frame-first_frame+1);
	sum+=st*st;

	LoadImages();  
	sh=GetNH(mincol,maxcol,MAXROW,last_frame-first_frame+1);
	sum+=sh*sh;

	LoadImages();  
	sv=GetNV(mincol,maxcol,MAXROW,last_frame-first_frame+1);
	sum+=sv*sv;

	LoadImages();  
	stv=GetNTV(mincol,maxcol,MAXROW,last_frame-first_frame+1);
	sum+=stv*stv;

	LoadImages();  
	sth=GetNTH(mincol,maxcol,MAXROW,last_frame-first_frame+1);
	sum+=sth*sth;

	LoadImages();  
	svh=GetNVH(mincol,maxcol,MAXROW,last_frame-first_frame+1);
	sum+=svh*svh;

	LoadImages();  
	stvh=GetNTVH(mincol,maxcol,MAXROW,last_frame-first_frame+1);
	sum+=stvh*stvh;

	printf("\nImage average S       : %7.2f (gl)\n",s);
	printf("Total system noise    : %7.3f (gl)\n\n",sqrt(sum));
	printf("Fixed/spatial noise (gl) | Temporal noise (gl) | Variation effect \n");
	printf("-------------------------|---------------------|---------------------------\n");
	printf("Nh    : %7.3f          | Nth   : %7.3f     | Column - scan effects\n",sh,sth);
	printf("Nv    : %7.3f          | Ntv   : %7.3f     | Row - non-uniformity & dcs\n",sv,stv);
	printf("Nvh   : %7.3f          | Ntvh  : %7.3f     | Pixel -  detector noise\n",svh,stvh);
	printf("                         | Nt    : %7.3f     | Frame \n",st);

	//printf("                                         Fixed/spatial noise         Temporal noise \n");
	//printf("Column variations [scan effects]         Nh    : %7.3f (gl)        Nth   : %7.3f (gl)\n",sh,sth);
	//printf("Row variations [non-uniformity and dcs]  Nv    : %7.3f (gl)        Ntv   : %7.3f (gl)\n",sv,stv);
	//printf("Pixel variations [detector noise]        Nvh   : %7.3f (gl)        Ntvh  : %7.3f (gl)\n",svh,stvh);
	//printf("Frame variations                                                     Nt    : %7.3f (gl)\n",st);


	//printf("Nvh/Ntvh: %f \n",svh/stvh);
	//printf("Ntv/Ntvh: %f \n",stv/stvh);
	//printf("Nv/Ntvh : %f \n",sv/stvh);
	//printf("Nth/Ntvh: %f \n",sth/stvh);
	//printf("Nh/Ntvh : %f \n",sh/stvh);
	//printf("Nt/Ntvh : %f \n\n",st/stvh);

	return(0);

}

/******************************************************************/
/*                                                                */
/* GetS                                                           */
/*                                                                */
/* s average over all the images                                  */
/*                                                                */
/******************************************************************/
double GetS(int mincol, int maxcol, int maxrow, int maxframe)
{
	double s;

	oprDH(mincol, maxcol,maxrow,maxframe); // results in one column for each frame
	oprDV(mincol, mincol,maxrow,maxframe); // results in one pixel for each frame
	oprDT(mincol, mincol,1,maxframe);      // results in one pixel

	s = Image[0][0][mincol];

	return s;
}

/******************************************************************/
/*                                                                */
/* GetNT                                                          */
/*                                                                */
/* average for each pixel over all frames                         */
/*                                                                */
/******************************************************************/
double GetNT(int mincol, int maxcol, int maxrow, int maxframe)
{
	int    frame;
	double sum, mean, stddev;
	double pixval;

	oprSDT(mincol, maxcol,maxrow,maxframe);  // results in one matrix for each frame
	oprDH(mincol, maxcol,maxrow,maxframe);   // results in one column for each frame
	oprDV(mincol, mincol,maxrow,maxframe);   // results in one pixel for each frame

	// calc the variance 

	for (sum=0.0,sum=0.0,frame=0;frame<maxframe;frame++)
		sum+=(Image[frame][0][mincol]);

	mean=sum/maxframe;  

	for (sum=0.0,frame=0;frame<maxframe;frame++)
	{
		pixval=(Image[frame][0][mincol]-mean);
		sum+=pixval*pixval;
	}

	stddev=sqrt(sum/(double)(maxframe-1));

	return stddev;
}

/******************************************************************/
/*                                                                */
/* GetNVH                                                         */
/*                                                                */
/* average for each frame over all pixels                         */
/*                                                                */
/******************************************************************/
double GetNVH(int mincol, int maxcol, int maxrow, int maxframe)
{
	int    row,col;
	double sum, mean, stddev;
	double pixval;
	double count;

	oprSDH(mincol, maxcol,maxrow,maxframe);  // results in one matrix for each frame
	oprSDV(mincol, maxcol,maxrow,maxframe);  // results in one matrix for each frame
	oprDT(mincol, maxcol,maxrow,maxframe);   // results in one matrix 

	// calc the variance 

	for (sum=0.0,row=0;row<maxrow;row++)
		for (col=mincol ;col<=maxcol;col++)
			sum+=(Image[0][row][col]);

	mean=(sum/maxrow)/(maxcol-mincol+1);  

	for (sum=0.0,row=0;row<maxrow;row++)
	{
		for (col=mincol;col<=maxcol;col++)
		{
			pixval=(Image[0][row][col]-mean);
			sum+=pixval*pixval;
		}
	}

	count=(double)maxrow*(double)(maxcol-mincol+1)-1;
	stddev=sqrt(sum/count);

	return stddev;
}

/******************************************************************/
/*                                                                */
/* GetNTV                                                         */
/*                                                                */
/* average for each row and frame over all columns                */
/*                                                                */
/******************************************************************/
double GetNTV(int mincol, int maxcol, int maxrow, int maxframe)
{
	int frame,row;
	double sum, mean, stddev;
	double pixval;
	double count;

	oprSDV(mincol,maxcol,maxrow,maxframe);  // results in one matrix for each frame
	oprSDT(mincol,maxcol,maxrow,maxframe);  // results in one matrix for each frame
	oprDH(mincol,maxcol,maxrow,maxframe);   // results in one column for each frame

	// calc the variance 

	for (sum=0.0,row=0;row<maxrow;row++)
		for (frame=0;frame<maxframe;frame++)
			sum+=(Image[frame][row][mincol]);

	mean=(sum/maxrow)/maxframe;  

	for (sum=0.0,row=0;row<maxrow;row++)
	{
		for (frame=0;frame<maxframe;frame++)
		{
			pixval=(Image[frame][row][mincol]-mean);
			sum+=pixval*pixval;
		}
	}

	count=(double)maxrow*(double)maxframe-1;
	stddev=sqrt(sum/count);

	return stddev;
}

/******************************************************************/
/*                                                                */
/* GetNTH                                                         */
/*                                                                */
/* average for each column and frame over all rows                */
/*                                                                */
/******************************************************************/
double GetNTH(int mincol, int maxcol, int maxrow, int maxframe)
{

	int frame,col;
	double sum, mean, stddev;
	double count;
	double pixval;

	oprSDH(mincol,maxcol,maxrow,maxframe);  // results in one matrix for each frame
	oprSDT(mincol,maxcol,maxrow,maxframe);  // results in one matrix for each frame
	oprDV(mincol,maxcol,maxrow,maxframe);   // results in one row for each frame

	// calc the variance 

	for (sum=0.0,col=mincol;col<=maxcol;col++)
		for (frame=0;frame<maxframe;frame++)
			sum+=(Image[frame][0][col]);

	mean=(sum/(maxcol-mincol+1))/maxframe;  

	for (sum=0.0,col=mincol;col<=maxcol;col++)
	{
		for (frame=0;frame<maxframe;frame++)
		{
			pixval=(Image[frame][0][col]-mean);
			sum+=pixval*pixval;
		}
	}

	count = (double)(maxcol-mincol+1)*(double)(maxframe)-1;
	stddev=sqrt(sum/count);

	return stddev;
}

/******************************************************************/
/*                                                                */
/* GetNV                                                          */
/*                                                                */
/* average for each column over all frames and rows               */
/*                                                                */
/******************************************************************/
double GetNV(int mincol, int maxcol, int maxrow, int maxframe)
{

	int row;
	double sum, mean, stddev;
	double pixval;
	double count;

	oprSDV(mincol, maxcol,maxrow,maxframe);  // results in one matrix for each frame
	oprDH(mincol, maxcol,maxrow,maxframe);   // results in one column for each frame
	oprDT(mincol, mincol,maxrow,maxframe);   // results in one column

	// calc the variance 

	for (sum=0.0,row=0;row<maxrow;row++)
		sum+=(Image[0][row][mincol]);

	mean=(sum/maxrow);  

	for (sum=0.0,row=0;row<maxrow;row++)
	{
		pixval=(Image[0][row][mincol]-mean);
		sum+=pixval*pixval;
	}

	count=(double)maxrow-1;
	stddev=sqrt(sum/count);

	return stddev;
}


/******************************************************************/
/*                                                                */
/* GetNH                                                          */
/*                                                                */
/* average for each row over all frames and cols                  */
/*                                                                */
/******************************************************************/
double GetNH(int mincol, int maxcol, int maxrow, int maxframe)
{

	int col;
	double count;
	double pixval;
	double sum, mean, stddev;

	oprSDH(mincol, maxcol,maxrow,maxframe);  // results in one matrix for each frame
	oprDV(mincol, maxcol,maxrow,maxframe);   // results in one row for each frame
	oprDT(mincol, maxcol,1,maxframe);        // results in one row 

	// calc the variance 

	for (sum=0.0,col=mincol;col<=maxcol;col++)
		sum+=(Image[0][0][col]);

	mean=(sum/(maxcol-mincol+1));  

	for (sum=0.0,col=mincol;col<=maxcol;col++)
	{
		pixval=(Image[0][0][col]-mean);
		sum+=pixval*pixval;
	}

	count=(double)(maxcol-mincol+1)-1;
	stddev=sqrt(sum/count);
	return stddev;
}

/******************************************************************/
/*                                                                */
/* GetNTVH                                                        */
/*                                                                */
/* noise for each row ,   frame & column                          */
/*                                                                */
/******************************************************************/
double GetNTVH(int mincol, int maxcol, int maxrow, int maxframe)
{

	int    frame,row,col;
	double sum, mean, stddev;
	double pixval;
	double count;

	oprSDH(mincol, maxcol,maxrow,maxframe); // results in one matrix for each frame
	oprSDV(mincol, maxcol,maxrow,maxframe); // results in one matrix for each frame
	oprSDT(mincol, maxcol,maxrow,maxframe); // results in one matrix for each frame

	// calc the variance 

	for (sum=0.0,row=0;row<maxrow;row++)
		for (col=mincol;col<=maxcol;col++)
			for (frame=0;frame<maxframe;frame++)
				sum+=(Image[frame][row][col]);

	mean=((sum/maxrow)/maxframe)/(maxcol-mincol+1);  

	for (sum=0.0,row=0;row<maxrow;row++)
	{
		for (col=mincol;col<=maxcol;col++)
		{
			for (frame=0;frame<maxframe;frame++)
			{
				pixval=(Image[frame][row][col]-mean);
				sum+=pixval*pixval;
			}
		}
	}

	count=(double)maxrow*(double)(maxcol-mincol+1)*(double)maxframe-1;
	stddev=sqrt(sum/count);

	return stddev;
}

/******************************************************************/
/*                                                                */
/* oprDV                                                          */
/*                                                                */
/* this operation takes place between 0 and  < the                */
/* upper limits specified as function parameters.                 */
/* operator DH averages over rows, for each pixel                 */
/* the result is left in the first (0th) row of the               */
/* Images vector,                                                 */
/*                                                                */
/******************************************************************/
int oprDV(int mincol, int maxcol, int maxrow, int maxframe)
{

	int row,col,frame;
	double sum;

	if (maxrow  > MAXROW) maxrow=MAXROW;
	if (maxcol  > MAXCOL) maxcol=MAXCOL;
	if (maxframe > MAX_FRAMES) maxframe = MAX_FRAMES;

	for (frame=0;frame<maxframe;frame++){
		for (col= mincol;col<=maxcol;col++){
			sum=0.0;
			for (row=0;row<maxrow;row++){
				sum+= Image[frame][row][col];
			}
			Image[frame][0][col]=(float)(sum/maxrow);
		}
	}

	return 0;
}

/******************************************************************/
/*                                                                */
/* oprDH                                                          */
/*                                                                */
/* this operation takes place between 0 and  < the                */
/* upper limits specified as function parameters.                 */
/* operator DV averages over columns, for each pixel              */
/* the result is left in the user specified "mincol" col of the   */
/* Images vector,                                                 */
/*                                                                */
/******************************************************************/
int oprDH(int mincol, int maxcol, int maxrow, int maxframe)
{

	int row,col,frame;
	double sum;

	if (maxrow  > MAXROW) maxrow=MAXROW;
	if (maxcol  > MAXCOL) maxcol=MAXCOL;
	if (maxframe > MAX_FRAMES) maxframe = MAX_FRAMES;

	for (frame=0;frame<maxframe;frame++){
		for (row=0;row<maxrow;row++){
			sum=0.0;
			for (col=mincol;col<=maxcol;col++){
				sum+= Image[frame][row][col] ;
			}
			Image[frame][row][mincol] = (float)(sum/(maxcol-mincol+1));
		}
	}
	return 0;
}

/******************************************************************/
/*                                                                */
/* oprDT                                                          */
/*                                                                */
/* this operation takes place between 0 and  < the                */
/* upper limits specified as function parameters.                 */
/* operator DT averages over frame, for each pixel                */
/* the result is left in the first (0th) image of the             */
/* Images vector,                                                 */
/*                                                                */
/******************************************************************/
int oprDT(int mincol, int maxcol, int maxrow, int maxframe)
{

	int row,col,frame;
	double sum;

	if (maxrow  > MAXROW) maxrow=MAXROW;
	if (maxcol  > MAXCOL) maxcol=MAXCOL;

	if (maxframe > MAX_FRAMES) maxframe = MAX_FRAMES;

	for (row=0;row<maxrow;row++)
	{
		for (col=mincol;col<=maxcol;col++)
		{
			sum = 0.0;
			for (frame=0;frame<maxframe;frame++)
			{
				sum += Image[frame][row][col];
			}
			Image[0][row][col]= (float)(sum/maxframe);
		}
	}

	return 0;
}

/******************************************************************/
/*                                                                */
/* oprSDT                                                         */
/*                                                                */
/* this operation takes place between 0 and  < the                */
/* upper limits specified as function parameters.                 */
/* operator SDT averages over frames, for each pixel              */
/* the result is  subtracted from all other images                */
/*                                                                */
/******************************************************************/
int oprSDT(int mincol, int maxcol, int maxrow, int maxframe)
{
	int    row,col,frame;
	double sum;

	if (maxrow  > MAXROW) maxrow=MAXROW;
	if (maxcol  > MAXCOL) maxcol=MAXCOL;
	if (maxframe > MAX_FRAMES) maxframe = MAX_FRAMES;

	for (row=0;row<maxrow;row++)
	{
		for (col=mincol;col<=maxcol;col++)
		{
			sum=0.0;
			for (frame=0;frame<maxframe;frame++)
			{
				sum+= Image[frame][row][col] ;
			}
			tmpImage[0][row][col]= (float)(sum/maxframe);
		}
	}

	/* now subtract the average */

	for (row=0;row<maxrow;row++)
		for (col=mincol;col<=maxcol;col++)
			for (frame=0;frame<maxframe;frame++)
				Image[frame][row][col]-= tmpImage[0][row][col];

	return 0;
}

/******************************************************************/
/*                                                                */
/* oprSDH                                                         */
/*                                                                */
/* this operation takes place between 0 and  < the                */
/* upper limits specified as function parameters.                 */
/* operator SDV averages over columns, for each pixel             */
/* the result is  subtracted from all other images                */
/*                                                                */
/******************************************************************/
int oprSDH(int mincol, int maxcol, int maxrow, int maxframe)
{
	int row,col,frame;
	double sum;

	if (maxrow  > MAXROW) maxrow=MAXROW;
	if (maxcol  > MAXCOL) maxcol=MAXCOL;
	if (maxframe > MAX_FRAMES) maxframe = MAX_FRAMES;

	for (frame=0;frame<maxframe;frame++)
	{
		for (row=0;row<maxrow;row++)
		{
			sum=0.0;
			for (col= mincol ;col<=maxcol;col++)
			{
				sum+=  (Image[frame][row][col]);
			}
			tmpImage[frame][row][0]= (float)(sum/(maxcol-mincol+1));
		}
	}

	/* now subtract the average */

	for (frame=0;frame<maxframe;frame++)
		for (row=0;row<maxrow;row++)
			for (col= mincol;col<=maxcol;col++)
				Image[frame][row][col]-=tmpImage[frame][row][0];

	return 0;
}

/******************************************************************/
/*                                                                */
/* oprSDV                                                         */
/*                                                                */
/* this operation takes place between 0 and  < the                */
/* upper limits specified as function parameters.                 */
/* operator SDH averages over rows, for each pixel                */
/* the result is  subtracted from all other images                */
/*                                                                */
/******************************************************************/
int oprSDV(int mincol, int maxcol, int maxrow, int maxframe)
{

	int row,col,frame;
	double sum;

	if (maxrow  > MAXROW) maxrow=MAXROW;
	if (maxcol  > MAXCOL) maxcol=MAXCOL;
	if (maxframe > MAX_FRAMES) maxframe = MAX_FRAMES;

	for (frame=0;frame<maxframe;frame++)
	{
		for (col=mincol;col<=maxcol;col++)
		{
			sum=0.0;
			for (row=0;row<maxrow;row++)
			{
				sum+=Image[frame][row][col];
			}
			tmpImage[frame][0][col] = (float)(sum/maxrow);

		}
	}

	/* now subtract the average */

	for (frame=0;frame<maxframe;frame++)
		for (col=mincol;col<=maxcol;col++)
			for (row=0;row<maxrow;row++)
				Image[frame][row][col]-=tmpImage[frame][0][col];

	return 0;
}

/******************************************************************/
/*                                                                */
/* matrix3d                                                       */
/*                                                                */
/******************************************************************/
float ***matrix3d(int nsl,int nsh,int nrl,int nrh,int ncl,int nch)
{
	int i,j;
	float ***m;

	m=(float ***) malloc((unsigned) (nsh-nsl+1)*sizeof(float**));
	if (!m)  {
		printf("memory allocation failure in matrix3d()\n");
		exit(-1);
	}
	m -= nsl;

	for(i=nsl;i<=nsh;i++) {  
		m[i]=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
		if (!m[i])  {
			printf("memory allocation failure in matrix3d()\n");
			exit(-1);
		}
		m[i] -= nrl;

		for(j=nrl;j<=nrh;j++) {
			m[i][j]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
			if (!m[i][j])  {
				printf("memory allocation failure in matrix3d()\n");
				exit(-1);
			}
			m[i][j] -= ncl;
		}

	}

	return m;
}

/******************************************************************/
/*                                                                */
/* matrix3dUL                                                     */
/*                                                                */
/******************************************************************/
unsigned long ***matrix3dUL(int nsl,int nsh,int nrl,int nrh,int ncl,int nch)
{
	int i,j;
	unsigned long ***m;

	m=(unsigned long ***) malloc((unsigned) (nsh-nsl+1)*sizeof(unsigned long**));
	if (!m)  {
		printf("memory allocation failure in matrix3dUL()\n");
		exit(-1);
	}
	m -= nsl;

	for(i=nsl;i<=nsh;i++) {  
		m[i]=(unsigned long **) malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned long*));
		if (!m[i])  {
			printf("memory allocation failure in matrix3dUL()\n");
			exit(-1);
		}
		m[i] -= nrl;

		for(j=nrl;j<=nrh;j++) {
			m[i][j]=(unsigned long *) malloc((unsigned) (nch-ncl+1)*sizeof(unsigned long));
			if (!m[i][j])  {
				printf("memory allocation failure in matrix3dUL()\n");
				exit(-1);
			}
			m[i][j] -= ncl;
		}

	}

	return m;
}
