#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include<unistd.h>
#include<pthread.h>
#include<string.h>
#include<sys/stat.h>
#include"sgy_reader.cpp"
#define angle_r 80

char filename[]=R"(data/model_data.sgy)";
int trace_win=28;//can be determined by find_window_size
int span=trace_win+1;
int time_win=56;//can be determined by find_window_size

#ifdef _WIN32
#include<windows.h>
LARGE_INTEGER frequency;
LARGE_INTEGER start_time;
void tic() 
{
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&start_time);
}
double toc() 
{
    LARGE_INTEGER end_time;
    QueryPerformanceCounter(&end_time);
    double elapsed_time = (double)(end_time.QuadPart - start_time.QuadPart) / frequency.QuadPart;
    return elapsed_time;
}

int TH4=atoi(getenv("NUMBER_OF_PROCESSORS"));
#elif __linux__
int TH4=sysconf(_SC_NPROCESSORS_CONF);
struct timespec start_time;
void tic() 
{
    clock_gettime(CLOCK_MONOTONIC, &start_time);
}
double toc() 
{
    struct timespec end_time;
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    double elapsed_time = (end_time.tv_sec - start_time.tv_sec) +
                          (end_time.tv_nsec - start_time.tv_nsec) / 1e9;
    return elapsed_time;
}
#endif
int finish;
void padarray2(double** input,double** output,int padx,int pady,int padder,int nrow,int ncol)
{
	//0 stands for pading with itp
	//1 stands for padding with mean value
	//2 stands for 'replicate' in matlab
	//3 stands for 'symmetric' in matlab
	/*//call:
	float** A=(float**)malloc(4*sizeof(float*));
	for(int i=0;i<4;i++)A[i]=(float*)malloc(5*sizeof(float));
	for(int i=0;i<4;i++)for(int j=0;j<5;j++)A[i][j]=3*i+2;
	float** B=(float**)malloc((4+2)*sizeof(float*));
	for(int i=0;i<4+2;i++)B[i]=(float*)malloc((5+3)*sizeof(float));*/
	float itp=0;
	int i,j;

	for(i=0; i<nrow+2*padx; i++)for(j=0; j<ncol+2*pady; j++)output[i][j]=0;
	if((padx>nrow||pady>ncol)&&padder==3)printf("warning: pad size exceeds the input size\n");
	switch(padder)
	{
	case 0:
		for(i=0; i<2*padx+nrow; i++)for(j=0; j<2*pady+ncol; j++)output[i][j]=itp;
		for(i=padx; i<nrow+padx; i++)for(j=pady; j<ncol+pady; j++)output[i][j]=input[i-padx][j-pady];
		break;
	case 1:
		for(i=0; i<nrow; i++)for(j=0; j<ncol; j++)itp+=input[i][j];
		itp/=(nrow*ncol);
		for(i=0; i<2*padx+nrow; i++)for(j=0; j<2*pady+ncol; j++)output[i][j]=itp;
		for(i=padx; i<nrow+padx; i++)for(j=pady; j<ncol+pady; j++)output[i][j]=input[i-padx][j-pady];
		break;
	case 2:
		for(i=padx; i<nrow+padx; i++)for(j=pady; j<ncol+pady; j++)output[i][j]=input[i-padx][j-pady];
		for(i=pady; i<pady+ncol; i++)
		{
			for(j=0; j<padx; j++)output[j][i]=output[padx][i];
			for(j=padx+nrow; j<2*padx+nrow; j++)output[j][i]=output[padx+nrow-1][i];
		}
		for(i=0; i<2*padx+nrow; i++)
		{
			for(j=0; j<pady; j++)output[i][j]=output[i][pady];
			for(j=pady+ncol; j<2*pady+ncol; j++)output[i][j]=output[i][pady+ncol-1];
		}
		break;
	case 3:
		for(i=padx; i<nrow+padx; i++)for(j=pady; j<ncol+pady; j++)output[i][j]=input[i-padx][j-pady];
		for(i=padx; i<nrow+padx; i++)
		{
			for(j=0; j<pady; j++)output[i][j]=input[i-padx][pady-j-1];
			for(j=pady+ncol; j<2*pady+ncol; j++)output[i][j]=input[i-padx][2*ncol+pady-j-1];
		}
		for(j=pady; j<ncol+pady; j++)
		{
			for(i=0; i<padx; i++)output[i][j]=input[padx-i-1][j-pady];
			for(i=padx+nrow; i<2*padx+nrow; i++)output[i][j]=input[2*nrow+padx-i-1][j-pady];
		}
		for(i=0; i<padx; i++)
		{
			for(j=0; j<pady; j++)output[i][j]=input[padx-i-1][pady-j-1];
			for(j=ncol+pady; j<ncol+2*pady; j++)output[i][j]=input[padx-i-1][2*ncol+pady-j-1];
		}
		for(i=nrow+padx; i<nrow+2*padx; i++)
		{
			for(j=ncol+pady; j<ncol+2*pady; j++)output[i][j]=input[2*nrow+padx-i-1][2*ncol+pady-j-1];
			for(j=0; j<pady; j++)output[i][j]=input[2*nrow+padx-i-1][pady-j-1];
		}
	}
}
void padarray2f(float** input,float** output,int padx,int pady,int padder,int nrow,int ncol)
{
	//0 stands for pading with itp
	//1 stands for padding with mean value
	//2 stands for 'replicate' in matlab
	//3 stands for 'symmetric' in matlab
	/*//call:
	float** A=(float**)malloc(4*sizeof(float*));
	for(int i=0;i<4;i++)A[i]=(float*)malloc(5*sizeof(float));
	for(int i=0;i<4;i++)for(int j=0;j<5;j++)A[i][j]=3*i+2;
	float** B=(float**)malloc((4+2)*sizeof(float*));
	for(int i=0;i<4+2;i++)B[i]=(float*)malloc((5+3)*sizeof(float));*/

	float itp=0;
	int i,j;

	for(i=0; i<nrow+2*padx; i++)for(j=0; j<ncol+2*pady; j++)output[i][j]=0;
	if((padx>nrow||pady>ncol)&&padder==3)printf("warning: pad size exceeds the input size\n");
	switch(padder)
	{
	case 0:
		for(i=0; i<2*padx+nrow; i++)for(j=0; j<2*pady+ncol; j++)output[i][j]=itp;
		for(i=padx; i<nrow+padx; i++)for(j=pady; j<ncol+pady; j++)output[i][j]=input[i-padx][j-pady];
		break;
	case 1:
		for(i=0; i<nrow; i++)for(j=0; j<ncol; j++)itp+=input[i][j];
		itp/=(nrow*ncol);
		for(i=0; i<2*padx+nrow; i++)for(j=0; j<2*pady+ncol; j++)output[i][j]=itp;
		for(i=padx; i<nrow+padx; i++)for(j=pady; j<ncol+pady; j++)output[i][j]=input[i-padx][j-pady];
		break;
	case 2:
		for(i=padx; i<nrow+padx; i++)for(j=pady; j<ncol+pady; j++)output[i][j]=input[i-padx][j-pady];
		for(i=pady; i<pady+ncol; i++)
		{
			for(j=0; j<padx; j++)output[j][i]=output[padx][i];
			for(j=padx+nrow; j<2*padx+nrow; j++)output[j][i]=output[padx+nrow-1][i];
		}
		for(i=0; i<2*padx+nrow; i++)
		{
			for(j=0; j<pady; j++)output[i][j]=output[i][pady];
			for(j=pady+ncol; j<2*pady+ncol; j++)output[i][j]=output[i][pady+ncol-1];
		}
		break;
	case 3:
		for(i=padx; i<nrow+padx; i++)for(j=pady; j<ncol+pady; j++)output[i][j]=input[i-padx][j-pady];
		for(i=padx; i<nrow+padx; i++)
		{
			for(j=0; j<pady; j++)output[i][j]=input[i-padx][pady-j-1];
			for(j=pady+ncol; j<2*pady+ncol; j++)output[i][j]=input[i-padx][2*ncol+pady-j-1];
		}
		for(j=pady; j<ncol+pady; j++)
		{
			for(i=0; i<padx; i++)output[i][j]=input[padx-i-1][j-pady];
			for(i=padx+nrow; i<2*padx+nrow; i++)output[i][j]=input[2*nrow+padx-i-1][j-pady];
		}
		for(i=0; i<padx; i++)
		{
			for(j=0; j<pady; j++)output[i][j]=input[padx-i-1][pady-j-1];
			for(j=ncol+pady; j<ncol+2*pady; j++)output[i][j]=input[padx-i-1][2*ncol+pady-j-1];
		}
		for(i=nrow+padx; i<nrow+2*padx; i++)
		{
			for(j=ncol+pady; j<ncol+2*pady; j++)output[i][j]=input[2*nrow+padx-i-1][2*ncol+pady-j-1];
			for(j=0; j<pady; j++)output[i][j]=input[2*nrow+padx-i-1][pady-j-1];
		}
	}
}
void normalization(float **data,long long nrow,int ncol)
{
	int i,j;
	float data_min=FLT_MAX,data_max=-FLT_MAX;
	for(i=0;i<nrow;i++)
	{
		for(j=0;j<ncol;j++)
		{
			if(data_max<data[i][j])data_max=data[i][j];
			if(data_min>data[i][j])data_min=data[i][j];
		}
	}
	for(i=0;i<nrow;i++)for(j=0;j<ncol;j++)data[i][j]=2*(data[i][j]-data_min)/(data_max-data_min)-1;
}
void gaussfilt2(double** input,double** output,double sigma,int nrow,int ncol)
{	//paras↓
	int padder=2;//pad mode
	float L;
	L=2*ceil(2*sigma)+1;//kernal size, odd number required
	//paras↑
	/*//call:
	float** A=(float**)malloc(sizeof(float*)*3);
	for(i=0;i<3;i++)A[i]=(float*)malloc(sizeof(float)*4);
	float** B=(float**)malloc(sizeof(float*)*Nx);
	for(i=0;i<3;i++)B[i]=(float*)malloc(sizeof(float)*4);
	for(i=0;i<3;i++)for(j=0;j<4;j++)A[i][j]=2*i+j;*/
	int i,j,kx,ky;
	double sum=0;
	int l=(L-1)/2;

	double** gaussfilter=(double**)malloc(sizeof(double*)*L);
	for(i=0; i<L; i++)gaussfilter[i]=(double*)malloc(sizeof(double)*L);
	for(i=0; i<L; i++)
	{
		for(j=0; j<L; j++)
		{
			gaussfilter[i][j]=exp(-((i-l)*(i-l)+(j-l)*(j-l))/(2.0*sigma*sigma));
			sum+=gaussfilter[i][j];
		}
	}
	for(i=0; i<L; i++)for(j=0; j<L; j++)gaussfilter[i][j]=gaussfilter[i][j]/sum;//kernal normalization

	double** pad=(double**)calloc(2*l+nrow,sizeof(double*));
	for(i=0; i<nrow+2*l; i++)pad[i]=(double*)calloc(ncol+2*l,sizeof(double));
	padarray2(input,pad,l,l,padder,nrow,ncol);
	for(i=0; i<nrow; i++)for(j=0; j<ncol; j++)output[i][j]=0;
	for(i=0; i<nrow; i++)for(j=0; j<ncol; j++)for(kx=0; kx<L; kx++)for(ky=0; ky<L; ky++)output[i][j]+=gaussfilter[kx][ky]*pad[2*l+i-kx][2*l+j-ky];
	free(pad);
	free(gaussfilter);
}
float corhrow(float** input,int nrow,int ncol)
{
	int i,j,cnt_zero=0,main_row;
	float C=1.0,dot,sqt1,sqt2=0;
	main_row=(int)(nrow/2)+nrow%2-1;
	for(j=0; j<ncol; j++)sqt2+=input[main_row][j]*input[main_row][j];
	if(sqt2==0)return 0;
	for(i=0; i<nrow; i++)
	{
		dot=0;
		sqt1=0;
		for(j=0; j<ncol; j++)
		{
			dot+=input[i][j]*input[main_row][j];
			sqt1+=input[i][j]*input[i][j];
		}
		if(sqt1!=0)C+=dot/(sqrt(sqt1)*sqrt(sqt2));// 
		else
		{
			cnt_zero++;
			if(cnt_zero>round(main_row/2))return 0;
		}
	}
	if(C<trace_win/5)return 0;
	return C;
}

float median(float* arr,int l)
{
	int i,j,p;
	float t,mid,* x;
	x=(float*)malloc(l*sizeof(float));
	for(i=0; i<l; i++)x[i]=*(arr++);
	for(i=0; i<l-1; i++)for(j=0; j<l-1; j++)if(x[j]>x[j+1])
	{
		t=x[j+1];
		x[j+1]=x[j];
		x[j]=t;
	}
	if(l%2==0)
	{
		p=l/2;
		mid=(x[p]+x[p-1])/2.0;
	}
	else mid=x[(l-1)/2];
	free(x);
	return mid;
}
float mean(float* arr,int l)
{
	int i;
	float sum=0;
	for(i=0;i<l;i++)sum+=arr[i];
	return sum/l;
}
struct paras
{
	float* data;
	int tracesL;
	int samplesL;
	int thread;
	int pad_x;
	int pad_t;
	int* outcome;
};
struct manynames
{
	char basename[50];
	char path1[150];
	char path2[150];
	char file[20];
	char suf[20];
};
void findnames(char* filename,manynames* outfile)
{
	char* pbasename,* pfile,* psuf=NULL;
	strcpy(outfile->path1,"");
	strcpy(outfile->path2,"");
	strcpy(outfile->suf,"");
	struct stat st;
	if(stat(filename,&st)!=0){ printf("file not found\n");return; }
	printf("\n");
	char* sepa_sym;
	if(sepa_sym=strrchr(filename,'/'))
	{
		*sepa_sym='\0';
		pbasename=sepa_sym+1;
		strcpy(outfile->path1,filename);
		strcpy(outfile->path2,filename);
		strcat(outfile->path2,"/");
	}
	else if(sepa_sym=strrchr(filename,'\\'))
	{
		*sepa_sym='\0';
		pbasename=sepa_sym+1;
		strcpy(outfile->path1,filename);
		strcpy(outfile->path2,filename);
		strcat(outfile->path2,"\\");
		if(outfile->path1[strlen(outfile->path1)-1]=='\\')outfile->path1[strlen(outfile->path1)-1]='\0';
	}
	else { pbasename=filename; }
	strcpy(outfile->basename,pbasename);
	char tmpbase[50];
	strcpy(tmpbase,outfile->basename);
	char* dot;
	pfile=tmpbase;
	if(dot=strrchr(tmpbase,'.')){ *dot='\0';psuf=dot+1;strcpy(outfile->suf,psuf); }
	strcpy(outfile->file,pfile);
}
void* calc_angles_pt(void* args)
{
	float* AL;
	int* theta_win_c1;
	int thread,i,j,k,kk,tracesL,samplesL,pad_x,pad_t,traces,samples,temp_theta,angle,start;
	paras* par=(paras*)args;
	AL=par->data;
	thread=par->thread;
	tracesL=par->tracesL;
	samplesL=par->samplesL;
	pad_x=par->pad_x;
	pad_t=par->pad_t;
	theta_win_c1=par->outcome;

	traces=tracesL-2*pad_x;
	samples=samplesL-2*pad_t;

	float corr,corr1;
	float** corh_mat=(float**)calloc(trace_win,sizeof(float*)); for(int ii=0; ii<trace_win; ii++)corh_mat[ii]=(float*)calloc(time_win,sizeof(float));
	for(i=pad_x+1+(int)ceil(traces*thread/TH4); i<pad_x+1+traces*(thread+1)/TH4; i++)
	{
		for(j=pad_t+1; j<pad_t+1+samples; j++)
		{
			corr=0;
			corr1=0;
			temp_theta=0;
			for(angle=-angle_r; angle<=angle_r; angle++)
			{
				for(k=i-trace_win/2; k<i+trace_win/2; k++)
				{
					start=roundf((j-time_win/2+(k-i)*tanf(angle*3.1415926/180)));
					for(kk=0;kk<time_win; kk++)corh_mat[k-(i-trace_win/2)][kk]=AL[(k-1)*samplesL+start-1+kk];
				}
				corr1=corhrow(corh_mat,trace_win,time_win);
				//if(i==pad_x+__&&j==pad_t+__)printf("%d\t%f\t%d\n",angle,corr,corr1 > corr);
				if(corr1>corr)
				{
					temp_theta=angle;
					corr=corr1;
				}
			}
			theta_win_c1[(i-1)*samplesL+j-1]=temp_theta;
		}
		finish++;
		if(finish%10==0)printf("%d/%d\n",finish,traces);
	}
	free(corh_mat);
	return args;
}
void interpolation_near(double** mtx,int nrow,int ncol,int idx)
{
	int l0=20,L0=50,dl=5,dL=5;
	int i,j,ii,jj,l,L,cnt;
	double mean,tmp;
	if(idx==1)
		for(i=0;i<nrow;i++)
		{
			for(j=0;j<ncol;j++)
			{
				l=l0;L=L0;mean=0;cnt=0;//Initial value of l and L
				while(mtx[i][j]<=0)
				{
					for(ii=0;ii<2*l+1;ii++)
					{
						for(jj=0;jj<2*L+1;jj++)
						{
							if(0<=i+ii-l&&i+ii-l<nrow&&0<=j+jj-L&&j+jj-L<ncol)
							{
								if((tmp=mtx[i+ii-l][j+jj-L])>0)
								{
									mean+=tmp;
									cnt++;
								}
							}
						}
					}
					if(cnt==0)
					{
						if(2*l<nrow)l+=dl;
						else L+=dL;
						continue;
					}
					mtx[i][j]=mean/cnt;
				}
				
			}
		}
	else if(idx==0)
		for(i=0;i<nrow;i++)
		{
			for(j=0;j<ncol;j++)
			{
				l=l0;L=L0;mean=0,cnt=0;
				while(mtx[i][j]>=0)
				{
					for(ii=0;ii<2*l+1;ii++)//遍历周围元素
					{
						for(jj=0;jj<2*L+1;jj++)
						{
							if(0<=i+ii-l&&i+ii-l<nrow&&0<=j+jj-L&&j+jj-L<ncol)//保证索引值为自然数
							{
								if((tmp=mtx[i+ii-l][j+jj-L])<0)
								{
									mean+=tmp;
									cnt++;
								}
							}
						}
					}
					if(cnt==0)
					{
						if(2*l<nrow)l+=dl;
						else L+=dL;
						continue;
					}
					mtx[i][j]=mean/cnt;
				}
			}
		}
	else
		printf("unknow instructions of interpolation\n");
}
void interpolation_easy(double** mtx,int nrow,int ncol,int idx)
{
	int i,j,cnt=0;
	double mean=0;
	if(idx==1)
	{
		for(i=0; i<nrow; i++)for(j=0; j<ncol; j++)if(mtx[i][j]>0) { mean+=mtx[i][j]; cnt++; }
		mean/=cnt;
		for(i=0; i<nrow; i++)for(j=0; j<ncol; j++)if(mtx[i][j]<=0)mtx[i][j]=mean;
	}
	else if(idx==0)
	{
		for(i=0; i<nrow; i++)for(j=0; j<ncol; j++)if(mtx[i][j]<0) { mean+=mtx[i][j]; cnt++; }
		mean/=cnt;
		for(i=0; i<nrow; i++)for(j=0; j<ncol; j++)if(mtx[i][j]>=0)mtx[i][j]=mean;
	}
}
int call_times=0;
void wavefield_sep(float** seis,float** U,float** D,int traces,int samples,char* anglefile)
{
	int i,j,k;
	float pi=3.1415926;
	int pad_x,pad_t,tracesL,samplesL;

	pad_x=ceil(trace_win/2)>ceil((span-1)/2)?ceil(trace_win/2):ceil((span-1)/2);
	pad_t=ceil(time_win/2+tan(angle_r*pi/180)*trace_win/2);
	printf("pad_x=%d,pad_t=%d\n",pad_x,pad_t);
	tracesL=2*pad_x+traces;
	samplesL=2*pad_t+samples;

	float** AL=(float**)malloc(tracesL*sizeof(float*));for(i=0; i<tracesL; i++)AL[i]=(float*)malloc(samplesL*sizeof(float));
	padarray2f(seis,AL,pad_x,pad_t,3,traces,samples);

	float* A1=(float*)malloc(4*tracesL*samplesL);
	for(i=0; i<tracesL; i++)for(j=0; j<samplesL; j++)A1[i*samplesL+j]=AL[i][j];

	tic();
	int* theta_win_c1=(int*)calloc(tracesL*samplesL,sizeof(int));
	pthread_t th[TH4];
	for(i=0; i<TH4; i++)
	{
		paras* para=new paras;
		para->data=A1;
		para->thread=i;
		para->tracesL=tracesL;
		para->samplesL=samplesL;
		para->pad_x=pad_x;
		para->pad_t=pad_t;
		para->outcome=theta_win_c1;
		pthread_create(&th[i],NULL,calc_angles_pt,(void*)para);
	}
	for(i=0; i<TH4; i++)pthread_join(th[i],NULL);
	printf("consumed time=%lf\n",toc());

	int** theta_win_c1_temp=(int**)malloc(traces*sizeof(int*)); for(i=0; i<traces; i++)theta_win_c1_temp[i]=(int*)malloc(samples*sizeof(int));
	for(i=0; i<traces; i++)for(j=0; j<samples; j++)theta_win_c1_temp[i][j]=(float)theta_win_c1[(i+pad_x)*samplesL+j+pad_t];
	free(theta_win_c1);/**/

	FILE* fp; fp=fopen(anglefile,"w");
	for(i=0; i<samples; i++)
	{
		for(j=0; j<traces; j++)
		{
			fprintf(fp,"%d,",theta_win_c1_temp[j][i]);
		}fprintf(fp,"\n");
	}fclose(fp);

	double** theta_win_c1_pos=(double**)malloc(traces*sizeof(double*)); for(i=0; i<traces; i++)theta_win_c1_pos[i]=(double*)malloc(samples*sizeof(double));
	double** pos=(double**)malloc(traces*sizeof(float*)); for(i=0; i<traces; i++)pos[i]=(double*)malloc(samples*sizeof(double));
	double** posL=(double**)malloc(tracesL*sizeof(float*)); for(i=0; i<tracesL; i++)posL[i]=(double*)calloc(samplesL,sizeof(double));
	for(i=0; i<traces; i++)for(j=0; j<samples; j++)
	{
		theta_win_c1_pos[i][j]=(double)theta_win_c1_temp[i][j];
		//pos[i][j] = (double)theta_win_c1_temp[i][j];
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//方法1，直接取平均值
	//interpolation_easy(theta_win_c1_pos,traces,samples,1);
	//方法2，取附近的点
	interpolation_near(theta_win_c1_pos,traces,samples,1);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	gaussfilt2(theta_win_c1_pos,pos,10,traces,samples);
	free(theta_win_c1_pos);
	////
	double** theta_win_c1_neg=(double**)malloc(traces*sizeof(double*)); for(i=0; i<traces; i++)theta_win_c1_neg[i]=(double*)malloc(samples*sizeof(double));
	double** neg=(double**)malloc(traces*sizeof(double*)); for(i=0; i<traces; i++)neg[i]=(double*)malloc(samples*sizeof(double));
	double** negL=(double**)malloc(tracesL*sizeof(double*)); for(i=0; i<tracesL; i++)negL[i]=(double*)calloc(samplesL,sizeof(double));
	for(i=0; i<traces; i++)for(j=0; j<samples; j++)
	{
		theta_win_c1_neg[i][j]=(double)theta_win_c1_temp[i][j];
		//neg[i][j] = (double)theta_win_c1_temp[i][j];
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//方法1，直接取平均值
	//interpolation_easy(theta_win_c1_neg,traces,samples,0);
	//方法2，取附近的点
	interpolation_near(theta_win_c1_neg,traces,samples,0);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	gaussfilt2(theta_win_c1_neg,neg,10,traces,samples);
	free(theta_win_c1_neg);
	free(theta_win_c1_temp);
	for(i=0; i<traces; i++)for(j=0; j<samples; j++)
	{
		posL[i+pad_x][j+pad_t]=pos[i][j];
		negL[i+pad_x][j+pad_t]=neg[i][j];
	}

	free(pos);
	free(neg);

	int x_cor,t_cor;
	float* med_arr;
	float** UL=(float**)malloc(tracesL*sizeof(float*)); for(i=0; i<tracesL; i++)UL[i]=(float*)calloc(samplesL,sizeof(float));
	float** DL=(float**)malloc(tracesL*sizeof(float*)); for(i=0; i<tracesL; i++)DL[i]=(float*)calloc(samplesL,sizeof(float));
	med_arr=(float*)calloc(span,sizeof(float));
	for(i=pad_x; i<pad_x+traces; i++)
	{
		for(j=pad_t; j<pad_t+samples; j++)
		{
			double k_pos=tan(posL[i][j]*pi/180);

			for(k=1; k<span+1; k++)
			{
				x_cor=round(i+1+k-(1+span)/2);
				t_cor=round(j+1+k_pos*(k-(1+span)/2));
				if(t_cor>0&&t_cor<samplesL&&x_cor>0&&x_cor<tracesL)med_arr[k-1]=AL[x_cor-1][t_cor-1];
			}
			DL[i][j]=median(med_arr,span);

			double k_neg=tan(negL[i][j]*pi/180);
			for(k=1; k<span+1; k++)
			{
				x_cor=round(i+1+k-(1+span)/2);
				t_cor=round(j+1+k_neg*(k-(1+span)/2));
				if(t_cor>0&&t_cor<samplesL&&x_cor>0&&x_cor<tracesL)med_arr[k-1]=AL[x_cor-1][t_cor-1];
			}
			UL[i][j]=median(med_arr,span);
		}
	}
	free(med_arr);
	free(AL);
	free(posL);
	free(negL);
	for(i=0; i<traces; i++)for(j=0; j<samples; j++)
	{
		U[i][j]=UL[i+pad_x][j+pad_t];
		D[i][j]=DL[i+pad_x][j+pad_t];
	}
	free(UL);
	free(DL);
}
int main()
{
	printf("%d threads are used\n",TH4);
	long long traces,samples,ts[2];
	traces_samples(filename,ts);
	traces=ts[0];
	samples=ts[1];
	int i,j;
	float** A=(float**)malloc(traces*sizeof(float*)); for(i=0; i<traces; i++)A[i]=(float*)malloc(samples*sizeof(float));
	unsigned char* hdr=(unsigned char*)malloc((3600+240*traces)*sizeof(unsigned char));
	readsgy(filename,A,hdr,traces,samples);
	normalization(A,traces,samples);
	//A[i][j]=amplitude(ith trace, jth sample), in traces, samples
	/////////////////////////////////////////////////////////////////////////////////////////////

	struct stat st;
	struct manynames thisfile;
	findnames(filename,&thisfile);
	char outfile01[100],outfile02[100],anglefile[100];
#ifdef _WIN32
	if(stat("ups",&st)!=0)mkdir("ups");
	if(stat("downs",&st)!=0)mkdir("downs");
	if(stat("angles",&st)!=0)mkdir("angles");
#elif __linux__
	if(stat("ups",&st)!=0)mkdir("ups",0755);
	if(stat("downs",&st)!=0)mkdir("downs",0755);
	if(stat("angles",&st)!=0)mkdir("angles",0755);
#endif

	float** U=(float**)malloc(traces*sizeof(float*)); for(i=0; i<traces; i++)U[i]=(float*)malloc(samples*sizeof(float));
	float** D=(float**)malloc(traces*sizeof(float*)); for(i=0; i<traces; i++)D[i]=(float*)malloc(samples*sizeof(float));
	float** RU=(float**)malloc(traces*sizeof(float*)); for(i=0; i<traces; i++)RU[i]=(float*)malloc(samples*sizeof(float));
	float** RD=(float**)malloc(traces*sizeof(float*)); for(i=0; i<traces; i++)RD[i]=(float*)malloc(samples*sizeof(float));
	float** R0=(float**)malloc(traces*sizeof(float*)); for(i=0; i<traces; i++)R0[i]=(float*)malloc(samples*sizeof(float));

	strcpy(anglefile,"./angles/angle_raw.csv");
	wavefield_sep(A,U,D,traces,samples,anglefile);
	finish=0;

	sprintf(outfile01,"./ups/%s_up0.sgy",thisfile.file);
	sprintf(outfile02,"./downs/%s_down0.sgy",thisfile.file);
	revisegy(U,outfile01,hdr,traces,samples);
	revisegy(D,outfile02,hdr,traces,samples);

	for(i=0; i<traces; i++)for(j=0; j<samples; j++)
	{
		RD[i][j]=A[i][j]-U[i][j];
		RU[i][j]=A[i][j]-D[i][j];
	}


	for(int iter=0; iter<2; iter++)
	{
		printf("iter=%d\n",iter+1);
		char outfile1[100],outfile2[100];
		sprintf(anglefile,"angles/angle_down%d.csv",iter+1);
		wavefield_sep(RD,R0,D,traces,samples,anglefile); finish=0;
		for(i=0; i<traces; i++)for(j=0; j<samples; j++)RD[i][j]=RD[i][j]-R0[i][j];
		sprintf(anglefile,"angles/angle_up%d.csv",iter+1);
		wavefield_sep(RU,U,R0,traces,samples,anglefile); finish=0;
		for(i=0; i<traces; i++)for(j=0; j<samples; j++)RU[i][j]=RU[i][j]-R0[i][j];

		sprintf(outfile1,"ups/%s_up%d.sgy",thisfile.file,iter+1);
		sprintf(outfile2,"downs/%s_down%d.sgy",thisfile.file,iter+1);
		revisegy(U,outfile1,hdr,traces,samples);
		revisegy(D,outfile2,hdr,traces,samples);
	}

	free(A);
	free(RU);
	free(RD);
	free(R0);
	free(U);
	free(D);
	free(hdr);
	return 0;
}

