#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
void float2ibm(float** x,int nrow,int ncol)
{
	int i,j;
	char sgn;
	unsigned char E,* p;
	for(i=0; i<nrow; i++)
	{
		for(j=0; j<ncol; j++)
		{
			if((x[i][j])!=0)
			{
				p=(unsigned char*)&x[i][j];
				sgn=1-2*(p[0]>>7);
				E=p[0]<<1;
				E>>=1;
				x[i][j]=sgn*(*(p+1)*65536+*(p+2)*256+*(p+3))*pow(16,(int)E-70);//
			}
		}
	}
}
void ibm2float(float** x,int nrow,int ncol)
{
	int i,j,k;
	unsigned char* p,* pfn,E;
	unsigned fn;
	float ft;
	char sgn;
	for(i=0; i<nrow; i++)
	{
		for(j=0; j<ncol; j++)
		{
			if((ft=x[i][j])!=0)
			{
				if(!(ft<=DBL_MAX&&ft>=-DBL_MAX))break;
				E=64;
				p=(unsigned char*)&x[i][j];
				pfn=(unsigned char*)&fn;
				for(k=0; k<4; k++)pfn[k]=p[k];
				sgn=(pfn[3]>>7);
				if(fabs(ft)>=1)
				{
					while((int)ft!=0)
					{
						E++;
						ft/=16;
					}
				}
				else
				{
					E++;
					while((int)ft*16==0)
					{
						E--;
						ft*=16;
					}
				}
				fn=((sgn<<31)|(E<<24))|(unsigned)((1-2*sgn)*x[i][j]*pow(16,70-(int)E));
				for(k=0; k<4; k++)p[k]=pfn[3-k];
			}
		}
	}
}
long long getFileSize(char* filePath)
{//https://blog.csdn.net/codears/article/details/111405309
	FILE* f;
	if(NULL==(f=fopen(filePath,"rb"))) { printf("fopen error\n"); return -1; }
	if(0!=fseek(f,0,SEEK_END)) { printf("getFileSize fseek error\n"); return -1; }
	long long fileSize=ftell(f);
	if(fileSize<0) { printf("ftell error\n"); }
	fclose(f);
	return fileSize;
}
void traces_samples(char* filename,long long* x0y1)
{
	FILE* fsgy;
	if((fsgy=fopen(filename,"rb"))==NULL) { printf("file not found\n");exit(-1); }
	int i;
	long long l;
	fseek(fsgy,0,2);
	l=getFileSize(filename);
	// printf("file size=%fMb\n",(float)l/(1024*1024));
	unsigned char fmt[2],sample[2];
	fseek(fsgy,3224,0);
	for(i=0; i<2; i++)fmt[i]=fgetc(fsgy);
	if((fmt[0]+fmt[1]!=1)&&(fmt[0]+fmt[1]!=5))
		printf("warning: format(byte3224-3225,start from 0) incorrect(%d,%d), check files first, traces and samples needed to define\n",fmt[0],fmt[1]);
	fseek(fsgy,3220,0);
	for(i=0; i<2; i++)sample[i]=fgetc(fsgy);
	if(fmt[0]==0)
	{
		x0y1[1]=sample[0]*256+sample[1];
		// printf("little endian\n");
	}
	else if(fmt[1]==0)
	{
		x0y1[1]=sample[1]*256+sample[0];
		// printf("big endian\n");
	}
	x0y1[0]=(l-3600)/(240+4*x0y1[1]);
}
void readsgy(char* filename,float** data,unsigned char* header,int nrow,int ncol)
{//header has 3600bytes volume header followed by 240*traces trace header 
	FILE* fsgy;
	int i,j;
	if((fsgy=fopen(filename,"rb"))==NULL) { printf("file not found\n"); exit(-1); }
	// printf("trace=%d\nsamples=%d\n",nrow,ncol);
	unsigned char fmt[2];
	fseek(fsgy,3224,0);
	for(i=0; i<2; i++)fmt[i]=fgetc(fsgy);
	size_t check;
	if(fmt[0]+fmt[1]==1)//sun
	{
		// printf("sun\n");
		fseek(fsgy,0,0);
		check=fread(header,3600,1,fsgy);
		for(i=0; i<nrow; i++)
		{
			check=fread(header+3600+i*240,1L,240L,fsgy);
			check=fread(data[i],4L,ncol,fsgy);
		}
		if(fmt[1]==0)
		{
			unsigned char* p,tmp;
			for(i=0; i<nrow; i++)
			{
				for(j=0; j<ncol; j++)
				{
					p=(unsigned char*)&data[i][j];
					tmp=p[0]; p[0]=p[3]; p[3]=tmp;
					tmp=p[1]; p[1]=p[2]; p[2]=tmp;
				}
			}
		}
		float2ibm(data,nrow,ncol);
	}
	else if(fmt[0]+fmt[1]==5)//pc
	{
		// printf("pc\n");
		fseek(fsgy,0,0);
		check=fread(header,3600,1,fsgy);
		for(i=0; i<nrow; i++)
		{
			check=fread(header+3600+i*240,1L,240L,fsgy);
			check=fread(data[i],4L,ncol,fsgy);
		}
		if(fmt[1]==5)
		{
			unsigned char* p,tmp;
			for(i=0; i<nrow; i++)
			{
				for(j=0; j<ncol; j++)
				{
					p=(unsigned char*)&data[i][j];
					tmp=p[0]; p[0]=p[3]; p[3]=tmp;
					tmp=p[1]; p[1]=p[2]; p[2]=tmp;
				}
			}
		}
	}
	else printf("format(byte3224-3225,start from 0) incorrect(%d,%d), check files first, traces and samples needed to define\n",fmt[0],fmt[1]);
	// printf("finished reading\n");
	fclose(fsgy);//data[i][j]=amplitude(ith trace, jth sample), in traces, samples
}
void revisegy(float** A,char* filename,unsigned char* header,int ntraces,int nsamples)
{
	int i,j;
	FILE* fsgy;
	fsgy=fopen(filename,"wb");
	fwrite(header,3600L,1,fsgy);
	float** Atmp=(float**)malloc(ntraces*sizeof(float*)); for(i=0; i<ntraces; i++)Atmp[i]=(float*)malloc(nsamples*sizeof(float));
	for(i=0; i<ntraces; i++)for(j=0; j<nsamples; j++)Atmp[i][j]=A[i][j];
	if(header[3224]+header[3225]==1)//sun,ibm
	{
		if(header[3225]==0)
		{
			unsigned char* p,tmp;
			for(i=0; i<ntraces; i++)
			{
				for(j=0; j<nsamples; j++)
				{
					p=(unsigned char*)&Atmp[i][j];
					tmp=p[0]; p[0]=p[3]; p[3]=tmp;
					tmp=p[1]; p[1]=p[2]; p[2]=tmp;
				}
			}
		}

		ibm2float(Atmp,ntraces,nsamples);
		for(i=0; i<ntraces; i++)
		{
			fwrite(header+3600+240*i,240L,1,fsgy);
			fwrite(Atmp[i],nsamples*4L,1,fsgy);
		}
		free(Atmp);
	}
	else if(header[3224]+header[3225]==5)//default=pc,ieee
	{
		if(header[3225]==5)
		{
			unsigned char* p,tmp;
			for(i=0; i<ntraces; i++)
			{
				for(j=0; j<nsamples; j++)
				{
					p=(unsigned char*)&Atmp[i][j];
					tmp=p[0]; p[0]=p[3]; p[3]=tmp;
					tmp=p[1]; p[1]=p[2]; p[2]=tmp;
				}
			}
		}
		for(i=0; i<ntraces; i++)
		{
			fwrite(header+3600+240*i,240L,1,fsgy);
			fwrite(Atmp[i],nsamples*4L,1,fsgy);
		}
	}
	else printf("unknown type!\n");
	fclose(fsgy);
}
