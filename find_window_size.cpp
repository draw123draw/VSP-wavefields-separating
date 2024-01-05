#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"sgy_reader.cpp"
char filename[]="./data/model_data.sgy";
char fbfilename[]="./data/model_data_fb.txt";
struct cmpx
{
    float re;
    float im;
};
void fft(cmpx* x,cmpx* f,int N)
{
    int i,ii;
    float sr,si,pi=3.1415926;
    for(i=0; i<N; i++)
    {
        sr=0;
        si=0;
        for(ii=0; ii<N; ii++)
        {
            sr+=x[ii].re*cos(2*pi*ii*i/N)+x[ii].im*sin(2*pi*ii*i/N);
            si+=-x[ii].re*sin(2*pi*ii*i/N)+x[ii].im*cos(2*pi*ii*i/N);
        }
        f[i].re=sr;
        f[i].im=si;
    }
}
int findmaxidx(float* arr,int n)
{
    int maxidx=0;
    for(int i=1;i<n;i++)if(arr[i]>arr[maxidx])maxidx=i;
    return maxidx;
}
int find_delta_t(char* fbfilename,int nums,float dt)
{
    char split='\n';
    int i,delta_t=0;
    float sum=0;
    FILE* fp;
    fp=fopen(fbfilename,"r");
    fseek(fp,0,SEEK_SET);
    float* arr=(float*)malloc(sizeof(float)*nums);
    for(i=0;i<nums;i++)fscanf(fp,"%f ",&arr[i]);
    for(i=0;i<nums-1;i++)sum+=arr[i+1]-arr[i];
    sum=dt*ceil(abs(sum)/(float)(nums-1));

    free(arr);
    fclose(fp);
    delta_t=ceil(sum);
    return delta_t;
}
int main()
{
    //read the segy file
    long long traces,samples,ts[2];
    traces_samples(filename,ts);
    traces=ts[0];
    samples=ts[1];
    int i,j,fm;
    float dt;
    float** A=(float**)malloc(traces*sizeof(float*)); for(i=0; i<traces; i++)A[i]=(float*)malloc(samples*sizeof(float));
    unsigned char* hdr=(unsigned char*)malloc((3600+240*traces)*sizeof(unsigned char));
    readsgy(filename,A,hdr,traces,samples);
    dt=(float)(hdr[3216]*256+hdr[3217])/1e3;
    float* A_f_mean=(float*)calloc(samples,sizeof(float));
    //mean trace value
    cmpx* in=(cmpx*)malloc(samples*sizeof(cmpx));
    cmpx* out=(cmpx*)malloc(samples*sizeof(cmpx));
    for(i=0;i<traces;i++)
    {
        for(j=0;j<samples;j++)
        {
            in[j].re=A[i][j];
            in[j].im=0;
        }
        fft(in,out,samples);
        for(j=0;j<samples;j++)A_f_mean[j]+=out[j].re*out[j].re+out[j].im*out[j].im;
    }
    fm=findmaxidx(A_f_mean,samples);
    int best_time_win;
    best_time_win=2*ceil(2000/fm);
    int delta_t,best_trace_win;

    //determining the delta_t
    // delta_t=2;
    delta_t=find_delta_t(fbfilename,samples,dt);

    best_trace_win=ceil(best_time_win/delta_t);
    printf("The best time window is %d\nThe best trace window is %d\nSlightly larger than the reference value may be more stable\n",best_time_win,best_trace_win);
    free(A);
    free(A_f_mean);
    free(hdr);
    free(in);
    free(out);
    return 0;
}

