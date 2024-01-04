/*There is only 'fftw3-3.lib' in the lib folder, as a result, the compiler may not find it in the system of Linux.
However, this can be solved using 'find_window_size.cpp'*/
#include<stdio.h>
#include<math.h>
#include"fftw3.h"
#include"sgy_reader.cpp"
char filename[]=R"(data\model_data.sgy)";
char fbfilename[]=R"(data\model_data_fb.txt)";
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
	FILE *fp;
	fp=fopen(fbfilename,"r");
	fseek(fp, 0, SEEK_SET);
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
    fftw_plan p;
    fftw_complex* in=(fftw_complex*)fftw_malloc(samples*sizeof(fftw_complex));
    fftw_complex* out=(fftw_complex*)fftw_malloc(samples*sizeof(fftw_complex));
    for(i=0;i<traces;i++)
    {
        for(j=0;j<samples;j++)
        {
            in[j][0]=A[i][j];
            in[j][1]=0;
        }
        p=fftw_plan_dft_1d(samples,in,out,FFTW_BACKWARD,FFTW_MEASURE);
        fftw_execute(p);
        for(j=0;j<samples;j++)A_f_mean[j]+=out[j][0]*out[j][0]+out[j][1]*out[j][1];
    }
    fm=findmaxidx(A_f_mean,samples);
    int best_time_win;
    best_time_win=2*ceil(2000/fm);
    int delta_t,best_trace_win;

    //determining the delta_t
    // delta_t=2;
    delta_t=find_delta_t(fbfilename,samples,dt);

    best_trace_win=ceil(best_time_win/delta_t);
    printf("The best time window may be %d\nThe best trace window may be %d\nSlightly larger than the reference value may be more stable",best_time_win,best_trace_win);
    free(A);
    free(A_f_mean);
    free(hdr);
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    return 0;
}

