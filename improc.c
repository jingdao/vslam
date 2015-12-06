#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "lapack.h"

typedef struct {
	int width,height;
	unsigned char* data;
} Bitmap;

void eigenvalue(int N,double* A,double* lambda_real,double* lambda_imag,double* v) {
	int info,ldvl=1,ldvr=N,lwork=15*N;	
	double *work = malloc(lwork * sizeof(double));
	char jobvl = 'N', jobvr = 'V';
	dgeev_(&jobvl,&jobvr, &N, A, &N, lambda_real, lambda_imag,
	    NULL,&ldvl, v, &ldvr ,work, &lwork, &info);
	if (info!=0) {
		printf("Error in subroutine dgeev_ (info=%d)\n",info);
	}
	free(work);
}

void pca_transform(unsigned char *src, unsigned char* dst,int width,int height) {
	double C[9];
	C[0] = 0;
	C[1] = 0;
	C[2] = 0;
	C[4] = 0;
	C[5] = 0;
	C[8] = 0;
	double rmean=0,gmean=0,bmean=0;
	unsigned char *p = src;
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			rmean += *p++;
			gmean += *p++;
			bmean += *p++;
		}
	}
	rmean /= width * height;
	gmean /= width * height;
	bmean /= width * height;
	p = src;
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			double r = *p++;
			double g = *p++;
			double b = *p++;
			C[0] += (r-rmean) * (r-rmean);
			C[1] += (r-rmean) * (g-gmean);
			C[2] += (r-rmean) * (b-bmean);
			C[4] += (g-gmean) * (g-gmean);
			C[5] += (g-gmean) * (b-bmean);
			C[8] += (b-bmean) * (b-bmean);
		}
	}
	C[0] /= width * height;
	C[1] /= width * height;
	C[2] /= width * height;
	C[4] /= width * height;
	C[5] /= width * height;
	C[8] /= width * height;
	C[3] = C[1];
	C[6] = C[2];
	C[7] = C[5];
	double lambda_real[3];
	double lambda_imag[3];
	double V[9];
	eigenvalue(3,C,lambda_real,lambda_imag,V);
	int maxL = 0;
	maxL = lambda_real[1] > lambda_real[maxL] ? 1 : maxL;
	maxL = lambda_real[2] > lambda_real[maxL] ? 2 : maxL;
	if (V[maxL*3] < 0) {
		V[maxL*3] = -V[maxL*3];
		V[maxL*3+1] = -V[maxL*3+1];
		V[maxL*3+2] = -V[maxL*3+2];
	}
	double *T = malloc(width*height*sizeof(double));
	p = src;
	double *q = T;
	double ymin=DBL_MAX,ymax=-DBL_MAX;
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			double r = *p++;
			double g = *p++;
			double b = *p++;
			double y = (r-rmean)*V[maxL*3] + (g-gmean)*V[maxL*3+1] + (b-bmean)*V[maxL*3+2];
			*q++ = y;
			ymin = y < ymin ? y : ymin;
			ymax = y > ymax ? y : ymax;
		}
	}
	q = T;
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			*dst++ = (*q++ - ymin) / (ymax - ymin) * 255;
		}
	}
	free(T);
}

void histeq(unsigned char* src,unsigned char* dst,int width,int height) {
	int count[256];
	memset(count,0,256*sizeof(int));
	unsigned char *p = src;
	for (int i=0;i<height;i++) 
		for (int j=0;j<width;j++) 
			count[*p++]++;
	for (int i=1;i<256;i++)
		count[i] += count[i-1];
	for (int i=0;i<height;i++) 
		for (int j=0;j<width;j++)
			*dst++ = count[*src++] * 255 / width / height;

}

void rgb2gray(unsigned char* src,unsigned char* dst,int width,int height) {
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			*dst++ = (src[0] + src[1] + src[2]) / 3;
			src += 3;
		}
	}
}

void writePGM(char* filename,unsigned char* imageBuffer,int width,int height) {
	FILE* f = fopen(filename,"w");
	if (!f)
		return;
	fprintf(f,"P5\n%d %d\n255\n",width,height);
	fwrite(imageBuffer,1,width*height,f);
	printf("Saved to %s\n",filename);
	fclose(f);
}

Bitmap* readPPM(char* filename) {
	FILE* f = fopen(filename,"r");
	if (!f)
		return NULL;
	Bitmap* b = malloc(sizeof(Bitmap));
	char buffer[64];
	fgets(buffer,64,f);
	fgets(buffer,64,f);
	if (sscanf(buffer,"%d %d",&b->width,&b->height)==2) {
		b->data = malloc(b->width * b->height * 3);
		fread(b->data,1,b->width*b->height*3,f);
		return b;
	} else {
		free(b);
		return NULL;
	}
}

int main(int argc,char* argv[]) {
	Bitmap* b = readPPM(argv[1]);
	if (b) {
		unsigned char* pgm = malloc(b->width*b->height);
		unsigned char* pgm2 = malloc(b->width*b->height);
//		rgb2gray(b->data,pgm,b->width,b->height);
//		histeq(pgm,pgm2,b->width,b->height);
		pca_transform(b->data,pgm2,b->width,b->height);
		writePGM(argv[2],pgm2,b->width,b->height);
		free(pgm);
		free(pgm2);
		free(b->data);
		free(b);
	}	
}
