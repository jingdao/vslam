#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "geometry.h"

struct Blob {
	Coord coord;
	int value;
};

int compareBlob(const void* v1,const void* v2) {
	Blob* b1 = (Blob*)v1;
	Blob* b2 = (Blob*)v2;
	return abs(b2->value) - abs(b1->value);
}

void drawPoint(int* imageBuffer,int width,int height,Coord point,int pointSize,int color);

void conv2(unsigned char* src,int* dest,int* kernel,int width,int height,int kernel_size) {
	for (int i=0;i<height-kernel_size;i++) {
		for (int j=0;j<width-kernel_size;j++) {

		}
	}
}

void laplacian(int* src,int width,int height,Blob* blob) {
	for (int i=1;i<height-1;i++) {
		for (int j=1;j<width-1;j++) {
			int acc = 0;
			acc += -8 * (src[i*width+j] & 0xFF);
			acc += src[(i-1)*width+j-1] & 0xFF;
			acc += src[(i-1)*width+j] & 0xFF;
			acc += src[(i-1)*width+j+1] & 0xFF;
			acc += src[i*width+j-1] & 0xFF;
			acc += src[i*width+j+1] & 0xFF;
			acc += src[(i+1)*width+j-1] & 0xFF;
			acc += src[(i+1)*width+j] & 0xFF;
			acc += src[(i+1)*width+j+1] & 0xFF;
			blob->coord.x = j;
			blob->coord.y = i;
			blob->value = acc;
			blob++;
		}
	}
}

void sobel(int* src,int width,int height, Blob* blob) {
	for (int i=1;i<height-1;i++) {
		for (int j=1;j<width-1;j++) {
			int gx=0,gy=0;
			gx -= 2 * src[i*width+j-1];
			gx += 2 * src[i*width+j+1];
			gx -= src[(i-1)*width+j-1];
			gx += src[(i-1)*width+j+1];
			gx -= src[(i+1)*width+j-1];
			gx += src[(i+1)*width+j+1];
			gy -= 2 * src[(i-1)*width+j];
			gy += 2 * src[(i+1)*width+j];
			gy -= src[(i-1)*width+j-1];
			gy += src[(i+1)*width+j-1];
			gy -= src[(i-1)*width+j+1];
			gy += src[(i+1)*width+j+1];
			blob->coord.x = j;
			blob->coord.y = i;
			blob->value = abs(gx) + abs(gy);
			blob++;
		}
	}
}

void pointDetection(int* src,int width,int height) {
	int numBlobs = (width-2)*(height-2);
	Blob* blob = new Blob[numBlobs];
//	laplacian(src,width,height,blob);
	sobel(src,width,height,blob);
	qsort(blob,numBlobs,sizeof(Blob),compareBlob);
	for (size_t i=0;i<200;i++)
		drawPoint(src,width,height,blob[i].coord,2,255<<16);
	delete[] blob;
}

//int main(int argc,char* argv[]) {
//	return 0;
//}
