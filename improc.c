#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>
#include "lapack.h"
#define POINT_SIZE 5
#define MATCH_THRESHOLD 0.3
#define SAVE_KEYPOINTS 1

typedef struct {
	int width,height;
	unsigned char* data;
} Bitmap;

typedef struct {
	int size;
	int dimension;
	float* coord;
	void* data;
} Descriptor;

typedef struct {
	int numImages;
	int* image_index;
	int* desc_index;
	float* xi;
	float* yi;
} Match;

typedef struct {
	int numMatches;
	Match* matches;
} MatchArray;

FILE* matchFile;
char* matchFileName;

void colormap (float f,unsigned char *r,unsigned char *g,unsigned char *b) {
	*r=0;
	*g=0;
	*b=0;
	if (f<=0) {
		*b = 128;
	} else if (f <= 0.25) {
		*g = (unsigned char) f / 0.25 * 255;
		*b = (unsigned char) 128 * (1 - f / 0.25);
	} else if (f <= 0.5) {
		*g = 255;
		*r = (unsigned char) (f - 0.25) / 0.25 * 255;
	} else if (f <= 0.75) {
		*r = 255;
		*g = (unsigned char) 255 + (0.5 - f) / 0.25 * 127;
	} else if (f <= 1) {
		*r = 255;
		*g = (unsigned char) 128 * (1 - f) / 0.25;
	} else {
		*r = 255;
	}
}

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

void descriptor_transform(Descriptor** desc, int numDescriptors,int numKeypoints,unsigned char* pixels) {
	int dimension = desc[0]->dimension;
	double* C = calloc(dimension * dimension , sizeof(double));
	double* mean = malloc(dimension * sizeof(double));
	double* lambda_real = malloc(dimension * sizeof(double));
	double* lambda_imag = malloc(dimension * sizeof(double));
	double* V = malloc(dimension * dimension * sizeof(double));
	float* component = calloc(numKeypoints * 3 ,sizeof(float));
	for (int i=0;i<dimension;i++)
		mean[i] = 0.0f;
	for (int m=0;m<numDescriptors;m++) {
		unsigned char* c = desc[m]->data;
		for (int n=0;n<desc[m]->size;n++) {
			for (int i=0;i<dimension;i++) {
				mean[i] += *c++;
			}
		}
	}
	for (int i=0;i<dimension;i++)
		mean[i] /= numKeypoints;
	for (int m=0;m<numDescriptors;m++) {
		unsigned char* val = desc[m]->data;
		for (int n=0;n<desc[m]->size;n++) {
			unsigned char* val_n = val + n * dimension;
			for (int i=0;i<dimension;i++) {
				for (int j=i;j<dimension;j++) {
					C[i*dimension+j] += ((double)val_n[i] - mean[i]) * ((double)val_n[j] - mean[j]);
				}
			}
		}
	}
	for (int i=0;i<dimension;i++) {
		for (int j=i;j<dimension;j++) {
			C[i*dimension+j] /= numKeypoints * numKeypoints;
			C[j*dimension+i] = C[i*dimension+j];
		}
	}
	eigenvalue(dimension,C,lambda_real,lambda_imag,V);
	int tmp,l1=0,l2=1,l3=2;
	if (lambda_real[l2] > lambda_real[l1]) {
		tmp = l1; l1 = l2; l2 = tmp;
	}
	if (lambda_real[l3] > lambda_real[l1]) {
		tmp = l1; l1 = l3; l3 = tmp;
	}
	if (lambda_real[l3] > lambda_real[l2]) {
		tmp = l2; l2 = l3; l3 = tmp;
	}
	for (int i=3;i<dimension;i++) {
		if (lambda_real[i] > lambda_real[l1]) {
			l3 = l2; l2 = l1; l1 = i;			
		} else if (lambda_real[i] > lambda_real[l2]) {
			l3 = l2; l2 = i;
		} else if (lambda_real[i] > lambda_real[l3]) {
			l3 = i;
		}
	}
	for (int i=0;i<numKeypoints*3;i++)
		component[i] = 0.0f;
	float* rgb = component;
	for (int m=0;m<numDescriptors;m++) {
		unsigned char* val = desc[m]->data;
		for (int n=0;n<desc[m]->size;n++) {
			unsigned char* val_n = val + n * dimension;
			for (int i=0;i<dimension;i++) {
				rgb[0] += V[l1*dimension+i] * ((double)val_n[i] - mean[i]);
				rgb[1] += V[l2*dimension+i] * ((double)val_n[i] - mean[i]);
				rgb[2] += V[l3*dimension+i] * ((double)val_n[i] - mean[i]);
			}
			rgb += 3;
		}
	}
	float rmin=component[0],rmax=component[0];
	float gmin=component[1],gmax=component[1];
	float bmin=component[2],bmax=component[2];
	for (int i=1;i<numKeypoints;i++) {
		if (component[i*3] < rmin) rmin = component[i*3];
		else if (component[i*3] > rmax) rmax = component[i*3];
		if (component[i*3 + 1] < gmin) gmin = component[i*3 + 1];
		else if (component[i*3 + 1] > gmax) gmax = component[i*3 + 1];
		if (component[i*3 + 2] < bmin) bmin = component[i*3 + 2];
		else if (component[i*3 + 2] > bmax) bmax = component[i*3 + 2];
	}
	for (int i=0;i<numKeypoints;i++) {
		pixels[i*3] = (component[i*3] - rmin) / (rmax - rmin) * 255;
		pixels[i*3 + 1] = (component[i*3 + 1] - gmin) / (gmax - gmin) * 255;
		pixels[i*3 + 2] = (component[i*3 + 2] - bmin) / (bmax - bmin) * 255;
	}
	free(C);
	free(mean);
	free(lambda_real);
	free(lambda_imag);
	free(V);
	free(component);
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

void gray2rgb(unsigned char* src,unsigned char* dst,int width,int height) {
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			*dst++ = *src;
			*dst++ = *src;
			*dst++ = *src;
			src++;
		}
	}
}

void drawKeyPoint(unsigned char* dst,Descriptor* desc,int k,int width,int height,unsigned char r,unsigned char g,unsigned char b) {
	int x = (int) desc->coord[k*2];
	int y = (int) desc->coord[k*2+1];
	int top = y - POINT_SIZE < 0 ? 0 : y - POINT_SIZE;
	int bottom = y + POINT_SIZE > height ? height : y + POINT_SIZE;
	int left = x - POINT_SIZE < 0 ? 0 : x - POINT_SIZE;
	int right = x + POINT_SIZE > width ? width : x + POINT_SIZE;
	//TOP
	int i=left;
	unsigned char* c = dst + (top*width+left)*3;
	while (i++ < right) {
		*c++ = r;
		*c++ = g;
		*c++ = b;
	}	
	//BOTTOM
	i=left;
	c = dst + (bottom*width+left)*3;
	while (i++ < right) {
		*c++ = r;
		*c++ = g;
		*c++ = b;
	}	
	//LEFT
	i=top;
	c = dst + (top*width+left)*3;
	while (i++ < bottom) {
		c[0] = r;
		c[1] = g;
		c[2] = b;
		c += width * 3;
	}	
	//RIGHT
	i=top;
	c = dst + (top*width+right)*3;
	while (i++ < bottom) {
		c[0] = r;
		c[1] = g;
		c[2] = b;
		c += width * 3;
	}
}

double desc_dist(unsigned char* d1,unsigned char* d2,int dimension) {
	double res = 0;
	for (int i=0;i<dimension;i++) {
		res += (*d1 - *d2) * (*d1 - *d2);
		d1++;
		d2++;
	}
	return res;
}

MatchArray matchDescriptors(Descriptor** desc,int numDescriptors) {
	int* matches = malloc(numDescriptors * sizeof(int));
	bool** visited = malloc(numDescriptors * sizeof(bool*));
	MatchArray matcharray = {0,NULL};
	for (int m=0;m<numDescriptors;m++) {
		Descriptor* d = desc[m];
		visited[m] = calloc(d->size,sizeof(bool));
	}
	for (int m=0;m<numDescriptors;m++) {
		Descriptor* d1 = desc[m];
		for (int i=0;i<d1->size;i++) {
			if (visited[m][i])
				continue;
			memset(matches,0,numDescriptors * sizeof(int));
			bool foundMatch = false;
			for (int n=m+1;n<numDescriptors;n++) {
				Descriptor* d2 = desc[n];
				unsigned char* p1 = (unsigned char*)d1->data + d1->dimension * i;
				unsigned char* p2 = (unsigned char*)d2->data;
				double s1 = desc_dist(p1,p2,d1->dimension);
				p2 += d1->dimension;
				double s2 = desc_dist(p1,p2,d1->dimension);
				int i1 = 0;
				if (s2 < s1) {
					double tmp = s1;
					s1 = s2;
					s2 = tmp;
					i1 = 1;
				}
				p2 += d1->dimension;
				for (int j=2;j<d2->size;j++,p2+=d1->dimension) {
					if (visited[n][j])
						continue;
					double score = desc_dist(p1,p2,d1->dimension);
					if (score < s1) {
						s2 = s1;
						i1 = j;
						s1 = score;
					} else if (score < s2) {
						s2 = score;
					}
				}
				if (s1 < MATCH_THRESHOLD * MATCH_THRESHOLD * s2) {
					visited[m][i] = true;
					visited[n][i1] = true;
					foundMatch = true;
					matches[n] = i1 + 1;
				}
			}
			if (foundMatch) {
				matcharray.numMatches++;
				matcharray.matches = realloc(matcharray.matches,matcharray.numMatches * sizeof(Match));
				if (!matcharray.matches) {
					matcharray.numMatches = 0;
					return matcharray;
				}
				Match* mt = matcharray.matches + matcharray.numMatches - 1;
				mt->numImages = 1;
				mt->image_index = malloc(sizeof(int));
				mt->desc_index = malloc(sizeof(int));
				mt->xi = malloc(sizeof(float));
				mt->yi = malloc(sizeof(float));
				if (!mt->image_index || !mt->desc_index || !mt->xi || !mt->yi) {
					mt->numImages = 0;
					continue;
				}
				mt->image_index[0] = m;
				mt->desc_index[0] = i;
				mt->xi[0] = d1->coord[i*2];
				mt->yi[0] = d1->coord[i*2+1];
				for (int n=m+1;n<numDescriptors;n++) {
					if (matches[n]) {
						Descriptor* d2 = desc[n];
						int i1 = matches[n] - 1;
						mt->numImages++;
						mt->image_index = realloc(mt->image_index,mt->numImages * sizeof(int));
						mt->desc_index = realloc(mt->desc_index,mt->numImages * sizeof(int));
						mt->xi = realloc(mt->xi,mt->numImages * sizeof(float));
						mt->yi = realloc(mt->yi,mt->numImages * sizeof(float));
						if (!mt->image_index || !mt->desc_index || !mt->xi || !mt->yi) {
							mt->numImages = 0;
							break;
						}
						mt->image_index[mt->numImages - 1] = n;
						mt->desc_index[mt->numImages - 1] = i1;
						mt->xi[mt->numImages - 1] = d2->coord[i1*2];
						mt->yi[mt->numImages - 1] = d2->coord[i1*2+1];
					}
				}
			}
		}
	}
	for (int m=0;m<numDescriptors;m++)
		free(visited[m]);
	free(visited);
	free(matches);
	return matcharray;
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

void writePPM(char* filename,unsigned char* imageBuffer,int width,int height) {
	FILE* f = fopen(filename,"w");
	if (!f)
		return;
	fprintf(f,"P6\n%d %d\n255\n",width,height);
	fwrite(imageBuffer,1,width*height*3,f);
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
		fclose(f);
		return b;
	} else {
		free(b);
		fclose(f);
		return NULL;
	}
}

Descriptor* readKeyFile(char* filename) {
	FILE* f = fopen(filename,"r");
	if (!f)
		return NULL;
	Descriptor* desc = malloc(sizeof(Descriptor));
	if (fscanf(f,"%d %d",&desc->size,&desc->dimension)==2) {
		desc->data = malloc(desc->size * desc->dimension);
		desc->coord = malloc(desc->size * 2 * sizeof(float));
		unsigned char* c = desc->data;
		float* d = desc->coord;
		float x,y,scale,orientation;
		for (int i=0;i<desc->size;i++) {
			fscanf(f,"%f %f %f %f",&y,&x,&scale,&orientation);
			*d++ = x;
			*d++ = y;
			for (int j=0;j<desc->dimension;j++) {
				fscanf(f,"%hhu",c++);
			}
		}
		fclose(f);
		return desc;
	} else {
		printf("Error reading %s\n",filename);
		free(desc);
		fclose(f);
		return NULL;
	}
}

void displayHelp() {
	printf("./improc matchFile *.key\n");
}

int main(int argc,char* argv[]) {
	if (argc < 3) {
		displayHelp();
		return 1;
	}
	Descriptor** desc = malloc((argc-2) * sizeof(Descriptor*));
	matchFile = fopen(argv[1],"w");
	if (!matchFile) {
		printf("Cannot open %s\n",argv[1]);
		return 1;
	}
	matchFileName = argv[1];
	for (int i=2;i<argc;i++)
		desc[i-2] = readKeyFile(argv[i]);
	MatchArray matcharray = matchDescriptors(desc,argc-2);
	for (int i=0;i<matcharray.numMatches;i++) {
		Match m = matcharray.matches[i];
		fprintf(matchFile,"%d %f %f",m.image_index[0],m.xi[0],m.yi[0]);
		for (int j=1;j<m.numImages;j++) {
			fprintf(matchFile," %d %f %f",m.image_index[j],m.xi[j],m.yi[j]);
		}
		fprintf(matchFile,"\n");
	}
	printf("Wrote %d matches to %s\n",matcharray.numMatches,matchFileName);
#if SAVE_KEYPOINTS
	int numKeypoints = 0;
	for (int i=0;i<argc-2;i++)
		numKeypoints += desc[i]->size;
	unsigned char* pixels = malloc(numKeypoints * 3);
	descriptor_transform(desc,argc-2,numKeypoints,pixels);
	char buffer[256];
	unsigned char* rgb = pixels;
	for (int i=2;i<argc;i++) {
		int l = strlen(argv[i]);
		int j;
		for (j=l-1;j>=0;j--)
			if (argv[i][j] == '.')
				break;
		memcpy(buffer,argv[i],j+1);
		sprintf(buffer+j+1,"ppm");
		Bitmap* bmp = readPPM(buffer);
		if (bmp) {
			//draw descriptor at each keypoint
			unsigned char* gray = malloc(bmp->width * bmp->height);
			rgb2gray(bmp->data,gray,bmp->width,bmp->height);
			gray2rgb(gray,bmp->data,bmp->width,bmp->height);
			sprintf(buffer+j,"-desc.ppm");
			for (int k=0;k<desc[i-2]->size;k++) {
				drawKeyPoint(bmp->data,desc[i-2],k,bmp->width,bmp->height,rgb[0],rgb[1],rgb[2]);
				rgb += 3;
			}
			writePPM(buffer,bmp->data,bmp->width,bmp->height);
			free(bmp->data);
			free(bmp);
			//draw matched keypoints
			sprintf(buffer+j,".ppm");
			bmp = readPPM(buffer);
			rgb2gray(bmp->data,gray,bmp->width,bmp->height);
			gray2rgb(gray,bmp->data,bmp->width,bmp->height);
			sprintf(buffer+j,"-match.ppm");
			for (int k=0;k<matcharray.numMatches;k++) {
				Match mt = matcharray.matches[k];
				for (int l=0;l<mt.numImages;l++) {
					if (mt.image_index[l] == i-2) {
						unsigned char r,g,b;
						colormap(1.0 * k / matcharray.numMatches, &r, &g, &b);
						drawKeyPoint(bmp->data,desc[i-2],mt.desc_index[l],bmp->width,bmp->height,r,g,b);
						break;
					}
				}
			}
			writePPM(buffer,bmp->data,bmp->width,bmp->height);
			free(bmp->data);
			free(bmp);
			free(gray);
		}
	}
	free(pixels);
#endif
	fclose(matchFile);
	for (int i=0;i<matcharray.numMatches;i++) {
		Match m = matcharray.matches[i];
		if (m.image_index) free(m.image_index);
		if (m.desc_index) free(m.desc_index);
		if (m.xi) free(m.xi);
		if (m.yi) free(m.yi);
	}
	if (matcharray.matches)
		free(matcharray.matches);
	for (int i=0;i<argc-2;i++) {
		free(desc[i]->data);
		free(desc[i]->coord);
		free(desc[i]);
	}
	free(desc);
	return 0;
}
