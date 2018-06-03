#include <stdio.h>
#include <vector>
#include <string>
#include <string.h>
#include "OpenCV.h"
#include "Detect3.h"
using namespace std;

int main(int argc, char** argv)
{
	if (argc != 3)
	{
		printf("usage:./bin image_list outfile\n");
		return -1;
	}
	char *in_file = argv[1];
	char *out_file = argv[2];
	FILE *fp = fopen(in_file, "rt");
	FILE *fo = fopen(out_file, "wt");
	if (fp == NULL) {
		printf("Cannot open file %s to read!\n", in_file);
	}
	if (fo == NULL) {
		printf("Cannot open file %s to read!\n", out_file);
	}
	int cnt = 0;
	char content[1024];
	while (!feof(fp)) {
		memset(content, '\0', 1024);
		fgets(content, 1024, fp);
		char *temp = strtok(content, "\t\r\n");
		if (temp == NULL || temp[0] == '\0') {
			continue;
		}
		string image_path = temp;
		temp = strtok(NULL, "\t\r\n");
		if (temp == NULL || temp[0] == '\0') {
			continue;
		}
		int label = atoi(temp);
		IplImage* iplimg = cvLoadImage(image_path.c_str(), -1);
		if (iplimg == 0) {
			printf("load image error [%s]\n", image_path.c_str());
			continue;
		}
		IplImage* src = cvCreateImage(cvSize(iplimg->width,iplimg->height), 8, 1);
		if (iplimg->nChannels == 3) {
			cvCvtColor(iplimg, src, CV_BGR2GRAY);
		} else {
			cvCopy(iplimg, src);
		}
		if (0) {
			uchar *p = (uchar*)(src->imageData);
			for (int ii = 0; ii < src->height; ii++) {
				for (int jj = 0; jj < src->width; jj++) {
					printf("%d ", p[ii * src->width + jj]);
				}
			}
			printf("\n");
		}
		CvMat *tmp_mat = cvCreateMat(src->height, src->width, CV_8U);
		cvConvert(src, tmp_mat);
		CtRectF detect_rect;
		detect_rect.x = 0;
		detect_rect.y = 0;
		detect_rect.width = src->width;
		detect_rect.height = src->height;
		extract_feat(detect_rect, (CtMat*)tmp_mat);

		int feat_cnt = 0;
		string feat_str;
		fprintf(fo, "%d", label);
		for (int ii = 0; ii < _tmpl->rows; ii++) {
			float* p_data = (float*)(_tmpl->data.ptr + ii * _tmpl->step);
			for (int jj = 0; jj < _tmpl->cols; jj++) {
				//feat_cnt++;
				//if (fabs(p_data[jj]) < 0.000001f) {
				//	continue;
				//}
				fprintf(fo, " %d:%.6f", feat_cnt, p_data[jj]);
				feat_cnt++;
			}
		}
		//if (feat_cnt != 18252) {
		if (feat_cnt != 1728) {
			printf("feat size not match: [%d]\n", feat_cnt);
			exit(-1);
		}
		fprintf(fo, "\n");
		if (cnt % 100 == 0) {
			printf("deal num: [%d]\r", cnt);
			fflush(stdout);
		}
		cnt++;
		cvReleaseImage(&iplimg);
		cvReleaseImage(&src);
		extract_feat_release();
	}
	fclose(fp);
	fclose(fo);
	return 0;
}
