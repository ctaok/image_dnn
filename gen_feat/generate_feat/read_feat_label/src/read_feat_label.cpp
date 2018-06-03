#include "stdio.h"
#include "stdlib.h"
#include "lex.h"
#include "base_func.h"

vector<double> _global_mean_vec;
vector<double> _global_var_vec;
long long _global_frame_cnt = 0;

void swap32(char *c) {
    char t0 = c[0];
    char t1 = c[1];
    c[0] = c[3];
    c[1] = c[2];
    c[2] = t1;
    c[3] = t0;
}
void swap16(char *c) {
    char t0 = c[0];
    c[0] = c[1];
    c[1] = t0;
}

int main(int argc, char* argv[])
{
    if (argc != 2) {
        printf("usage: /bin in_feat_file!\n");
        exit(0);
    }
    string in_feat_file = argv[1];
	int pos = in_feat_file.rfind(".");
	string pre_str = in_feat_file.substr(0, pos);
	string tail_str = in_feat_file.substr(pos + 1, in_feat_file.length() - pos - 1);
	if (strcmp(tail_str.c_str(), "label") != 0 && strcmp(tail_str.c_str(), "feat") != 0) {
		printf("not support type of %s\n", tail_str.c_str());
		exit(0);
	}

	char feat_path_str[1024];
	char label_path_str[1024];
	sprintf(feat_path_str, "%s.feat", pre_str.c_str());
	sprintf(label_path_str, "%s.label", pre_str.c_str());

    FILE* fi_feat = fopen(feat_path_str, "rb");
    FILE* fi_label = fopen(label_path_str, "rb");
    int sample_num = 0;
    int sample_period = 0;
    int sample_size = 0;
    int sample_kind = 0;
    int pre_len = 3;
    fread(&sample_num,    sizeof(int), 1, fi_feat);
    fread(&sample_period, sizeof(int), 1, fi_feat);
    fread(&sample_size,   sizeof(int), 1, fi_feat);
    fread(&sample_kind,   sizeof(short), 1, fi_feat);
    swap32((char*)&sample_num);
    swap32((char*)&sample_period);
    swap16((char*)&sample_size);
    swap16((char*)&sample_kind);

    vector<float> feat_vec;
    int label;
    //int mean_feat_dim = sample_size / sizeof(float);
    int mean_feat_dim = 18260;
    //int mean_feat_dim = 1728;

    feat_vec.resize(mean_feat_dim);
    for (int ii = 0; ii < pre_len; ii++) {
        fread(&feat_vec[0], sizeof(float), mean_feat_dim, fi_feat);
        fread(&label, sizeof(int), 1, fi_label);
    }
    vector<float> feat_vec_256;
    feat_vec_256.resize(256);
    for (int ii = 0; ii < sample_num - pre_len; ii++) {
		printf("[%d]\t", ii);
        fread(&feat_vec[0], sizeof(float), mean_feat_dim, fi_feat);
        fread(&label, sizeof(int), 1, fi_label);
        for (int jj = 0; jj < mean_feat_dim; jj++) {
            swap32((char*)&feat_vec[jj]);
            printf("%d[%.6f] ", jj, feat_vec[jj]);
        }
		printf("\tlabel: [%d]\n", label);
    }
    fclose(fi_feat);
    fclose(fi_label);

    return 0;
}
