#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <math.h>
#include "generate_feat.h"

typedef double Elem_t;

int create_feat();
int main(int argc, char *argv[])
{
    if (argc != 1) {
        printf("cmd: ./gen_feat\n");
        return 1;
    }

    create_feat();
    return 0;
}

int create_feat()
{
	//int feat_dim     = 18260;
	int feat_dim     = 1728;
    int label_dim    = 1;
    int max_sent_num = 1000;
    int max_sent_len = 1;
    int global_mean  = 1;
    int is_ftrl      = 0;
    int is_dnn       = 1;

	const int max_text_len = 300000;
    char *file_name = "out.feat";
    char *dir_name = "feat";

    if (max_sent_len != 1) {
        printf("The input data length must be equal to 1!\n");
        return 1;
    }
    if (access(dir_name, 0) != 0) {
        printf("the %s directory does not exists!\n", dir_name);
        return 1;
    }

    vector<Elem_t> mean_vec;
    vector<Elem_t> var_vec;
    int sent_cnt = 0;
    FILE *fp_out = NULL;

    char content[max_text_len];
    if (global_mean) {
        mean_vec.resize(feat_dim, 0.0f);
        var_vec.resize (feat_dim, 0.0f);
        sprintf(content, "%s/%s", dir_name, "global_mean_var");
        fp_out = fopen(content, "wt");
        if (fp_out == NULL) {
            printf("Cannot open file %s to write!\n", content);
            return 1;
        }
    }

    GenerateFeat gen_feat;
	gen_feat.init(feat_dim, label_dim, max_sent_num,
                max_sent_len, dir_name, true);

    vector< vector<float> > feat_mat;
    vector< int > label_vec;
    feat_mat.resize(max_sent_len);
    for (int ii = 0; ii < max_sent_len; ii++) {
        feat_mat[ii].resize(feat_dim);
    }
    label_vec.resize(max_sent_len);

    FILE *fp = fopen(file_name, "rt");
    if (fp == NULL) {
        printf("Cannot open file %s to read!\n", file_name);
    }
    int line_id = 0;
    while (!feof(fp)) {
        memset(content, '\0', max_text_len);
        memset(&feat_mat[0][0], 0, feat_dim * sizeof(float));
        fgets(content, max_text_len, fp);
        char *temp = strtok(content, "\r\n");
        if (temp == NULL) {
            continue;
        }
        line_id++;
		if (line_id % 100 == 0) {
			printf("%d\r", line_id);
			fflush(stdout);
		}

        string str(temp);
        int pos = str.find(" ");

        int pre_pos = 0;
        int cur_pos = 0;

        label_vec[0] = atoi(str.substr(pos-1, 1).c_str());
        label_vec[0] = (label_vec[0] == 1) ? 1 : 0;
        pre_pos = str.find(":", 0) + 1;
	    cur_pos = str.find(" ", 0);
	    if (pre_pos < cur_pos){
            float labelvalue = atof(str.substr(pre_pos, cur_pos - pre_pos).c_str());
            if (label_vec[0] == 1)
                memcpy(&label_vec[0], &labelvalue, sizeof(int));
        }


        pre_pos = 0;
        cur_pos = 0;
	    string sub_str = str.substr(pos + 1, str.size() - pos - 1);
        /*for (int ii = 0; ii < feat_dim; ii++) {
            sprintf(content, " %d:", ii);
            pre_pos = sub_str.find(content, 0);
            if (pre_pos == -1) {
                continue;
            }
            pre_pos = sub_str.find(":", pre_pos) + 1;
            cur_pos = sub_str.find(" ", pre_pos);
            //
            float value = atof(sub_str.substr(pre_pos, cur_pos - pre_pos).c_str());
            //sprintf(content, "%d:%s ", ii, sub_str.substr(pre_pos, cur_pos - pre_pos).c_str());

            feat_mat[0][ii] = value;
        }*/
		int pos_st = 0;
		int pos_ed = 0;
		while (pos_st < sub_str.length()) {
			pos_ed = sub_str.find(" ", pos_st);
			if (pos_ed == -1) {
				pos_ed = sub_str.length();
			}
			//string tmp_str = sub_str.substr(pos_st, pos_ed - pos_st);
			int pos_mid = sub_str.find(":", pos_st);
			int feat_id = atoi((sub_str.substr(pos_st, pos_mid - pos_st)).c_str());
			float feat_val = atof((sub_str.substr(pos_mid + 1, pos_ed - pos_mid - 1)).c_str());
			feat_mat[0][feat_id] = feat_val;
			pos_st = pos_ed + 1;
		}
        
        if (is_dnn) {
            gen_feat.add_one_sent(feat_mat, label_vec, 1);
            if (gen_feat.sent_num() == gen_feat.max_sent_num()) {
                gen_feat.write();
            }
        }
        if (global_mean) {
            for (int ii = 0; ii < 1; ii++) {
                for (int jj = 0; jj < feat_dim; jj++) {
                    mean_vec[jj] += feat_mat[ii][jj];
                    var_vec [jj] += feat_mat[ii][jj] * feat_mat[ii][jj];
                }
            }
            sent_cnt++;
        }
    }
    fclose(fp);
    if (is_dnn && gen_feat.sent_num() != 0) {
        gen_feat.write();
    }

    if (global_mean) {
        printf("%d\n", sent_cnt);
        for (int ii = 0; ii < feat_dim; ii++) {
            mean_vec[ii] = mean_vec[ii] / (Elem_t) sent_cnt;
            var_vec[ii]  = var_vec[ii]  / (Elem_t) sent_cnt;
            var_vec[ii]  = var_vec[ii] - mean_vec[ii] * mean_vec[ii];
            fprintf(fp_out, "%f %f\n", mean_vec[ii], var_vec[ii]);
        }
        fclose(fp_out);
    }
    return 0;
}

