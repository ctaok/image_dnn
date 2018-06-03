#include "stdio.h"
#include "stdlib.h"
#include "lex.h"
#include "base_func.h"

int main(int argc, char* argv[])
{
    if (argc != 3) {
        printf("usage: /bin in_file out_file!\n");
    }
    char *in_file = argv[1];
    char *out_file = argv[2];

    printf("begin to load file!\n");
    int cnt_sent = 0;
    FILE *fp_in = fopen(in_file, "rt");
    FILE *fp_out = fopen(out_file, "wt");
    char content[MAX_TEXT_LEN];
    vector<string> sent_str_vec;
    while (!feof(fp_in)) {
        memset(content, '\0', MAX_TEXT_LEN);
        fgets(content, MAX_TEXT_LEN, fp_in);
        char *temp = strtok(content, "\r\n");
        if (temp == NULL || temp[0] == '#') {
            continue;
        }
        cnt_sent++;
        sent_str_vec.push_back(temp);
    }
    printf("load sentence [%d]!\n", cnt_sent);

    printf("begin to join!\n");
    int idx = -1;
    int sent_st = 1;
    srand(12345);
    for (int ii = 0; ii < sent_str_vec.size(); ii++) {
        if (idx == 0) {
            fprintf(fp_out, "\n");
            sent_st = 1;
        }
        if (idx <= 0) {
            int tmp = rand() % 20;
            if (tmp <10) {
                idx = 1;
            } else if (tmp < 14) {
                idx = 2;
            } else if (tmp < 17) {
                idx = 3;
            } else if (tmp < 19) {
                idx = 4;
            } else {
                idx = 5;
            }
        }
        if (sent_str_vec[ii].length() > 100) {
            if (sent_st == 0) {
                fprintf(fp_out, "\n");
            }
            fprintf(fp_out, "%s", sent_str_vec[ii].c_str());
            idx = 0;
        } else {
            fprintf(fp_out, "%s", sent_str_vec[ii].c_str());
            sent_st = 0;
            idx--;
        }
    }
    printf("join over!\n");
    return 0;
}
