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

    printf("begin to sort!\n");
    int cnt_single = 0;
    int item_num = sent_str_vec.size();
    sort(sent_str_vec.begin(), sent_str_vec.end());
    for (int ii = 0; ii < item_num; ii++) {
        int jj = 0;
        for (jj = ii + 1; jj < item_num; jj++) {
            if (sent_str_vec[ii] != sent_str_vec[jj]) {
                break;
            }
        }
        cnt_single++;
        fprintf(fp_out, "%s\n", sent_str_vec[ii].c_str());
        ii = jj - 1;
    }
    printf("single sentence [%d]!\n", cnt_single);
    return 0;
}
