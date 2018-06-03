#include "stdio.h"
#include "stdlib.h"
#include "lex.h"
#include "base_func.h"

int main(int argc, char* argv[])
{
    if (argc != 4) {
        printf("usage: /bin in_file out_file lex_dict!\n");
    }
    char *in_file = argv[1];
    char *out_file = argv[2];
    char *punc_dict = argv[3];

	//过滤句末无标点的句子
    CLex punc_lex;
    punc_lex.load_lexicon(punc_dict);
    int punc_unk = punc_lex.find_word("<unk>");

    FILE *fp_in = fopen(in_file, "rt");
    FILE *fp_out = fopen(out_file, "wt");
    char content[MAX_TEXT_LEN];
    int sent_id[MAX_TEXT_LEN];
    while (!feof(fp_in)) {
        fgets(content, MAX_TEXT_LEN, fp_in);
        char *temp = strtok(content, "\r\n");
        if (temp == NULL || temp[0] == '#') {
            continue;
        }
        string word_str;
        seg_by_char(temp, word_str);
        int sent_num = punc_lex.analysis_sent(word_str.c_str(), sent_id);
        if (sent_id[sent_num - 1] < punc_unk && sent_id[sent_num - 1] > 0) {
            fprintf(fp_out, "%s\n", temp);
        }
    }
    return 0;
}
