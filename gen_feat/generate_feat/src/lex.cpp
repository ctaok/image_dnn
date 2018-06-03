#include "lex.h"
#include <time.h>

CLex::CLex() {
    _up_limit  = 3000000;
    _pst_voc = NULL;
    _ps_voc  = NULL;
    _add_voc = 5;
    _sen_max = 256;
    _unk_num  = 0;
    _base_unk = 0;
}

CLex::~CLex() {
    if (_pst_voc != NULL) {
        delete [] _pst_voc;
        _pst_voc = NULL;
    }

    if (_ps_voc  != NULL) {
        delete [] _ps_voc;
        _ps_voc  = NULL;
    }
}

bool CLex::generate_lexicon(char* corpus_name, char* lexicon_name) {
    if (corpus_name == NULL || lexicon_name == NULL) {
        printf("CLex::generate_lexicon\tThe file name should not be NULL!\n");
        return false;
    }

    FILE* fpin   = NULL;
    FILE* fpout  = NULL;
    FILE* tempfp = NULL;
    fpin = fopen(corpus_name, "rt");

    if (fpin == NULL) {
        printf("CLex::generate_lexicon\tCannot open file %s to read!\n", corpus_name);
        return false;
    }

    //
    int ii = 0;
    int jj = 0;
    int part_file  = 0;
    int line_num   = 0;
    int item_num   = 0;
    int vocab_size = 20000000;
    char sz_content[4096];
    char temp_file[1024];
    char* cur_voc = NULL;
    char* nxt_voc = NULL;
    char* temp = NULL;
    std::vector<ConstVoc> vocab_vec;
    vocab_vec.resize(vocab_size);
    std::vector<int> seq_vec;
    ConstVoc voc_item;

    while (!feof(fpin)) {
        memset(sz_content, '\0', 4096);
        fgets(sz_content, 4096, fpin);

        if (sz_content == NULL || strlen(sz_content) < 2) {
            continue;
        }

        temp = strtok(sz_content, "\r\n");

        if (temp == NULL || strlen(temp) == 0) {
            continue;
        }

        //
        temp = strtok(temp, " ");

        while (temp != NULL) {
            if (strlen(temp) > 250) {
                printf("CLex::generate_lexicon\tThe word length is too long! %s\n", temp);
                temp = strtok(NULL, " ");
                continue;
            }

            voc_item.init(temp);
            vocab_vec[item_num] = voc_item;
            item_num++;
            temp = strtok(NULL, " ");
        }

        line_num++;

        if (line_num % 100000 == 0) {
            printf("File %d, line %d\n", part_file, line_num);
        }

        if (item_num > vocab_size - 4096) {
            memset(temp_file, '\0', 1024);
            snprintf(temp_file, 1024, "%s_Lex_%d.txt", corpus_name, part_file);
            tempfp = fopen(temp_file, "wt");

            if (tempfp == NULL) {
                printf("CLex::generate_lexicon\tCannot open file %s to write!\n", temp_file);
                return false;
            }

            //item_num = (int)vocab_vec.size();
            if (item_num == 0) {
                continue;
            }

            seq_vec.resize(item_num);

            for (ii = 0; ii < item_num; ii++) {
                seq_vec[ii] = ii;
            }

            QuickSort<ConstVoc> sort_class(&vocab_vec[0], (unsigned int*)(&seq_vec[0]),
                    item_num, 0, item_num - 1, 0);
            printf("\nBegin to sort lexicon in %s\n", temp_file);
            sort_class.sort_with_idx();
            printf("end sort\n\n");

            for (ii = 0; ii < item_num; ii++) {
                cur_voc = vocab_vec[seq_vec[ii]]._voc;

                for (jj = ii + 1; jj < item_num; jj++) {
                    nxt_voc = vocab_vec[seq_vec[jj]]._voc;

                    if (strcmp(cur_voc, nxt_voc) != 0) {
                        break;
                    }
                }

                fprintf(tempfp, "%s\n", cur_voc);
                ii = jj - 1;
            }

            fclose(tempfp);
            //
            vocab_vec.clear();
            vocab_vec.resize(vocab_size);
            item_num = 0;
            part_file++;
        }
    }

    memset(temp_file, '\0', 1024);
    snprintf(temp_file, 1024, "%s_Lex_%d.txt", corpus_name, part_file);
    tempfp = fopen(temp_file, "wt");

    if (tempfp == NULL) {
        printf("CLex::generate_lexicon\tCannot open file %s to write!\n", temp_file);
        return false;
    }

    //item_num = (int)vocab_vec.size();
    if (item_num != 0) {
        seq_vec.resize(item_num);

        for (ii = 0; ii < item_num; ii++) {
            seq_vec[ii] = ii;
        }

        QuickSort<ConstVoc> sort_class(&vocab_vec[0], (unsigned int*)(&seq_vec[0]), item_num, 0,
                item_num - 1, 0);
        printf("\nBegin to sort lexicon in %s\n", temp_file);
        sort_class.sort_with_idx();
        printf("end sort\n\n");

        for (ii = 0; ii < item_num; ii++) {
            cur_voc = vocab_vec[seq_vec[ii]]._voc;

            for (jj = ii + 1; jj < item_num; jj++) {
                nxt_voc = vocab_vec[seq_vec[jj]]._voc;

                if (strcmp(cur_voc, nxt_voc) != 0) {
                    break;
                }
            }

            fprintf(tempfp, "%s\n", cur_voc);
            ii = jj - 1;
        }

        fclose(tempfp);
        //
        vocab_vec.clear();
        part_file++;
    }

    //
    fpout = fopen(lexicon_name, "wt");

    if (fpout == NULL) {
        printf("CLex::generate_lexicon\tCannot open file %s to write!\n", lexicon_name);
        return false;
    }

    for (ii = 0; ii < part_file; ii++) {
        memset(temp_file, '\0', 1024);
        snprintf(temp_file, 1024, "%s_Lex_%d.txt", corpus_name, ii);
        tempfp = fopen(temp_file, "rt");

        if (tempfp == NULL) {
            printf("CLex::generate_lexicon\tCannot open file %s to read!\n", temp_file);
            return false;
        }

        while (!feof(tempfp)) {
            memset(sz_content, '\0', 4096);
            fgets(sz_content, 4096, tempfp);

            if (sz_content == NULL || strlen(sz_content) < 2) {
                continue;
            }

            temp = strtok(sz_content, "\r\n");

            if (temp == NULL) {
                continue;
            }

            voc_item.init(temp);
            vocab_vec.push_back(voc_item);
        }

        fclose(tempfp);
        remove(temp_file);
    }

    item_num = (int)vocab_vec.size();
    line_num = 0;

    if (item_num != 0) {
        seq_vec.resize(item_num);

        //fprintf(fpout, "%d\n", item_num);
        for (ii = 0; ii < item_num; ii++) {
            seq_vec[ii] = ii;
        }

        QuickSort<ConstVoc> sort_class(&vocab_vec[0], (unsigned int*)(&seq_vec[0]), item_num, 0,
                item_num - 1, 0);
        sort_class.sort_with_idx();

        for (ii = 0; ii < item_num; ii++) {
            cur_voc = vocab_vec[seq_vec[ii]]._voc;

            for (jj = ii + 1; jj < item_num; jj++) {
                nxt_voc = vocab_vec[seq_vec[jj]]._voc;

                if (strcmp(cur_voc, nxt_voc) != 0) {
                    break;
                }
            }

            fprintf(fpout, "%s\n", cur_voc);
            line_num++;
            ii = jj - 1;
        }

        fclose(fpout);
    }

    //
    memset(temp_file, '\0', 1024);
    snprintf(temp_file, 1024, "%s_temp.txt", lexicon_name);
    tempfp = fopen(temp_file, "wt");

    if (tempfp == NULL) {
        printf("CLex::generate_lexicon\tCannot open file %s to write!\n", temp_file);
        return false;
    }

    fpin = fopen(lexicon_name, "rt");

    if (fpin == NULL) {
        printf("CLex::generate_lexicon\tCannot open file %s to read!\n", lexicon_name);
        return false;
    }

    fprintf(tempfp, "%d\n", line_num);
    line_num = 1;

    while (!feof(fpin)) {
        memset(sz_content, '\0', 4096);
        fgets(sz_content, 4096, fpin);

        if (sz_content == NULL || strlen(sz_content) < 2) {
            continue;
        }

        temp = strtok(sz_content, "\r\n");

        if (temp == NULL || strlen(temp) == 0) {
            continue;
        }

        fprintf(tempfp, "%d\t%s\n", line_num, temp);
        line_num++;
    }

    fclose(tempfp);
    fclose(fpin);
    //
    remove(lexicon_name);
    rename(temp_file, lexicon_name);
    //
    printf("CLex::generate_lexicon\tGenerate Lexicon file %s successfully!\n\n", lexicon_name);
    //
    return true;
}

bool CLex::sort_lexicon(char* corpus_name, char* SortedLex, char* lexicon_name) {
    if (lexicon_name == NULL) {
        char templexicon_name[1024];
        snprintf(templexicon_name, 1024, "%s_lexicon.txt", corpus_name);
        generate_lexicon(corpus_name, templexicon_name);
        load_lexicon(templexicon_name);
    } else {
        load_lexicon(lexicon_name);
    }

    FILE* fpin  = NULL;
    FILE* fpout = NULL;
    int ii = 0;
    int line_num = 0;
    int word_id = 0;
    int voc_size = 0;
    int item_num = 0;
    std::vector<int> voc_count_vec;
    std::vector<int> seq_vec;
    char sz_content[4096];
    char* temp = NULL;

    fpin = fopen(corpus_name, "rt");

    if (fpin == NULL) {
        printf("CLex::SortLexicon\tCannot open file %s to read!\n", corpus_name);
        return false;
    }

    fpout = fopen(SortedLex, "wt");

    if (fpout == NULL) {
        printf("CLex::SortLexicon\tCannot open file %s to write!\n", SortedLex);
        return false;
    }

    voc_size  = get_voc_size();
    voc_count_vec.resize(voc_size, 0);

    //
    while (!feof(fpin)) {
        memset(sz_content, '\0', 4096);
        fgets(sz_content, 4096, fpin);

        if (sz_content == NULL || strlen(sz_content) < 2) {
            continue;
        }

        temp = strtok(sz_content, "\r\n");

        if (temp == NULL || strlen(temp) == 0) {
            continue;
        }

        temp = strtok(temp, " ");

        while (temp != NULL) {
            if (strlen(temp) > 250) {
                printf("CLex::SortLexicon\tThe word length is too long! %s\n", temp);
                temp = strtok(NULL, " ");
                continue;
            }

            word_id = find_word(temp);

            //word_id = 0;
            if (word_id == -1 || word_id > voc_size) {
                printf("CLex::SortLexicon\tCannot find word %s(%d) in Lexicon!\n", temp, word_id);
            } else {
                voc_count_vec[word_id - 1]++;
            }

            temp = strtok(NULL, " ");
        }

        line_num++;

        if (line_num % 100000 == 0) {
            printf("line %d\n", line_num);
        }
    }

    item_num = voc_size;

    if (item_num != 0) {
        seq_vec.resize(item_num);

        for (ii = 0; ii < item_num; ii++) {
            seq_vec[ii] = ii;
        }

        QuickSort<int> sort_class(&voc_count_vec[0], (unsigned int*)(&seq_vec[0]),
                item_num, 0, item_num - 1, 0);
        printf("\nBegin to sort\n");
        sort_class.sort_with_idx();
        printf("end sort\n\n");

        //
        for (ii = 0; ii < item_num; ii++) {
            fprintf(fpout, "%s\t%d\n", get_lexicon_by_id(seq_vec[item_num - ii - 1] + 1),
                    voc_count_vec[seq_vec[item_num - ii - 1]]);
        }
    }

    fclose(fpin);
    fclose(fpout);
    return true;
}

bool CLex::combine_lexicon(char* lexicon1, char* lexicon2, char* comb_lexicon) {
    FILE* fpin1 = NULL;
    FILE* fpin2 = NULL;
    FILE* fpout = NULL;
    int ii = 0;
    int jj = 0;
    int line_num = 0;
    int item_num = 0;
    std::vector<ConstVoc> vocab_vec;
    std::vector<int>        seq_vec;
    ConstVoc voc_item;
    char sz_content[4096];
    char* cur_voc = NULL;
    char* nxt_voc = NULL;
    char* temp = NULL;
    //
    fpin1 = fopen(lexicon1, "rt");

    if (fpin1 == NULL) {
        printf("CLex::CombineLexicon\tCannot open file %s to read!\n", lexicon1);
        return false;
    }

    fpin2 = fopen(lexicon2, "rt");

    if (fpin2 == NULL) {
        printf("CLex::CombineLexicon\tCannot open file %s to read!\n", lexicon2);
        return false;
    }

    memset(sz_content, '\0', 4096);
    fgets(sz_content, 4096, fpin1);
    temp = strtok(sz_content, "\r\n");
    item_num += atoi(temp);
    //
    memset(sz_content, '\0', 4096);
    fgets(sz_content, 4096, fpin2);
    temp = strtok(sz_content, "\r\n");
    item_num += atoi(temp);
    //
    vocab_vec.resize(item_num);
    //
    line_num = 0;

    while (!feof(fpin1)) {
        memset(sz_content, '\0', 4096);
        fgets(sz_content, 4096, fpin1);

        if (sz_content == NULL || strlen(sz_content) < 2) {
            continue;
        }

        temp = strtok(sz_content, "\r\n");

        if (temp == NULL || strlen(temp) == 0) {
            continue;
        }

        temp = strtok(temp, "\t");
        temp = strtok(NULL, "\t");

        if (strlen(temp) > 250) {
            printf("CLex::CombineLexicon\tThe word length is too long! %s\n", temp);
            continue;
        }

        voc_item.init(temp);
        vocab_vec[line_num] = voc_item;
        line_num++;

        if (line_num % 100000 == 0) {
            printf("line %d\n", line_num);
        }
    }

    while (!feof(fpin2)) {
        memset(sz_content, '\0', 4096);
        fgets(sz_content, 4096, fpin2);

        if (sz_content == NULL || strlen(sz_content) < 2) {
            continue;
        }

        temp = strtok(sz_content, "\r\n");

        if (temp == NULL || strlen(temp) == 0) {
            continue;
        }

        temp = strtok(temp, "\t");
        temp = strtok(NULL, "\t");

        if (strlen(temp) > 250) {
            printf("CLex::CombineLexicon\tThe word length is too long! %s\n", temp);
            continue;
        }

        voc_item.init(temp);
        vocab_vec[line_num] = voc_item;
        line_num++;

        if (line_num % 100000 == 0) {
            printf("line %d\n", line_num);
        }
    }

    fclose(fpin1);
    fclose(fpin2);
    //
    item_num = line_num;
    line_num = 0;

    if (item_num != 0) {
        fpout = fopen(comb_lexicon, "wt");

        if (fpout == NULL) {
            printf("CLex::CombineLexicon\tCannot open file %s to write!\n", comb_lexicon);
            return false;
        }

        seq_vec.resize(item_num);

        for (ii = 0; ii < item_num; ii++) {
            seq_vec[ii] = ii;
        }

        QuickSort<ConstVoc> sort_class(&vocab_vec[0], (unsigned int*)(&seq_vec[0]),
                item_num, 0, item_num - 1, 0);
        printf("\nBegin to sort\n");
        sort_class.sort_with_idx();
        printf("end sort\n\n");

        //
        //fprintf(fpout, "%d\n", item_num);
        for (ii = 0; ii < item_num; ii++) {
            cur_voc = vocab_vec[seq_vec[ii]]._voc;

            for (jj = ii + 1; jj < item_num; jj++) {
                nxt_voc = vocab_vec[seq_vec[jj]]._voc;

                if (strcmp(cur_voc, nxt_voc) != 0) {
                    break;
                }
            }

            fprintf(fpout, "%s\n", cur_voc);
            line_num++;
            ii = jj - 1;
        }

        fclose(fpout);
        //
        char temp_file[1024];
        memset(temp_file, '\0', 1024);
        snprintf(temp_file, 1024, "%s_temp.txt", comb_lexicon);
        FILE* tempfp = fopen(temp_file, "wt");

        if (tempfp == NULL) {
            printf("CLex::CombineLexicon\tCannot open file %s to write!\n", temp_file);
            return false;
        }

        fpin1 = fopen(comb_lexicon, "rt");

        if (fpin1 == NULL) {
            printf("CLex::CombineLexicon\tCannot open file %s to read!\n", comb_lexicon);
            return false;
        }

        fprintf(tempfp, "%d\n", line_num);
        line_num = 1;

        while (!feof(fpin1)) {
            memset(sz_content, '\0', 4096);
            fgets(sz_content, 4096, fpin1);

            if (sz_content == NULL || strlen(sz_content) < 2) {
                continue;
            }

            temp = strtok(sz_content, "\r\n");

            if (temp == NULL || strlen(temp) == 0) {
                continue;
            }

            fprintf(tempfp, "%d\t%s\n", line_num, temp);
            line_num++;
        }

        fclose(tempfp);
        fclose(fpin1);
        //
        remove(comb_lexicon);
        rename(temp_file, comb_lexicon);
        //
        printf("CLex::CombineLexicon\tCombine Lexicon file %s successfully!\n\n", comb_lexicon);
    } else {
        return false;
    }

    return true;
}

bool CLex::map_sri_lexicon(char* lexicon1, char* lexicon2) {
    FILE* fpin = NULL;
    FILE* fpout = NULL;
    int ii = 0;
    int word_num = 0;
    int tm_word_id = 0;
    int lm_word_id = 0;
    char sz_out_file[1024];
    char sz_content[1024];
    char* temp = NULL;
    std::string str;
    fpin = fopen(lexicon1, "rt");

    if (fpin == NULL) {
        printf("CLex::MapLexicon\tCannot open %s file to read!\n", lexicon1);
        return false;
    }

    str.append(lexicon1);
    ii = str.find_last_of("/");
    snprintf(sz_out_file, 1024, "%s/Lexicons_Map.txt", str.substr(0, ii).c_str());
    fpout = fopen(sz_out_file, "wt");

    if (fpin == NULL) {
        printf("CLex::MapLexicon\tCannot open %s file to write!\n", sz_out_file);
        return false;
    }

    load_lexicon(lexicon2);
    //
    memset(sz_content, '\0', 1024);
    fgets(sz_content, 1024, fpin);
    temp     = strtok(sz_content, "\r\n");
    word_num = atoi(temp);
    fprintf(fpout, "%d\n", word_num);
    ii = 0;

    while (!feof(fpin)) {
        memset(sz_content, '\0', 1024);
        fgets(sz_content, 1024, fpin);

        if (sz_content == NULL || strlen(sz_content) < 2) {
            continue;
        }

        temp = strtok(sz_content, "\r\n");

        if (temp == NULL || strlen(temp) == 0) {
            continue;
        }

        temp      = strtok(temp, "\t");
        tm_word_id = atoi(temp);
        temp      = strtok(NULL, "\t");
        lm_word_id = find_word(temp);

        if (lm_word_id != -1) {
            fprintf(fpout, "%d\t%d\n", tm_word_id, lm_word_id);
        } else {
            fprintf(fpout, "%d\t0\n", tm_word_id);
        }

        ii++;

        if (ii % 10000 == 0) {
            printf("ii = %d\n", ii);
        }
    }

    fclose(fpin);

    if (ii != word_num) {
        printf("The TM Lexicon number (%d) is not equal to the real number (%d)!\n", word_num, ii);
        fclose(fpout);
        remove(sz_out_file);
        return false;
    }

    for (ii = 0; ii < 3; ii++) {
        fprintf(fpout, "%d\t%d\n", word_num + ii + 1, get_voc_size() + ii + 1);
    }

    return true;
}

bool CLex::load_lexicon(const char* filename) {
    FILE* fp = fopen(filename, "rt");

    if (fp == NULL) {
        printf("CLex::load_lexicon\tCannot open %s file to read!\n", filename);
        exit(1);
    }

    _voc_size = 0;
    fscanf(fp, "%d", &_voc_size);
    //load lexicon without <s>, </s>, UNK
    _pst_voc = new s_vocabulary[_voc_size + _add_voc + _sen_max];

    if (_pst_voc == NULL) {
        printf("CLex::load_lexicon\t_pst_vocChi can't allocate enough memory!\n");
        exit(1);
    }

    _ps_voc = new char[_voc_size + _add_voc + _sen_max + 1][256];

    if (_ps_voc == NULL) {
        printf("CLex::load_lexicon\t_ps_vocChi can't allocate enough memory!\n");
        exit(1);
    }

    _ps_voc[0][0] = '\0';

    for (int i = 0; i < int(_voc_size); i++) {
        fscanf(fp, "%d\t%s\n", &_pst_voc[i].word_id, _pst_voc[i].pc_word);
        strcpy(_ps_voc[_pst_voc[i].word_id], _pst_voc[i].pc_word);

        if (ferror(fp)) {
            return false;
        }
    }

    fclose(fp);

    //Sort lexicon
    sort_lexicon(0, _voc_size - 1);
    //
    //Add the special mark to lexicon : <s>, </s>,<unk>,_x1,_x2
    _pst_voc[_voc_size].word_id = _voc_size + 1;
    strcpy(_pst_voc[_voc_size].pc_word, "<s>");
    strcpy(_ps_voc[_voc_size + 1], "<s>");
    _pst_voc[_voc_size + 1].word_id = _voc_size + 2;
    strcpy(_pst_voc[_voc_size + 1].pc_word, "</s>");
    strcpy(_ps_voc[_voc_size + 2], "</s>");
    _pst_voc[_voc_size + 2].word_id = _voc_size + 3;
    strcpy(_pst_voc[_voc_size + 2].pc_word, "<unk>");
    strcpy(_ps_voc[_voc_size + 3], "<unk>");

    _pst_voc[_voc_size + 3].word_id = _voc_size + 4;
    strcpy(_pst_voc[_voc_size + 3].pc_word, "_x1");
    strcpy(_ps_voc[_voc_size + 4], "_x1");
    _pst_voc[_voc_size + 4].word_id = _voc_size + 5;
    strcpy(_pst_voc[_voc_size + 4].pc_word, "_x2");
    strcpy(_ps_voc[_voc_size + 5], "_x2");
    //
    _unk_id = _voc_size + 3;
    _x1_id  = _voc_size + 4;
    _x2_id  = _voc_size + 5;
    //
    //_voc_size += _add_voc;
    return true;
}

int CLex::add_unk_word(char* word) {
    _pst_voc[_voc_size + _add_voc + _unk_num].word_id = _voc_size + _add_voc + _unk_num + 1;
    strcpy(_pst_voc[_voc_size + _add_voc + _unk_num].pc_word, word);
    strcpy(_ps_voc[_voc_size + _add_voc + _unk_num + 1], word);
    _unk_num++;

    if (_unk_num + 10 == _sen_max) {
        printf("CLex::add_unk_word\tThere are too many unk words! _unk_num = %d\n", _unk_num);
        exit(1);
    }

    return _voc_size + _add_voc + _unk_num; //UnkNum++, so here just return _voc_size+UnkNum
}

void CLex::sort_lexicon(int left, int right) {
    char word1[256];
    char word2[256];
    int id1 = 0;
    int i = 0;
    int j = 0;
    int middle = 0;

    i = left;
    j = right;
    middle = (left + right) / 2 ;
    strcpy(word1, _pst_voc[middle].pc_word);

    //
    do {
        while (strcmp(_pst_voc[i].pc_word, word1) < 0) {
            i++;
        }

        while (strcmp(_pst_voc[j].pc_word, word1) > 0) {
            j--;
        }

        if (i <= j) {
            id1 =  _pst_voc[i].word_id;
            strcpy(word2, _pst_voc[i].pc_word);
            _pst_voc[i].word_id = _pst_voc[j].word_id;
            strcpy(_pst_voc[i].pc_word, _pst_voc[j].pc_word);
            _pst_voc[j].word_id = id1;
            strcpy(_pst_voc[j].pc_word, word2);
            i++;
            j--;
        }

    } while (i <= j);

    if (left < j) {
        sort_lexicon(left, j);
    }

    if (right > i) {
        sort_lexicon(i, right);
    }
}

char ret_str_lex_get_lexicon_by_id[] = "**";
char* CLex::get_lexicon_by_id(int wordid) {
    if ((wordid < 0) || (wordid > int(_voc_size + _add_voc + _unk_num))) {
        return ret_str_lex_get_lexicon_by_id;
    } else {
        return _ps_voc[wordid];
    }
}

int CLex::find_word(const char* word) {
    int sign = -1;
    int ul_top = 0;
    int ul_middle = 0;
    int ul_bottom = _voc_size;

    while (ul_top <= ul_bottom) {   //*1
        ul_middle = (ul_top + ul_bottom) / 2;

        if (strcmp(word, _pst_voc[ul_middle].pc_word) == 0) {      //*2
            return _pst_voc[ul_middle].word_id;
        }//*2

        if (strcmp(word, _pst_voc[ul_middle].pc_word) > 0) {
            ul_top = ul_middle + 1;
        }

        if (strcmp(word, _pst_voc[ul_middle].pc_word) < 0) {
            if (ul_middle > 0) {
                ul_bottom = ul_middle - 1;
            } else {
                break;
            }
        }
    }//*1

    for (int ii = _voc_size; ii < _voc_size + _add_voc + _unk_num; ii++) {
        if (strcmp(word, _pst_voc[ii].pc_word) == 0) {
            return _pst_voc[ii].word_id;
        }
    }

    return sign;
}

int CLex::analysis_sentence(char* sen, char** word, int* wordid) {
    int sign = 0;
    int len = 0;
    int idx = 0;
    char sep[] = " \r\t\n";
    len = strlen(sen);
    char sz_lex[200];
    int newlen = 0;

    while (strlen(sen) > 1) {
        len = strlen(sen);
        sen = strtok(sen, sep);
        newlen = strlen(sen);
        strcpy(sz_lex, sen);
        strcpy(word[sign], sz_lex);
        idx = find_word(word[sign]);

        if (idx != -1) {
            wordid[sign] = idx;
        } else {
            if (strcmp("<s>", word[sign]) == 0) {
                wordid[sign] = _voc_size + 2;
            } else if (strcmp("</s>", word[sign]) == 0) {
                wordid[sign] = _voc_size + 3;
            } else {
                wordid[sign] = _voc_size + 1; //unk word id = _voc_size+1
            }
        }

        sign++;

        if ((len - newlen) < 2) {
            break;
        }

        sen += strlen(sen) + 1;

        while (sen[0] == ' ') {
            sen++;
        }
    }

    return sign;
}

