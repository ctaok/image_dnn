/*
* lex.h
* Author: fuxiaoyin
* Created on: 2016-11-05
* Copyright (c) Ainirobot.com, Inc. All Rights Reserved
*/
#ifndef _LEX_H
#define _LEX_H

#include "commons.h"
#include "error.h"
#include "quick_sort.h"

class ConstVoc {
public:
    ConstVoc() {
        memset(_voc, '\0', 256);
    }
    ~ConstVoc() {};
    //
    void init(char* in_str) {
        memset(_voc, '\0', 256);
        strcpy(_voc, in_str);
    }
    ConstVoc& operator = (const ConstVoc& item) {
        strcpy(_voc, item._voc);
        return (*this);
    }
    bool operator < (const class ConstVoc& item) {
        return (strcmp(_voc, item._voc) < 0) ? true : false;
    }
    bool operator > (const class ConstVoc& item) {
        return (strcmp(_voc, item._voc) > 0) ? true : false;
    }
public:
    char _voc[256];
};

class CLex {
public:
    CLex();
    ~CLex();
    //
    bool  load_lexicon_null(char* filename);
    bool  load_lexicon(const char* filename);
    int   analysis_sent(const char* sen, int* wordid);
    int   analysis_sent(const char* sen, int* wordid, vector<string> &unk_str_vec);
    int   remove_unk_words(int *word_ids, int word_num);
    char* get_lexicon_by_id(int wordid);
	string get_lexicon_by_id(int *word_ids, int word_num, vector<string> unk_str_vec, int &idx_unk);
    bool  generate_lexicon(char* CorpusName, char* LexiconName);
    bool  sort_lexicon(char* CorpusName, char* SortedLex, char* LexiconName);
    bool  combine_lexicon(char* Lexicon1,   char* Lexicon2,  char* CombLexicon);
    bool  map_sri_lexicon(char* Lexicon1,   char* Lexicon2);
    //
    int  get_voc_size() {
        return _voc_size;
    };
    int  get_lex_x1() {
        return _x1_id;
    };
    int  get_lex_x2() {
        return _x2_id;
    };
    int  get_lex_unk() {
        return _unk_id;
    };
    void reset_unk_word() {
        _unk_num  = _base_unk;
    }
    void set_unk_word() {
        _base_unk = _unk_num;
    }
    void clear_unk_word() {
        _unk_num  = 0;
    }
    //
    int  find_word(const char* word);
    int  add_unk_word(char* Word);
    //
private:
    void sort_lexicon(int left, int right);
private:
    int _voc_size;
    struct s_vocabulary {
        int word_id;
        char pc_word[256];
    }* _pst_voc;
    char(*_ps_voc)[256];
    int _unk_num, _base_unk;
    //
    int _up_limit;
    int _add_voc;
    int _x1_id;
    int _x2_id;
    int _unk_id;
    int _sen_max;
};

#endif
