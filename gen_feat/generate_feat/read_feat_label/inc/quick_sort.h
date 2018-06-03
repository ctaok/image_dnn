/*
* quick_sort.h
* Author: fuxiaoyin
* Created on: 2016-10-24
* Copyright (c) Ainirobot.com, Inc. All Rights Reserved
*/
#ifndef QUICK_SORT_H
#define QUICK_SORT_H

#include <iostream>
#include <assert.h>

//ANS_BEG

//------------------------------
//  QuickSort
//------------------------------

template <class T>
class QuickSort {
protected:
    T* _array;
    unsigned int* _sort_idx;
    unsigned int _size;
    unsigned int _low;
    unsigned int _high;
    unsigned int _owner;

public:
    QuickSort(T* num, unsigned int n, unsigned int low, unsigned int high, int owner) {
        _size    = n;
        _low   = low;
        _high  = high;
        _owner  = owner;
        _sort_idx = NULL;

        if (_owner == 1) {
            _array = new T[_size];

            for (unsigned int ii = low; ii <= high; ii++) {
                _array[ii] = num[ii];
            }
        } else {
            _array = num;
        }
    }

    QuickSort(T* num, unsigned int* sort_idx, unsigned int n,
            unsigned int low, unsigned int high, int owner) {

        unsigned int ii;
        _size    = n;
        _low   = low;
        _high  = high;
        _owner  = owner;

        if (_owner == 1) {
            _array = new T[_size];

            for (ii = low; ii <= high; ii++) {
                _array[ii] = num[ii];
            }

            sort_idx = new unsigned int[_size];

            for (ii = low; ii <= high; ii++) {
                sort_idx[ii] = ii;
            }
        } else {
            _array = num;
            _sort_idx = sort_idx;
        }
    }

    virtual ~QuickSort() {
        if (_array && _owner) {
            delete []_array;
        }

        if (_sort_idx && _owner) {
            delete []_sort_idx;
        }
    }
    void sort() {
        partion_sort(_low, _high);
    }

    void sort_with_idx() {
        begin_sort_with_idx(_low, _high);
    }

    void out_result() {
        if (_sort_idx == NULL) {
            for (unsigned int i = 0; i < _size; i++) {
                std::cout << _array[i] << " " ;
            }
        } else {
            for (unsigned int i = 0; i < _size; i++) {
                std::cout << _array[_sort_idx[i]] << " " ;
            }
        }

        std::cout << std::endl;
    }

private:
    void partion_sort(int left, int right) {
        int i, j, mid;
        i = left;
        j = right;

        mid = (left + right) / 2;
        T pivotkey = _array[mid];

        do {
            while ((_array[i] < pivotkey) && (i < right)) {
                i++;
            }

            while ((_array[j] > pivotkey) && (j > left)) {
                j--;
            }

            if (i <= j) {
                T tempval  = _array[i];
                _array[i] = _array[j];
                _array[j] = tempval;
                i++;
                j--;
            }
        } while (i <= j);

        if (left < j) {
            partion_sort(left, j);
        }

        if (right > i) {
            partion_sort(i, right);
        }
    }
    void begin_sort_with_idx(unsigned int low, unsigned int high) {
        partion_with_idx(low, high);
    }

    void partion_with_idx(int left, int right) {
        int i, j, middle, tempidx;
        i = left;
        j = right;
        //
        middle = (left + right) / 2;
        T itemp = _array[_sort_idx[middle]];

        //
        do {
            while ((_array[_sort_idx[i]] < itemp) && (i < right)) {
                i++;
            }

            while ((_array[_sort_idx[j]] > itemp) && (j > left)) {
                j--;
            }

            if (i <= j) {
                tempidx =  _sort_idx[i];
                _sort_idx[i] = _sort_idx[j];
                _sort_idx[j] = tempidx;
                i++;
                j--;
            }
        } while (i <= j);

        if (left < j) {
            partion_with_idx(left, j);
        }

        if (right > i) {
            partion_with_idx(i, right);
        }
    }
};

//ANS_END

#endif

