#ifndef INC_BASE_FUNC_H_
#define INC_BASE_FUNC_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

#define MAX_TEXT_LEN 1024
using namespace std;

int seg_by_char(string in_str, string &out_str)
{
	out_str = in_str;
	for (int ii = 0; ii < out_str.size(); ii++) {
		if (out_str[ii] < 0) {
			ii += 2;
			if (ii + 1 < out_str.size()) {
				out_str.insert(ii, " ");
			}
		} else if (ii + 1 < out_str.size() && out_str[ii + 1] <0){
			ii++;
			out_str.insert(ii, " ");
		}
	}
	return 0;
}

#endif  // INC_BASE_FUNC_H_

