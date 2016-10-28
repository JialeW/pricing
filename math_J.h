//
//  math.h
//  matrix
//
//  Created by jiale on 3/6/16.
//  Copyright © 2016 jiale. All rights reserved.
//

#ifndef math_h
#define math_h
#include <vector>
using std::vector;

template <typename T>
T mean(vector<T> myvector)
{
    T sum = 0;
    for(int i = 0; i < myvector.size() ; ++ i)
    {
        sum += myvector[i];
    }
    return (T)(sum / myvector.size());
}

template <typename T>
T var(vector<T> myvector)
{
    T Mean = mean(myvector);
    T sum = 0;
    for(int i = 0; i < myvector.size() ; ++ i)
    {
        sum += (myvector[i] - Mean)*(myvector[i] - Mean);
    }
    return (T)(sum/(myvector.size()-1));
};

template <typename T>
T vector_times(vector<T> _row_t,vector<T> _col_t)
throw(std::runtime_error)
{
    if (_row_t.size() != _col_t.size()) {
        throw std::runtime_error("matrix<T>::vector_times::col num not match row");
    }
    T result = 0;
    for(int i = 0; i < _row_t.size(); ++ i)
    {
        result += _row_t[i]*_col_t[i];
    }
    return result;
};






#endif /* math_h */
