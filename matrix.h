//
//  matrix.h
//  matrix
//
//  Created by jiale on 3/5/16.
//  Copyright © 2016 jiale. All rights reserved.
//

#ifndef matrix_h
#define matrix_h


#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include "math_J.h"
using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::ostringstream;
using std::istringstream;
using std::runtime_error;

template <class T>
class matrix {
public:
    typedef vector<vector<T> > matrix_t;
    typedef vector<T> matrix_row_t;
    typedef vector<T> matrix_col_t;
    
public:
    matrix_t m;
    int row, col;// row,col means length, _row means index,from 0
    std::map<int, string> tittle, index;
public:
    matrix()
    :row(0),col(0)
    {
    };
    matrix(int _row, int _col)
    : row(_row),col(_col)
    {
        m.resize(row);
        for(int i = 0; i < row; i++) {
            m[i].resize(col);
        }
    }
    virtual ~matrix()
    {
        
    }
    matrix(const matrix<T> &obj)
    {
        this->m   = obj.m;
        this->row = obj.row;
        this->col = obj.col;
    }
    matrix( vector<T>& myvector, int _row , int  _col)
    throw(runtime_error)
    {
        if(_col*_row != (int)myvector.size()) {
            throw runtime_error("matrix<T>::constructor wrong row/col num");
        }
        row = _row;
        col = _col;
        m.resize(row);
        for(int i = 0; i < row ; ++ i)
        {
            m[i].resize(col);
            for(int j = 0; j < col; ++j)
            {
                m[i][j] = myvector[(i*col)+j];
            }
        }

    }
    void clear()
    {
        row = 0;
        col = 0;
        m.resize(0);
    }
    
    matrix<T> &operator=(const matrix<T> &obj)
    {
        if(this != &obj) {
            this->m   = obj.m;
            this->row = obj.row;
            this->col = obj.col;
        }
            return *this;
    }
    
    // Manipulators.
    void resize(int _row, int _col)
    {
        row = _row;
        col = _col;
        m.resize(row);
        for(int i = 0; i < row; i++) {
            m[i].resize(col);
        }
    }
    void add_row_InTheEnd(matrix_row_t &_row_t)
    {
        m.push_back(_row_t);
        row ++;
    }
    
    void add_row(int _row, matrix_row_t &_row_t)
    throw(runtime_error)
    {
        if(_row >= (row + 1) || _row < 0) {
            throw runtime_error("matrix<T>::add_row:: row out-of-bounds");
        }
        m.resize(row + 1);
        for(int i = row; i > _row; -- i){
            m[i] = m[i-1];
        }
        m[_row] = _row_t;
        row = row + 1;
    }
    void del_row(int _row)
    {
       
        while (_row < row - 1) {
            m[_row] = m[++_row];
        }
        
        m.resize(_row);
        row = _row;
    }
    //
    void add_col(int _col, matrix_col_t &_col_t)
    throw(runtime_error)
    {
        if(_col > col || _col < 0) {
            throw runtime_error("matrix<T>::add_col:: col out-of-bounds");
        }
        for(int i = 0; i < row; i++) {
            m[i].resize(col + 1);
        }
        for(int i = 0; i < row; i++) {
            for(int j = col; j > _col; -- j) {
                m[i][j] = m[i][j-1];
            }
            m[i][col] = _col_t[i];
        }
        col = col + 1;
    }
    void del_col(int _col)
    {
        col = col - 1;
        for(int i = 0; i < row; i++) {
            for(int j = _col; j < col ; ++ j) {
                m[i][j] = m[i][j+1];
            }
            m[i].resize(col);
        }
    }
    void show_row(int _rol){
        for (int j = 0; j < col; ++ j) {
            cout << " " << m[_rol][j] ;
        }
        cout <<" ; ";
    }
    void show(){
        cout << "{\n";
        for (int i = 0; i < row; ++ i) {
            show_row(i);
            cout <<"\n";
        }
        cout <<"}\n";
    }
    // Input/Output.
    friend
    ostream &operator<<(ostream &out, matrix &obj)
    {
        out << " { " << endl;
        for(int i = 0; i < obj.row; i++) {
            for(int j = 0; j < obj.col; j++) {
                if(j == (obj.col-1)) {
                    out << obj.m[i][j] << endl;
                } else if(j == 0) {
                    out << "   " << obj.m[i][j] << " , ";
                } else {
                    out << obj.m[i][j] << " , ";
                }
            }
        }
        out << " } " << endl;
        return out;
    }
    string stream_out() const
    {
        ostringstream out;
        out << row << " " << col << " ";
        for(int i = 0; i < row; i++) {
            for(int j = 0; j < col; j++) {
                out << m[i][j] << " ";
            }
        }
        return out.str();
    }
    void stream_in(string _m)
    {
        istringstream in(_m,istringstream::in);
        in >> row;
        in >> col;
        m.resize(row);
        for(int i = 0; i < row; i++) {
            m[i].resize(col);
        }
        for(int i = 0; i < row; i++) {
            for(int j = 0; j < col; j++) {
                in >> m[i][j];
            }
        }
    }
    // dot know
    string to_html() const
    {
        ostringstream out;
        out << "<table border=\"1\">" << endl;
        for(int i = 0; i < row; i++) {
            out << "<tr>" << endl;
            for(int j = 0; j < col; j++) {
                out << "<td>" << m[i][j] << "</td>" << endl;
            }
            out << "</tr>" << endl;
        }
        out << "</table>" << endl;
        return out.str();
    }
    // Accessors.
    matrix_t &get_data()
    {
        return m;
    }
    int get_row() const
    {
        return row;
    }
    int get_col() const
    {
        return col;
    }
    matrix_col_t get_row(int _row) const
    {
        return m[_row];
    }
    matrix_row_t get_col(int _col) const
    {
        matrix_col_t _col_t;
        _col_t.resize(row);
        for(int i = 0; i < row; i++) {
            _col_t[i] = m[i][_col];
        }
        return _col_t;
    }
    void head(){
        cout << "====== Head ======= \n{ \n";
        for (int i = 0; i <= 4 && i < row; ++i) {
                show_row(i);
                cout << "\n";
            }
        cout << "}\n";
    }
    void tail(){
        cout << "====== Tail ======= \n{ \n";
        if (row >= 5) {
            for (int i = row - 5; i < row; ++i) {
                show_row(i);
                cout << "\n";
            }
        }
        
        else if (row < 5)
        {
            for (int i = 0; i < row; ++i) {
                show_row(i);
                cout << "\n";
            }
        }
        cout << "}\n";
    }
    
// save one line of stock data to the  class
    matrix_row_t save_stock_line(string aLine)
    {
        matrix_row_t temp_row;
        temp_row.resize(7);
        char* pc_line = (char*)aLine.c_str();
        char* pchar;	// one piece of the line
        
        pchar = strtok(pc_line, ",");
        temp_row[0] = pchar; //Date
        temp_row[1] = strtok(NULL, ",");//Open
        temp_row[2] = strtok(NULL, ",");//High
        temp_row[3] = strtok(NULL, ",");//Low
        temp_row[4] = strtok(NULL, ",");//Close
        temp_row[5] = strtok(NULL, ",");//Volume
        temp_row[6] = strtok(NULL, ",");//Adj Close
        return temp_row;
    }
// save tittle if read historical data prices from yahoo
    void save_stock_tittle()
    {
        tittle[0] = "Date";
        tittle[1] = "Open";
        tittle[2] = "High";
        tittle[3] = "Low";
        tittle[4] = "Volume";
        tittle[5] = "Open";
        tittle[6] = "Adj Close";
    }
// read stock data from a csv file on my Desktop
    void read_yahoo_csv(string stock_name)
    {
        std::string file_name = "/Users/jiale/Desktop/" + stock_name +".csv";
        std::fstream myfile(file_name);
        if (!myfile.is_open())      // check if the file exists
        {
            cout << "Can not find file . please check stock name and date!" << endl;
            system("Pause");// return with error code
        }
        string oneLine;
        getline(myfile, oneLine);//jump tittle,read data from the second line
        save_stock_tittle();
        for (int i = 0; i < 2000; ++i)
        {
            if (myfile.eof())
            {
                cout << "End of file, stop reading." << endl;
                break;
            }
            getline(myfile, oneLine);
            if (oneLine.size() == 0)
                break;
            matrix_row_t tmpStockLine= save_stock_line(oneLine);
            add_row_InTheEnd(tmpStockLine);
        }
    }
//read txt data, The accompanying data file data.txt contains daily closing prices for five securities for a
    void read_data_txt(string file_name)
    {
        std::fstream myfile(file_name);
        if (!myfile.is_open())      // check if the file exists
        {
            cout << "Can not find file . please check stock name and date!" << endl;
            system("Pause");// return with error code
        }
        string oneLine;
        for (int i = 0; i < 1000; ++i)
        {
            if (myfile.eof())
            {
                break;
            }
            getline(myfile, oneLine);
            if (oneLine.size() == 0)
                break;
            matrix_row_t tmpStockLine= save_data_line(oneLine);
            add_row_InTheEnd(tmpStockLine);
        }
    }
// save one line from data.txt
    matrix_row_t save_data_line(string aLine)//five securities
    {
        matrix_row_t temp_row;
        temp_row.resize(5);
        char* pc_line = (char*)aLine.c_str();
        char* pchar;	// one piece of the line
        
        pchar = strtok(pc_line, " ");
        temp_row[0] = atof(pchar); //Date
        temp_row[1] = atof(strtok(NULL, " "));//Open
        temp_row[2] = atof(strtok(NULL, " "));//High
        temp_row[3] = atof(strtok(NULL, " "));//Low
        temp_row[4] = atof(strtok(NULL, " "));//Clos
        return temp_row;
    }
// time of two matrix (constructor)
    matrix (matrix m1,matrix m2)
    throw(runtime_error)
    {
        if (m1.col != m2.row) {
            throw runtime_error("matrix<T>::times::col num not match row");
        }
        else
        {   //initialize matrix m
            row = m1.row;
            col = m2.col;
            m.resize(row);
            for(int i = 0; i < row; ++i) {
                m[i].resize(col);
                for(int j = 0; j < col; ++j)
                {
                    vector<T> mycol;
                    for(int k = 0; k < m2.row; ++k)
                    {
                        mycol.push_back(m2.m[k][j]);
                    }
                    m[i][j] = vector_times( m1.m[i], mycol);
                }
            }
        }
    }
    // time of two matrix (function)
    matrix operator *(matrix m2)
    throw(runtime_error)
    {
        matrix m1(0,m2.get_col());
        if (col != m2.row) {
            throw runtime_error("matrix<T>::times::col num not match row");
        }
        else
        {
            for(int i = 0; i < row; ++i)
            {
                vector<T> m1_row;
                for(int j = 0; j < col; ++j)
                {
                    vector<T> mycol;
                    for(int k = 0; k < m2.row; ++k)
                    {
                        mycol.push_back(m2.m[k][j]);
                    }
                    m1_row.push_back( vector_times( m[i], mycol));
                }
                m1.add_row_InTheEnd(m1_row);
            }
        }
        return m1;
    }

// matirx times a constant number
    void times_number(T num)
    {
        for (int i = 0; i < row; ++ i)
        {
            for (int j = 0; j < col; ++ j)
            {
                m[i][j] *= num;
            }
        }
    }
    
  //determinant
    double det()
    {
        int n = row;
        double d = 0;
        int c, subi, i, j, subj;
        matrix<double> submat(n -1,n-1);
        if (n == 2)
        {
            return( (m[0][0] * m[1][1]) - (m[1][0] * m[0][1]));
        }
        else
        {
            for(c = 0; c < n; c++)
            {
                subi = 0;
                for(i = 1; i < n; i++)
                {
                    subj = 0;
                    for(j = 0; j < n; j++)
                    {
                        if (j == c)
                        {
                            continue;
                        }
                        submat.m[subi][subj] = m[i][j];
                        subj++;
                    }
                    subi++;
                }
                d += (pow(-1 ,c) * m[0][c] * submat.det( ));
            }
        }
        return d;
    }
    
// get this matrix's algebraic_complement
    matrix algebraic_complement()
    {
        matrix m2(row,col);
        m2.m = this->m;
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < col; ++j)
                {
                    matrix temp_m(row,col);
                    temp_m.m = this->m;
                    temp_m.del_row(i);
                    temp_m.del_col(j);
                    if((i+j)%2 == 0)
                    {
                        m2.m[j][i] = (T)temp_m.det();
                    }
                    else
                    {
                        m2.m[j][i] = (T)0 - temp_m.det();
                    }
                }
            }
            return m2;
    }
// get inverse of this matrix
    matrix inverse()
    {
        matrix * temp_m = this;
        matrix m2 = this->algebraic_complement();
        m2.times_number(1/(* temp_m).det());
        return  m2;
    }
     
// get transpose of this matrix
    matrix transpose()
    {
        matrix m2(0,row);
        for (int j = 0; j < col; ++ j)
        {
            vector<T> tempvector;
            for (int i = 0; i < row; ++ i)
            {
                tempvector.push_back(m[i][j]);
            }
            m2.add_row_InTheEnd(tempvector);
        }
        return m2;
    }
//use this matrix minue another matrix m2 ,and return the result
    matrix minus(matrix m2)
    {
         matrix m3(0,col);
        if (row == m2.get_col() && col == m2.get_row()) {
            for (int i = 0; i < row; ++ i)
            {
                vector<T> tempvector;
                for (int j = 0; j < col; ++ j)
                {
                    tempvector.push_back(m[i][j] - m2.m[i][j]);
                }
                m3.add_row_InTheEnd(tempvector);
            }
        }
        return  m3;
    }
};



#endif /* matrix_h */
