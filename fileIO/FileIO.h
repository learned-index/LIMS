/*
* FileIO.h
*
* Created on 04/17/2021
* Author: ty_yan
*/

#ifndef FILEIO_H
#define FILEIO_H

#include <iostream>
#include <fstream>
#include <string>
#include "../entities/Point.h"
#include "../entities/Reference.h"

using namespace std;

// load data set from filename.txt
class InputReader
{
public:
    ifstream fin;
    vector<Point> points;

    InputReader();
    InputReader(string);
    vector<Point> LoadRangeQuery(string);
    vector<Point> LoadPointQuery(string);
    void parse();
    Clu_Point getCluster();
};

class OutputPrinter{
public:
    ofstream fout;
    Clu_Point cluster;
    vector<Point> result;
    string outputFileName;
    mainRef_Point mainRef_point;

    // 输出一维值
    OutputPrinter();
    OutputPrinter(string, Clu_Point);
    // 输出结果点
    OutputPrinter(string, vector<Point>);
    // 输出距离环号对应值
    OutputPrinter(string, mainRef_Point &);
    OutputPrinter(string, vector<double> &);

    // 按照一维值顺序以二进制存储数据文件
    OutputPrinter(int, mainRef_Point &);

    void print();
    void print_result();
    void print_refCircle();
    void print_errResult(string, vector<double> &);
    void print_ivalue(string, vector<Point> &);
};

#endif