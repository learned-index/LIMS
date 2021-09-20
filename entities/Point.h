/*
* Point.h
*
* Author: ty_yan
*/

#ifndef Point_H
#define Point_H

#include <vector>
#include <string>

using namespace std;

// each point data
class Point{
public:
    // coordinate of each point of n-dimentional data
    vector<double> coordinate;

    // id of each point
    unsigned id;

    // 1-d value of point
    unsigned long long i_value;

    Point();
    Point(vector<double> &);
    Point(vector<double> &, unsigned);
    void setIValue(unsigned long long);
};

class InsertPt{
public:
    // coordinate of each point of n-dimentional data
    vector<double> coordinate;

    // id of each point
    unsigned id;

    // 1-d value of point
    double i_value;

    InsertPt();
    InsertPt(vector<double> &, unsigned);
    void setIValue(double);
};


// each cluster of data
class Clu_Point{
public:
    vector<Point> clu_point;

    Clu_Point();
    Clu_Point(vector<Point> &);
};

// all data point
class All_Point{
public:
    vector<Clu_Point> all_point;

    All_Point();
    All_Point(vector<Clu_Point> &);
};

#endif