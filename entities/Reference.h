/*
* Reference.h
*
* Author: ty_yan
*/

#ifndef Reference_H
#define Reference_H

#include "Point.h"

using namespace std;

// each reference point defination
class Ref_Point
{
public:
    Point point;
    double r;
    double r_low;

    vector<double> dis;

    // TODO : delete circle bound (it seems like not used)
    vector<vector<double> > dict_circle;

    // polynomial regression model coefficients of this reference point
    vector<double> coeffs;

    Ref_Point(Point &, double, double);
    Ref_Point();
    void setDisArr(vector<double> &);
    void setDict(vector<vector<double> > &);
    void setCoeffs(vector<double> &);
};


// main reference point of each cluster
class mainRef_Point
{
public:
    Point point;
    double r;
    double r_low;

    vector<Ref_Point> ref_points;
    vector<double> dis;

    // point sorted by i value
    vector<Point> iValuePts;

    // polynomial regression model coefficients of this reference point
    vector<double> coeffs;

    // list for insert point
    vector<InsertPt> insert_list;

    // store the function to fit distribute of i value list
    // y = ax + b
    double a;
    double b;
    int err_min;
    int err_max;

    // lower and upper of each circle
    //2021/5/15 TODO : delete circle bound (it seems like not used)
    vector<vector<double> > dict_circle;

    mainRef_Point(Point &, double, double, vector<Ref_Point> &);
    mainRef_Point();
    void setMainRefDisArr(vector<double> &);
    void setIValuePts(vector<Point> &);
    void setLinear(double, double, int, int);
    void setDict(vector<vector<double> > &);
    void setCoeffs(vector<double> &);
    void setInsertPt(vector<InsertPt> &);
};

// all reference point of whole data
class Ref_Set
{
public:
    vector<mainRef_Point> ref_set;

    Ref_Set(vector<mainRef_Point> &);
    Ref_Set();
};

#endif