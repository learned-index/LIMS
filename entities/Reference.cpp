/*
* Reference.cpp
*
* Author: ty_yan
*/

#include "Reference.h"
#include "Point.h"

Ref_Point::Ref_Point(Point &point, double r, double r_low)
{
    this->point = point;
    this->r = r;
    this->r_low = r_low;
}

Ref_Point::Ref_Point()
{
}

void Ref_Point::setDisArr(vector<double> &dis){
    this->dis = dis;
}

void Ref_Point::setDict(vector<vector<double> > &dict){
    this->dict_circle = dict;
}

void Ref_Point::setCoeffs(vector<double> &coeffs){
    this->coeffs = coeffs;
}

mainRef_Point::mainRef_Point(Point &point, double r, double r_low, vector<Ref_Point> &ref_points)
{
    this->point = point;
    this->r = r;
    this->r_low = r_low;
    this->ref_points = ref_points;
}

mainRef_Point::mainRef_Point()
{
}

void mainRef_Point::setMainRefDisArr(vector<double> &dis){
    this->dis = dis;
}

void mainRef_Point::setDict(vector<vector<double> > &dict){
    this->dict_circle = dict;
}


void mainRef_Point::setIValuePts(vector<Point> &iValuePts){
    this->iValuePts = iValuePts;
}

void mainRef_Point::setInsertPt(vector<InsertPt> &insert_list){
    this->insert_list = insert_list;
}

void mainRef_Point::setLinear(double a, double b, int err_min, int err_max){
    this->a = a;
    this->b = b;
    this->err_min = err_min;
    this->err_max = err_max;
}

void mainRef_Point::setCoeffs(vector<double> &coeffs){
    this->coeffs = coeffs;
}

Ref_Set::Ref_Set(vector<mainRef_Point> &ref_set)
{
    this->ref_set = ref_set;
}

Ref_Set::Ref_Set()
{
}