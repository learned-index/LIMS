/*
* Point.cpp
*
* Author: ty_yan
*/

#include "Point.h"
#include <string>

using namespace std;

Point::Point(){
}

Point::Point(vector<double> &coordinate){
    this->coordinate = coordinate;
}

Point::Point(vector<double> &coordinate, unsigned id){
    this->coordinate = coordinate;
    this->id = id;
}

InsertPt::InsertPt(){
}

InsertPt::InsertPt(vector<double> &coordinate, unsigned id){
    this->coordinate = coordinate;
    this->id = id;
}

void InsertPt::setIValue(double dis){
    this->i_value = dis;
}

void Point::setIValue(unsigned long long i_value){
    this->i_value = i_value;
}

Clu_Point::Clu_Point(){
}

Clu_Point::Clu_Point(vector<Point> &clu_point){
    this->clu_point = clu_point;
}

All_Point::All_Point(){
}

All_Point::All_Point(vector<Clu_Point> &all_point){
    this->all_point = all_point;
}