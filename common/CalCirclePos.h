/*
* CalCirclePos.h
*
* Author: ty_yan
*/

#ifndef CalCirclePos_H
#define CalCirclePos_H

#include "../entities/Point.h"

#include <math.h>

using namespace std;

/*
* Input : 
*   Point refPt
*   double radius_refPt
*   Point queryPt
*   double radius_queryPt
* Output:
*   relation of these two points
*/
class CalCirclePos
{
private:
    unsigned dim;
public:
    double dis_upper;
    double dis_lower;
    /* 
    *  1 ：外切+外离
    *  2 ：相交
    *  3 ：内含
    */
    unsigned label;

    ~CalCirclePos();
    CalCirclePos(Point &, double, Point &, double);
    double CalDisOfPt(Point &, Point &);
};

CalCirclePos::~CalCirclePos()
{
}

CalCirclePos::CalCirclePos(Point &refPt, double radius_refPt, Point &queryPt, double radius_queryPt)
{
    this->dim = refPt.coordinate.size();
    double distance = CalDisOfPt(refPt, queryPt);
    // 外切+外离
    if(distance >= radius_refPt + radius_queryPt){
        this->label = 1;
        this->dis_lower = 0x3f3f3f;
        this->dis_upper = 0x3f3f3f;
    }else if(distance >= fabs(radius_refPt - radius_queryPt)){
        this->label = 2;
        // 没有吧参考点圆心包进去
        if(distance > radius_queryPt){
            this->dis_upper = radius_refPt;
            this->dis_lower = distance - radius_queryPt;
        // 吧参考点圆心包进去了
        }else{
            this->dis_lower = 0.0;
            this->dis_upper = radius_refPt;
        }
    }else{
        this->label = 3;
        if(distance > radius_queryPt){
            this->dis_lower = distance - radius_queryPt;
            this->dis_upper = distance + radius_queryPt;
        }else{
            this->dis_upper = distance + radius_queryPt;
            this->dis_lower = 0.0;
        }
    }
}

double CalCirclePos::CalDisOfPt(Point &point_a, Point &point_b){
    double total = 0.0;
    for(unsigned i = 0; i < dim; i++){
        total += pow(point_a.coordinate[i] - point_b.coordinate[i], 2);
    }
    return sqrt(total);
}

#endif