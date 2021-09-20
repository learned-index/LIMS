/*
* CaculateIValue.h
*
* Author: ty_yan
*/

#ifndef CaculateIValue_H
#define CaculateIValue_H

#include "ChooseRef.h"
#include "../entities/Reference.h"
#include "../entities/Point.h"
#include "../common/Constants.h"

#include <string>
#include <algorithm>
#include <math.h>

using namespace std;

// This function is to caculate the 1 dimentional value of each point by multi-reference point
class CaculateIValue
{
private:
    unsigned num_refSet;
    unsigned dim;
    int num_data;
    unsigned n;

    Clu_Point cluster;
    mainRef_Point mainRef_point;
public:
    CaculateIValue();
    CaculateIValue(Clu_Point &, mainRef_Point);

    double CaculateEuclideanDis(Point &, Point &);

    int BinarySearch(vector<double> &, double);

    void RadixSort();

    void CaculLinearFun_A();
    void CaculLinearFun_B();
    
    Clu_Point getCluster();
    mainRef_Point getMainRef_Point();
};

CaculateIValue::CaculateIValue()
{
}

CaculateIValue::CaculateIValue(Clu_Point &clu, mainRef_Point mainRefPoint)
{
    this->cluster = clu;
    this->mainRef_point = mainRefPoint;
    this->num_refSet = mainRef_point.ref_points.size();
    this->dim = cluster.clu_point[0].coordinate.size();
    this->num_data = cluster.clu_point.size();

    int n = 1;
    while(n <= Constants::NUM_CIRCLE - 1)
        n *= 10;
    this->n = n;

    
    int split = ceil((double)num_data / Constants::NUM_CIRCLE);
    for(int i = 0; i < num_data; i++){
        unsigned long long i_value = 0;
        double dis_pt_mainRefPt = CaculateEuclideanDis(cluster.clu_point[i], mainRef_point.point);
        int rank_cir = BinarySearch(mainRef_point.dis, dis_pt_mainRefPt) / split;
        i_value += rank_cir;
        for(unsigned j = 0; j < num_refSet; j++){
            double dis_pt_refPt = CaculateEuclideanDis(cluster.clu_point[i], mainRef_point.ref_points[j].point);
            int rank_refCir = BinarySearch(mainRef_point.ref_points[j].dis, dis_pt_refPt) / split;
            i_value = i_value * n + rank_refCir;
        }
        cluster.clu_point[i].setIValue(i_value);
    }

    // radix sort i value list
    RadixSort();

    // caculate y = ax + b to fit i value list distribute
    // CaculLinearFun_A();
    CaculLinearFun_B();
}

// This function search postion(rank) of point in reference point dis list
// which means the dis between ref point & point must exit
// The function of search rank of query point in ref point dis list need to be realized in another function
int CaculateIValue::BinarySearch(vector<double> &disList, double target){
    int low = 0;
    int high = num_data - 1;
    int middle = 0;
    while(low <= high) {
		middle = (high - low) / 2 + low;
		if(target < disList[middle]) {
			high = middle - 1;
		} else if(target > disList[middle]) {
			low = middle + 1;
		}else{
            if((middle > 0 && disList[middle - 1] < target) || (middle == 0))//找到开始位置
            {
                break;
            }
            else
                high = middle - 1;
        }
    }
    return middle;
}

double CaculateIValue::CaculateEuclideanDis(Point &point_a, Point &point_b){
    double total = 0.0;
    for(unsigned i = 0; i < dim; i++){
        total += pow(point_a.coordinate[i] - point_b.coordinate[i], 2);
    }
    return sqrt(total);
}

/* Sort point by i_value using radix sort
* Store sorted array in this format:
* Point position : array index
* array content : point coordinate, point i_value
*/ 
void CaculateIValue::RadixSort(){
    vector<vector<Point> > buckets(num_data);
    vector<Point> cache = cluster.clu_point;
    for(unsigned pos = 0; pos < num_refSet + 1; pos++){
        for(int i = 0; i < num_data; i++){
            unsigned long long  value = cache[i].i_value;
            int index = value % (int)pow(n, pos + 1) / pow(n, pos);
            buckets[index].push_back(cache[i]);
        }
        int j = 0;
        for(auto & theBucket : buckets){
            for(Point & pt : theBucket){
                cache[j++] = pt;
            }
            theBucket.clear();
        }
    }
    mainRef_point.setIValuePts(cache);
}

// This function is to caculate a , b in y = ax + b of i value list
void CaculateIValue::CaculLinearFun_A(){
    double sum_xy = 0.0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_pow2x = 0.0;
    vector<double> iValueList;
    double num_pow = pow(10, num_refSet + 1);
    for(int i = 0; i < num_data; i++){
        unsigned long long i_value = mainRef_point.iValuePts[i].i_value;
  
        double x = i_value / num_pow;

        iValueList.push_back(x);
        sum_x += x;
        sum_y += i;
        sum_pow2x += x * x;
        sum_xy += x * i;
    }
    double avg_x = sum_x / num_data;
    double avg_y = sum_y / num_data;
    double avg_xy = sum_xy / num_data;
    double avg_pow2x = sum_pow2x / num_data;
    double a = (avg_xy - avg_x * avg_y) / (avg_pow2x - avg_x * avg_x);
    double b = avg_y - a * avg_x;
    int err_min = 0;
    int err_max = 0;
    for(int i = 0; i < num_data; i++){
        int CaculY = (int)(a * iValueList[i] + b);
        if(CaculY > i){
            err_min = max(err_min, CaculY - i);
        }else if(CaculY < i){
            err_max = max(err_max, i - CaculY);
        }
    }
    mainRef_point.setLinear(a, b, err_min, err_max);
}

void CaculateIValue::CaculLinearFun_B(){
    double sum_xy = 0.0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_pow2x = 0.0;
    vector<double> iValueList;
    vector<int> pos;
    // 环号为两位数时适用，环号若不是两位数要修改下方式子
    double num_pow = pow(10, num_refSet + 1);

    unsigned long long before = 0;
    int total_num = 0;
    for(int i = 0; i < num_data; i++){
        unsigned long long i_value = mainRef_point.iValuePts[i].i_value;

        if(before == i_value){
            continue;
        }

        before = i_value;
        total_num++;

        double x = i_value / num_pow;

        iValueList.push_back(x);
        pos.push_back(i);
        sum_x += x;
        sum_y += i;
        sum_pow2x += x * x;
        sum_xy += x * i;
    }
    pos.push_back(mainRef_point.iValuePts.size());

    // cout << total_num << endl;

    double avg_x = sum_x / total_num;
    double avg_y = sum_y / total_num;
    double avg_xy = sum_xy / total_num;
    double avg_pow2x = sum_pow2x / total_num;
    double a = (avg_xy - avg_x * avg_y) / (avg_pow2x - avg_x * avg_x);
    double b = avg_y - a * avg_x;
    int err_min = 0;
    int err_max = 0;
    for(int i = 0; i < total_num; i++){
        int CaculY = (int)(a * iValueList[i] + b);
        if(CaculY > pos[i]){
            err_min = max(err_min, CaculY - pos[i]);
        }else if(CaculY < pos[i]){
            err_max = max(err_max, pos[i+1]-1 - CaculY);
        }
    }
    mainRef_point.setLinear(a, b, err_min, err_max);
    // cout << "err bound is : \"min : " << err_min << " max : " << err_max << "\"" << endl;
}

Clu_Point CaculateIValue::getCluster(){
    return cluster;
}

mainRef_Point CaculateIValue::getMainRef_Point(){
    return mainRef_point;
}

#endif