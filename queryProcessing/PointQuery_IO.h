/*
* PointQuery_IO.h
* Author: ty_yan
*/

#ifndef PointQuery_IO_H
#define PointQuery_IO_H

#include "../entities/Point.h"
#include "../entities/Reference.h"
#include "../common/Constants.h"

#include "math.h"
#include <iostream>
#include <fstream>

#define FLOAT_ERROR 1e-8

using namespace std;

class PointQuery_IO{
public:

PointQuery_IO();
PointQuery_IO(Point &, mainRef_Point &, int &, int,double, int &);
double CaculateEuclideanDis(Point &, Point &);
int DoublingSearch(vector<double> &, bool, double, int);
int BinarySearchForDoubling(vector<double> &, double, int, int);
};

PointQuery_IO::PointQuery_IO(){}

PointQuery_IO::PointQuery_IO(Point &queryPt, mainRef_Point &mainRef_point, int &query_res, int clu_id,double dist_qt_mainRef, int &page){
    // double dist_qt_mainRef = CaculateEuclideanDis(queryPt, mainRef_point.point);

    vector<int> rank_list;
    unsigned refPtSize = mainRef_point.ref_points.size();
    rank_list.resize(refPtSize + 1);

    int len_dis = mainRef_point.dis.size();
    int split = ceil((double)len_dis / Constants::NUM_CIRCLE);
    int num_pt = mainRef_point.iValuePts.size();

    double pred_pos = mainRef_point.coeffs[0];
    for(int i = 1; i <= Constants::COEFFS; ++i)
        pred_pos += mainRef_point.coeffs[i] * pow(dist_qt_mainRef, i);
    int mainPred_pos = pred_pos;

    if(mainPred_pos < 0)
        mainPred_pos = 0;
    if(mainPred_pos >= len_dis)
        mainPred_pos = len_dis - 1;

    int mainAcc_pos = 0;
    if(fabs(mainRef_point.dis[mainPred_pos]-dist_qt_mainRef)<FLOAT_ERROR)
        mainAcc_pos = mainPred_pos;
    else if(mainRef_point.dis[mainPred_pos] > dist_qt_mainRef)
        mainAcc_pos = DoublingSearch(mainRef_point.dis, true, dist_qt_mainRef,mainPred_pos);
    else
        mainAcc_pos = DoublingSearch(mainRef_point.dis, false, dist_qt_mainRef, mainPred_pos);

    if(fabs(mainRef_point.dis[mainAcc_pos]-dist_qt_mainRef)>FLOAT_ERROR){
        return;
    }else{
        rank_list[0] = mainAcc_pos / split;
    }

    // config pos in other reference points
    for(unsigned i = 0; i < refPtSize; ++i){
        double dist_qt_Ref = CaculateEuclideanDis(queryPt, mainRef_point.ref_points[i].point);

        int len_Coeff = mainRef_point.ref_points[i].coeffs.size();
        double refPred_pos = mainRef_point.ref_points[i].coeffs[0];
        for(int j = 1; j <= len_Coeff; ++j)
            refPred_pos += mainRef_point.ref_points[i].coeffs[j] * pow(dist_qt_Ref, j);
        int ref_pos = refPred_pos;
        if(ref_pos < 0)
            ref_pos = 0;
        if(ref_pos >= len_dis)
            ref_pos = len_dis - 1;

        int refAcc_pos = 0;
        if(fabs(mainRef_point.ref_points[i].dis[ref_pos]-dist_qt_Ref)<FLOAT_ERROR)
            refAcc_pos = ref_pos;
        else if(mainRef_point.ref_points[i].dis[ref_pos] > dist_qt_Ref)
            refAcc_pos = DoublingSearch(mainRef_point.ref_points[i].dis, true, dist_qt_Ref, ref_pos);
        else
            refAcc_pos = DoublingSearch(mainRef_point.ref_points[i].dis, false, dist_qt_Ref, ref_pos);

        if(fabs(mainRef_point.ref_points[i].dis[refAcc_pos]-dist_qt_Ref)>FLOAT_ERROR){
            return;
        }else{
            rank_list[i+1] = refAcc_pos / split;
        }   
        
    }
    int n = 1;
    while(n <= Constants::NUM_CIRCLE - 1)
        n *= 10;

    unsigned long long iValue_qrPt = rank_list[0];
    for(unsigned i = 1; i <= refPtSize; i++){
        iValue_qrPt = iValue_qrPt * n + rank_list[i];
    }
    

    double num_pow = pow(10, refPtSize + 1);
    double X_qrPt = iValue_qrPt / num_pow;
    int caculY_qrPt = (int)(mainRef_point.a * X_qrPt + mainRef_point.b);  

    int low = caculY_qrPt - mainRef_point.err_min < 0 ? 0 : caculY_qrPt - mainRef_point.err_min;
    low = low > num_pt - 1 ? num_pt - 1 : low;
    int high = caculY_qrPt + mainRef_point.err_max < num_pt ? caculY_qrPt + mainRef_point.err_max : num_pt - 1;
    int middle = 0;

    while(low <= high)
    {
        middle = (high - low) / 2 + low;
        if(iValue_qrPt < mainRef_point.iValuePts[middle].i_value)
            high = middle - 1;
        else if(iValue_qrPt > mainRef_point.iValuePts[middle].i_value)
            low = middle + 1;
        else  //mid=target
        {
            if((middle > 0 && mainRef_point.iValuePts[middle - 1].i_value < iValue_qrPt) ||(middle == 0))//找到开始位置
            {
                break;
            }
            else
                high = middle - 1;
        }
    }


    if(mainRef_point.iValuePts[middle].i_value == iValue_qrPt)//find point
    {
        int start_page = middle / Constants::PAGE_SIZE;
        int start_pos = middle;
        while(mainRef_point.iValuePts[++middle].i_value == iValue_qrPt);
        middle--;
        int end_page = middle / Constants::PAGE_SIZE;
        int end_pos = middle;
        ifstream fin;
        for(int s = start_page; s <= end_page; ++s){
            ++page;
            vector<vector<double> > pt_page;
            string fileName = "./data/cluster_" + to_string(clu_id) + "_" + to_string(s);
            fin.open(fileName, ios::binary);
            int size_page;
            fin.read((char*)&size_page,4);
            pt_page.reserve(size_page);
            int size_dim;
            fin.read((char*)&size_dim,4);
            for(int q = 0 ; q < size_page; q++){
                vector<double> coordinate;
                coordinate.resize(size_dim);
                fin.read(reinterpret_cast<char*>(coordinate.data()),size_dim * sizeof(coordinate.front()));
                pt_page.push_back(coordinate);
            }
            fin.close();
            if(start_page == end_page){
                for(int x = start_pos % Constants::PAGE_SIZE; x <= end_pos % Constants::PAGE_SIZE; ++x){
                    bool flag = false;
                    for(int i = 0; i < size_dim; i++){
                        if(fabs(queryPt.coordinate[i] - pt_page[x][i]) < FLOAT_ERROR)
                            continue;
                        else{
                            flag = true;
                            break;
                        }
                    }
                    if(!flag){
                        query_res = mainRef_point.iValuePts[s*Constants::PAGE_SIZE+x].id;
                        return;
                    }
                }
            }else if(s == start_page){
                for(int x = start_pos % Constants::PAGE_SIZE; x < Constants::PAGE_SIZE; ++x){
                    bool flag = false;
                    for(int i = 0; i < size_dim; i++){
                        if(fabs(queryPt.coordinate[i] - pt_page[x][i]) < FLOAT_ERROR)
                            continue;
                        else{
                            flag = true;
                            break;
                        }
                    }
                    if(!flag){
                        query_res = mainRef_point.iValuePts[s*Constants::PAGE_SIZE+x].id;
                        return;
                    }                
                }
            }else if(s == end_page){
                for(int x = 0; x <= end_pos % Constants::PAGE_SIZE; ++x){
                    bool flag = false;
                    for(int i = 0; i < size_dim; i++){
                        if(fabs(queryPt.coordinate[i] - pt_page[x][i]) < FLOAT_ERROR)
                            continue;
                        else{
                            flag = true;
                            break;
                        }
                    }
                    if(!flag){
                        query_res = mainRef_point.iValuePts[s*Constants::PAGE_SIZE+x].id;
                        return;
                    }  
                }
            }else{
                for(int x = 0; x < Constants::PAGE_SIZE; ++x){
                    bool flag = false;
                    for(int i = 0; i < size_dim; i++){
                        if(fabs(queryPt.coordinate[i] - pt_page[x][i]) < FLOAT_ERROR)
                            continue;
                        else{
                            flag = true;
                            break;
                        }
                    }
                    if(!flag){
                        query_res = mainRef_point.iValuePts[s*Constants::PAGE_SIZE+x].id;
                        return;
                    } 
                }
            }
        }
    }
}

int PointQuery_IO::DoublingSearch(vector<double> &disList, bool flag, double target, int start){
    int len = disList.size();
    if(flag){
        if(start == 0)
            return 0;
        int k = 1;
        int end = start;
        while(disList[end] > target){
            start = end;
            end -= k;
            k = k << 1;
            if(end < 0){
                end = 0;
                break;
            }
        }
        if(fabs(disList[end] - target) < FLOAT_ERROR)
            return end;
        if(start == end && start == 0)
            return 0;
        // vector<double> search_range(disList.begin() + end, disList.begin() + start);
        return BinarySearchForDoubling(disList, target, end, start);
    }else{
        if(start == len - 1)
            return len - 1;
        int k = 1;
        int end = start;
        while(disList[end] < target){
            start = end;
            end += k;
            k = k << 1;
            if(end > len - 1){
                end = len - 1;
                break;
            }
        }
        if(fabs(disList[end] - target) < FLOAT_ERROR)
            return end;
        if(start == end && start == len - 1)
            return len - 1;
        // vector<double> search_range(disList.begin() + start, disList.begin() + end);
        return BinarySearchForDoubling(disList, target, start, end);
    }
}

int PointQuery_IO::BinarySearchForDoubling(vector<double> &disList, double target, int low, int high){
    int middle = 0;
    while(low <= high) {
		middle = (high - low) / 2 + low;
        if(fabs(target-disList[middle])<FLOAT_ERROR){
            if((middle > 0 && disList[middle - 1] < target) || (middle == 0))//找到开始位置
            {
                break;
            }
            else
                high = middle - 1;
		}else if(target > disList[middle]) {
			low = middle + 1;
		}else{
            high = middle - 1;
		}
	}
    return middle;
}

double PointQuery_IO::CaculateEuclideanDis(Point &point_a, Point &point_b){
    double total = 0.0;
    unsigned dim = point_a.coordinate.size();
    for(unsigned i = 0; i < dim; ++i){
        total += pow(point_a.coordinate[i] - point_b.coordinate[i], 2);
    }
    return sqrt(total);
}

#endif