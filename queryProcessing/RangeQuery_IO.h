/*
* RangeQuery.h
*
* Author: ty_yan
*/

#ifndef RangeQuery_IO_H
#define RangeQuery_IO_H

#include "../entities/Point.h"
#include "../entities/Reference.h"
#include "../common/CalCirclePos.h"
#include "../common/Constants.h"

#include <vector>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>

using namespace std;

class RangeQuery_IO
{
public:
    // unsigned digits;
    unsigned N;
    Point queryPoint;
    double radius;
    double num_pow;
    unsigned clu_id;
    int num_pt;
    int last_page;

    RangeQuery_IO();
    RangeQuery_IO(Point &, double, mainRef_Point &, vector<int> &, CalCirclePos &, vector<CalCirclePos> &, unsigned, int &,double &);
    int DoublingSearch(vector<double> &, bool, double, int);
    // int BinarySearch(vector<double> &, double);
    int BinarySearchForDoubling(vector<double> &, double, int, int);
    void DoRangeQuery(unsigned, vector<vector<int> > &, unsigned long long, mainRef_Point &,vector<int> &, int &,double &);

    double CaculateEuclideanDis(vector<double> &, Point &);
};

RangeQuery_IO::RangeQuery_IO()
{
}

RangeQuery_IO::RangeQuery_IO(Point &queryPoint, double radius, mainRef_Point &mainRef_point, vector<int> &rangeQueryRes, CalCirclePos &mainRefPtCircle, vector<CalCirclePos> &ref_query, unsigned clu_id, int &page,double &io_time)
{
    // chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    // find out query rank range
    // store as vector<vector<unsigned> >
    // TODO BinarySearch -> learned model
    this->clu_id = clu_id;
    vector<vector<int> > rank_list;
    vector<int> mainRankRange;

    this->queryPoint = queryPoint;
    this->radius = radius;
    this->num_pow = pow(10, mainRef_point.ref_points.size() + 1);
    this->num_pt = mainRef_point.iValuePts.size();

    // Plan B
    int len_dis = mainRef_point.dis.size();
    int split = ceil((double)len_dis / Constants::NUM_CIRCLE);

    unsigned refPtSize = mainRef_point.ref_points.size();

    rank_list.resize(refPtSize + 1);
    mainRankRange.resize(2);

    // use polynomial regression model
    // 2021/5/18
    double mainRefPre_disLower = mainRef_point.coeffs[0];
    double mainRefPre_disUpper = mainRef_point.coeffs[0];
    for(int i = 1; i <= Constants::COEFFS; ++i){
        mainRefPre_disLower += mainRef_point.coeffs[i] * pow(mainRefPtCircle.dis_lower, i);
        mainRefPre_disUpper += mainRef_point.coeffs[i] * pow(mainRefPtCircle.dis_upper, i);
    }
    // 取整
    int mainRefPrePos_lower = mainRefPre_disLower;
    int mainRefPrePos_upper = mainRefPre_disUpper;

    if(mainRefPrePos_lower < 0)
        mainRefPrePos_lower = 0;
    if(mainRefPrePos_lower >= len_dis)
        mainRefPrePos_lower = len_dis - 1;

    int pred_lower = 0;

    if(mainRef_point.dis[mainRefPrePos_lower] == mainRefPtCircle.dis_lower){
        pred_lower = mainRefPrePos_lower;
    }else if(mainRef_point.dis[mainRefPrePos_lower] > mainRefPtCircle.dis_lower){
        pred_lower = DoublingSearch(mainRef_point.dis, true, mainRefPtCircle.dis_lower, mainRefPrePos_lower);
    }else{
        pred_lower = DoublingSearch(mainRef_point.dis, false, mainRefPtCircle.dis_lower, mainRefPrePos_lower);
    }

    if(mainRef_point.dis[pred_lower] > mainRefPtCircle.dis_lower)
        mainRankRange[0] = pred_lower / split;
    else
        mainRankRange[0] = (pred_lower + 1) / split;

    if(mainRefPrePos_upper < 0)
        mainRefPrePos_upper = 0;
    if(mainRefPrePos_upper >= len_dis)
        mainRefPrePos_upper = len_dis - 1;

    int pred_upper = 0;
    
    if(mainRef_point.dis[mainRefPrePos_upper] == mainRefPtCircle.dis_upper){
        pred_upper = mainRefPrePos_upper;
    }else if(mainRef_point.dis[mainRefPrePos_upper] > mainRefPtCircle.dis_upper){
        pred_upper = DoublingSearch(mainRef_point.dis, true, mainRefPtCircle.dis_upper, mainRefPrePos_upper);
    }else{
        pred_upper = DoublingSearch(mainRef_point.dis, false, mainRefPtCircle.dis_upper, mainRefPrePos_upper);
    }

    if(mainRef_point.dis[pred_upper] > mainRefPtCircle.dis_upper)
        mainRankRange[1] = (pred_upper - 1) / split;
    else
        mainRankRange[1] = pred_upper / split;

    rank_list[0] = mainRankRange;

    for(unsigned i = 0; i < refPtSize; ++i){

        vector<int> refRankRange;
        refRankRange.resize(2);

        // 2021/5/18
        double refPre_disLower = mainRef_point.ref_points[i].coeffs[0];
        double refPre_disUpper = mainRef_point.ref_points[i].coeffs[0];
        for(int j = 1; j <= Constants::COEFFS; ++j){
            refPre_disLower += mainRef_point.ref_points[i].coeffs[j] * pow(ref_query[i].dis_lower, j);
            refPre_disUpper += mainRef_point.ref_points[i].coeffs[j] * pow(ref_query[i].dis_upper, j);
        }
        // 取整
        int refPrePos_lower = refPre_disLower;
        int refPrePos_upper = refPre_disUpper;

        if(refPrePos_lower < 0)
            refPrePos_lower = 0;
        if(refPrePos_lower >= len_dis)
            refPrePos_lower = len_dis - 1;

        int refPred_lower = 0;

        if(mainRef_point.ref_points[i].dis[refPrePos_lower] == ref_query[i].dis_lower){
            refPred_lower = refPrePos_lower;
        }else if(mainRef_point.ref_points[i].dis[refPrePos_lower] > ref_query[i].dis_lower){
            refPred_lower = DoublingSearch(mainRef_point.ref_points[i].dis, true, ref_query[i].dis_lower, refPrePos_lower);
        }else{
            refPred_lower = DoublingSearch(mainRef_point.ref_points[i].dis, false, ref_query[i].dis_lower, refPrePos_lower);
        }

        if(mainRef_point.ref_points[i].dis[refPred_lower] > ref_query[i].dis_lower)
            refRankRange[0] = refPred_lower / split;
        else
            refRankRange[0] = (refPred_lower + 1) / split;

        if(refPrePos_upper < 0)
            refPrePos_upper = 0;
        if(refPrePos_upper >= len_dis)
            refPrePos_upper = len_dis - 1;

        int refPred_upper = 0;
        
        if(mainRef_point.ref_points[i].dis[refPrePos_upper] == ref_query[i].dis_upper){
            refPred_upper = refPrePos_upper;
        }else if(mainRef_point.ref_points[i].dis[refPrePos_upper] > ref_query[i].dis_upper){
            refPred_upper = DoublingSearch(mainRef_point.ref_points[i].dis, true, ref_query[i].dis_upper, refPrePos_upper);
        }else{
            refPred_upper = DoublingSearch(mainRef_point.ref_points[i].dis, false, ref_query[i].dis_upper, refPrePos_upper);
        }

        if(mainRef_point.ref_points[i].dis[refPred_upper] > ref_query[i].dis_upper)
            refRankRange[1] = (refPred_upper - 1) / split;
        else
            refRankRange[1] = refPred_upper / split;
        
        rank_list[i + 1] = refRankRange;
    }

    int n = 1;
    while(n <= Constants::NUM_CIRCLE - 1)
        n *= 10;
    this->N = n;

    unsigned i = 0;
	unsigned long long queryRange = 0;

    this->last_page = -1;




    DoRangeQuery(i, rank_list, queryRange, mainRef_point, rangeQueryRes, page,io_time);


}

/* DoublingSearch
* param 1 : list
* param 2 : search forward(true) / search backward(false)
* param 3 : target value
*/
int RangeQuery_IO::DoublingSearch(vector<double> &disList, bool flag, double target, int start){
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
        if(disList[end] == target)
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
        if(disList[end] == target)
            return end;
        if(start == end && start == len - 1)
            return len - 1;
        // vector<double> search_range(disList.begin() + start, disList.begin() + end);
        return BinarySearchForDoubling(disList, target, start, end);
    }
}

int RangeQuery_IO::BinarySearchForDoubling(vector<double> &disList, double target, int low, int high){
    int middle = 0;
    while(low <= high) {
		middle = (high - low) / 2 + low;
		if(target < disList[middle]) {
			high = middle - 1;
		}else if(target > disList[middle]) {
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

void RangeQuery_IO::DoRangeQuery(unsigned i, vector<vector<int> > &rank_list, unsigned long long queryRange, mainRef_Point &mainRef_point, vector<int> &rangeQueryRes, int &page,double &io_time){
    if (i == rank_list.size() - 1) {
		unsigned long long query_lower = queryRange * N + rank_list[i][0];
		unsigned long long query_upper = queryRange * N + rank_list[i][1];

        double x_lower = query_lower / num_pow;
        double x_upper = query_upper / num_pow;
        int CaculY_lower = (int)(mainRef_point.a * x_lower + mainRef_point.b);
        int CaculY_upper = (int)(mainRef_point.a * x_upper + mainRef_point.b);
        unsigned long long find_queryLower;
        int pos_lower;
        unsigned long long find_queryUpper;
        int pos_upper;

        // Binary search for the first one which >= query_lower
        int low = CaculY_lower - mainRef_point.err_min < 0 ? 0 : CaculY_lower - mainRef_point.err_min;
        low = low > num_pt - 1 ? num_pt - 1 : low;
        int high = CaculY_lower + mainRef_point.err_max < num_pt ? CaculY_lower + mainRef_point.err_max : num_pt - 1;
        int middle = 0;


        while(low <= high)
        {
            middle = (high - low) / 2 + low;
            if(query_lower < mainRef_point.iValuePts[middle].i_value)
                high = middle - 1;
            else if(query_lower > mainRef_point.iValuePts[middle].i_value)
                low = middle + 1;
            else  //mid=target
            {
                if((middle > 0 && mainRef_point.iValuePts[middle - 1].i_value < query_lower) ||(middle == 0))//找到开始位置
                {
                    break;
                }
                else
                    high = middle - 1;
            }
        }

        if(middle < low && low < num_pt){
            pos_lower = low;
            find_queryLower = mainRef_point.iValuePts[low].i_value;
        }else{
            pos_lower = middle;
            find_queryLower = mainRef_point.iValuePts[middle].i_value;
        }


        // Binary search for the last one which <= query_upper
        low = CaculY_upper - mainRef_point.err_min < 0 ? 0 : CaculY_upper - mainRef_point.err_min;
        low = low > num_pt - 1 ? num_pt - 1 : low;
        high = CaculY_upper + mainRef_point.err_max < num_pt ? CaculY_upper + mainRef_point.err_max : num_pt - 1;
        middle = 0;

        while(low <= high)
        {
            middle = (high - low) / 2 + low;
            if(mainRef_point.iValuePts[middle].i_value > query_upper)
                high = middle - 1;
            else if(mainRef_point.iValuePts[middle].i_value < query_upper)
                low = middle + 1;
            else  //mid=target
            {
                if((middle < num_pt - 1 && mainRef_point.iValuePts[middle+1].i_value > query_upper) || (middle == num_pt - 1))//找到结束位置
                {
                    break;
                }
                else
                    low = middle + 1;
            }
        }
        if(query_upper == mainRef_point.iValuePts[middle].i_value){
            pos_upper = middle;
            find_queryUpper = query_upper;
        }else if(middle < low){
            pos_upper = low - 1;
            find_queryUpper = mainRef_point.iValuePts[low - 1].i_value;
        }else{
            if(middle == 0){
                pos_upper = 0;
                find_queryUpper = mainRef_point.iValuePts[0].i_value;
            }else{
                pos_upper = middle - 1;
                find_queryUpper = mainRef_point.iValuePts[middle - 1].i_value;
            }
        }


        if(find_queryLower>=find_queryUpper){
            if(find_queryLower == find_queryUpper && find_queryLower >= query_lower && find_queryUpper <= query_upper){
                // for(int x = pos_lower; x <= pos_upper; ++x){
                //     if(CaculateEuclideanDis(mainRef_point.iValuePts[x], queryPoint) <= radius)
                //         rangeQueryRes.push_back(mainRef_point.iValuePts[x]);
                // }  
                ifstream fin;
                int start_page = pos_lower / Constants::PAGE_SIZE;
                int end_page = pos_upper / Constants::PAGE_SIZE;

                // cout << pos_lower << " " << pos_upper << endl;
                // cout << mainRef_point.iValuePts[pos_lower].i_value << " " << mainRef_point.iValuePts[pos_upper].i_value << endl;

                for(int s = start_page; s <= end_page; ++s){
                    if(s == last_page)
                        continue;
                    ++page;
                    vector<vector<double> > pt_page;
                    chrono::steady_clock::time_point start = chrono::steady_clock::now();
                    string fileName = "./data/cluster_" + to_string(clu_id) + "_" + to_string(s);
                    // cout << fileName << endl;
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
                    chrono::steady_clock::time_point end = chrono::steady_clock::now();
                    io_time+=chrono::duration_cast<chrono::microseconds>(end - start).count();
                    for(int x = 0; x < size_page; ++x){
                        if(CaculateEuclideanDis(pt_page[x], queryPoint) <= radius){
                            rangeQueryRes.push_back(mainRef_point.iValuePts[s*Constants::PAGE_SIZE+x].id);
                        }
                    }       
                }
                this->last_page = end_page;
            } 
        }else{
            // for(int x = pos_lower; x <= pos_upper; ++x){
            //     if(CaculateEuclideanDis(mainRef_point.iValuePts[x], queryPoint) <= radius)
            //         rangeQueryRes.push_back(mainRef_point.iValuePts[x]);
            // }
            ifstream fin;
            int start_page = pos_lower / Constants::PAGE_SIZE;
            int end_page = pos_upper / Constants::PAGE_SIZE;

            // cout << pos_lower << " " << pos_upper << endl;
            // cout << mainRef_point.iValuePts[pos_lower].i_value << " " << mainRef_point.iValuePts[pos_upper].i_value << endl;



            for(int s = start_page; s <= end_page; ++s){
                if(s == last_page)
                    continue;
                ++page;
                chrono::steady_clock::time_point start = chrono::steady_clock::now();
                vector<vector<double> > pt_page;
                string fileName = "./data/cluster_" + to_string(clu_id) + "_" + to_string(s);
                // cout << fileName << endl;
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
                chrono::steady_clock::time_point end = chrono::steady_clock::now();
                io_time+=chrono::duration_cast<chrono::microseconds>(end - start).count();
                for(int x = 0; x < size_page; ++x){
                    if(CaculateEuclideanDis(pt_page[x], queryPoint) <= radius){
                        rangeQueryRes.push_back(mainRef_point.iValuePts[s*Constants::PAGE_SIZE+x].id);
                    }
                }       
            }
            this->last_page = end_page;
        }

	}
	else {
		for (int j = rank_list[i][0]; j <= rank_list[i][1]; ++j) {
			queryRange = queryRange * N + j;

			DoRangeQuery(++i, rank_list, queryRange, mainRef_point, rangeQueryRes, page,io_time);
			i--;

			queryRange = (int)queryRange / N;
		}
	}
}


double RangeQuery_IO::CaculateEuclideanDis(vector<double> &coordinate_a, Point &point_b){
    double total = 0.0;
    unsigned dim = coordinate_a.size();
    for(unsigned i = 0; i < dim; ++i){
        total += pow(coordinate_a[i] - point_b.coordinate[i], 2);
    }
    return sqrt(total);
}

#endif