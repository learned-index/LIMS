/*
* KNN.h
* Author: ty_yan
*/

#ifndef KNN_H
#define KNN_H

#include "../entities/Point.h"
#include "../entities/Reference.h"
#include "../common/CalCirclePos.h"
#include "../common/Constants.h"

#include <queue>

using namespace std;

class KNN{
public:
    Point queryPoint;
    double num_pow;
    unsigned N;
    unsigned k;
    unsigned clu_id;
    int num_pt;

    KNN();
    KNN(Point &, mainRef_Point &, priority_queue<pair<double, unsigned int> > &,CalCirclePos &, vector<CalCirclePos> &, unsigned, int8_t (*)[100000], unsigned);

    double CaculateEuclideanDis(Point &, Point &);
    int DoublingSearch(vector<double> &, bool, double, int);
    int BinarySearch(vector<double> &, double, int, int);
    void DoRangeQueryForKNN(unsigned, vector<vector<int> > &, unsigned long long, mainRef_Point &, priority_queue<pair<double, unsigned int> > &, int8_t (*)[100000]);
};

KNN::KNN(){}

KNN::KNN(Point &KNN_pt, mainRef_Point &mainRef_point, priority_queue<pair<double, unsigned int> > &KNNRes_queue, CalCirclePos &mainRefPtCircle, vector<CalCirclePos> &ref_query, unsigned k, int8_t (*arr)[100000], unsigned clu_id)
{
    unsigned refPtSize = mainRef_point.ref_points.size();
    this->queryPoint = KNN_pt;
    this->num_pow = pow(10, refPtSize + 1);
    this->k = k;
    this->clu_id = clu_id;
    this->num_pt = mainRef_point.iValuePts.size();
    int n = 1;
    while(n <= Constants::NUM_CIRCLE - 1)
        n *= 10;
    this->N = n;

    int len_dis = mainRef_point.dis.size();
    int split = ceil((double)len_dis / Constants::NUM_CIRCLE);

    vector<vector<int> > rank_list;
    vector<int> mainRankRange;

    rank_list.resize(refPtSize + 1);
    mainRankRange.resize(2);

    double mainRefPre_disLower = mainRef_point.coeffs[0];
    double mainRefPre_disUpper = mainRef_point.coeffs[0];
    for(int i = 1; i <= Constants::COEFFS; ++i){
        mainRefPre_disLower += mainRef_point.coeffs[i] * pow(mainRefPtCircle.dis_lower, i);
        mainRefPre_disUpper += mainRef_point.coeffs[i] * pow(mainRefPtCircle.dis_upper, i);
    }

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

    if(mainRef_point.dis[pred_lower] >= mainRefPtCircle.dis_lower)
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

    // 计算其余参考点的多项式函数及环号预测
    for(unsigned i = 0; i < refPtSize; ++i){
        vector<int> refRankRange;
        refRankRange.resize(2);

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

        if(mainRef_point.ref_points[i].dis[refPred_lower] >= ref_query[i].dis_lower)
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

    // 拼接一维值，计算每个区间
    unsigned i = 0;
	unsigned long long queryRange = 0;

    DoRangeQueryForKNN(i, rank_list, queryRange, mainRef_point, KNNRes_queue, arr);
}

void KNN::DoRangeQueryForKNN(unsigned i, vector<vector<int> > &rank_list, unsigned long long queryRange, mainRef_Point &mainRef_point, priority_queue<pair<double, unsigned int> > &KNNRes_queue, int8_t (*arr)[100000]){
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
                for(int x = pos_lower; x <= pos_upper; ++x){
                    if(arr[clu_id][x] == 0){
                        arr[clu_id][x] = 1;
                        double dis_pt_qryPt = CaculateEuclideanDis(mainRef_point.iValuePts[x], queryPoint);
                        if(KNNRes_queue.size() < k){
                            pair<double, unsigned int> elem(dis_pt_qryPt, mainRef_point.iValuePts[x].id);
                            KNNRes_queue.push(elem);
                        }else{
                            if(KNNRes_queue.top().first > dis_pt_qryPt){
                                KNNRes_queue.pop();
                                pair<double, unsigned int> elem(dis_pt_qryPt, mainRef_point.iValuePts[x].id);
                                KNNRes_queue.push(elem);
                            }
                        }
                    }
                }   
            } 
        }else{
            for(int x = pos_lower; x <= pos_upper; ++x){
                if(arr[clu_id][x] == 0){
                    arr[clu_id][x] = 1;
                    double dis_pt_qryPt = CaculateEuclideanDis(mainRef_point.iValuePts[x], queryPoint);
                    if(KNNRes_queue.size() < k){
                        pair<double, unsigned int> elem(dis_pt_qryPt, mainRef_point.iValuePts[x].id);
                        KNNRes_queue.push(elem);
                    }else{
                        if(KNNRes_queue.top().first > dis_pt_qryPt){
                            KNNRes_queue.pop();
                            pair<double, unsigned int> elem(dis_pt_qryPt, mainRef_point.iValuePts[x].id);
                            KNNRes_queue.push(elem);
                        }
                    }
                }
            }
        }

	}
	else {
		for (int j = rank_list[i][0]; j <= rank_list[i][1]; ++j) {
			queryRange = queryRange * N + j;

			DoRangeQueryForKNN(++i, rank_list, queryRange, mainRef_point, KNNRes_queue, arr);
			i--;

			queryRange = (int)queryRange / N;
		}
	}
}

int KNN::DoublingSearch(vector<double> &disList, bool flag, double target, int start){
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
        return BinarySearch(disList, target, end, start);
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
        return BinarySearch(disList, target, start, end);
    }
}

int KNN::BinarySearch(vector<double> &disList, double target, int low, int high){
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

double KNN::CaculateEuclideanDis(Point &point_a, Point &point_b){
    double total = 0.0;
    unsigned dim = point_a.coordinate.size();
    for(unsigned i = 0; i < dim; ++i){
        total += pow(point_a.coordinate[i] - point_b.coordinate[i], 2);
    }
    return sqrt(total);
}

#endif