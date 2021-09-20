/*
* ChooseRef.h
*
* Author: ty_yan
*/
#ifndef ChooseRef_H
#define ChooseRef_H

#include "../entities/Reference.h"
#include "../entities/Point.h"
#include "../model/PolynomialRegression.h"
#include "../common/Constants.h"

#include <limits.h>
#include <math.h>
#include <algorithm>

using namespace std;

// caculate reference of each cluster
class ChooseRef
{
private:
    // number of reference point
    unsigned num_ref;
    // dimentional of data point
    unsigned dim;
    // total number of data set
    int num_data;

    Point mainRefPoint;
    vector<Ref_Point> RefPoint_Set;
    // double r_mainRefPoint;

    vector<double> mainRefDisArr;
public:
    mainRef_Point main_pointSet;
    vector<double> Pos;
    ChooseRef();
    ChooseRef(unsigned, Clu_Point &, Point &,vector<Point> &,int);


    // Point ChooseMainRefPoint(Clu_Point &);
    vector<Ref_Point> ChooseRefPoint(Clu_Point &);
    vector<Ref_Point> ChooseRefPoint_Input(Clu_Point &,vector<Point> &);
    double CaculateR(Clu_Point &, Point &);

    double CaculateEuclideanDis(Point &, Point &);

    vector<double> CaculateDisArr(Clu_Point &, Point &);

    vector<vector<double> > CaculateCircleBound(vector<double> &);

    mainRef_Point getMainRefPoint();

    vector<double> TrainModel(vector<double> &, vector<double> &);
    double CaculR2(vector<double> &,vector<double> &);
};

ChooseRef::ChooseRef()
{
}

ChooseRef::ChooseRef(unsigned num_ref, Clu_Point &cluster, Point &main_pivot,vector<Point> &oth_pivot,int type)
{
    this->num_ref = num_ref;
    this->dim = cluster.clu_point[0].coordinate.size();
    this->num_data = cluster.clu_point.size();

    for(int i = 0; i < num_data; i++)
        Pos.push_back(i);

    // this->mainRefPoint = ChooseMainRefPoint(cluster);
    this->mainRefPoint = main_pivot;
    // corner
    if(type == 0)
        this->RefPoint_Set = ChooseRefPoint(cluster);
    else
        // by input
        this->RefPoint_Set = ChooseRefPoint_Input(cluster,oth_pivot);
    this->mainRefDisArr = CaculateDisArr(cluster, mainRefPoint);
    // this->r_mainRefPoint = CaculateR(cluster, mainRefPoint);
    this->main_pointSet = mainRef_Point(mainRefPoint, mainRefDisArr[num_data - 1], mainRefDisArr[0], RefPoint_Set);
    this->main_pointSet.setMainRefDisArr(mainRefDisArr);
    vector<double> mainRefPt_coeffs = TrainModel(main_pointSet.dis, Pos);
    this->main_pointSet.setCoeffs(mainRefPt_coeffs);

}




vector<Ref_Point> ChooseRef::ChooseRefPoint(Clu_Point &cluster){
    vector<Ref_Point> ref_points;
    for(unsigned i = 0; i < num_ref; i++){
        double max = cluster.clu_point[0].coordinate[i];
        Point assRef_point = cluster.clu_point[0];
        for(int j = 1; j < num_data; j++){
            if(cluster.clu_point[j].coordinate[i] > max){
                max = cluster.clu_point[j].coordinate[i];
                assRef_point = cluster.clu_point[j];
            }
        }

        // double r = CaculateR(cluster, assRef_point);
        vector<double> dis = CaculateDisArr(cluster, assRef_point);
        Ref_Point ref_point = Ref_Point(assRef_point, dis[num_data - 1], dis[0]);

        ref_point.setDisArr(dis);
        // ref_point.setDict(CaculateCircleBound(ref_point.dis));

        vector<double> coeffs = TrainModel(ref_point.dis, Pos);
        ref_point.setCoeffs(coeffs);
        // cout << "ref " << i << " RMSE : " << CaculR2(ref_point.dis, coeffs) << ",";

        ref_points.push_back(ref_point);       
    }
    return ref_points;
}


vector<Ref_Point> ChooseRef::ChooseRefPoint_Input(Clu_Point &cluster,vector<Point> &oth_pivot){
    vector<Ref_Point> ref_points;
    for(unsigned i = 0; i < num_ref; i++){

        Point assRef_point = oth_pivot[i];


        // double r = CaculateR(cluster, assRef_point);
        vector<double> dis = CaculateDisArr(cluster, assRef_point);
        Ref_Point ref_point = Ref_Point(assRef_point, dis[num_data - 1], dis[0]);

        ref_point.setDisArr(dis);
        // ref_point.setDict(CaculateCircleBound(ref_point.dis));

        vector<double> coeffs = TrainModel(ref_point.dis, Pos);
        ref_point.setCoeffs(coeffs);
        // cout << "ref " << i << " RMSE : " << CaculR2(ref_point.dis, coeffs) << ",";

        ref_points.push_back(ref_point);       
    }
    return ref_points;
}

vector<double> ChooseRef::CaculateDisArr(Clu_Point &cluster, Point &refPoint){
    vector<double> dis;
    for(int i = 0; i < num_data; i++){
        double distance = CaculateEuclideanDis(refPoint, cluster.clu_point[i]);
        // Keep two decimal places
        // double a = (int)(distance * 100);
        // distance = a / 100;

        dis.push_back(distance);
    }

    /* sort distance vector
    * PlanA we get 1-d value by rank of distance list so we need to sort distance list
    * PlanB we get 1-d value by the point of which circle it is so there is no need for sort this array
    */
    // 2021/5/15 update : dist arr need to be sorted cuz we get rank of circle by (pos of arr)/(arr.len / num_cir)
    sort(dis.begin(), dis.end());
    return dis;
}

vector<vector<double> > ChooseRef::CaculateCircleBound(vector<double> &dis){
    /* 5/11
    * we use model to caculate number of circle of a distance between query point & reference point
    * so here we sort distance list and devide list into circle_num(50) pieces 
    * and record lower & upper of each piece and piece number
    */
    sort(dis.begin(), dis.end());
    int circle_num = 30;
    int split = dis.size() / circle_num + 1;
    vector<vector<double> > dict_circle;
    for(int i = 0; i < circle_num; i++){
        vector<double> circle_bound;
        // lower bound
        circle_bound.push_back(dis[i*split]);
        // upper bound
        if(i < circle_num - 1){
            circle_bound.push_back(dis[(i+1)*split-1]);
        }else{
            circle_bound.push_back(dis[dis.size()-1]);
        }
        dict_circle.push_back(circle_bound);
    }
    return dict_circle;
}

double ChooseRef::CaculateR(Clu_Point &cluster, Point &point){
    double r = 0.0;

    for(int i = 0; i < num_data; i++){
        double dis = CaculateEuclideanDis(point, cluster.clu_point[i]);
        r = r > dis ? r : dis;
    }
    return r;
}  

double ChooseRef::CaculateEuclideanDis(Point &point_a, Point &point_b){
    double total = 0.0;
    for(unsigned i = 0; i < dim; i++){
        total += pow(point_a.coordinate[i] - point_b.coordinate[i], 2);
    }
    return sqrt(total);
}

mainRef_Point ChooseRef::getMainRefPoint(){
    // main_pointSet.setDict(CaculateCircleBound(mainRefDisArr));
    return main_pointSet;
}

vector<double> ChooseRef::TrainModel(vector<double> &X, vector<double> &Y){
    PolynomialRegression<double> polyreg;
    vector<double> coeffs;
    polyreg.fitIt(X, Y, Constants::COEFFS, coeffs);
    return coeffs;
}

double ChooseRef::CaculR2(vector<double> &X, vector<double> &coeff){
    double SSE = 0.0;
    double SST = 0.0;
    double avg = (0 + num_data - 1) / 2;
    for(int i = 0; i < num_data; ++i){
        double pre = coeff[0];
        for(int j = 1; j <= Constants::COEFFS; ++j)
            pre += coeff[j] * pow(X[i], j);
        SSE += pow(pre - i, 2);
        SST += pow(i - avg, 2);
    }
    // return 1 - SSE/SST;
    return sqrt(SSE/num_data);
}
#endif