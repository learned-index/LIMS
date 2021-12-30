#include "fileIO/FileIO.h"
#include "caculateRef/ChooseRef.h"
#include "caculateRef/CaculateIValue.h"
#include "entities/Reference.h"
#include "entities/Point.h"
#include "queryProcessing/RangeQuery.h"
#include "queryProcessing/PointQuery.h"
#include "queryProcessing/RangeQuery_IO.h"
#include "queryProcessing/PointQuery_IO.h"
#include "queryProcessing/KNN.h"
#include "queryProcessing/KNN_IO.h"
#include "common/CalCirclePos.h"

#include "dis.h"
// #include "dis_string.h"
#include <sys/stat.h>
#include <sys/types.h>

#include <string>
#include <chrono>
#include <ctime>
#include <sstream>
#include <queue>
#include <algorithm>

#define random(a,b) (rand() % (b-a+1))+a
#define rand() ((rand()%10000)*10000+rand()%10000)

using namespace std;

double CaculateEuclideanDis(Point &point_a, Point &point_b){
    double total = 0.0;
    for(unsigned i = 0; i < point_a.coordinate.size(); i++){
        total += pow(point_a.coordinate[i] - point_b.coordinate[i], 2);
    }
    return sqrt(total);
}

double CaculateEuclideanDis2(vector<double> &point_a, Point &point_b){
    double total = 0.0;
    for(unsigned i = 0; i < point_a.size(); i++){
        total += pow(point_a[i] - point_b.coordinate[i], 2);
    }
    return sqrt(total);
}

vector<Point> LoadPointForQuery(string filename){
    ifstream fin;
    vector<Point> queryPoints;
    fin.open(filename.c_str());
    if(!fin){
        cout << filename << " file could not be opened\n";
        exit(0);
    }
    string line;
    while(getline(fin,line)){
        stringstream ss(line); 
        string value;
        vector<double> coordinate;
        while(getline(ss, value, ','))
            coordinate.push_back(atof(value.c_str()));

        Point point = Point(coordinate);
        queryPoints.push_back(point);
    }
    fin.close();
    return queryPoints;
}

static bool AscendingSort(const InsertPt& point_a, const InsertPt& point_b){
    return point_a.i_value < point_b.i_value;
}

void Insert(string filename,vector<mainRef_Point> &all_ref, unsigned num_clu){
    ifstream fin;
    vector<InsertPt> insertPt;
    fin.open(filename);
    if(!fin){
        cout << filename << " file could not be opened\n";
        exit(0);
    }
    string line;
    unsigned id = 0;
    while(getline(fin,line)){
        stringstream ss(line); 
        string value;
        vector<double> coordinate;
        while(getline(ss, value, ','))
            coordinate.push_back(atof(value.c_str()));

        InsertPt point = InsertPt(coordinate,id++);
        insertPt.push_back(point);
    }
    fin.close();


    int insert_size = insertPt.size();
    vector<vector<InsertPt> > all_insertPt;
    all_insertPt.resize(num_clu);
    for(int index = 0; index < insert_size; ++index){
        double dis = CaculateEuclideanDis2(insertPt[index].coordinate,all_ref[0].point);
        int clu_id = 0;
        for(unsigned i = 1; i < num_clu;++i){
            double dis_oth = CaculateEuclideanDis2(insertPt[index].coordinate,all_ref[i].point);
            if(dis_oth < dis){
                dis = dis_oth;
                clu_id = i;
            }
        }
        insertPt[index].setIValue(dis);
        all_insertPt[clu_id].push_back(insertPt[index]);
    }
    for(unsigned index = 0; index < num_clu; ++index){
        sort(all_insertPt[index].begin(),all_insertPt[index].end(),AscendingSort);
        all_ref[index].setInsertPt(all_insertPt[index]);


        ofstream fout;
        unsigned page_num = ceil((double)all_ref[index].insert_list.size() / Constants::PAGE_SIZE);
        unsigned count = 0;
        int size_dim = all_ref[index].insert_list[0].coordinate.size();
        for(unsigned i = 0; i < page_num; i++){
            string outputFileName = "./data/cluster_" + to_string(index) + "_insert_" + to_string(i);
            fout.open(outputFileName, ios::binary);
            int size_page = count + Constants::PAGE_SIZE > all_ref[index].insert_list.size() ? all_ref[index].insert_list.size()-count : Constants::PAGE_SIZE;
            fout.write((const char*)&size_page,4);
            fout.write((const char*)&size_dim,4);
            for(int j = 0; j < size_page; ++j)
                fout.write(reinterpret_cast<char *>(all_ref[index].insert_list[count+j].coordinate.data()),size_dim*sizeof(all_ref[index].insert_list[count+j].coordinate.front()));
            fout.close();
            count += Constants::PAGE_SIZE;

        }

    }
}

int BinarySearch(vector<InsertPt> &disList, double target, int low, int high){
    int middle = 0;
    while(low <= high) {
		middle = (high - low) / 2 + low;
		if(target < disList[middle].i_value) {
			high = middle - 1;
		}else if(target > disList[middle].i_value) {
			low = middle + 1;
		}else{
            if((middle > 0 && disList[middle - 1].i_value < target) || (middle == 0))//找到开始位置
            {
                break;
            }
            else
                high = middle - 1;
		}
	}
    return middle;
}





// cluster evaluate method
double dist[NUUM] ;
int fgg[NUUM];
int stack1[NUUM],top=0;

double sm(double a,double b){ 
	if(a>b) return b;
	return a;
}

void work(int dim){
	top=1;
	stack1[0]=rand()%num;
	fgg[stack1[0]]=1;
	for(int j=0;j<num;j++) dist[j]=dis(j,stack1[0],dim);
	for(int i=1;i<pnum;i++){
		double max=0;
		int k=-1;
		for(int j=0;j<num;j++)if(fgg[j]==0){
			if(dist[j]>max){
				max=dist[j];
				k=j;
			}
		}
		stack1[top++]=k;	
		fgg[k]=1;
		for(int j=0;j<num;j++)dist[j]=sm(dist[j],dis(j,k,dim));
	}
}

void clu(int dim,string path){
	vector<vector<int> > cluster;
	cluster.resize(pnum);
	for(int index = 0; index < num; ++index){
		double dist = dis(index,stack1[0],dim);
		int pivot_index = 0;
		for(int i = 1; i < pnum; ++i){
			double dis_o = dis(index,stack1[i],dim);
			if(dis_o < dist){
				dist = dis_o;
				pivot_index = i;
			}
		}
		cluster[pivot_index].push_back(index);
	}
	ofstream fout;
	fout.setf(ios::fixed);
	fout.precision(8);

	for(int i = 0; i < pnum; ++i){
		int len = cluster[i].size();
		// cout << len << endl; 
		if(len < 3){
			continue;
		}
		string filename = path + "/clu_" + to_string(i) + ".txt";
		fout.open(filename);
		for(int j = 0; j < len; ++j){
			fout << cluster[i][j];
			for(int d = 0; d < dim; ++d)
				fout << "," << loc[cluster[i][j]][d];
			fout << endl;
		}

		fout.close();
	}
}


void op(string filename,int dim){	
	ofstream out1;
    out1.open(filename);
	out1.setf(ios::fixed);
	out1.precision(8);


	for(int i=0;i<top;i++){
		out1<<loc[stack1[i]][0];	
		for(int d = 1; d < dim; ++d)
			out1<< "," << loc[stack1[i]][d];
		out1 << endl;	
	}	
	out1.close(); 
}


double CalErrAvg(vector<Clu_Point> &data, vector<vector<Point> > &ref_data, vector<double> &dist_max, int k, int num_ref, int num_alldata){
    double avg_err = 0;
    // double max_err = 0, min_err = 0x3f3f3f;
    for(int i = 0; i < k; ++i){
        int len = data[i].clu_point.size();
        double clu_err = 0;
        for(int r = 0; r < num_ref; ++r){
            double err = 0;
            vector<double> arr_dis;
            arr_dis.reserve(len);
            for(int l = 0; l < len; ++l)
                arr_dis.push_back(CaculateEuclideanDis(ref_data[i][r],data[i].clu_point[l]));
            sort(arr_dis.begin(),arr_dis.end());
            // cout << arr_dis[0] << " "<< arr_dis[len-1] <<  "       ";
            if(!r)
                dist_max[i] = arr_dis[len-1];

            double sum_xy = 0, sum_x = 0, sum_y = 0, sum_pow2x = 0;
            for(int y = 0; y < len; ++y){
                sum_x += arr_dis[y];
                sum_y += y;
                sum_pow2x += arr_dis[y] * arr_dis[y];
                sum_xy += arr_dis[y] * y;
            }
            double avg_x = sum_x / len;
            double avg_y = sum_y / len;
            double avg_xy = sum_xy / len;
            double avg_pow2x = sum_pow2x / len;
            double a = (avg_xy - avg_x * avg_y) / (avg_pow2x - avg_x * avg_x);
            double b = avg_y - a * avg_x;
            // cout << a << " " << b << " ";
            // double max_err = 0, min_err = 0x3f3f3f;
            // vector<int> err_arr;
            // err_arr.resize(len);
            for(int c = 0; c < len; ++c){
                int cal_y = (int)(a * arr_dis[c] + b);
                // err_arr[c] = abs(cal_y - c);
                // max_err = max_err > err_arr[c] ? max_err : err_arr[c];
                // min_err = min_err < err_arr[c] ? min_err : err_arr[c];
                err += abs(cal_y - c);
            }
            // cout << max_err << " " << min_err << endl;
            // for(int c = 0; c < len; ++c){
            //     err += (err_arr[c]-min_err) / (max_err - min_err);
            // }
            // cout << err / len << endl;
            
            // clu_err = clu_err + err / len;
            clu_err = clu_err + err;
        }
        // clu_err /= num_ref;
        // max_err = max_err > clu_err ? max_err : clu_err;
        // min_err = min_err < clu_err ? min_err : clu_err;
        // avg_err += clu_err * len / num_alldata;
        avg_err += clu_err;
    }
    // avg_err /= k;
    // cout << min_err << " " << max_err  << " " << avg_err << endl;
    // return (avg_err - min_err)/(max_err- min_err);
    // return avg_err;
    
    return avg_err / (num_ref * num_alldata);
}


double CalSihCoe(vector<Clu_Point> &data, int k){
    double s = 0;
    double num_data = 0;
    for(int i = 0; i < k; ++i){
        // cout << i << " ";
        int len = data[i].clu_point.size();
        num_data += len;
        
        for(int j = 0; j < len; ++j){
            double a = 0, b = 0x3f3f3f;


            for(int k = 0; k < len; ++k)
                a += CaculateEuclideanDis(data[i].clu_point[j],data[i].clu_point[k]);

            a = a / (len - 1);
            

            for(int l = 0; l < k; ++l){
                if(l == i)
                    continue;
                double bi = 0;
                int oth_len = data[l].clu_point.size();
                for(int o = 0; o < oth_len; ++o)
                    bi += CaculateEuclideanDis(data[i].clu_point[j],data[l].clu_point[o]);
                bi /= oth_len;
                b = b < bi ? b : bi;
            }
            s = s + (b-a) / (a < b ? b : a);
            
        }
    }
    return s/num_data; 
}

double CalOverlapRate(vector<vector<Point> > &ref_data, vector<double> &dist_max, int k){
    double overlap_rate = 0;
    for(int i = 0; i < k; ++i){
        double rate = 0;
        for(int j = 0; j < k; ++j){
            if(i == j)
                continue;
            double dist_cl1_cl2 = CaculateEuclideanDis(ref_data[i][0], ref_data[j][0]);

            if(dist_cl1_cl2 > dist_max[i] + dist_max[j])
                rate += 0;
            else
                rate = rate + (min(dist_cl1_cl2 + dist_max[j], dist_max[i]) - max(dist_cl1_cl2 - dist_max[j], 0.0))/dist_max[i];
            
        }
        overlap_rate = overlap_rate + rate / (k-1);
    }
    return overlap_rate/k;
}

double CalR2(vector<Clu_Point> &data, vector<vector<Point> > &ref_data, int k, int num_ref){
    double r2 = 0;
    for(int i = 0; i < k; ++i){
        int len = data[i].clu_point.size();
        double clu_r2 = 0;
        for(int r = 0; r < num_ref; ++r){
            vector<double> arr_dis;
            arr_dis.reserve(len);
            for(int l = 0; l < len; ++l)
                arr_dis.push_back(CaculateEuclideanDis(ref_data[i][r],data[i].clu_point[l]));
            sort(arr_dis.begin(),arr_dis.end());

            // cout << arr_dis[0] << "   "<< arr_dis[len-1] << endl;

            double sum_xy = 0, sum_x = 0, sum_y = 0, sum_pow2x = 0;
            for(int y = 0; y < len; ++y){
                sum_x += arr_dis[y];
                sum_y += y;
                sum_pow2x += arr_dis[y] * arr_dis[y];
                sum_xy += arr_dis[y] * y;
            }
            double avg_x = sum_x / len;
            double avg_y = sum_y / len;
            double avg_xy = sum_xy / len;
            double avg_pow2x = sum_pow2x / len;
            double a = (avg_xy - avg_x * avg_y) / (avg_pow2x - avg_x * avg_x);
            double b = avg_y - a * avg_x;


            double sum = 0,sum_var = 0;

            for(int c = 0; c < len; ++c){
                int cal_y = (int)(a * arr_dis[c] + b);
                sum += pow(c-cal_y,2);
                sum_var += pow(c-avg_y,2);
            }
            sum /= len;
            sum_var /= len;
            clu_r2 = clu_r2 + 1 - sum / sum_var;
        }
        clu_r2 /= num_ref;
        // cout << clu_r2 << endl;

        r2 += clu_r2;
    }
    r2 /= k;



    return r2;
}


int main(int argc, const char* argv[]){
    srand((int)time(0));
    // begining
    if(argc < 5){
        cout << "Plz follow format of : ./main [number of reference point] [number of cluster] [dimention of point] [evaluate the clustering or not(0 /1 : [min number of cluster] [max number of cluster] [step size])]" << endl;
        return 0;
    }

    unsigned num_ref = stoi(argv[1]);
    unsigned num_clu = stoi(argv[2]);
    unsigned dim = stoi(argv[3]);
    unsigned flag_eva = stoi(argv[4]);


    if(num_ref > dim + 1){
        cout << "The number of reference point should be smaller than dimentional of point data. Plz reinput parameter" << endl;
        return 0;
    }





    // Pretreatment : clustering and cluster evaluation

    if(flag_eva){
        unsigned min_clu = stoi(argv[5]);
        unsigned max_clu = stoi(argv[6]);
        unsigned step = stoi(argv[7]);

        string data_path = "/home/qrstu/yty/dataset/color/color_32.txt";
        for(unsigned i = min_clu; i < max_clu + 1; i += step){
            // cout << "k is : " << i << endl;

            // pnum = i;
            readm(data_path, dim, 0);
            int num_alldata = num;
            cout << num_alldata << endl;
            // for(int k=0;k<=num;k++)fgg[k]=0;
            // work(dim);

            string clu_data_path = "/home/qrstu/yty/LIMS/outputFiles/K_" + to_string(i);
            // int isCreate = mkdir(clu_data_path.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
            // if( !isCreate )
            //     cout << "create path success : " << clu_data_path << endl;
            // else{
            //     cout << "create path failed ! error code :  " << isCreate << clu_data_path <<endl;
            //     return 0;
            // }

            // op(clu_data_path+"/ref.txt",dim);
            // clu(dim, clu_data_path);


            pnum = num_ref-1;
            for(unsigned j = 0; j < i; ++j){
                string clu_path = clu_data_path + "/clu_" + to_string(j) + ".txt";
                readm(clu_path, dim, 1);
                for(int k=0;k<=num;k++)fgg[k]=0;
                work(dim);
                op(clu_data_path+"/ref_"+to_string(j)+".txt",dim);
            }


            // reread cluster result data file and caculate average err and silhouette coefficient
            vector<Clu_Point> data;
            data.reserve(i);
            vector<vector<Point> > ref_data;
            ref_data.resize(i);

            string ref_path = clu_data_path+"/ref.txt";
            ifstream fin;
            fin.open(ref_path);
            if(!fin){
                cout << ref_path << " file could not be opened\n";
                exit(0);
            }
            string line;
            int index = 0;
            while(getline(fin, line)){
                stringstream ss(line);
                string value;
                vector<double> coor;
                while(getline(ss, value, ',')){
                    coor.push_back(atof(value.c_str()));
                }
                Point main_ref = Point(coor);
                ref_data[index++].push_back(main_ref);
            }
            fin.close();


            for(unsigned j = 0; j < i; ++j){
                InputReader inputReader(clu_data_path + "/clu_" + to_string(j) + ".txt");
                data.push_back(inputReader.getCluster());

                string oth_ref_path = clu_data_path + "/ref_" + to_string(j) + ".txt";
                fin.open(oth_ref_path);
                if(!fin){
                    cout << oth_ref_path << " file could not be opened\n";
                    exit(0);
                }
                while(getline(fin, line)){
                    stringstream ss(line);
                    string value;
                    vector<double> coor;
                    while(getline(ss, value, ',')){
                        coor.push_back(atof(value.c_str()));
                    }
                    Point oth_ref = Point(coor);
                    ref_data[j].push_back(oth_ref);
                }
                fin.close();
            }


            vector<double> dist_max;
            dist_max.resize(i);

            // double s = CalSihCoe(data,i);
            double avg_err = CalErrAvg(data,ref_data,dist_max,i,num_ref,num_alldata);
            double rate = CalOverlapRate(ref_data,dist_max,i);
            double r2 = CalR2(data,ref_data,i,num_ref);


            // cout << "K is : " << i << " Silhouette Coefficient : " << s 
            // << " average err of RP model : " << avg_err << " overlap rate is : " << rate << endl;

            cout << "K is : " << i  << " average err of RP model : " 
            << avg_err << " overlap rate is : " << rate << " r2 is :" << r2 << endl;

        }

        return 0;
    }







    // all data
    vector<Clu_Point> all_data;
    vector<mainRef_Point> all_refSet;


    vector<Point> pivots;
    pivots.reserve(num_clu);
    vector<vector<Point> > oth_pivots;
    oth_pivots.reserve(num_clu);
    string filename = "/home/qrstu/yty/LIMS/inputFiles/ref/ref.txt";
    ifstream fin;
    fin.open(filename);
    if(!fin){
        cout << filename << " file could not be opened\n";
        exit(0);
    }
    string line;
    while(getline(fin,line)){
        stringstream ss(line); 
        string value;
        vector<double> coordinate;
        while(getline(ss, value, ',')){
            coordinate.push_back(atof(value.c_str()));
        }
        Point pivot_pt = Point(coordinate);
        pivots.push_back(pivot_pt);
    }
    fin.close();





    for(unsigned p = 0; p < num_clu; ++p){
        string filename = "/home/qrstu/yty/LIMS/inputFiles/ref/ref_" + to_string(p) +".txt";
        ifstream fin;
        fin.open(filename);
        if(!fin){
            cout << filename << " file could not be opened\n";
            exit(0);
        }
        string line;

        // getline(fin,line);
        // stringstream ss(line); 
        // string value;
        // vector<double> coordinate_main;
        // while(getline(ss, value, ',')){
        //     coordinate_main.push_back(atof(value.c_str()));
        // }
        // Point pivot_pt = Point(coordinate_main);
        // pivots.push_back(pivot_pt);


        vector<Point> other_pivot;
        while(getline(fin,line)){
            stringstream ss(line); 
            string value;
            vector<double> coordinate;
            while(getline(ss, value, ',')){
                coordinate.push_back(atof(value.c_str()));
            }
            Point pivot_pt = Point(coordinate);
            other_pivot.push_back(pivot_pt);
        }
        oth_pivots.push_back(other_pivot);
        fin.close();
    }


    double build_time = 0.0;
    for(unsigned i = 0; i < num_clu; i++){
        InputReader inputReader("/home/qrstu/yty/LIMS/inputFiles/clu/8d_" + to_string(i) + ".txt");
        all_data.push_back(inputReader.getCluster());

        if(all_data[i].clu_point.empty()){
            cout << "Plz do not load a null file" << endl;
            return 0;
        }

        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        ChooseRef ref_point(num_ref - 1, all_data[i], pivots[i],oth_pivots[i],1);
        CaculateIValue calIValue(all_data[i], ref_point.getMainRefPoint());
        chrono::steady_clock::time_point end = chrono::steady_clock::now();

        build_time += chrono::duration_cast<chrono::milliseconds>(end - begin).count();

        all_data[i] = calIValue.getCluster();
        all_refSet.push_back(calIValue.getMainRef_Point());
    }
    // cout << "build time is : " << build_time << endl;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    for(unsigned i = 0; i < num_clu; i++){
        OutputPrinter output(i, all_refSet[i]);
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    build_time += chrono::duration_cast<chrono::milliseconds>(end - begin).count();

    cout << "build time is : " << build_time << endl;



    unsigned K = 5;
    double time = 0.0;
    int page = 0;

    filename = "/home/qrstu/yty/LIMS/inputFiles/KNN.txt";
    vector<Point> list_KNN = LoadPointForQuery(filename);
    
    for(unsigned m = 0; m < list_KNN.size(); ++m){


        // 初始半径
        double init_r = 0.05;
        // 递增半径
        double delta_r = 0.05;



        // 标志位数组，判断每个点检查过没有
        int8_t arr[num_clu][20000] = {0};
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();

        priority_queue<pair<double, unsigned int> > KNNRes_queue;
        while(KNNRes_queue.size() < K || KNNRes_queue.top().first > init_r + delta_r){

            for(unsigned i = 0; i < num_clu; ++i){
                CalCirclePos mainRefPtCircle(all_refSet[i].point, all_refSet[i].r, list_KNN[m], init_r);
                if(mainRefPtCircle.label == 1){
                    continue;
                }
                if(mainRefPtCircle.label == 3 && mainRefPtCircle.dis_upper < all_refSet[i].r_low){
                    continue;
                }
                if(mainRefPtCircle.label == 3 && init_r > all_refSet[i].r){
                    for(unsigned j = 0; j < all_refSet[i].iValuePts.size(); ++j){
                        double dis_pt_qrpt = CaculateEuclideanDis(all_refSet[i].iValuePts[j], list_KNN[m]);
                        if(KNNRes_queue.size() < K){
                            arr[i][j] = 1;
                            pair<double, unsigned int> elem(dis_pt_qrpt,all_refSet[i].iValuePts[j].id);
                            KNNRes_queue.push(elem);
                        }else{
                            if(KNNRes_queue.top().first > dis_pt_qrpt && arr[i][j] == 0){
                                arr[i][j] = 1;
                                KNNRes_queue.pop();
                                pair<double, unsigned int> elem(dis_pt_qrpt,all_refSet[i].iValuePts[j].id);
                                KNNRes_queue.push(elem);
                            }
                        }
                    }
                    continue;
                }


                bool flag = true;
                vector<CalCirclePos> ref_query;
                ref_query.reserve(num_ref - 1);
                for(unsigned j = 0; j < num_ref - 1; ++j){
                    CalCirclePos RefPtCircle(all_refSet[i].ref_points[j].point, all_refSet[i].ref_points[j].r, list_KNN[m], init_r);
                    if(RefPtCircle.label == 1){
                        flag = false;
                        break;
                    }
                    if(RefPtCircle.label == 3 && RefPtCircle.dis_upper < all_refSet[i].ref_points[j].r_low){
                        flag = false;
                        break;
                    }
                    if(RefPtCircle.label == 3 && init_r > all_refSet[i].ref_points[j].r){
                        for(unsigned l = 0; l < all_refSet[i].iValuePts.size(); ++l){
                            double dis_pt_qrpt = CaculateEuclideanDis(all_refSet[i].iValuePts[l], list_KNN[m]);
                            if(KNNRes_queue.size() < K){
                                arr[i][l] = 1;
                                pair<double, unsigned int> elem(dis_pt_qrpt,all_refSet[i].iValuePts[l].id);
                                KNNRes_queue.push(elem);
                            }else{
                                if(KNNRes_queue.top().first > dis_pt_qrpt && arr[i][l] == 0){
                                    arr[i][l] = 1;
                                    KNNRes_queue.pop();
                                    pair<double, unsigned int> elem(dis_pt_qrpt,all_refSet[i].iValuePts[l].id);
                                    KNNRes_queue.push(elem);
                                }
                            }
                        }
                        break;
                    }
                    ref_query.push_back(RefPtCircle);
                }

                if(!flag)
                    continue;

                // KNN KNNQuery(list_KNNQry[m], all_refSet[i], KNNRes_queue, mainRefPtCircle, ref_query, K, arr, i);
                KNN_IO KNNQuery(list_KNN[m], all_refSet[i], KNNRes_queue, mainRefPtCircle, ref_query, K, arr, i, page);
                // OutputPrinter outputPrinter2;
                // outputPrinter2.print_ivalue(outputFileName + "ivalue.txt", all_refSet[i].iValuePts);
            }

            init_r += delta_r;
        }

        // 最后一轮准确范围查询
        for(unsigned i = 0; i < num_clu; ++i){
            CalCirclePos mainRefPtCircle(all_refSet[i].point, all_refSet[i].r, list_KNN[m], init_r);
            if(mainRefPtCircle.label == 1){
                continue;
            }
            if(mainRefPtCircle.label == 3 && mainRefPtCircle.dis_upper < all_refSet[i].r_low){
                continue;
            }
            if(mainRefPtCircle.label == 3 && init_r > all_refSet[i].r){
                for(unsigned j = 0; j < all_refSet[i].iValuePts.size(); ++j){
                    double dis_pt_qrpt = CaculateEuclideanDis(all_refSet[i].iValuePts[j], list_KNN[m]);
                    if(KNNRes_queue.size() < K){
                        arr[i][j] = 1;
                        pair<double, unsigned int> elem(dis_pt_qrpt,all_refSet[i].iValuePts[j].id);
                        KNNRes_queue.push(elem);
                    }else{
                        if(KNNRes_queue.top().first > dis_pt_qrpt && arr[i][j] == 0){
                            KNNRes_queue.pop();
                            pair<double, unsigned int> elem(dis_pt_qrpt,all_refSet[i].iValuePts[j].id);
                            KNNRes_queue.push(elem);
                        }
                    }
                }
                continue;
            }


            bool flag = true;
            vector<CalCirclePos> ref_query;
            ref_query.reserve(num_ref - 1);
            for(unsigned j = 0; j < num_ref - 1; ++j){
                CalCirclePos RefPtCircle(all_refSet[i].ref_points[j].point, all_refSet[i].ref_points[j].r, list_KNN[m], init_r);
                if(RefPtCircle.label == 1){
                    flag = false;
                    break;
                }
                if(RefPtCircle.label == 3 && RefPtCircle.dis_upper < all_refSet[i].ref_points[j].r_low){
                    flag = false;
                    break;
                }
                if(RefPtCircle.label == 3 && init_r > all_refSet[i].ref_points[j].r){
                    for(unsigned l = 0; l < all_refSet[i].iValuePts.size(); ++l){
                        double dis_pt_qrpt = CaculateEuclideanDis(all_refSet[i].iValuePts[l], list_KNN[m]);
                        if(KNNRes_queue.size() < K){
                            arr[i][l] = 1;
                            pair<double, unsigned int> elem(dis_pt_qrpt,all_refSet[i].iValuePts[l].id);
                            KNNRes_queue.push(elem);
                        }else{
                            if(KNNRes_queue.top().first > dis_pt_qrpt && arr[i][l] == 0){
                                arr[i][l] = 1;
                                KNNRes_queue.pop();
                                pair<double, unsigned int> elem(dis_pt_qrpt,all_refSet[i].iValuePts[l].id);
                                KNNRes_queue.push(elem);
                            }
                        }
                    }
                    break;
                }
                ref_query.push_back(RefPtCircle);
            }

            if(!flag)
                continue;

            // KNN KNNQuery(list_KNNQry[m], all_refSet[i], KNNRes_queue, mainRefPtCircle, ref_query, K, arr, i);
            KNN_IO KNNQuery(list_KNN[m], all_refSet[i], KNNRes_queue, mainRefPtCircle, ref_query, K, arr, i, page);
        
        }



        chrono::steady_clock::time_point end = chrono::steady_clock::now();

        time += chrono::duration_cast<chrono::microseconds>(end - begin).count();

    }

    time /= list_KNN.size();
    page /= list_KNN.size();
    cout << " KNN LIMS index time is : " << time << " [µs]" <<  endl;
    cout << "avg page : " << page << endl;


    // range query
    string rangeQuery_filename = "/home/qrstu/yty/LIMS/inputFiles/range_0.2.txt";
    vector<Point> list_rangeQry = LoadPointForQuery(rangeQuery_filename);
    double r = 0.2;
    time = 0.0;
    double io_time = 0.0;
    int queryPt_num = 0;
    page = 0;
    for(unsigned m = 0; m < list_rangeQry.size(); ++m){


        vector<int> rangeQueryRes;
        rangeQueryRes.reserve(10000);

        chrono::steady_clock::time_point begin = chrono::steady_clock::now();

        for(unsigned i = 0; i < num_clu; ++i){
            CalCirclePos mainRefPtCircle(all_refSet[i].point, all_refSet[i].r, list_rangeQry[m], r);
            if(mainRefPtCircle.label == 1){
                continue;
            }
            if(mainRefPtCircle.label == 3 && mainRefPtCircle.dis_upper < all_refSet[i].r_low){
                continue;
            }
            if(mainRefPtCircle.label == 3 && r > all_refSet[i].r){

                ifstream fin;
                int start_page = 0;
                int end_page = all_refSet[i].iValuePts.size() / Constants::PAGE_SIZE;

                for(int s = start_page; s <= end_page; ++s){
                    vector<vector<double> > pt_page;
                    string fileName = "./data/cluster_" + to_string(i) + "_" + to_string(s);
                    ++page;
                    fin.open(fileName, ios::binary);
                    int size_page;
                    fin.read((char*)&size_page,4);
                    pt_page.reserve(size_page);
                    int size_dim;
                    fin.read((char*)&size_dim,4);
                    for(int q = 0 ; q < size_page; q++){
                        vector<double> coordinate;
                        coordinate.resize(size_dim);
                        fin.read(reinterpret_cast<char*>(&coordinate), sizeof(coordinate));
                        pt_page.push_back(coordinate);
                    } 
                    fin.close();
                    for(int l = 0+s*Constants::PAGE_SIZE; l < 0+s*Constants::PAGE_SIZE+size_page; ++l)
                        rangeQueryRes.push_back(all_refSet[i].iValuePts[l].id);
                    // rangeQueryRes.insert(rangeQueryRes.end(), all_refSet[i].iValuePts.begin()+s*Constants::PAGE_SIZE, all_refSet[i].iValuePts.begin()+s*Constants::PAGE_SIZE+size_page);
                }
                continue;
            }

            bool flag = true;
            vector<CalCirclePos> ref_query;
            ref_query.reserve(num_ref - 1);
            for(unsigned j = 0; j < num_ref - 1; ++j){
                CalCirclePos RefPtCircle(all_refSet[i].ref_points[j].point, all_refSet[i].ref_points[j].r, list_rangeQry[m], r);
                if(RefPtCircle.label == 1){
                    flag = false;
                    break;
                }
                ref_query.push_back(RefPtCircle);
            }

            if(!flag)
                continue;


            RangeQuery_IO rangeQuery( list_rangeQry[m], r, all_refSet[i], rangeQueryRes, mainRefPtCircle, ref_query,i,page,io_time);
        }
        chrono::steady_clock::time_point end = chrono::steady_clock::now();

        time += chrono::duration_cast<chrono::microseconds>(end - begin).count();
        queryPt_num += rangeQueryRes.size();
    }
    time /= list_rangeQry.size();
    io_time /= list_rangeQry.size();
    queryPt_num /= list_rangeQry.size();
    page /= list_rangeQry.size();
    cout << "range query LIMS index time is : " << time << " [µs]" <<  endl;
    cout << "io time : " << io_time << " [µs]" << endl;
    cout << "query point number avg is : " << queryPt_num << endl;
    cout << "avg page : " << page << endl;




    // // test pointquery
    // double time = 0.0;
    // int page = 0;
    // int acc = 0;
    // // filename = "/home/qrstu/yty/dataset/syn/guass/range_0.01_2d/point.txt";
    // // vector<Point> list_point = LoadPointForQuery(filename);

    // for(unsigned o = 0; o < 200; ++o){
    //     int cl = random(0,49);
    //     int pt = random(0,100);
    //     vector<double> queryPt;
    //     for(int d = 0; d < dim; ++d){
    //         double coor = random(0,100000) / 100000.0;
    //         queryPt.push_back(coor);
    //     }
    //     Point point = Point(queryPt);

    //     chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    //     int query_res = -1;
    //     // double dis = CaculateEuclideanDis(all_data[cl].clu_point[pt], all_refSet[0].point);
    //     double dis = CaculateEuclideanDis(point, all_refSet[0].point);
    //     int clu_id = 0;
    //     for(unsigned i = 1; i < num_clu; ++i){
    //         // double dis_other = CaculateEuclideanDis(all_data[cl].clu_point[pt], all_refSet[i].point);
    //         double dis_other = CaculateEuclideanDis(point, all_refSet[i].point);
    //         if(dis_other < dis){
    //             dis = dis_other;
    //             clu_id = i;
    //         }
    //     }
    //     // PointQuery_IO pointQuery_IO(all_data[cl].clu_point[pt], all_refSet[clu_id], query_res,clu_id,dis,page);
    //     PointQuery_IO pointQuery_IO(point, all_refSet[clu_id], query_res,clu_id,dis,page);
    //     chrono::steady_clock::time_point end = chrono::steady_clock::now();
    //     time += chrono::duration_cast<chrono::nanoseconds>(end - begin).count();
    //     if(query_res!= -1)
    //         acc++;
    // }
    // time /= 200;
    // cout << "point query LIMS index avg time : " << time << "[ns]" << endl;
    // cout << "page avg : " << page<< endl;
    // cout << "acc " <<  acc << endl;






    // test insert and query

    // string insert_file = "/home/qrstu/yty/LIMS/inputFiles/insert/skew_0.001.txt";
    // Insert(insert_file,all_refSet,num_clu);


    // rangeQuery_filename = "/home/qrstu/yty/LIMS/inputFiles/range_0.2.txt";
    // list_rangeQry = LoadPointForQuery(rangeQuery_filename);
    // r = 0.2;
    // time = 0.0;
    // io_time = 0.0;
    // queryPt_num = 0;
    // page = 0;
    // for(unsigned m = 0; m < list_rangeQry.size(); ++m){

    //     vector<int> rangeQueryRes;
    //     rangeQueryRes.reserve(100000);

    //     chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    //     for(unsigned i = 0; i < num_clu; ++i){
    //         CalCirclePos mainRefPtCircle(all_refSet[i].point, all_refSet[i].r, list_rangeQry[m], r);
    //         if(mainRefPtCircle.label == 1){
    //             continue;
    //         }
    //         if(mainRefPtCircle.label == 3 && mainRefPtCircle.dis_upper < all_refSet[i].r_low){
    //             continue;
    //         }
    //         if(mainRefPtCircle.label == 3 && r > all_refSet[i].r){

    //             ifstream fin;
    //             int start_page = 0;
    //             int end_page = all_refSet[i].iValuePts.size() / Constants::PAGE_SIZE;

    //             for(int s = start_page; s <= end_page; ++s){
    //                 vector<vector<double> > pt_page;
    //                 string fileName = "./data/cluster_" + to_string(i) + "_" + to_string(s);
    //                 fin.open(fileName, ios::binary);
    //                 int size_page;
    //                 fin.read((char*)&size_page,4);
    //                 pt_page.reserve(size_page);
    //                 int size_dim;
    //                 fin.read((char*)&size_dim,4);
    //                 for(int q = 0 ; q < size_page; q++){
    //                     vector<double> coordinate;
    //                     coordinate.resize(size_dim);
    //                     fin.read(reinterpret_cast<char*>(&coordinate), sizeof(coordinate));
    //                     pt_page.push_back(coordinate);
    //                 } 
    //                 fin.close();
    //                 for(int l = 0+s*Constants::PAGE_SIZE; l < 0+s*Constants::PAGE_SIZE+size_page; ++l)
    //                     rangeQueryRes.push_back(all_refSet[i].iValuePts[l].id);
    //             }
    //             continue;
    //         }


    //         double dis_qrpt_mainRef = CaculateEuclideanDis(all_refSet[i].point,list_rangeQry[m]);
    //         int insert_pt_size = all_refSet[i].insert_list.size();
    //         double lower = 0.0;
    //         double upper = 0.0;
    //         int low_pos = -1;
    //         int up_pos = -1;
    //         if(dis_qrpt_mainRef - r < all_refSet[i].insert_list[insert_pt_size-1].i_value){
    //             if(dis_qrpt_mainRef - r < 0){
    //                 lower = 0;
    //                 low_pos = 0;
    //             }else{
    //                 lower = dis_qrpt_mainRef - r;
    //             }

    //             if(dis_qrpt_mainRef + r < all_refSet[i].insert_list[insert_pt_size-1].i_value){
    //                 upper = dis_qrpt_mainRef + r;
    //             }else{
    //                 upper = all_refSet[i].insert_list[insert_pt_size-1].i_value;
    //                 up_pos = insert_pt_size-1;
    //             }

    //             if(low_pos == -1)
    //                 low_pos = BinarySearch(all_refSet[i].insert_list,lower,0,insert_pt_size);
    //             if(up_pos == -1)
    //                 up_pos = BinarySearch(all_refSet[i].insert_list,upper,0,insert_pt_size);
                
                
                
    //             int start_page = low_pos / Constants::PAGE_SIZE;
    //             int end_page = up_pos / Constants::PAGE_SIZE;

                

    //             for(int s = start_page; s <= end_page; ++s){
    //                 ++page;
    //                 vector<vector<double> > pt_page;
    //                 chrono::steady_clock::time_point start = chrono::steady_clock::now();
    //                 string fileName = "./data/cluster_" + to_string(i) + "_insert_" + to_string(s);
    //                 // cout << fileName << endl;
    //                 fin.open(fileName, ios::binary);
    //                 int size_page;
    //                 fin.read((char*)&size_page,4);
    //                 pt_page.reserve(size_page);
    //                 int size_dim;
    //                 fin.read((char*)&size_dim,4);
    //                 for(int q = 0 ; q < size_page; q++){
    //                     vector<double> coordinate;
    //                     coordinate.resize(size_dim);
    //                     fin.read(reinterpret_cast<char*>(coordinate.data()),size_dim * sizeof(coordinate.front()));
    //                     pt_page.push_back(coordinate);
    //                 }
    //                 fin.close();
    //                 chrono::steady_clock::time_point end = chrono::steady_clock::now();
    //                 io_time+=chrono::duration_cast<chrono::microseconds>(end - start).count();
    //                 for(int x = 0; x < size_page; ++x){
    //                     if(CaculateEuclideanDis2(pt_page[x], list_rangeQry[m]) <= r){
    //                         rangeQueryRes.push_back(all_refSet[i].insert_list[s*Constants::PAGE_SIZE+x].id);
    //                     }
    //                 }       
    //             }


                
    //             // for(; low_pos <= up_pos; ++low_pos){
    //             //     if(CaculateEuclideanDis2(all_refSet[i].insert_list[low_pos].coordinate,list_rangeQry[m]) <= r)
    //             //         rangeQueryRes.push_back(all_refSet[i].insert_list[low_pos].id);
    //             // }
    //         }



    //         bool flag = true;
    //         vector<CalCirclePos> ref_query;
    //         ref_query.reserve(num_ref - 1);
    //         for(unsigned j = 0; j < num_ref - 1; ++j){
    //             CalCirclePos RefPtCircle(all_refSet[i].ref_points[j].point, all_refSet[i].ref_points[j].r, list_rangeQry[m], r);
    //             if(RefPtCircle.label == 1){
    //                 flag = false;
    //                 break;
    //             }
    //             ref_query.push_back(RefPtCircle);
    //         }

    //         if(!flag)
    //             continue;



    //         RangeQuery_IO rangeQuery( list_rangeQry[m], r, all_refSet[i], rangeQueryRes, mainRefPtCircle, ref_query,i,page,io_time);
    //     }
    //     chrono::steady_clock::time_point end = chrono::steady_clock::now();

    //     time += chrono::duration_cast<chrono::microseconds>(end - begin).count();
    //     queryPt_num += rangeQueryRes.size();
    // }
    // time /= list_rangeQry.size();
    // io_time /= list_rangeQry.size();
    // queryPt_num /= list_rangeQry.size();
    // page /= list_rangeQry.size();
    // cout << "range query LIMS index time is : " << time << " [µs]" <<  endl;
    // cout << "io time : " << io_time << " [µs]" << endl;
    // cout << "query point number avg is : " << queryPt_num << endl;
    // cout << "avg page : " << page << endl;

    return 0;
}
