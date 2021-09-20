#include <string>
#include <cstring>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/time.h> 
#define NUUM 10000000

using namespace std;

float* loc[NUUM];


int pnum,num,PairA=0; // pnum : number of cluster / pivots
// ptype : type of distance calculate function
int ptype=2;
double mm;
double dis(int i,int j,int dim){
	if(ptype==2){
		double sum=0;
		for(int k=0;k<dim;k++) {
			sum+=pow(loc[i][k]-loc[j][k],2);	
		//	cout<<loc[i][k]<<" "<<loc[j][k]<<endl;
		}
		return pow(sum,0.5);	
	}
	if(ptype==1){
		double sum=0;
		for(int k=0;k<dim;k++) sum+=fabs(loc[i][k]-loc[j][k]);
		return sum;	
	}

	if(ptype==5){
		double sum=0;
		for(int k=0;k<dim;k++) sum+=pow(loc[i][k]-loc[j][k],5);
		return pow(sum,0.2);	
	}

	if(ptype==0){
		double max=0,p;
		for(int k=0;k<dim;k++) 
		{			
			p=abs(loc[i][k]-loc[j][k]);
			if(p>max) max=p;
		}
		return max;	
	}	
}
void readm(string filename, int dim, int type){ //re-read the file
	// string filename = "/home/qrstu/yty/LIMS/inputFiles/clu/8d_" + to_string(p) + ".txt";
	// string filename = "/home/qrstu/yty/LIMS/inputFiles/skew_8d.txt";
	ifstream ain1;
//	cout<<"Input pnum : ";
//	cin>>pnum;
	ain1.open(filename, ifstream::in);	
	// ain1>>dim>>num>>ptype;

	num = 2000000;
	ptype = 2;
	//cout<<dim<<" "<<num<<" "<<ptype<<endl; // dim = 20 num = 1000 dist_metric = L0
	//PairA=num/100; 
	PairA = 100000; // 10000对于dynamic太慢了
	//num=100000;
	for(int i=0;i<num;i++)loc[i]=new float[dim+3];


	string line;
	int i = 0;
	// type == 1, load each cluster data file
	if(type){
		while(getline(ain1,line)){
			// if(i%300000==0) cout<<"read finish "<<i<<endl;
			stringstream ss(line);
			string value;
			int j = 0;
			// pass the id of each point
			getline(ss,value,',');
			while(getline(ss,value,',')){
				loc[i][j] = atof(value.c_str());
				++j;
			}
			++i;
		}
	}else{
		// type == 0, load whole data file
		while(getline(ain1,line)){
			// if(i%300000==0) cout<<"read finish "<<i<<endl;
			stringstream ss(line);
			string value;
			int j = 0;
			while(getline(ss,value,',')){
				loc[i][j] = atof(value.c_str());
				++j;
			}
			++i;
		}
	}
	num = i;


	// for(int i=0;i<num;i++){
	// 	if(i%30000==0) cout<<"read finish "<<i<<endl;
	// 	for(int j=0;j<dim;j++)
	// 		ain1>>loc[i][j];
	// }
	ain1.close();	
	int j=1;
	mm=0;
	double mm1=0;
	for(int i=1;i<num;i++) if(dis(0,i,dim)>mm1){
		j=i;
		mm1=dis(0,i,dim); // mm1 距离loc[0]最远的点
	}	
	for(int i=1;i<num;i++) if(dis(j,i,dim)>mm) mm=dis(j,i,dim); // mm 距离loc[j]最远的点

}
