#include <string>
#include <cstring>
#include <cmath>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <sys/time.h> 
#define NUUM 200000

#define CandA 300
// #define FILENAME "../sample.txt" 
using namespace std;

string loc[NUUM];
int dim;
int pnum=2,num,PairA=0;
int ptype=0;
double mm=0;
int ** table1;
// double dis(int ii,int jj){
// 	string s1=loc[ii] ;
// 	string s2=loc[jj];
// 	double dist = 0.0;
// 	if (table1 == NULL) {	// init. table if NULL
// 		table1 = new int*[300];
// 		for (int i = 0; i<300; i++)
// 			table1[i] = new int[300];
// 	}


// 	int n = s1.length(), m = s2.length();
// 	if (n == 0) return m;
// 	if (m == 0) return n;
// 	//printf("%d %d\n",m,n);
// 	for (int i = 0; i <= n; i++) table1[i][0] = i;
// 	for (int j = 0; j <= m; j++) table1[0][j] = j;

// 	for (int i = 1; i <= n; i++) {
// 		for (int j = 1; j <= m; j++) {
// 			int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;	// ith character of s, jth character of t
// 			table1[i][j] = 1 + min(table1[i - 1][j], table1[i][j - 1]);
// 			table1[i][j] = min(table1[i - 1][j - 1] + cost, table1[i][j]);
// 		}
// 	}
// 	return (double)(table1[n][m]);
// }


int minDistance(string word1, string word2) {
	int l1 = word1.length() + 1;
	int l2 = word2.length() + 1;
	int md[l1][l2];
	for (int j = 0; j < l2; j++) { md[0][j] = j; } // 插入
	for (int i = 0; i < l1; i++) { md[i][0] = i; } // 删除
	for (int i = 1; i < l1; i++) {
		for (int j = 1; j < l2; j++) {
			md[i][j] = min(min(md[i-1][j]+1, md[i][j-1]+1), md[i-1][j-1] + (word1[i-1] == word2[j-1] ? 0 : 1));
		}
	}
	return md[l1-1][l2-1];
}
void readm(int p){ //re-read the file
	ifstream ain1;
	cout<<"Pnum : 5";
	//cin>>pnum;
	string filename = "/home/qrstu/yty/dataset/word/fft_clu/word_" + to_string(p) + ".txt";
	// string filename = "/home/qrstu/yty/dataset/word/signature.txt";
	ain1.open(filename, ifstream::in);	
	num=0;
	string str;
	while(getline(ain1, str)){
		stringstream ss(str);
		string value;
		getline(ss, value, ',');
		string id = value;

        getline(ss, value, ',');
		loc[num]=value;
		num++;
	}
	cout<<"num "<<num<<endl;
	PairA=num/100;
	ain1.close();	
	mm=0;
	for(int i=0;i<num;i++) if(loc[i].length()>mm) mm=loc[i].length();
	cout<<"maxdis "<<mm<<endl;
}
