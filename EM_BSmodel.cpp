#include<iostream>
#include<random>
#include<cmath>
#include<math.h>
#include<limits>
#include<algorithm>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace std;

//r: 利子率
//sigma: ボラティリティ
//T: 満期
//K: 行使価格
//S_0: 現在の株価

//EM近似で価格計算
//パスと分割数を変えたときの変化を見たい


double Error(int n,int N){ //
  double r=0.1,S_0=62.0,K=60.0,sigma=0.2;
  //int N=10000; //パス
  int T=1;
  //int n=1000; //分割数
  double h=double(T)/double(n); //刻み幅

  random_device seed;
  mt19937 engine(seed()); //メルセンヌ・ツイスターで擬似乱数を生成
  normal_distribution<> dist(0.0,1.0);
  vector<double> z(n);
  vector<double> S(n);
  vector<double> v(N);
  double sum=0.0;
  S.at(0)=S_0;
  vector<double> error;

  //理論値の計算
  double dplus=0.0,dminus=0.0,price=0.0;
  dplus=1.0/(sigma*pow(T,0.5))*(log(S_0/K)+(r+0.5*pow(sigma,2))*T);
  dminus=1.0/(sigma*pow(T,0.5))*(log(S_0/K)+(r-0.5*pow(sigma,2))*T);
  price=S_0*0.5*(1+erf(dplus/pow(2,0.5)))-K*exp(-r*T)*0.5*(1+erf(dminus/pow(2,0.5)));

  //標本近似
  for(int j=0;j<N;j++){ //パスの本数
    for(int i=0;i<n-1;i++){ //1本のpathに対しEM近似(time_step固定)
      z.at(i)=dist(engine);
      S.at(i+1)=S.at(i)+r*S.at(i)*h+S.at(i)*sigma*pow(h,0.5)*z.at(i);
    }
    v.at(j)=exp(-r*T)*max(S.at(n-1)-K,0.0);
    sum+=v.at(j);
    //error.push_back(abs(sum/(j+1)-price));
  }

  return abs(sum/double(N)-price);
}

void visualization(vector<double> x, vector<double> y){ //グラフをかく
  plt::backend("Agg");
  plt::plot(x,y);
  plt::xlabel("h");
  plt::ylabel("error");
  plt::save("Graph.png");
  return;
}

int main(){

  int N=5000; //発生させるパス
  int n=10000; //time step =>S(0)~S(T)のstep、S(0)~S(t1)~S(T)のstepの2つ（timestepは2つから)の上限
  //vector<vector<double>> error(n,vector<double>(N)); //行を分割数、列をパスにしたときの誤差をみる
  vector<double> time_steps,error;

  for(int i=1;i<=n;i++){ //time_stepを変化
    time_steps.push_back(-log10(1.0/double(i)));
    error.push_back(log10(Error(i,N)));
    //cout<<time_steps.at(i-1)<<endl;
  }

  visualization(time_steps,error);
}
