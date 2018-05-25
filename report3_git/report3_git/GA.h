#pragma once

#include<vector>
#include <string>
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>

//#define __ENABLE_SINGLE_POINT_MUTATION__

class GA
{
private:
	double individualMutationRate = 0.05;//個体突然変異率
	double genomMutationRate = 0.1;
	double alpha = 0.3;
	std::vector<double> varMax, varMin;//変数の最小値・最大値
	bool isChanged = false;
public:
	double resultSumValue;//評価関数の合計

	class Data//データ格納用クラス
	{
	public:
		std::vector<int> num;//座標
		double functionValue;//与えられた関数の値
		double result;

		Data(int _var_num)//コンストラクタ
		{
			num.resize(_var_num);//isIncludedの配列の長さの設定
		}
	};

	class CityData
	{
	public:
		int cityNum;
		std::vector<double> point;

		CityData(int _point_dim)
		{
			point.resize(_point_dim);
		}
	};

	std::vector<Data> data, prev_data;//操作前後で値を保持するために2個
	Data eliteData;
	std::vector<CityData> model;
	GA(int _max_genom_list, int _var_num, std::vector<GA::CityData> model);	//コンストラクタ
	bool selection();//選択
	void blxAlphaCrossover();
	void pmxCrossover();
	void mutation();//突然変異
	void calc(bool enableDisplay, bool enableOneLine = false);//評価関数の計算
private:
	void calcResult(bool enableSort = false);
	int random(int min, int max);
	double random(int min, double max);
	double random(double min, int max);
	double random(double min, double max);
	void displayValues(bool enableOneLine);
	Data searchRank(int num);
public:
	~GA();//デコンストラクタ
};