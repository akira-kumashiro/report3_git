#pragma once

#include<vector>
#include <string>
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <numeric>

//#define __ENABLE_SINGLE_POINT_MUTATION__
//#define __ENABLE_DOUBLE_POINT_MUTATION__
#define __ENABLE_SEGMENTAL_MUTATION__

class GA
{
private:
	double individualMutationRate = 0.4;//個体突然変異率
	double genomMutationRate = 0.3;
	double crossoverRate = 0.8;
	bool isChanged = false;
	double alpha = 1;
	int localMinNum = 0;
	std::vector<int> cityTemp;
	int genNum = 1;
public:
	double resultSumValue;//評価関数の合計

	class Data//データ格納用クラス
	{
	public:
		std::vector<int> num;//座標
		double functionValue;//与えられた関数の値
		double result;

		Data(int _var_num) ://コンストラクタ
			num(std::vector<int>(_var_num, -1))
		{

		}
	};

	class PointXY
	{
	public:
		std::vector<double> point;

		PointXY(double _x, double _y) :
			point(std::vector<double>{_x, _y})
		{

		}
	};

	std::vector<Data> data, prev_data;//操作前後で値を保持するために2個
	Data eliteData;
	std::vector<PointXY> model;
	GA(int _max_genom_list, int _var_num, std::vector<GA::PointXY> _model);	//コンストラクタ
	bool selection();//選択
	void pmxCrossover();
	void mutation();//突然変異
	void calc(bool enableDisplay, bool enableOneLine = false);//評価関数の計算
	void displayValues(bool enableOneLine);
	void calcResult(bool enableSort = false);
private:
	int random(int min, int max);
	double random(int min, double max);
	double random(double min, int max);
	double random(double min, double max);
	Data searchRank(int num);
	void setEmptyNum(void);
public:
	~GA();//デコンストラクタ
};
