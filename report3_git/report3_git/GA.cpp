#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _var_num, std::vector<GA::CityData> _model) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//dataの初期化
	eliteData(_var_num)
{
	//もらった変数をクラス内変数に格納
	model = _model;

	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].x.size(); j++)
		{
			data[i].x[j] = random(varMin[j], varMax[j]);//遺伝子の初期設定
		}
	}
	prev_data = data;
	calcResult();

	displayValues(false);
}

bool GA::selection()
{
	int max_num = 0;//最も評価の良い個体の番号
	bool ret = isChanged;//最も評価の良い個体の変化の監視(デバッグ用)
	isChanged = false;

	eliteData = searchRank(0);//最も評価の良い個体を保持
	prev_data = data;
	for (int i = 0; i < data.size(); i++)
	{
		double selector = random(0.0, 1.0);//乱数を生成
		double needle = 0;//ルーレットの針を生成
		int j = 0;
		for (; ; j++)
		{
			needle += (prev_data[j].result / resultSumValue);//ルーレットの針を乱数の値まで進める
			if (needle > selector)
				break;
			if (j == (data.size() - 1))
				break;
		}
		data[i] = prev_data[j];
	}
	return ret;
}

void GA::blxAlphaCrossover()
{
	prev_data = data;

	for (int i = 0; i < data.size(); i += 2)//2個ずつ交叉
	{
		for (int j = 0; j < data[i].x.size(); j++)
		{
			double ave = (data[i].x[j] + data[i + 1].x[j]) / 2;
			double length = std::abs((data[i].x[j] - data[i + 1].x[j]));

			data[i].x[j] = random(ave - length * (1 / 2 + alpha), ave + length * (1 / 2 + alpha));
			data[i + 1].x[j] = random(ave - length * (1 / 2 + alpha), ave + length * (1 / 2 + alpha));
		}
	}
}

void GA::pmxCrossover()
{
	prev_data = data;

	for (int i = 0; i < data.size(); i += 2)//2個ずつ交叉
	{
		int del1 = random(0, (int)data[i].x.size() - 1);
		int del2 = random(del1, (int)data[i].x.size());

		for (int j = del1; j < del2; j++)
		{
			data[i + 1].x[j] = prev_data[i].x[j];
			data[i].x[j] = prev_data[i + 1].x[j];
		}
	}
}

void GA::mutation()
{
	for (int i = 0; i < data.size(); i++)
	{
		if (random(0.0, 1.0) <= individualMutationRate)//個体突然変異率の計算
		{
#ifdef __ENABLE_SINGLE_POINT_MUTATION__
			int pos = random(0, (int)data[i].x.size() - 1);
			data[i].x[pos] += random(varMin[pos] * random(0, genomMutationRate), varMax[pos] * random(0, genomMutationRate));
#else
			for (int j = 0; j < data[i].x.size(); j++)
			{
				if (random(0.0, 1.0) <= genomMutationRate)
					data[i].x[j] = random(varMin[j], varMax[j]);
			}
#endif
		}
	}
}

void GA::calc(bool enableDisplay, bool enableOneLine)
{
	int minNum = 0;
	calcResult();
	for (int i = 0; i < data.size(); i++)//評価関数が最小の奴を検索
	{
		if (data[i].result < data[minNum].result)
			minNum = i;
	}
	if (searchRank(0).functionValue - eliteData.functionValue > 0.00001)
		isChanged = true;
	//評価関数が最もいいやつを保存
	data[minNum] = eliteData;

	calcResult();

	if (enableDisplay)
		displayValues(enableOneLine);
}

void GA::calcResult(bool enableSort)
{
	int maxNum = 0;

	for (int i = 0; i < data.size(); i++)
	{
		data[i].functionValue = std::sin(data[i].x[0] + data[i].x[1]) + std::pow((data[i].x[0] - data[i].x[1]), 2.0) - 1.5*data[i].x[0] + 2.5*data[i].x[1] + 1;//与えられた関数

		if (data[maxNum].functionValue < data[i].functionValue)//座標の中で最も関数が大きいやつを検索
			maxNum = i;
	}
	double seg = data[maxNum].functionValue;//評価関数の切片を与えられた関数が最も大きいやつにセット
	double seg2 = searchRank(data[0].x.size() - 2).functionValue - seg;
	resultSumValue = 0;
	double coefficient = 0.001 / data[0].x.size();//評価関数用の定数

	for (int i = 0; i < data.size(); i++)
	{
		bool flag = true;

		for (int j = 0; j < data[i].x.size(); j++)
		{
			if (data[i].x[j] > varMax[j] || data[i].x[j] < varMin[j])//座標が場外にいるやつの処理
				flag = false;
		}
		data[i].result = seg2 == 0 ? 0 : (data[i].functionValue - seg) / seg2 / coefficient;//与えられた関数の値から切片で設定した値を引いて2乗する→与えられた関数の値が小さいやつが強くなる
																							//data[i].result = std::abs(data[i].functionValue - seg);

		if (!flag)//場外に出たやつの処理
			data[i].result *= coefficient;
		resultSumValue += data[i].result;
	}
	if (enableSort)
		std::sort(data.begin(), data.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });
}

int GA::random(int min, int max)
{
	//乱数の設定
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(engine);
}

double GA::random(int min, double max) { return random((double)min, max); }
double GA::random(double min, int max) { return random(min, (double)max); }
double GA::random(double min, double max)
{
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(engine);
}

void GA::displayValues(bool enableOneLine)
{
	std::vector<Data> data_temp;
	if (enableOneLine)
	{
		data_temp.push_back(eliteData);
	}
	else
	{
		data_temp = data;
		std::sort(data_temp.begin(), data_temp.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });
	}

	for (int i = 0; i < data_temp.size(); i++)
	{
		for (int j = 0; j < data_temp[i].x.size(); j++)
		{
			printf_s("%10.7lf,", data_temp[i].x[j]);//デバッグ用
		}
		printf_s(" \t f(x,y)=%10.7lf\t Result=%10.7lf\n", data_temp[i].functionValue, data_temp[i].result);
	}
}

GA::Data GA::searchRank(int num)//評価がいい順
{
	std::vector<Data> data_temp = data;
	std::sort(data_temp.begin(), data_temp.end(), [](const Data& x, const Data& y) { return x.functionValue < y.functionValue; });
	return data_temp[num];
}

GA::~GA()
{

}