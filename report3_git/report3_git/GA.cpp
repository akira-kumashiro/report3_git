#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _var_num, std::vector<GA::PointXY> _model) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//dataの初期化
	eliteData(_var_num),
	cityTemp(_var_num)
{
	//もらった変数をクラス内変数に格納
	model = _model;
	std::iota(cityTemp.begin(), cityTemp.end(), 0);

	//data[0].num = std::vector<int>{ 22, 7, 26, 15,12, 23, 0, 27, 5, 11, 8, 25, 2, 28, 4, 20, 1, 19, 9, 3, 14, 17, 13, 16, 21, 10, 18, 24, 6 };
	setEmptyNum();
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

void GA::pmxCrossover()
{
	prev_data = data;

	//Step1
	for (int i = 0; i < data.size(); i += 2)//2個ずつ交叉
	{
		if (random(0.0, 1.0) <= crossoverRate)
		{
			int del1 = random(0, (int)data[0].num.size() - 1);
			int del2 = random(del1, (int)data[0].num.size());
			for (int j = 0; j < data[i].num.size(); j++)
			{
				if (j<del1 || j>del2)
				{
					data[i].num[j] = -1;
					data[i + 1].num[j] = -1;
				}
				else
				{
					data[i + 1].num[j] = prev_data[i].num[j];
					data[i].num[j] = prev_data[i + 1].num[j];
				}
			}
		}
	}

	//Step2
	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].num.size(); j++)
		{
			if (data[i].num[j] == -1)
			{
				bool isIncluded = false;
				for (int k = 0; k < data[i].num.size(); k++)
				{
					if (prev_data[i].num[j] == data[i].num[k])
					{
						isIncluded = true;
						break;
					}
				}
				if (!isIncluded)
					data[i].num[j] = prev_data[i].num[j];
			}
		}
	}
	//Step3
	setEmptyNum();
}

void GA::mutation()
{
	for (int i = 0; i < data.size(); i++)
	{
		if (random(0.0, 1.0) <= individualMutationRate)//個体突然変異率の計算
		{
#ifdef __ENABLE_SINGLE_POINT_MUTATION__
			int pos = random(0, (int)data[i].x.size() - 1);
			data[i].num[pos] = -1;
#endif
#ifdef __ENABLE_DOUBLE_POINT_MUTATION__
			int del1 = random(0, (int)data[i].num.size() - 1);
			int del2 = random(0, (int)data[i].num.size() - 1);

			data[i].num[del1] = -1;
			data[i].num[del2] = -1;
#endif
#ifdef __ENABLE_SEGMENTAL_MUTATION__
			int a = random(0, 2);
			if (!a)
			{
				int del1 = random(0, (int)data[i].num.size() - 1);
				int del2 = random(del1, (int)data[i].num.size());

				for (int j = del1; j < del2; j++)
				{
					data[i].num[j] = -1;
				}
			}
			else if (a == 2)
			{
				int del1 = random(0, (int)data[i].num.size() - 1);
				int del2 = random(0, (int)data[i].num.size() - 2);

				int temp = data[i].num[del1];
				data[i].num.erase(data[i].num.begin() + del1);
				data[i].num.insert(data[i].num.begin() + del2, temp);
			}
			else
			{
				for (int j = 0; j < data[i].num.size(); j++)
				{
					if (random(0.0, 1.0) <= genomMutationRate)
						data[i].num[j] = -1;
				}
			}
#endif
#if !defined(__ENABLE_SINGLE_POINT_MUTATION__) && !defined(__ENABLE_DOUBLE_POINT_MUTATION__) && !defined(__ENABLE_SEGMENTAL_MUTATION__)
			for (int j = 0; j < data[i].num.size(); j++)
			{
				if (random(0.0, 1.0) <= genomMutationRate)
					data[i].num[j] = -1;
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
	double nowElite = searchRank(0).functionValue;
	if (eliteData.functionValue - nowElite > 0)
	{
		_RPT0(_CRT_WARN, "changed\n");
		isChanged = true;
		alpha = 70;
	}
	else
	{
		if (!isChanged)
		{
			localMinNum++;
			if (localMinNum > 5)
			{
				alpha = alpha / 1.001;
				localMinNum = 0;
			}
			//alpha += 0.1;
		}
	}
	//評価関数が最もいいやつを保存
	if (data[minNum].result < eliteData.result && data[minNum].functionValue > eliteData.functionValue)
		data[minNum] = eliteData;

	calcResult();

	if (enableDisplay)
		displayValues(enableOneLine);
}

void GA::calcResult(bool enableSort)
{
	int maxNum = 0;

	setEmptyNum();

	for (int i = 0; i < data.size(); i++)
	{
		data[i].functionValue = 0;
		for (int j = 0; j < data[i].num.size(); j++)
		{
			double temp = 0;
			for (int k = 0; k < model[data[i].num[j]].point.size(); k++)
			{
				temp += std::pow(model[data[i].num[j]].point[k] - model[data[i].num[j + 1 == data[i].num.size() ? 0 : j + 1]].point[k], 2.0);
			}
			data[i].functionValue += std::pow(temp, 0.5);
		}

		if (data[maxNum].functionValue < data[i].functionValue)//座標の中で最も関数が大きいやつを検索
			maxNum = i;
	}
	resultSumValue = 0;
	double coefficient = 0.001 / data[0].num.size();//評価関数用の定数

	for (int i = 0; i < data.size(); i++)
	{
		//data[i].result = std::exp(-data[i].functionValue*alpha*coefficient);
		data[i].result = 1 / std::pow(data[i].functionValue, alpha);
		//data[i].result = 1 / data[i].functionValue;
		//resultSumValue += data[i].result;
	}

	double normalization = searchRank(0).result;

	for (int i = 0; i < data.size(); i++)
	{
		data[i].result = data[i].result / normalization;
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
		for (int j = 0; j < data_temp[i].num.size(); j++)
		{
			printf_s("%2d,", data_temp[i].num[j]);//デバッグ用
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

void GA::setEmptyNum(void)
{
	for (int i = 0; i < data.size(); i++)
	{
		std::vector<int> noPlacedCity;
		noPlacedCity = cityTemp;

		for (int j = 0; j < data[i].num.size(); j++)
		{
			if (data[i].num[j] != -1)
				noPlacedCity[data[i].num[j]] = -1;
		}
		auto itr = noPlacedCity.begin();
		while (itr != noPlacedCity.end())
		{
			if ((*itr) == -1)
				itr = noPlacedCity.erase(itr);
			else
				itr++;
		}

		std::random_device rnd;
		std::mt19937 engine(rnd());
		std::shuffle(noPlacedCity.begin(), noPlacedCity.end(), engine);

		for (int j = 0, point = 0; j < data[i].num.size(); j++)
		{
			if (data[i].num[j] == -1)
				data[i].num[j] = noPlacedCity[point++];
		}
	}
}

GA::~GA()
{

}
