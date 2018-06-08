#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _var_num, std::vector<GA::PointXY> _model) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//data‚Ì‰Šú‰»
	eliteData(_var_num),
	cityTemp(_var_num)
{
	//‚à‚ç‚Á‚½•Ï”‚ğƒNƒ‰ƒX“à•Ï”‚ÉŠi”[
	model = _model;
	for (int i = 0; i < data[0].num.size(); i++)
	{
		cityTemp[i] = i;
	}
	setEmptyNum();
	prev_data = data;
	calcResult();

	displayValues(false);

}

GA::GA(std::vector<Data> _data, std::vector<GA::PointXY> _model) :
	data(_data),
	eliteData(_data[0].num.size()),
	cityTemp(_data[0].num.size())
{
	//‚à‚ç‚Á‚½•Ï”‚ğƒNƒ‰ƒX“à•Ï”‚ÉŠi”[
	model = _model;

	for (int i = 0; i < data[0].num.size(); i++)
	{
		cityTemp[i] = i;
	}
	prev_data = data;
	calcResult();

	displayValues(false);
}

bool GA::selection()
{
	int max_num = 0;//Å‚à•]‰¿‚Ì—Ç‚¢ŒÂ‘Ì‚Ì”Ô†
	bool ret = isChanged;//Å‚à•]‰¿‚Ì—Ç‚¢ŒÂ‘Ì‚Ì•Ï‰»‚ÌŠÄ‹(ƒfƒoƒbƒO—p)
	isChanged = false;

	eliteData = searchRank(0);//Å‚à•]‰¿‚Ì—Ç‚¢ŒÂ‘Ì‚ğ•Û
	prev_data = data;
	for (int i = 0; i < data.size(); i++)
	{
		double selector = random(0.0, 1.0);//—”‚ğ¶¬
		double needle = 0;//ƒ‹[ƒŒƒbƒg‚Ìj‚ğ¶¬
		int j = 0;
		for (; ; j++)
		{
			needle += (prev_data[j].result / resultSumValue);//ƒ‹[ƒŒƒbƒg‚Ìj‚ğ—”‚Ì’l‚Ü‚Åi‚ß‚é
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
	for (int i = 0; i < data.size(); i += 2)//2ŒÂ‚¸‚ÂŒğ³
	{
		int del1 = random(0, (int)data[0].num.size() - 1);
		int del2 = random(del1, (int)data[0].num.size());
		/*data[i].num[0] = 0;
		data[i + 1].num[0] = 0;*/
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
	genNum += 0.001;
	individualMutationRate = individualMutationRate / std::log(genNum);


	for (int i = 0; i < data.size(); i++)
	{
		if (random(0.0, 1.0) <= individualMutationRate)//ŒÂ‘Ì“Ë‘R•ÏˆÙ—¦‚ÌŒvZ
		{
#ifdef __ENABLE_SINGLE_POINT_MUTATION__
			int pos = random(0, (int)data[i].x.size() - 1);
			data[i].num[pos] = -1;
#else
			/*for (int j = 0; j < data[i].num.size(); j++)
			{
				if (random(0.0, 1.0) <= genomMutationRate)
					data[i].num[j] = -1;
			}*/
			int del1 = random(0, (int)data[i].num.size() - 1);
			int del2 = random(del1, (int)data[i].num.size());
			/*for (int j = del1; j < del2; j++)
			{
				data[i].num[j] = -1;
			}*/
			//int pos1 = data[i].num[del1];
			//int pos2 = data[i].num[del2];

			data[i].num[del1] = -1;
			data[i].num[del2] = -1;
#endif
		}
	}
}

void GA::calc(bool enableDisplay, bool enableOneLine)
{
	int minNum = 0;
	calcResult();
	for (int i = 0; i < data.size(); i++)//•]‰¿ŠÖ”‚ªÅ¬‚Ì“z‚ğŒŸõ
	{
		if (data[i].result < data[minNum].result)
			minNum = i;
	}
	double nowElite = searchRank(0).functionValue;
	if (eliteData.functionValue - nowElite > 1)
	{
		_RPT0(_CRT_WARN, "changed\n");
		isChanged = true;
		alpha = 1;
	}
	else
	{
		if (!isChanged)
		{
			localMinNum++;
			if (localMinNum > 20)
			{
				alpha += 1;
				localMinNum = 0;
			}
		}
	}
	//•]‰¿ŠÖ”‚ªÅ‚à‚¢‚¢‚â‚Â‚ğ•Û‘¶
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
				/*double diff = model[data[i].num[j]].point[k];
				if (j + 1 == data[i].num.size())
				{
					diff -= model[data[i].num[0]].point[k];
				}
				else
				{
					diff -= model[data[i].num[j + 1]].point[k];
				}
				temp += std::pow(diff, 2.0);*/
				temp += std::pow(model[data[i].num[j]].point[k] - model[data[i].num[j + 1 == data[i].num.size() ? 0 : j + 1]].point[k], 2.0);
			}
			data[i].functionValue += std::pow(temp, 0.5);
		}

		if (data[maxNum].functionValue < data[i].functionValue)//À•W‚Ì’†‚ÅÅ‚àŠÖ”‚ª‘å‚«‚¢‚â‚Â‚ğŒŸõ
			maxNum = i;
	}
	double seg = data[maxNum].functionValue;//•]‰¿ŠÖ”‚ÌØ•Ğ‚ğ—^‚¦‚ç‚ê‚½ŠÖ”‚ªÅ‚à‘å‚«‚¢‚â‚Â‚ÉƒZƒbƒg
	double seg2 = searchRank(data[0].num.size() - 2).functionValue - seg;
	resultSumValue = 0;
	double coefficient = 0.001 / data[0].num.size();//•]‰¿ŠÖ”—p‚Ì’è”

	for (int i = 0; i < data.size(); i++)
	{
		//data[i].result = std::exp(-data[i].functionValue*alpha*coefficient) / coefficient;
		//data[i].result = seg2 == 0 ? 0 : (data[i].functionValue - seg) / seg2 / coefficient;//—^‚¦‚ç‚ê‚½ŠÖ”‚Ì’l‚©‚çØ•Ğ‚Åİ’è‚µ‚½’l‚ğˆø‚¢‚Ä2æ‚·‚é¨—^‚¦‚ç‚ê‚½ŠÖ”‚Ì’l‚ª¬‚³‚¢‚â‚Â‚ª‹­‚­‚È‚é
		data[i].result = 1 / std::pow(data[i].functionValue, 0.5);
		//data[i].result = 1 / data[i].functionValue;
		/*if (model[data[i].num[1]].point[0] - model[data[i].num[0]].point[0] > 0)
		{
			data[i].result *= (1 + coefficient * 100);
			if (model[data[i].num[1]].point[1] - model[data[i].num[0]].point[1] > 0)
				data[i].result *= (1 + coefficient * 100);
		}*/
		resultSumValue += data[i].result;
	}
	if (enableSort)
		std::sort(data.begin(), data.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });
}

int GA::random(int min, int max)
{
	//—”‚Ìİ’è
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
			printf_s("%2d,", data_temp[i].num[j]);//ƒfƒoƒbƒO—p
		}
		printf_s(" \t f(x,y)=%10.7lf\t Result=%10.7lf\n", data_temp[i].functionValue, data_temp[i].result);
	}
}

GA::Data GA::searchRank(int num)//•]‰¿‚ª‚¢‚¢‡
{
	std::vector<Data> data_temp = data;
	std::sort(data_temp.begin(), data_temp.end(), [](const Data& x, const Data& y) { return x.functionValue < y.functionValue; });
	return data_temp[num];
}

void GA::setEmptyNum(void)
{
	for (int i = 0; i < data.size(); i++)
	{
		std::vector<int> noPlacedCity, noPlacedCityTemp;
		noPlacedCity = cityTemp;

		for (int j = 0; j < data[i].num.size(); j++)
		{
			if (data[i].num[j] != -1)
				noPlacedCity[data[i].num[j]] = -1;
		}

		for (int j = 0; j < noPlacedCity.size(); j++)
		{
			if (noPlacedCity[j] != -1)
				noPlacedCityTemp.push_back(noPlacedCity[j]);
		}
		noPlacedCity = noPlacedCityTemp;

		for (int j = 0; j < noPlacedCity.size(); j++)
		{
			int point = random(0, (int)noPlacedCityTemp.size() - 1);
			noPlacedCity[j] = noPlacedCityTemp[point];
			noPlacedCityTemp.erase(noPlacedCityTemp.begin() + point);
		}

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