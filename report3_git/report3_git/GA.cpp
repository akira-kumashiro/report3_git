#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _var_num, std::vector<GA::CityData> _model) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//data�̏�����
	eliteData(_var_num)
{
	//��������ϐ����N���X���ϐ��Ɋi�[
	model = _model;

	/*for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].num.size(); j++)
		{
			data[i].num[j] = random(varMin[j], varMax[j]);//��`�q�̏����ݒ�
		}
	}*/
	setEmptyNum();
	prev_data = data;
	calcResult();

	displayValues(false);
}

bool GA::selection()
{
	int max_num = 0;//�ł��]���̗ǂ��̂̔ԍ�
	bool ret = isChanged;//�ł��]���̗ǂ��̂̕ω��̊Ď�(�f�o�b�O�p)
	isChanged = false;

	eliteData = searchRank(0);//�ł��]���̗ǂ��̂�ێ�
	prev_data = data;
	for (int i = 0; i < data.size(); i++)
	{
		double selector = random(0.0, 1.0);//�����𐶐�
		double needle = 0;//���[���b�g�̐j�𐶐�
		int j = 0;
		for (; ; j++)
		{
			needle += (prev_data[j].result / resultSumValue);//���[���b�g�̐j�𗐐��̒l�܂Ői�߂�
			if (needle > selector)
				break;
			if (j == (data.size() - 1))
				break;
		}
		data[i] = prev_data[j];
	}
	return ret;
}

/*void GA::blxAlphaCrossover()
{
	prev_data = data;

	for (int i = 0; i < data.size(); i += 2)//2������
	{
		for (int j = 0; j < data[i].num.size(); j++)
		{
			double ave = (data[i].num[j] + data[i + 1].num[j]) / 2;
			double length = std::abs((data[i].num[j] - data[i + 1].num[j]));

			data[i].num[j] = random(ave - length * (1 / 2 + alpha), ave + length * (1 / 2 + alpha));
			data[i + 1].num[j] = random(ave - length * (1 / 2 + alpha), ave + length * (1 / 2 + alpha));
		}
	}
}*/

void GA::pmxCrossover()
{
	prev_data = data;

	int del1 = random(0, (int)data[0].num.size() - 1);
	int del2 = random(del1, (int)data[0].num.size());

	//Step1
	for (int i = 0; i < data.size(); i += 2)//2������
	{
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
				{
					data[i].num[j] = prev_data[i].num[j];
				}
			}
		}
	}
	//Steo3
	/*for (int i = 0; i < data.size(); i++)
	{
		std::vector<int> noPlacedCity, noPlacedCityTemp;
		for (int j = 0; j < data[i].num.size(); j++)
		{
			noPlacedCity[j] = j;
		}

		for (int j = 0; j < data[i].num.size(); j++)
		{
			if (data[i].num[j] != -1)
			{
				noPlacedCity[j] = -1;
			}
		}

		for (int j = 0; j < data[i].num.size(); j++)
		{
			if (noPlacedCity[j] == -1)
			{
				noPlacedCity.erase(noPlacedCity.begin() + j);
				j--;
			}
		}
		noPlacedCityTemp = noPlacedCity;

		for (int j = 0; j < noPlacedCity.size(); j++)
		{
			int point = random(0, (int)noPlacedCityTemp.size());
			noPlacedCity[j] = noPlacedCityTemp[point];
			noPlacedCityTemp.erase(noPlacedCityTemp.begin() + point);
		}

		for (int j = 0; j < data[i].num.size(); j++)
		{
			int point = 0;
			if (data[i].num[j] == -1)
			{
				data[i].num[j] = noPlacedCity[point];
			}
		}
	}*/
	setEmptyNum();
}

void GA::mutation()
{
	for (int i = 0; i < data.size(); i++)
	{
		if (random(0.0, 1.0) <= individualMutationRate)//�̓ˑR�ψٗ��̌v�Z
		{
#ifdef __ENABLE_SINGLE_POINT_MUTATION__
			int pos = random(0, (int)data[i].x.size() - 1);
			data[i].num[pos] = -1;
#else
			for (int j = 0; j < data[i].num.size(); j++)
			{
				if (random(0.0, 1.0) <= genomMutationRate)
					data[i].num[j] = -1;
			}
#endif
		}
	}
	setEmptyNum();
}

void GA::calc(bool enableDisplay, bool enableOneLine)
{
	int minNum = 0;
	calcResult();
	for (int i = 0; i < data.size(); i++)//�]���֐����ŏ��̓z������
	{
		if (data[i].result < data[minNum].result)
			minNum = i;
	}
	/*if (searchRank(0).functionValue - eliteData.functionValue > 1)
		isChanged = true;*/
		//�]���֐����ł��������ۑ�
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
		data[i].functionValue = 0;
		for (int j = 0; j < data[i].num.size(); j++)
		{
			double temp = 0;
			for (int k = 0; k < model[0].point.size(); k++)
			{
				temp += std::pow(model[data[i].num[j]].point[k] - model[data[i].num[j + 1 == data[i].num.size() ? 0 : j + 1]].point[k], 2.0);
			}
			data[i].functionValue += std::pow(temp, 0.5);
		}

		if (data[maxNum].functionValue < data[i].functionValue)//���W�̒��ōł��֐����傫���������
			maxNum = i;
	}
	double seg = data[maxNum].functionValue;//�]���֐��̐ؕЂ�^����ꂽ�֐����ł��傫����ɃZ�b�g
	double seg2 = searchRank(data[0].num.size() - 2).functionValue - seg;
	resultSumValue = 0;
	double coefficient = 0.001 / data[0].num.size();//�]���֐��p�̒萔

	for (int i = 0; i < data.size(); i++)
	{
		data[i].result = seg2 == 0 ? 0 : (data[i].functionValue - seg) / seg2 / coefficient;//�^����ꂽ�֐��̒l����ؕЂŐݒ肵���l��������2�悷�遨�^����ꂽ�֐��̒l����������������Ȃ�

		resultSumValue += data[i].result;
	}
	if (enableSort)
		std::sort(data.begin(), data.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });
}

int GA::random(int min, int max)
{
	//�����̐ݒ�
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
			printf_s("%2d,", data_temp[i].num[j]);//�f�o�b�O�p
		}
		printf_s(" \t f(x,y)=%10.7lf\t Result=%10.7lf\n", data_temp[i].functionValue, data_temp[i].result);
	}
}

GA::Data GA::searchRank(int num)//�]����������
{
	std::vector<Data> data_temp = data;
	std::sort(data_temp.begin(), data_temp.end(), [](const Data& x, const Data& y) { return x.functionValue < y.functionValue; });
	return data_temp[num];
}

void GA::setEmptyNum(void)
{
	std::vector<int> cityTemp(data[0].num.size());
	for (int i = 0; i < data[0].num.size(); i++)
	{
		cityTemp[i] = i;
	}

	for (int i = 0; i < data.size(); i++)
	{
		std::vector<int> noPlacedCity, noPlacedCityTemp;
		noPlacedCity = cityTemp;

		for (int j = 0; j < data[i].num.size(); j++)
		{
			if (data[i].num[j] != -1)
			{
				noPlacedCity[data[i].num[j]] = -1;
			}
		}

		for (int j = 0; j < data[i].num.size(); j++)
		{
			if (noPlacedCity.size() <= j)
			{
				break;
			}
			else if (noPlacedCity[j] == -1)
			{
				noPlacedCity.erase(noPlacedCity.begin() + j);
				j--;
			}
		}
		noPlacedCityTemp = noPlacedCity;

		for (int j = 0; j < noPlacedCity.size(); j++)
		{
			int point = random(0, (int)noPlacedCityTemp.size() - 1);
			noPlacedCity[j] = noPlacedCityTemp[point];
			noPlacedCityTemp.erase(noPlacedCityTemp.begin() + point);
		}

		for (int j = 0, point = 0; j < data[i].num.size(); j++)
		{
			if (data[i].num[j] == -1)
			{
				data[i].num[j] = noPlacedCity[point++];
			}
		}
	}
}

GA::~GA()
{

}