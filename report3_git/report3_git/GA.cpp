#include "stdafx.h"
#include "GA.h"

GA::GA(int _max_genom_list, int _var_num, std::vector<GA::CityData> _model) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//data�̏�����
	eliteData(_var_num)
{
	//��������ϐ����N���X���ϐ��Ɋi�[
	model = _model;

	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].x.size(); j++)
		{
			data[i].x[j] = random(varMin[j], varMax[j]);//��`�q�̏����ݒ�
		}
	}
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

void GA::blxAlphaCrossover()
{
	prev_data = data;

	for (int i = 0; i < data.size(); i += 2)//2������
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

	for (int i = 0; i < data.size(); i += 2)//2������
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
		if (random(0.0, 1.0) <= individualMutationRate)//�̓ˑR�ψٗ��̌v�Z
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
	for (int i = 0; i < data.size(); i++)//�]���֐����ŏ��̓z������
	{
		if (data[i].result < data[minNum].result)
			minNum = i;
	}
	if (searchRank(0).functionValue - eliteData.functionValue > 0.00001)
		isChanged = true;
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
		data[i].functionValue = std::sin(data[i].x[0] + data[i].x[1]) + std::pow((data[i].x[0] - data[i].x[1]), 2.0) - 1.5*data[i].x[0] + 2.5*data[i].x[1] + 1;//�^����ꂽ�֐�

		if (data[maxNum].functionValue < data[i].functionValue)//���W�̒��ōł��֐����傫���������
			maxNum = i;
	}
	double seg = data[maxNum].functionValue;//�]���֐��̐ؕЂ�^����ꂽ�֐����ł��傫����ɃZ�b�g
	double seg2 = searchRank(data[0].x.size() - 2).functionValue - seg;
	resultSumValue = 0;
	double coefficient = 0.001 / data[0].x.size();//�]���֐��p�̒萔

	for (int i = 0; i < data.size(); i++)
	{
		bool flag = true;

		for (int j = 0; j < data[i].x.size(); j++)
		{
			if (data[i].x[j] > varMax[j] || data[i].x[j] < varMin[j])//���W����O�ɂ����̏���
				flag = false;
		}
		data[i].result = seg2 == 0 ? 0 : (data[i].functionValue - seg) / seg2 / coefficient;//�^����ꂽ�֐��̒l����ؕЂŐݒ肵���l��������2�悷�遨�^����ꂽ�֐��̒l����������������Ȃ�
																							//data[i].result = std::abs(data[i].functionValue - seg);

		if (!flag)//��O�ɏo����̏���
			data[i].result *= coefficient;
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
		for (int j = 0; j < data_temp[i].x.size(); j++)
		{
			printf_s("%10.7lf,", data_temp[i].x[j]);//�f�o�b�O�p
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

GA::~GA()
{

}