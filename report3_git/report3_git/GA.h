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
	double individualMutationRate = 0.05;//�̓ˑR�ψٗ�
	double genomMutationRate = 0.1;
	double alpha = 0.3;
	std::vector<double> varMax, varMin;//�ϐ��̍ŏ��l�E�ő�l
	bool isChanged = false;
public:
	double resultSumValue;//�]���֐��̍��v

	class Data//�f�[�^�i�[�p�N���X
	{
	public:
		std::vector<int> num;//���W
		double functionValue;//�^����ꂽ�֐��̒l
		double result;

		Data(int _var_num)//�R���X�g���N�^
		{
			num.resize(_var_num);//isIncluded�̔z��̒����̐ݒ�
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

	std::vector<Data> data, prev_data;//����O��Œl��ێ����邽�߂�2��
	Data eliteData;
	std::vector<CityData> model;
	GA(int _max_genom_list, int _var_num, std::vector<GA::CityData> model);	//�R���X�g���N�^
	bool selection();//�I��
	void blxAlphaCrossover();
	void pmxCrossover();
	void mutation();//�ˑR�ψ�
	void calc(bool enableDisplay, bool enableOneLine = false);//�]���֐��̌v�Z
private:
	void calcResult(bool enableSort = false);
	int random(int min, int max);
	double random(int min, double max);
	double random(double min, int max);
	double random(double min, double max);
	void displayValues(bool enableOneLine);
	Data searchRank(int num);
public:
	~GA();//�f�R���X�g���N�^
};