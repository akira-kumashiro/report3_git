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
	double individualMutationRate = 0.8;//�̓ˑR�ψٗ�
	double genomMutationRate = 0.5;
	bool isChanged = false;
	double alpha = 1;
	int localMinNum = 0;
	std::vector<int> cityTemp;
	int genNum = 1;
public:
	double resultSumValue;//�]���֐��̍��v

	class Data//�f�[�^�i�[�p�N���X
	{
	public:
		std::vector<int> num;//���W
		double functionValue;//�^����ꂽ�֐��̒l
		double result;

		Data(int _var_num) ://�R���X�g���N�^
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

	std::vector<Data> data, prev_data;//����O��Œl��ێ����邽�߂�2��
	Data eliteData;
	std::vector<PointXY> model;
	GA(int _max_genom_list, int _var_num, std::vector<GA::PointXY> model);	//�R���X�g���N�^
	bool selection();//�I��
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
	void setEmptyNum(void);
public:
	~GA();//�f�R���X�g���N�^
};