//GitHub
//https://github.com/akira-kumashiro/report3_git

#include "stdafx.h"
#include "GA.h"
#include <conio.h>

#define MAX_GENERATION 30000
#define MAX_GENOM_LIST 50
#define VAR_NUM 28
#define __ENABLE_MUTATION__
/*
1    1150.0  1760.0
2     630.0  1660.0
3      40.0  2090.0
4     750.0  1100.0
5     750.0  2030.0
6    1030.0  2070.0
7    1650.0   650.0
8    1490.0  1630.0
9     790.0  2260.0
10     710.0  1310.0
11     840.0   550.0
12    1170.0  2300.0
13     970.0  1340.0
14     510.0   700.0
15     750.0   900.0
16    1280.0  1200.0
17     230.0   590.0
18     460.0   860.0
19    1040.0   950.0
20     590.0  1390.0
21     830.0  1770.0
22     490.0   500.0
23    1840.0  1240.0
24    1260.0  1500.0
25    1280.0   790.0
26     490.0  2130.0
27    1460.0  1420.0
28    1260.0  1910.0
29     360.0  1980.0
*/

int main()
{
	double varX[] = { 1150.0,630.0,40.0,750.0,750.0,1030.0,1650.0,1490.0,790.0,710.0,840.0,1170.0,970.0,510.0,750.0,1280.0,230.0,460.0,1040.0,590.0 , 830.0 , 490.0 ,1840.0 ,  1260.0 ,  1280.0 ,    490.0 , 1460.0 ,  1260.0 , 	  360.0 , };
	double varY[] = { 1760.0,1660.0,2090.0,1100.0,2030.0,2070.0, 650.0,1630.0,2260.0, 1310.0, 550.0, 2300.0,1340.0,  700.0,900.0,1200.0,  590.0,860.0, 950.0,1390.0, 1770.0,  500.0,1240.0,1500.0,790.0, 2130.0,1420.0,1910.0,1980.0 };
	std::vector<GA::CityData> cityData(VAR_NUM, 2);

	for (int i = 0; i < 29; i++)
	{
		cityData[i].cityNum = i;
		cityData[i].point[0] = varX[i];
		cityData[i].point[1] = varY[i];
	}

	//遺伝的アルゴリズム諸関数をまとめたクラスの宣言
	GA ga(MAX_GENOM_LIST, VAR_NUM, cityData);

	for (int i = 0; i <= MAX_GENERATION; i++)//メインのループ
	{
		bool change = ga.selection();//選択

		ga.blxAlphaCrossover();//交叉
#ifdef __ENABLE_MUTATION__
		ga.mutation();//突然変異
#endif
		if (i % (MAX_GENERATION / 10) == 0 || change)
		{
			std::cout << "i=" << std::to_string(i) << std::endl;
			ga.calc(true, i % (MAX_GENERATION / 10) != 0);//評価関数の計算
		}
		else
		{
			ga.calc(false);//評価関数の計算
		}
	}

	while (1)
	{
		if (_kbhit() && _getch() == 27)
			break;
	}
	return 0;
}
