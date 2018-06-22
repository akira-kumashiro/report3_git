//GitHub
//https://github.com/akira-kumashiro/report3_git

#include "stdafx.h"
#include "GA.h"
#include <conio.h>

#define MAX_GENERATION 15000
#define MAX_GENOM_LIST 50
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
	std::vector<GA::PointXY> cityData{
		GA::PointXY(1150.0,1760.0),
		GA::PointXY(630.0,1660.0),
		GA::PointXY(40.0,2090.0),
		GA::PointXY(750.0,1100.0),
		GA::PointXY(750.0,2030.0),
		GA::PointXY(1030.0,2070.0),
		GA::PointXY(1650.0,650.0),
		GA::PointXY(1490.0,1630.0),
		GA::PointXY(790.0,2260.0),
		GA::PointXY(710.0,1310.0),
		GA::PointXY(840.0,550.0),
		GA::PointXY(1170.0,2300.0),
		GA::PointXY(970.0,1340.0),
		GA::PointXY(510.0,700.0),
		GA::PointXY(750.0,900.0),
		GA::PointXY(1280.0,1200.0),
		GA::PointXY(230.0,590.0),
		GA::PointXY(460.0,860.0),
		GA::PointXY(1040.0,950.0),
		GA::PointXY(590.0,1390.0),
		GA::PointXY(830.0,1770.0),
		GA::PointXY(490.0,500.0),
		GA::PointXY(1840.0,1240.0),
		GA::PointXY(1260.0,1500.0),
		GA::PointXY(1280.0,790.0),
		GA::PointXY(490.0,2130.0),
		GA::PointXY(1460.0,1420.0),
		GA::PointXY(1260.0,1910.0),
		GA::PointXY(360.0,1980.0)
	};

	//遺伝的アルゴリズム諸関数をまとめたクラスの宣言
	std::vector<GA> ga(4, GA(MAX_GENOM_LIST, cityData.size(), cityData));
	GA gaJointed(MAX_GENOM_LIST * 2, cityData.size(), cityData);

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j <= MAX_GENERATION * 2; j++)//メインのループ
		{
			bool change = ga[i].selection();//選択

			ga[i].pmxCrossover();//交叉
#ifdef __ENABLE_MUTATION__
			ga[i].mutation();//突然変異
#endif
			if (j % (MAX_GENERATION / 5) == 0 || change)
			{
				std::cout << "i=" << std::to_string(j) << std::endl;
				ga[i].calc(true, j % (MAX_GENERATION / 5) != 0);//評価関数の計算
			}
			else
			{
				ga[i].calc(false);//評価関数の計算
			}
		}
		ga[i].calcResult(true);
		for (int j = 0; j < MAX_GENOM_LIST / 2; j++)//上位半分を結合クラスに
		{
			gaJointed.data[j + MAX_GENOM_LIST * i / 2] = ga[i].data[j + MAX_GENOM_LIST / 2];
		}
	}

	gaJointed.calc(true);

	for (int i = 0; i <= MAX_GENERATION; i++)//メインのループ
	{
		bool change = gaJointed.selection();//選択

		gaJointed.pmxCrossover();//交叉
#ifdef __ENABLE_MUTATION__
		gaJointed.mutation();//突然変異
#endif
		if (i % (MAX_GENERATION / 10) == 0 || change)
		{
			std::cout << "i=" << std::to_string(i) << std::endl;
			gaJointed.calc(true, i % (MAX_GENERATION / 10) != 0);//評価関数の計算
		}
		else
		{
			gaJointed.calc(false);//評価関数の計算
		}
	}

	while (1)
	{
		if (_kbhit() && _getch() == 27)
			break;
	}
	return 0;
}
