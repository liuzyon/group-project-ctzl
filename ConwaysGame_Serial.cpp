#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <vector>

using namespace std;

//Note that this is a serial implementation with a periodic grid
vector<vector<bool>> grid, new_grid;
int imax, jmax;
int max_steps = 100;

int num_neighbours(int ii, int jj)
{
	int ix, jx;
	int cnt = 0;
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			if (i != 0 || j != 0)
			{
				ix = (i + ii + imax) % imax;
				jx = (j + jj + jmax) % jmax;
				if (grid[ix][jx]) cnt++;
			}
	return cnt;
}

// 将每一次迭代后的状态输出
void grid_to_file(int it)
{
	stringstream fname;
	fstream f1;
	fname << "output" << "_" << it << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
			f1 << grid[i][j] << "\t";
		f1 << endl;
	}
	f1.close();
}

void do_iteration(void)
{
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
		{
			new_grid[i][j] = grid[i][j];
			int num_n = num_neighbours(i, j);   //  计算周围存活细胞数
			if (grid[i][j]) // 当细胞为存活状态时
			{
				if (num_n != 2 && num_n != 3)   // 如果周围的存活元胞数小于2或大于3，该元胞变成死亡状态
					new_grid[i][j] = false;
			}
			else if (num_n == 3) new_grid[i][j] = true; // 当前元胞为死亡状态时，如果周围有3个存活元胞时，该元胞变成存活状态
		}
	grid.swap(new_grid);    // 交换grid容器和new_gird容器内容
}

int main(int argc, char *argv[])
{
	srand(time(NULL));
	imax = 100;
	jmax = 100;
	// 改变容器大小和容量
	grid.resize(imax, vector<bool>(jmax));
	new_grid.resize(imax, vector<bool>(jmax));

	//set an initial random collection of points - You could set an initial pattern
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++) grid[i][j] = (rand() % 2);

	for (int n = 0; n < max_steps; n++)
	{
		cout << "it: " << n << endl;
		do_iteration(); // 执行游戏规则迭代
		grid_to_file(n);
	}

	return 0;
}
