#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <chrono>

using namespace std;

//Note that this is a serial implementation with a periodic grid
vector<vector<bool>> grid, new_grid;
int imax, jmax;
int max_steps = 100;

// how to calculate the numbers of neighbours.
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

// write in the status of each cell. 
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

void grid_to_ppm(int it, int mypit) {
	if (it % mypit != 0 && it != max_steps - 1) {

	}
	else {
		stringstream fname;
		fstream f1;
		fname << "output_image" << "_" << it << ".ppm";
		f1.open(fname.str().c_str(), ios_base::out);
		f1 << "P3" << endl;
		f1 << imax << " " << jmax << endl;
		f1 << "255" << endl;
		int r = 0;
		int g = 0;
		int b = 0;
		for (int i = 0; i < imax; i++) {
			for (int j = 0; j < jmax; j++) {
				g = grid[i][j] * 255;
				f1 << r << " " << g << " " << b << " ";
			}
			f1 << endl;
		}
		f1.close();
	}


}



//status of the cell
void do_iteration(void)
{
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
		{
			new_grid[i][j] = grid[i][j];
			int num_n = num_neighbours(i, j);
			if (grid[i][j])
			{
				if (num_n != 2 && num_n != 3)
					new_grid[i][j] = false;
			}
			else if (num_n == 3) new_grid[i][j] = true;
		}
	grid.swap(new_grid);
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	imax = 100;
	jmax = 100;
	grid.resize(imax, vector<bool>(jmax));
	new_grid.resize(imax, vector<bool>(jmax));

	clock_t start_serial = clock();

	//set an initial random collection of points - You could set an initial pattern
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++) grid[i][j] = (rand() % 2);

	for (int n = 0; n < max_steps; n++)
	{
		cout << "it: " << n << endl;
		do_iteration();
		// grid_to_file(n);
		grid_to_ppm(n, 20);
	}

	clock_t end_serial = clock();
	cerr << "Serial time of " << imax << " x " << jmax << " with " << max_steps << " generations(Seconds): " << (double)(end_serial - start_serial) / CLOCKS_PER_SEC << endl;

	return 0;
}