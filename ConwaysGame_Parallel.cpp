#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <chrono>
#include <omp.h>

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

void grid_to_ppm(int it) {

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

//	if (it % mypit != 0 && it != max_steps - 1) {
//        return;
//	}
//	else {
//		stringstream fname;
//		fstream f1;
//		fname << "output_image" << "_" << it << ".ppm";
//		f1.open(fname.str().c_str(), ios_base::out);
//		f1 << "P3" << endl;
//		f1 << imax << " " << jmax << endl;
//		f1 << "255" << endl;
//		int r = 0;
//		int g = 0;
//		int b = 0;
//		for (int i = 0; i < imax; i++) {
//			for (int j = 0; j < jmax; j++) {
//				g = grid[i][j] * 255;
//				f1 << r << " " << g << " " << b << " ";
//			}
//			f1 << endl;
//		}
//		f1.close();
//	}


}



//status of the cell
void do_iteration(void)
{
//    omp_set_num_threads(4);
    omp_set_nested(1);
//  Two iterations are collapsed into one large iteration
#pragma omp parallel for schedule(dynamic) collapse(2)
//#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < imax; i++)
//#pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < jmax; j++)
        {

// method 2: transfer to one iteration manually
//#pragma omp parallel for schedule(dynamic)
//        for (int ij=0; ij<imax*jmax; ++ij)
//		{
//            int i = ij / jmax;
//            int j = ij % jmax;

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

void iteration(int thread_number) {
    omp_set_num_threads(thread_number);
#pragma omp parallel
    {
#pragma omp for
        for (int n = 0; n < max_steps; n++)
        {
            cout << "iteration it: " << n << endl;
//            cout << "iteration it: " << n << " group numebrs: " <<  omp_get_num_threads() << endl;
#pragma omp critical    // only one thread can modify the array at one time
                do_iteration();
        }
    }
}

void print() {
//    omp_set_num_threads(thread_number);
//#pragma omp parallel
//    {
//#pragma omp for
        for (int m = 0; m < 10; m++)    // print m times， random print in all iterations
        {
            cout << "print it: " << m << endl;
//            cout << "print it: " << m << " group numebrs: " <<  omp_get_num_threads() << endl;
            grid_to_ppm(m);
        }
//    }
}

int main(int argc, char* argv[])
{
    std::cout << "MAX THREAD NUMBER: " << omp_get_max_threads() << endl;
    int MAX_THREADS = omp_get_max_threads();
    omp_set_nested(1);
    omp_set_num_threads(MAX_THREADS);
	srand(time(NULL));
	imax = 2000;
	jmax = 2000;
	grid.resize(imax, vector<bool>(jmax));
	new_grid.resize(imax, vector<bool>(jmax));

    double start_time = omp_get_wtime(); //start time - elapsed wall clock time in seconds

	// set an initial random collection of points - You could set an initial pattern
	// change it to a single for loop and parallelize it
	// increase the speed of consctruting larger random grids
#pragma omp parallel for schedule(dynamic)
	for (int ij = 0; ij < imax * jmax; ij++)
	{
		// get row index
		int i = ij / jmax;
		// get column index
        int j = ij % jmax;
		grid[i][j] = (rand() % 2);
	}

#pragma omp parallel
{
#pragma omp sections
    {
        // one thread to run iteration part
#pragma omp section
        {
//            std::cout << "max1: " << omp_get_max_threads();
            // gerente a thread group(max_size-2) to do iteration
            iteration(MAX_THREADS-2);
        }
        // one thread to run print
#pragma omp section
        {
//            std::cout << "max2: " << omp_get_max_threads();
            // one thread  to print
            print();
        }
	}
}

//#pragma omp parallel for
//	for (int n = 0; n < max_steps; n++)
//	{
//		cout << "it: " << n << endl;
//		do_iteration();
//		// grid_to_file(n);
//#pragma omp task
//		grid_to_ppm(n, 20);
//	}

    double end_time = omp_get_wtime();//end time
	cerr << "Parallel time of " << imax << " x " << jmax << " with " << max_steps << " generations(Seconds): " << (double)(end_time - start_time) << endl;

	return 0;
}