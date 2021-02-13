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
vector<vector<bool>> grid, new_grid, my_grid;
int imax, jmax;
int max_steps = 100;

// a new shared variable
int iter_num = 0;

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



//These codes use 9 thread to find the
//neighbors locations and judge the status of then
//return the number of live neighbors.
// According to the test these seem to become
//slower than serial codes.
//So we ban this method.
//int num_neighbours(int ii, int jj)
//{
//	int cnt = 0;
//#pragma omp parallel num_threads(9)
//	{
//		int i = (omp_get_thread_num() / 3) - 1;
//		int j = (omp_get_thread_num() % 3) - 1;
//
//		if (i != 0 || j != 0) {
//			int ix = (i + ii + imax) % imax;
//			int jx = (j + jj + jmax) % jmax;
//			if (grid[ix][jx])
//#pragma omp atomic
//				cnt++;
//		}
//	}
//
//	return cnt;
//}

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
    fname << "output_image" << "_" << int(max_steps / 2) << ".ppm";
    f1.open(fname.str().c_str(), ios_base::out);
    f1 << "P3" << endl;
    f1 << imax << " " << jmax << endl;
    f1 << "255" << endl;
    int r = 0;
    int g = 0;
    int b = 0;
    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {
            g = my_grid[i][j] * 255;
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
            //cout << "iteration it: " << n << endl;
//            cout << "iteration it: " << n << " group numebrs: " <<  omp_get_num_threads() << endl;
#pragma omp critical    // only one thread can modify the array at one time
                do_iteration();
// save the mid generation grid to a vector
#pragma omp atomic
				iter_num++;

				if (iter_num == int(max_steps / 2)) {
					my_grid.assign(grid.begin(), grid.end());
				}
        }
    }
}

void print(int print_wait_times) {
//    omp_set_num_threads(thread_number);
//#pragma omp parallel
//    {
//#pragma omp for


	// print m times， random print in all iterations(for iter_num<mid number, no use iteration)
	//when jmax and i max increase m should be larger
	// the value of m depends on the iteration time. when iter = mid number.
	for (int m = 0; m < print_wait_times; m++)
	{
		cout << "print it: " << iter_num << "$" << m << endl;
		//            cout << "print it: " << m << " group numebrs: " <<  omp_get_num_threads() << endl;

		grid_to_ppm(iter_num);

		if (iter_num > int(max_steps / 2)) {
			cout << "\n break with iter_num" << iter_num;
			break;

		}

	}
//    }
}

int main(int argc, char* argv[])
{
    std::cout << "MAX THREAD NUMBER: " << omp_get_max_threads() << endl;
    int MAX_THREADS = omp_get_max_threads();
    omp_set_nested(1);
    omp_set_num_threads(8);
	srand(time(NULL));
	imax = 2000;
	jmax = 2000;
	int print_wait_times = imax + jmax;//请仔细思考这个数字的选取，与do_iteration()速度有关。可大不可小。
	grid.resize(imax, vector<bool>(jmax));
	new_grid.resize(imax, vector<bool>(jmax));
	my_grid.resize(imax, vector<bool>(jmax));

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
            print(print_wait_times);
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