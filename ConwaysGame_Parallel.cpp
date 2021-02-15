#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <chrono>
#include <omp.h>

using namespace std;

// Note that this is a serial implementation with a periodic grid
vector<vector<bool>> grid, new_grid, my_grid_start, my_grid;
int imax, jmax;
int max_steps;

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

// print ppm file for start middle and final generations
void grid_to_ppm(int it) {
    stringstream fname;
    fstream f1;
	fname << "CTZL_" << imax << " x " << jmax << "_image_" << it << ".ppm";
    
    f1.open(fname.str().c_str(), ios_base::out);
    f1 << "P3" << endl;
    f1 << imax << " " << jmax << endl;
    f1 << "255" << endl;
    int r = 0;
    int g = 0;
    int b = 0;
    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {
			// in our case, this produces the mid point of all iterations
			if (it != 0 && it!= max_steps)
			{
				g = my_grid[i][j] * 255;
			}
			// for the initial case
            else if (it == 0)
			{
				g = my_grid_start[i][j] * 255;
			}
			// for the final status after all genenrations
			else if (it == max_steps)
			{
				g = grid[i][j] * 255;
			}
            f1 << r << " " << g << " " << b << " ";
        }
        f1 << endl;
    }
    f1.close();
}



// status of the cell
void do_iteration(void)
{
// Two iterations are collapsed into one large iteration
#pragma omp parallel for schedule(dynamic) collapse(2)
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
        {
// method 2: transfer to one iteration manually
//#pragma omp parallel for schedule(dynamic)
//    for (int ij=0; ij<imax*jmax; ++ij)
//	      {
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
// only one thread can modify the array at one time
#pragma omp critical    
                do_iteration();

// protect data.			
#pragma omp atomic
				iter_num++;
// save the mid generation grid to a vector
				if (iter_num == int(max_steps / 2)) {
					my_grid.assign(grid.begin(), grid.end());
				}
        }
    }
}

void print() {

	// begin by printing the initial state
	// this will print the saved my_grid_start which will not be affected as iteration goes
	//Put it inside parallel thread to save time, even though cost memory.
	grid_to_ppm(0);

	// using an empty while loop to wait for the iter_num to reach half of the max
	while (iter_num < int(max_steps / 2))
	{

	}
	// if max_steps / 2 is reached, print this saved status which is stored as my_grid
	grid_to_ppm(int(max_steps / 2));

	// then wait for all iterations completed and the last status will be printed outside at the end

}

int main(int argc, char* argv[])
{
    int MAX_THREADS = omp_get_max_threads();
    cout << "Max available threads number: " << MAX_THREADS << endl;

    int row_size = 1000;
    cout << "<Step 1>\nThe row size of Conways (Please input an integer value):" << endl;
    cin >> row_size;

    int column_size = 1000;
    cout << "<Step 2>\nThe column size of Conways (Please input an integer value):" << endl;
    cin >> column_size;

    int step_times = 100;
    cout << "<Step 3>\nThe iteration steps (Please input an integer value):" << endl;
    cin >> step_times;

    omp_set_nested(1);
    omp_set_num_threads(MAX_THREADS);
	srand(time(NULL));
	//initialise the size of grid and iteration times.
	imax = row_size;
	jmax = column_size;
	max_steps = step_times;
	grid.resize(imax, vector<bool>(jmax));
	new_grid.resize(imax, vector<bool>(jmax));
	my_grid_start.resize(imax, vector<bool>(jmax));
	my_grid.resize(imax, vector<bool>(jmax));

	std::cout << "<Processing>\nRunning... \nPlease wait for the result." << endl;

	// start time - elapsed wall clock time in seconds
    double start_time = omp_get_wtime(); 

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

    // store the current grid to my_grid_start and print it in parallel with the iteration
	my_grid_start.assign(grid.begin(), grid.end());

//use parallel threads in two sections.
#pragma omp parallel
{
#pragma omp sections nowait
    {
// one thread to run iteration part
#pragma omp section
        {
            // gerente a thread group(max_size-1) to do iteration
            iteration(MAX_THREADS-1);
        }
// one thread to run print part
#pragma omp section
        {
            // one thread  to print
            print();
        }
	}
}


	// when iteration ends print the last state of the grid
	grid_to_ppm(max_steps);

	//end time
    double end_time = omp_get_wtime();
	cerr << "<Result>\nParallel time of " << imax << " x " << jmax << " with " << max_steps << " generations(Seconds): " << (double)(end_time - start_time) << endl;

	return 0;
}