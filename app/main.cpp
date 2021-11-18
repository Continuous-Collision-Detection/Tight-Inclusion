#include <iostream>

#include <vector>

#include <array>
#include <tight_inclusion/inclusion_ccd.hpp>
#include <tight_inclusion/Rational.hpp>
#include <fstream>
#include <tight_inclusion/Timer.hpp>
#include <stdexcept>
#include <cstdlib>

using namespace inclusion_ccd;
void case_check()
{
#ifdef TIGHT_INCLUSION_DOUBLE
    std::cout << "using double precision values as inputs" << std::endl;
#else
    std::cout << "using single precision values as inputs" << std::endl;
#endif
    const Vector3d a0s(0.1, 0.1, 0.1);
    const Vector3d a1s(0, 0, 1);
    const Vector3d a0e(1, 0, 1);
    const Vector3d a1e(0, 1, 1);
    const Vector3d b0s(0.1, 0.1, 0.1);
    const Vector3d b1s(0, 0, 0);
    const Vector3d b0e(0, 1, 0);
    const Vector3d b1e(1, 0, 0);

    bool res;
    std::array<Scalar, 3> err = {{-1, -1, -1}};
    Scalar ms = 1e-8;
    Scalar toi;
    const Scalar tolerance = 1e-6;
    const Scalar t_max = 1;
    const int max_itr = 1e6;
    Scalar output_tolerance;
    const int CCD_TYPE = 1;
	res = inclusion_ccd::edgeEdgeCCD_double(a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, err, ms, toi, tolerance, t_max, max_itr, output_tolerance);
    /*res = inclusion_ccd::edgeEdgeCCD_double(
        a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, err, ms, toi, tolerance, t_max,
        max_itr, output_tolerance);*/
    if (!using_rational_method())
    {
        std::cout << "the double ccd result is " << res << std::endl;
    }
    else
    {
        std::cout << "the ccd using rational core result is " << res
                  << std::endl;
    }
#ifdef CHECK_QUEUE_SIZE
    //std::cout << "queue size max " << return_queue_size() << std::endl;
#endif
}

//# define TIDBG


#ifdef TIGHT_INCLUSION_RUN_EXAMPLES

#include "read_rational_csv.hpp"
std::string root_path(TICCD_EXAMPLE_QUERIES_DIR);
std::vector<std::string> simulation_folders = {
    {"chain", "cow-heads", "golf-ball", "mat-twist"}};
std::vector<std::string> handcrafted_folders = {
    {"erleben-sliding-spike", "erleben-spike-wedge", "erleben-sliding-wedge",
     "erleben-wedge-crack", "erleben-spike-crack", "erleben-wedges",
     "erleben-cube-cliff-edges", "erleben-spike-hole",
     "erleben-cube-internal-edges", "erleben-spikes", "unit-tests"}};
std::vector<std::string> fnames = {{"data_0_0.csv"}, {"data_0_1.csv"}};
void run_rational_data_single_method(
    const bool is_edge_edge,
    const bool is_simulation_data,
    const double minimum_seperation,
    const double tolerance,
    const long max_itr = 1e6)
{

    std::string sub_folder, filebase;
    Eigen::MatrixXd all_V;
    std::vector<bool> results;

    int total_number = -1;
    double new_timing = 0.0;
    int total_positives = 0;
    int new_false_positives = 0;
    int new_false_negatives = 0;

    int nbr_larger_tol = 0;
    int nbr_diff_tol = 0;
    Timer timer;

    int max_fnbr = is_simulation_data ? simulation_folders.size()
                                      : handcrafted_folders.size();
    const auto folders =
        is_simulation_data ? simulation_folders : handcrafted_folders;
    sub_folder = is_edge_edge ? "/edge-edge/" : "/vertex-face/";
    for (int fnbr = 0; fnbr < max_fnbr; fnbr++)
    {
        for (int ff = 0; ff < 2; ff++)
        {
			if (folders[fnbr] == "erleben-spike-hole"&&ff == 1)
			{
				continue;
			}
			all_V = read_rational_csv(
#ifdef TIDBG
				"D:\\vs\\collision\\interval\\Tight-Inclusion\\build\\edge-edge-0474.csv",
#else
				root_path + folders[fnbr] + sub_folder + fnames[ff],
#endif
                
				results);
            //assert(all_V.rows() % 8 == 0 && all_V.cols() == 3);
            if (all_V.rows() % 8 != 0 || all_V.cols() != 3)
            {
                std::cout << "wrong data happens" << std::endl;
                continue;
            }
            int v_size = all_V.rows() / 8;
            for (int i = 0; i < v_size; i++)
            {
                total_number += 1;
                Eigen::Matrix<double, 8, 3> V = all_V.middleRows<8>(8 * i);
                bool expected_result =
                    results[i * 8]; // args[k].result;// from mathematica

                bool new_result;

                const std::array<double, 3> err = {{-1, -1, -1}};

                double toi;
                const double t_max = 1;

                double output_tolerance = tolerance;

                int CCD_TYPE = 1;
                timer.start();
#ifdef TIDBG
				if (i != 6382) {
					continue;
				}
				new_result = edgeEdgeCCD_double(
					V.row(0), V.row(1), V.row(2), V.row(3), V.row(4),
					V.row(5), V.row(6), V.row(7), err, minimum_seperation,
					toi, tolerance, t_max, max_itr, output_tolerance,
					CCD_TYPE);
				if (i == 6382) {
					for (int row = 0; row < 8; row++) {
						std::cout << V(row, 0) << "," << V(row, 1) << "," << V(row, 2) << "," << toi << std::endl;
					}
				}
#else
				if (is_edge_edge)
				{
					new_result = edgeEdgeCCD_double(
						V.row(0), V.row(1), V.row(2), V.row(3), V.row(4),
						V.row(5), V.row(6), V.row(7), err, minimum_seperation,
						toi, tolerance, t_max, max_itr, output_tolerance,
						CCD_TYPE);
				}
				else
				{
					new_result = vertexFaceCCD_double(
						V.row(0), V.row(1), V.row(2), V.row(3), V.row(4),
						V.row(5), V.row(6), V.row(7), err, minimum_seperation,
						toi, tolerance, t_max, max_itr, output_tolerance,
						CCD_TYPE);
				}
#endif
               

                new_timing += timer.getElapsedTimeInMicroSec();
                std::cout << total_number << "\r" << std::flush;

                if (expected_result)
                {
                    total_positives++;
                }
                if (new_result != expected_result)
                {
                    if (new_result)
                    {
                        new_false_positives++;
                    }
                    else
                    {
                        new_false_negatives++;

                        std::cout << "false negative, "
                                  << root_path + folders[fnbr] + sub_folder
                                         + fnames[ff]
                                  << ", " << i << std::endl;
                        for (int j = 0; j < 8; j++)
                        {
                            std::cout << "v" << j << " " << V(j, 0) << ", "
                                      << V(j, 1) << ", " << V(j, 2)
                                      << std::endl;
                            if (j == 3)
                                std::cout << std::endl;
                        }

                        std::cout << "is edge? " << is_edge_edge << std::endl
                                  << std::endl;
                        throw "\n\n\nWRONG ANSWERS\n\n\n";
                        exit(0);
                    }
                }
            }
#ifdef TIDBG
			std::cout << "fps " << new_false_positives << " fns " << new_false_negatives << " total queries " << total_number + 1
				<<" total positives "<< total_positives << std::endl;

			exit(0);
#endif
        }
    }

    std::cout << "total number, " << total_number + 1 << std::endl;
    std::cout << "total positives, " << total_positives << std::endl;
    std::cout << "is_edge_edge? , " << is_edge_edge << std::endl;
    std::cout << "false_positives, " << new_false_positives << std::endl;
    std::cout << "false_negatives, " << new_false_negatives << std::endl;
    std::cout << "average time, " << new_timing / double(total_number + 1)
              << std::endl
              << std::endl;
    std::cout << "total time, " << new_timing << std::endl << std::endl;
    std::cout << "** finished " << std::endl;
}

void run_code()
{
    double ms = 0.0;
    double tolerance = 1e-6;

	std::cout << "\nRunning simulation Vertex-Face data" << std::endl;
	run_rational_data_single_method(
		/*is_edge_edge*/ false, /*is_simulation*/ true, ms, tolerance);
	std::cout << "finish simulation Vertex-Face data" << std::endl;
	std::cout << "\nRunning simulation Edge-Edge data" << std::endl;
	run_rational_data_single_method(
		/*is_edge_edge*/ true, /*is_simulation*/ true, ms, tolerance);
	std::cout << "finish simulation Edge-Edge data" << std::endl;

    std::cout << "\nRunning handcrafted Vertex-Face data" << std::endl;
    run_rational_data_single_method(
        /*is_edge_edge*/ false, /*is_simulation*/ false, ms, tolerance);
    std::cout << "finish handcrafted Vertex-Face data" << std::endl;
    std::cout << "\nRunning handcrafted Edge-Edge data" << std::endl;
    run_rational_data_single_method(
        /*is_edge_edge*/ true, /*is_simulation*/ false, ms, tolerance);
    std::cout << "finish handcrafted Edge-Edge data" << std::endl;


}
#endif

//void run_dbg() {
//	std::string file = "D:\\vs\\collision\\interval\\Tight-Inclusion\\build\\Release\\ee1simu1.csv";
//	std::vector<double> tois;
//	std::cout << "before read" << std::endl;
//	Eigen::MatrixXd all_V = read_csv(file, tois);
//	std::cout << "readed" << std::endl;
//	if (all_V.rows() != tois.size() || all_V.rows() % 8 != 0) {
//		std::cout << "wrong sizes, " << all_V.rows() << " " << tois.size() << std::endl;
//	}
//
//	int v_size = all_V.rows() / 8;
//	int counter = 0;
//	double tolerance = 1e-6;
//	double minimum_seperation = 0;
//	int max_itr = 1e6;
//	int total_positives = 0;
//	int no_zero = 0;
//	for (int i = 0; i < v_size; i++) {
//		double toi_float = tois[i * 8];
//		if (toi_float > 0) continue;
//		Eigen::Matrix<double, 8, 3> V = all_V.middleRows<8>(8 * i);
//		const std::array<double, 3> err = { {-1, -1, -1} };
//
//		double toi;
//		const double t_max = 1;
//
//		double output_tolerance = tolerance;
//
//		int CCD_TYPE = 1;
//		bool new_result = edgeEdgeCCD_double(
//			V.row(0), V.row(1), V.row(2), V.row(3), V.row(4),
//			V.row(5), V.row(6), V.row(7), err, minimum_seperation,
//			toi, tolerance, t_max, max_itr, output_tolerance,
//			CCD_TYPE);
//		if (new_result) {
//			total_positives += 1;
//		}
//		if (toi > 0) {
//			no_zero += 1;
//		}
//
//	}
//	std::cout << "total positives " << total_positives << std::endl;
//	std::cout << "no zero " << no_zero << std::endl;
//}

int main(int argc, char *argv[])
{
	//run_dbg();
#ifdef TIGHT_INCLUSION_RUN_EXAMPLES
    run_code();
#else
    case_check();
#endif

    std::cout << "using double precision? " << ticcd_using_double()
              << std::endl;
    std::cout << "done!" << std::endl;

    return 0;
}
