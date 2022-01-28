#include "read_rational_csv.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include <tight_inclusion/rational/rational.hpp>

namespace ticcd::rational {

    // for debug
    Eigen::MatrixXd
    read_csv(const std::string &inputFileName, std::vector<double> &results)
    {
        // be careful, there are n lines which means there are n/8 queries, but has
        // n results, which means results are duplicated
        results.clear();
        std::vector<std::array<double, 3>> vs;
        vs.clear();
        std::ifstream infile;
        infile.open(inputFileName);
        std::array<double, 3> v;
        if (!infile.is_open()) {
            std::cout << "Path Wrong!!!!" << std::endl;
            std::cout << "path, " << inputFileName << std::endl;
            return Eigen::MatrixXd(1, 1);
        }

        int l = 0;
        while (infile) // there is input overload classfile
        {
            l++;
            std::string s;
            if (!getline(infile, s))
                break;
            if (l == 1)
                continue; // skip the first row
            if (s[0] != '#') {
                std::istringstream ss(s);
                std::array<std::string, 5> record;
                int c = 0;
                while (ss) {
                    std::string line;
                    if (!getline(ss, line, ','))
                        break;
                    try {

                        record[c] = line;
                        c++;

                    } catch (const std::invalid_argument e) {
                        std::cout << "NaN found in file " << inputFileName
                                  << " line " << l << std::endl;
                        e.what();
                    }
                }

                double x = std::stod(record[0]);
                double y = std::stod(record[1]);
                double z = std::stod(record[2]);

                v[0] = x;
                v[1] = y;
                v[2] = z;
                vs.push_back(v);
                // if (record[6] != "1" && record[6] != "0") {
                //     std::cout
                //         << "ERROR:result position should be 1 or 0, but it is "
                //         << record[6] << std::endl;
                // }

                results.push_back(std::stod(record[4]));
            }
        }

        Eigen::MatrixXd all_v(vs.size(), 3);
        for (int i = 0; i < vs.size(); i++) {
            all_v(i, 0) = vs[i][0];
            all_v(i, 1) = vs[i][1];
            all_v(i, 2) = vs[i][2];
        }
        if (!infile.eof()) {
            std::cerr << "Could not read file " << inputFileName << "\n";
        }

        return all_v;
    }

    Eigen::MatrixXd read_rational_csv(
        const std::string &inputFileName, std::vector<bool> &results)
    {
        // be careful, there are n lines which means there are n/8 queries, but has
        // n results, which means results are duplicated
        results.clear();
        std::vector<std::array<double, 3>> vs;
        vs.clear();
        std::ifstream infile;
        infile.open(inputFileName);
        std::array<double, 3> v;
        if (!infile.is_open()) {
            std::cout << "Path Wrong!!!!" << std::endl;
            std::cout << "path, " << inputFileName << std::endl;
            Eigen::MatrixXd rst = Eigen::MatrixXd::Zero(1, 1);
            return rst;
        }

        int l = 0;
        while (infile) // there is input overload classfile
        {
            l++;
            std::string s;
            if (!getline(infile, s))
                break;
            if (s[0] != '#') {
                std::istringstream ss(s);
                std::array<std::string, 7>
                    record; // the first six are one vetex,
                            // the seventh is the result
                int c = 0;
                while (ss) {
                    std::string line;
                    if (!getline(ss, line, ','))
                        break;
                    try {
                        record[c] = line;
                        c++;
                    } catch (const std::invalid_argument e) {
                        std::cout << "NaN found in file " << inputFileName
                                  << " line " << l << std::endl;
                        e.what();
                    }
                }
                Rational rt;
                double x = rt.get_double(record[0], record[1]),
                       y = rt.get_double(record[2], record[3]),
                       z = rt.get_double(record[4], record[5]);
                v[0] = x;
                v[1] = y;
                v[2] = z;
                vs.push_back(v);
                // if (record[6] != "1" && record[6] != "0") {
                //     std::cout
                //         << "ERROR:result position should be 1 or 0, but it is "
                //         << record[6] << std::endl;
                // }
                results.push_back(bool(std::stoi(record[6])));
            }
        }
        Eigen::MatrixXd all_v(vs.size(), 3);
        for (int i = 0; i < vs.size(); i++) {
            all_v(i, 0) = vs[i][0];
            all_v(i, 1) = vs[i][1];
            all_v(i, 2) = vs[i][2];
        }
        if (!infile.eof()) {
            std::cerr << "Could not read file " << inputFileName << "\n";
        }

        return all_v;
    }

} // namespace ticcd::rational
