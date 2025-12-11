#include "../includes/Visual.hpp"
#include "gnuplot-iostream.h"
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <string>

void Histogram::plot(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& data, std::string filename)
{
    if (data.rows() != 2) {
        std::cerr << "Histogram::plot: data must be 2 x N\n";
        return;
    }

    Gnuplot gp;
    gp << "set terminal pngcairo size 800,600\n";
    gp << "set output '" << filename << "'\n";

    gp << "set style data histogram\n";
    gp << "set boxwidth 0.75\n";
    gp << "set style fill solid 1.0 border -1\n";
    gp << "set title '" << title_ << "'\n";
    gp << "set xlabel '" << axisXLabel_ << "'\n";
    gp << "set ylabel '" << axisYLabel_ << "'\n";

    std::vector<std::pair<int,int>> plot_data;
    for (int i = 0; i < data.cols(); ++i) {
        int x = data(0, i);
        int y = data(1, i);
        plot_data.emplace_back(x, y);
    }

    gp << "plot '-' using 1:2 with boxes lc rgb 'blue' notitle\n";
    gp.send1d(plot_data);

    std::cout << "Histogram saved to " << filename << "\n";
}


void Plot2D::plot(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& data, std::string filename)
{
    if (data.rows() != 2) {
        std::cerr << "Plot2D::plot: data must be 2 x N\n";
        return;
    }

    Gnuplot gp;
    gp << "set terminal pngcairo size 800,600\n";
    gp << "set output '" << filename << "'\n";

    gp << "set style data lines\n";
    gp << "set title '" << title_ << "'\n";
    gp << "set xlabel '" << axisXLabel_ << "'\n";
    gp << "set ylabel '" << axisYLabel_ << "'\n";

    std::vector<std::pair<int,int>> plot_data;
    for (int i = 0; i < data.cols(); ++i) {
        int x = data(0, i);
        int y = data(1, i);
        plot_data.emplace_back(x, y);
    }

    gp << "plot '-' using 1:2 with lines lc rgb 'red' notitle\n";
    gp.send1d(plot_data);

    std::cout << "Plot saved to " << filename << "\n";
}
