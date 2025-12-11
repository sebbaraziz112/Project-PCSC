/**
 * @file Visual.hpp
 * @brief Definition of classes for generating plots using gnuplot-iostream.
 */

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "gnuplot-iostream.h"

/**
 * @class Plot
 * @brief Abstract base class for all plot types.
 */
class Plot
{
public:
    Plot() {}
    virtual ~Plot() = default;

    /**
     * @brief Pure virtual method to plot data.
     * @param data Integer matrix of size 2xN representing the data.
     * @param filename Output file name for the plot.
     */
    virtual void plot(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& data,
                      std::string filename) = 0;

    void setTitle(const std::string& title)       { title_ = title; }
    void setAxisXLabel(const std::string& label)  { axisXLabel_ = label; }
    void setAxisYLabel(const std::string& label)  { axisYLabel_ = label; }

    std::string getTitle()  const { return title_; }
    std::string getAxisXLabel() const { return axisXLabel_; }
    std::string getAxisYLabel() const { return axisYLabel_; }

protected:
    std::string title_;
    std::string axisXLabel_;
    std::string axisYLabel_;
};

/**
 * @class Histogram
 */
class Histogram : public Plot
{
public:
    Histogram() : Plot() {}

    virtual void plot(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& data,
                      std::string filename) override;
};

/**
 * @class Plot2D
 */
class Plot2D : public Plot
{
public:
    Plot2D() : Plot() {}

    virtual void plot(const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& data,
                      std::string filename) override;
};
