/**
 * @file Visual.hpp
 * @brief Definition of classes for generating plots using gnuplot-iostream.
 *
 * This file contains abstract and concrete classes to generate
 * histograms and 2D plots from complex matrices.
 */

#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <string>
#include "gnuplot-iostream.h"  

/**
 * @class Plot
 * @brief Abstract base class for all plot types.
 *
 * Provides a common interface for all plots, including title
 * and axis labels.
 */
class Plot
{
public:
    /**
     * @brief Default constructor.
     */
    Plot() {}

    /**
     * @brief Virtual destructor.
     */
    virtual ~Plot() = default;

    /**
     * @brief Pure virtual method to plot data.
     * @param data Complex matrix of size 2xN representing the data.
     * @param filename Output file name for the plot.
     */
    virtual void plot(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& data,
                      std::string filename) = 0;

    /**
     * @brief Set the plot title.
     * @param title Title to display on the plot.
     */
    void setTitle(const std::string& title) { title_ = title; }

    /**
     * @brief Set the X-axis label.
     * @param label Label for the X-axis.
     */
    void setAxisXLabel(const std::string& label) { axisXLabel_ = label; }

    /**
     * @brief Set the Y-axis label.
     * @param label Label for the Y-axis.
     */
    void setAxisYLabel(const std::string& label) { axisYLabel_ = label; }

    /**
     * @brief Get the plot title.
     * @return The title of the plot.
     */
    std::string getTitle() const { return title_; }

    /**
     * @brief Get the X-axis label.
     * @return Label of the X-axis.
     */
    std::string getAxisXLabel() const { return axisXLabel_; }

    /**
     * @brief Get the Y-axis label.
     * @return Label of the Y-axis.
     */
    std::string getAxisYLabel() const { return axisYLabel_; }

protected:
    std::string title_;       /**< Plot title */
    std::string axisXLabel_;  /**< X-axis label */
    std::string axisYLabel_;  /**< Y-axis label */
};

/**
 * @class Histogram
 * @brief Class for generating histogram plots.
 *
 * Inherits from Plot and implements the plot() method to draw histograms.
 */
class Histogram : public Plot
{
public:
    /**
     * @brief Default constructor.
     */
    Histogram() : Plot() {}

    /**
     * @brief Plot data as a histogram.
     * @param data Complex matrix of size 2xN representing the histogram values.
     * @param filename Output file name for the histogram image.
     */
    virtual void plot(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& data,
                      std::string filename) override;
};

/**
 * @class Plot2D
 * @brief Class for generating 2D line plots.
 *
 * Inherits from Plot and implements the plot() method for 2D plotting.
 */
class Plot2D : public Plot
{
public:
    /**
     * @brief Default constructor.
     */
    Plot2D() : Plot() {}

    /**
     * @brief Plot data as a 2D line plot.
     * @param data Complex matrix of size 2xN representing the points to plot.
     * @param filename Output file name for the 2D plot image.
     */
    virtual void plot(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& data,
                      std::string filename) override;
};
