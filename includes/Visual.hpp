#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <string>
#include "gnuplot-iostream.h"  

class Plot
{
public:
    Plot() : {}
    virtual ~Plot() = default;

    virtual void plot(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& data, std::string filename) = 0;

    void setTitle(const std::string& title) { title_ = title; }
    void setAxisXLabel(const std::string& label) { axisXLabel_ = label; }
    void setAxisYLabel(const std::string& label) { axisYLabel_ = label; }


    std::string getTitle() const { return title_; }
    std::string getAxisXLabel() const { return axisXLabel_; }
    std::string getAxisYLabel() const { return axisYLabel_; }

protected:
    std::string title_;
    std::string axisXLabel_;
    std::string axisYLabel_;
};


class Histogram : public Plot
{
public:
    Histogram() : Plot() {}
    virtual void plot(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& data, std::string filename) override;
};

class Plot2D : public Plot
{
public:
    Plot2D() : Plot() {}
    virtual void plot(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& data, std::string filename) override;
};
