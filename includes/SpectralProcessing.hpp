#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <string>
#include "gnuplot-iostream.h"  

class Graph
{
public:
    Graph(const std::string& title = "") : title_(title) {}
    virtual ~Graph() = default;

    virtual void histogram(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& data) = 0;
    virtual void plot(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& data) = 0;

    void setTitle(const std::string& title) { title_ = title; }
    std::string getTitle() const { return title_; }

protected:
    std::string title_;
};


class Histogram : public Graph
{
public:
    Histogram(const std::string& title = "") : Graph(title) {}
    virtual void plot(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& data) override;
};

class Pdf_graph : public Graph
{
public:
    Pdf_graph(const std::string& title = "") : Graph(title) {}
    virtual void plot(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& data) override;
};
