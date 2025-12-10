#ifndef TRANSFORMS
#define TRANSFORMS

#include <string>
#include <Eigen/Dense>

template <typename T>
class AbstractTransform
{

public: 

    AbstractTransform(std::string name): name_(name){};

    virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> apply(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix) const = 0;

    virtual ~AbstractTransform();


private: 

    std::string name_;

};

template <typename T>
class IdentityTransform : public AbstractTransform<T>
{

public: 

    IdentityTransform(std::string name): AbstractTransform<T>(name){};

    virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> apply(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix) const override{
        return matrix;
    }

};

#endif