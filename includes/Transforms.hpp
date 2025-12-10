/**
 * @file Transforms.hpp
 * @brief Defines abstract and concrete matrix transformation classes.
 *
 * This file provides an abstract base class for matrix transformations
 * and a concrete identity transform that returns the input matrix unchanged.
 */

#ifndef TRANSFORMS
#define TRANSFORMS

#include <string>
#include <Eigen/Dense>
#include <cmath>

/**
 * @class AbstractTransform
 * @brief Abstract base class for matrix transformations.
 *
 * Defines a common interface for transformations that operate on
 * Eigen matrices. Derived classes must implement the apply() method.
 * @tparam T Type of the matrix elements (e.g., double, float, int).
 */
template <typename T>
class AbstractTransform
{
public:

    /**
     * @brief Constructor.
     * @param name Name of the transformation.
     */
    AbstractTransform(std::string name) : name_(name) {};

    /**
     * @brief Apply the transformation to a matrix.
     * @param matrix Input matrix to transform.
     * @return Transformed matrix.
     *
     * This is a pure virtual function that must be implemented
     * by derived classes.
     */
    virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> apply(
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix) const = 0;

    /**
     * @brief Virtual destructor.
     */
    virtual ~AbstractTransform() {};

private:
    std::string name_; /**< Name of the transformation */
};

/**
 * @class IdentityTransform
 * @brief Concrete class representing the identity transformation.
 *
 * Returns the input matrix unchanged.
 * @tparam T Type of the matrix elements (e.g., double, float, int).
 */
template <typename T>
class IdentityTransform : public AbstractTransform<T>
{
public:

    /**
     * @brief Constructor.
     * @param name Name of the transformation.
     */
    IdentityTransform(std::string name) : AbstractTransform<T>(name) {};

    /**
     * @brief Apply the identity transformation.
     * @param matrix Input matrix.
     * @return The same matrix as input (no changes applied).
     */
    virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> apply(
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix) const override
    {
        return matrix;
    }

};

#endif // TRANSFORMS
