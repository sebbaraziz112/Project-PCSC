#ifndef METHODS
#define METHODS

#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <cassert>




class AbstractMethod
{

public: 

    AbstractMethod(std::string name): name_(name){};
    ~AbstractMethod(){};

    template <typename T>
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> ComplexCast(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const;

    template <typename O>
    std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> OCast(std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const;

    template <typename T, typename O>
    std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> computeMethod(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const;
    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const = 0;

    int getNextPowerTwo(int M) const;

    std::string name_;

protected: 



};

class IdentityMethod: public AbstractMethod
{

public: 
    IdentityMethod(): AbstractMethod("Identity Method"){};
    ~IdentityMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

private: 

};

class NaiveDFTMethod: public AbstractMethod
{

public: 

    NaiveDFTMethod(): AbstractMethod("DFT Method"){};
    ~NaiveDFTMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

private: 

    std::complex<double> getResultElement(int u, int v, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const;


};

// 1D DFT Methods

class NaiveInvDFTMethod: public AbstractMethod
{

public: 

    NaiveInvDFTMethod(): AbstractMethod("Inverse DFT Method"){};
    ~NaiveInvDFTMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

private: 

    std::complex<double> getResultElement(int u, int v, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const;

};

class FFT1DMethod: public AbstractMethod
{

public: 

    FFT1DMethod(bool padd): AbstractMethod("FFT Method"), padd_(padd){};
    ~FFT1DMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

    std::vector<std::complex<double>> getCloserPowerTwo(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const;
    std::vector<std::complex<double>> fftCompute(std::vector<std::complex<double>> vector_c) const;
    std::vector<std::complex<double>> getEven(std::vector<std::complex<double>>& vector_c) const;
    std::vector<std::complex<double>> getOdd(std::vector<std::complex<double>>& vector_c) const;

private:

    

    bool padd_;

};

class Inv1DFFTMethod: public AbstractMethod
{

public: 

    Inv1DFFTMethod(bool padd): AbstractMethod("Inverse FFT Method"), padd_(padd){};
    ~Inv1DFFTMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

    std::vector<std::complex<double>> getCloserPowerTwo(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const;
    std::vector<std::complex<double>> fftCompute(std::vector<std::complex<double>> vector_c) const;
    std::vector<std::complex<double>> getEven(std::vector<std::complex<double>>& vector_c) const;
    std::vector<std::complex<double>> getOdd(std::vector<std::complex<double>>& vector_c) const;

private: 

    

    bool padd_;

};

class BlueStein1DMethod: public AbstractMethod
{
public: 

    BlueStein1DMethod(): AbstractMethod("BlueStein1D Method"){
        fftMethod_ = std::make_shared<FFT1DMethod>(false);
        invfftMethod_ = std::make_shared<Inv1DFFTMethod>(false);
    };
    ~BlueStein1DMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

private: 

    std::shared_ptr<FFT1DMethod> fftMethod_;
    std::shared_ptr<Inv1DFFTMethod> invfftMethod_;

};

class InvBlueStein1DMethod: public AbstractMethod
{

public: 

    InvBlueStein1DMethod(): AbstractMethod("Inverse BlueStein1D Method"){
        fftMethod_ = std::make_shared<FFT1DMethod>(false);
        invfftMethod_ = std::make_shared<Inv1DFFTMethod>(false);
    };
    ~InvBlueStein1DMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;


private: 

    std::shared_ptr<FFT1DMethod> fftMethod_;
    std::shared_ptr<Inv1DFFTMethod> invfftMethod_;



};

//////////////////////////////////////

// ND DFT Methods (1D and 2D)

class BlueSteinMethod: public AbstractMethod 
{

public: 

    BlueSteinMethod(): AbstractMethod("BlueStein Method"){
        bluestein1d_ = std::make_shared<BlueStein1DMethod>();
    }
    ~BlueSteinMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;


private: 

    std::shared_ptr<BlueStein1DMethod> bluestein1d_;

};

class InvBlueSteinMethod: public AbstractMethod 
{

public: 

    InvBlueSteinMethod(): AbstractMethod("Inverse BlueStein Method"){
        invbluestein1d_ = std::make_shared<InvBlueStein1DMethod>();
    }
    ~InvBlueSteinMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

private:

    std::shared_ptr<InvBlueStein1DMethod> invbluestein1d_;

};

//////////////////////////////////////

// Convolution Methods

template <int DIM>
class Convolution: public AbstractMethod
{

public: 

    
    Convolution(): AbstractMethod("Convolution Method"){};
    ~Convolution(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;
    void setKernel(std::complex<double> complexValue, int i, int j) {
        kernel_(i, j) = complexValue;
    }
    void normalizeKernel(){
        kernel_ = kernel_ / kernel_.sum();
    }
    Eigen::Matrix<std::complex<double>, DIM, DIM> kernel_;

protected: 

    std::complex<double> convolve(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix, int i, int j) const;

};

template <int DIM>
class MeanConv: public Convolution<DIM>
{

public: 

    MeanConv(): Convolution<DIM>() {
        for (int i = 0; i < DIM; i ++){
            for (int j = 0; j < DIM; j ++){
                this->setKernel(std::complex<double>(1.0/std::pow(DIM, 2), 0.0), i, j);
            }
        }
    }
    ~MeanConv(){};

};

template <int DIM>
class BlurrGaussianConv: public Convolution<DIM>
{

public: 

    BlurrGaussianConv(double sigmaX, double sigmaY): Convolution<DIM>(), sigmaX_(sigmaX), sigmaY_(sigmaY) {
        double xyc = DIM/2 + 1/2;
        for (int i = 0; i < DIM; i ++){
            for (int j = 0; j < DIM; j ++){
                this->setKernel(std::complex<double>(this->getGaussian(i, j, xyc, xyc), 0.0), i, j);
            }
        }
        this->normalizeKernel();
    }

    ~BlurrGaussianConv(){};

    double getGaussian(int x, int y, int xc, int yc){
        return (1.0 / (2 * M_PI * sigmaX_ * sigmaY_)) * std::exp((-std::pow((x - xc)/sigmaX_, 2) - std::pow((y - yc)/sigmaY_, 2))/2.0);
    }

    double sigmaX_;
    double sigmaY_;

};

class SobelX : public Convolution<3>
{
public:
    SobelX(): Convolution<3>() {
        kernel_ <<
            std::complex<double>(-1,0), std::complex<double>(0,0), std::complex<double>(1,0),
            std::complex<double>(-2,0), std::complex<double>(0,0), std::complex<double>(2,0),
            std::complex<double>(-1,0), std::complex<double>(0,0), std::complex<double>(1,0);
    }
};

class SobelY : public Convolution<3>
{
public:
    SobelY() : Convolution<3>() {
        kernel_ <<
            std::complex<double>( 1,0), std::complex<double>( 2,0), std::complex<double>( 1,0),
            std::complex<double>( 0,0), std::complex<double>( 0,0), std::complex<double>( 0,0),
            std::complex<double>(-1,0), std::complex<double>(-2,0), std::complex<double>(-1,0);
    }
};

class LaplacianClassical : public Convolution<3>
{
public:
    LaplacianClassical() : Convolution<3>() {
        kernel_ <<
            std::complex<double>(0,0),  std::complex<double>(-1,0), std::complex<double>(0,0),
            std::complex<double>(-1,0), std::complex<double>( 4,0), std::complex<double>(-1,0),
            std::complex<double>(0,0),  std::complex<double>(-1,0), std::complex<double>(0,0);
    }
};

class Laplacian8Connected : public Convolution<3>
{
public:
    Laplacian8Connected() : Convolution<3>() {
        kernel_ <<
            std::complex<double>(-1,0), std::complex<double>(-1,0), std::complex<double>(-1,0),
            std::complex<double>(-1,0), std::complex<double>( 8,0), std::complex<double>(-1,0),
            std::complex<double>(-1,0), std::complex<double>(-1,0), std::complex<double>(-1,0);
    }
};

class Laplacian4Connected : public Convolution<3>
{
public:
    Laplacian4Connected() : Convolution<3>() {
        kernel_ <<
            std::complex<double>(0,0),  std::complex<double>( 1,0), std::complex<double>(0,0),
            std::complex<double>(1,0),  std::complex<double>(-4,0), std::complex<double>(1,0),
            std::complex<double>(0,0),  std::complex<double>( 1,0), std::complex<double>(0,0);
    }
};

class LaplacianOfGaussian : public Convolution<5>
{
public:
    LaplacianOfGaussian() : Convolution<5>() {
        kernel_ <<
            std::complex<double>(0,0),  std::complex<double>(0,0),  std::complex<double>(-1,0), std::complex<double>(0,0),  std::complex<double>(0,0),
            std::complex<double>(0,0),  std::complex<double>(-1,0), std::complex<double>(-2,0), std::complex<double>(-1,0), std::complex<double>(0,0),
            std::complex<double>(-1,0), std::complex<double>(-2,0), std::complex<double>(16,0), std::complex<double>(-2,0), std::complex<double>(-1,0),
            std::complex<double>(0,0),  std::complex<double>(-1,0), std::complex<double>(-2,0), std::complex<double>(-1,0), std::complex<double>(0,0),
            std::complex<double>(0,0),  std::complex<double>(0,0),  std::complex<double>(-1,0), std::complex<double>(0,0),  std::complex<double>(0,0);
    }
};

//////////////////////////////////////

// Amplification and Attenuation

class BandFiltering: public AbstractMethod
{

public:
    
    BandFiltering(): AbstractMethod("Band Filtering Method"){
        amplificationFactor_ = 1.0;
        lowerFractionX_ = 0.0;
        lowerFractionY_ = 0.0;
        upperFractionX_ = 1.0;
        upperFractionY_ = 1.0;

    };
    ~BandFiltering(){};

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> getBandMatrix(int nRows, int nCols) const;

    template <typename T>
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> shiftCenter(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const;

    template <typename T>
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> shiftInverse(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const; 
    
    void setLowerFractionX(double lowerFraction){
        lowerFractionX_ = lowerFraction;
    }
    void setUpperFractionX(double upperFraction){
        upperFractionX_ = upperFraction;
    }
    void setLowerFractionY(double lowerFraction){
        lowerFractionY_ = lowerFraction;
    }
    void setUpperFractionY(double upperFraction){
        upperFractionY_ = upperFraction;
    }
    void setAmplificationFactor(double amplificationFactor){
        amplificationFactor_ = amplificationFactor;
    }

private: 

    double lowerFractionX_;
    double upperFractionX_;
    double lowerFractionY_;
    double upperFractionY_;

    double amplificationFactor_;

};

class LowPass: public BandFiltering
{

public: 

    LowPass(): BandFiltering(){
        this->setLowerFractionX(0.0);
        this->setLowerFractionY(0.0);
        this->setAmplificationFactor(1.0);
    }
    ~LowPass(){};

};

class HighPass: public BandFiltering
{

public: 

    HighPass(): BandFiltering(){
        this->setUpperFractionX(0.0);
        this->setUpperFractionY(0.0);
        this->setAmplificationFactor(1.0);
    }
    ~HighPass(){};

};

//////////////////////////////////////


// AbstractMethod

template <typename T> 
std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> AbstractMethod::ComplexCast(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const{
    if constexpr (std::is_same_v<T, std::complex<double>>) {
        std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> TVector;
        TVector.reserve(matrices.size());
        for (auto& mat : matrices) {
            TVector.push_back(mat);  
        }
        return TVector;
    }
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> complexVector;
    const int nRows = matrices[0].rows();
    const int nCols = matrices[0].cols();
    for (int i = 0; i < matrices.size(); i ++){
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> mat(nRows, nCols);
        complexVector.push_back(mat);
    }
    for (int chann = 0; chann < matrices.size(); chann ++){
        for (int i = 0; i < nRows; i ++){
            for (int j = 0; j < nCols; j ++){
                complexVector[chann](i, j) = static_cast<std::complex<double>>(matrices[chann](i, j));
            }
        }
    }
    return complexVector;
}

template <typename O>
std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> AbstractMethod::OCast(std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const{
    if constexpr (std::is_same_v<O, std::complex<double>>) {
        std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> OVector;
        for (auto& mat : matrices) {
            OVector.push_back(mat);  
        }
        return OVector;
    }
    std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> OVector;
    const int nRows = matrices[0].rows();
    const int nCols = matrices[0].cols();
    for (int i = 0; i < matrices.size(); i ++){
        Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic> mat(nRows, nCols);
        OVector.push_back(mat);
    }
    for (int chann = 0; chann < matrices.size(); chann ++){
        for (int i = 0; i < nRows; i ++){
            for (int j = 0; j < nCols; j ++){
                OVector[chann](i, j) = static_cast<O>(std::round(std::abs(matrices[chann](i, j))));
            }
        }
    }
    return OVector;
}

template <typename T, typename O>
std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> AbstractMethod::computeMethod(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const{
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultVector;
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> compMatrices = this->ComplexCast(matrices);
    for (int i = 0; i < matrices.size(); i ++){
        resultVector.push_back(this->computeMatrixMethod(compMatrices[i]));
    }
    std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> finalVector = this->OCast<O>(resultVector);
    return finalVector;
}

//////////////////////////////////////

// Convolution

template <int DIM>
std::complex<double> Convolution<DIM>::convolve(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix, int i, int j) const {
    std::complex<double> result = std::complex<double>(0.0, 0.0);
    const int nRows = matrix.rows(); const int nCols = matrix.cols();
    const int gap = std::ceil((double)(DIM)/2) - 1;
    int indexX = 0;
    int indexY = 0;
    for (int kx = -gap; kx <= gap; kx ++){
        for (int ky = -gap; ky <= gap; ky ++){
            indexX = i + kx; 
            indexY = j + ky;
            if ((indexX >= 0) && (indexX < nRows) && (indexY >= 0) && (indexY < nCols)){
                result += matrix(i + kx, j + ky) * kernel_(kx + gap, ky + gap);
            } 
        }
    }
    return result;
}

template <int DIM>
Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> Convolution<DIM>::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    const int nRows = matrix.rows(); const int nCols = matrix.cols();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> resultMatrix(nRows, nCols);
    for (int i = 0; i < nRows; i ++){
        for (int j = 0; j < nCols; j ++){
            resultMatrix(i, j) = this->convolve(matrix, i, j);
        }
    }
    return resultMatrix;

}

//////////////////////////////////////

// BandFiltering

template <typename T>
std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> BandFiltering::shiftCenter(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const {
    const int nChann = matrices.size();
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> centeredMatrices;
    centeredMatrices.reserve(nChann);

    for (int chann = 0; chann < nChann; chann++) {
        const auto& mat = matrices[chann];
        const int nRows = mat.rows();
        const int nCols = mat.cols();

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> newMatrix(nRows, nCols);

        int midRow = nRows / 2;
        int midCol = nCols / 2;
        newMatrix.topLeftCorner(midRow, midCol) = mat.bottomRightCorner(nRows - midRow, nCols - midCol);
        newMatrix.topRightCorner(midRow, nCols - midCol) = mat.bottomLeftCorner(nRows - midRow, midCol);
        newMatrix.bottomLeftCorner(nRows - midRow, midCol) = mat.topRightCorner(midRow, nCols - midCol);
        newMatrix.bottomRightCorner(nRows - midRow, nCols - midCol) = mat.topLeftCorner(midRow, midCol);

        centeredMatrices.push_back(newMatrix);
    }

    return centeredMatrices;
}

template <typename T>
std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> BandFiltering::shiftInverse(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const {
    const int nChann = matrices.size();
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> invMatrices;
    invMatrices.reserve(nChann);

    for (int chann = 0; chann < nChann; chann++) {
        const auto& mat = matrices[chann];
        const int nRows = mat.rows();
        const int nCols = mat.cols();

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> newMatrix(nRows, nCols);

        int midRow = nRows / 2;
        int midCol = nCols / 2;

        newMatrix.topLeftCorner(midRow, midCol) = mat.bottomRightCorner(nRows - midRow, nCols - midCol);
        newMatrix.topRightCorner(midRow, nCols - midCol) = mat.bottomLeftCorner(nRows - midRow, midCol);
        newMatrix.bottomLeftCorner(nRows - midRow, midCol) = mat.topRightCorner(midRow, nCols - midCol);
        newMatrix.bottomRightCorner(nRows - midRow, nCols - midCol) = mat.topLeftCorner(midRow, midCol);

        invMatrices.push_back(newMatrix);
    }

    return invMatrices;
}

#endif