
#include <Eigen/Dense>
#include "../includes/methods.hpp"

// Abstract Method

int AbstractMethod::getNextPowerTwo(int M) const{
    return std::pow(2, std::ceil(std::log2(M)));
}

std::vector<int> AbstractMethod::getFactors(int N) const{
    std::vector<int> factors;
    double root = std::sqrt(N);
    for (int a = root; a >= 1; a --){
        if (N % a == 0){
            factors.push_back(a);
            factors.push_back(N / a);
            
            return factors;
        }
    }
    factors.push_back(N);
    factors.push_back(1);
    return factors;
}

///////////////////////////////////////////////

// Identity Method

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> IdentityMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const{
    return matrix;
}

///////////////////////////////////////////////

// Naive DFT Method

std::complex<double> NaiveDFTMethod::getResultElement(int u, int v, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    std::complex<double> result(0, 0);
    const int nRows = matrix.rows();
    const int nCols = matrix.cols();
    for (int i = 0; i < nRows; i ++){
        for (int j = 0; j < nCols; j ++){
            result += matrix(i, j) * std::polar(1.0, -2*M_PI*(((double)(u*i)/nRows)+((double)(v*j)/nCols)));
        }
    }
    return result;
}

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> NaiveDFTMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const{
    
    const int nRows = matrix.rows();
    const int nCols = matrix.cols();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> resultMatrix(nRows, nCols);
    for (int i = 0; i < nRows; i ++){
        for (int j = 0; j < nCols; j ++){
            resultMatrix(i, j) = this->getResultElement(i, j, matrix);
        }
    }
    return resultMatrix;
}

///////////////////////////////////////////////

// Naive Inv DFT Method

std::complex<double> NaiveInvDFTMethod::getResultElement(int u, int v, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    std::complex<double> result(0, 0);
    const int nRows = matrix.rows();
    const int nCols = matrix.cols();
    for (int i = 0; i < nRows; i ++){
        for (int j = 0; j < nCols; j ++){
            result += matrix(i, j) * std::polar(1.0 / (nRows * nCols), 2*M_PI*(((double)(u*i)/nRows)+((double)(v*j)/nCols)));
        }
    }
    return result;
}

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> NaiveInvDFTMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const{
    
    const int nRows = matrix.rows();
    const int nCols = matrix.cols();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> resultMatrix(nRows, nCols);
    for (int i = 0; i < nRows; i ++){
        for (int j = 0; j < nCols; j ++){
            resultMatrix(i, j) = this->getResultElement(i, j, matrix);
        }
    }
    return resultMatrix;
}


////////////////////////////////////////////////

// FFT Method

std::vector<std::complex<double>> FFT1DMethod::getCloserPowerTwo(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    assert(matrix.rows() == 1);
    const int nCols = matrix.cols();
    if (!padd_){
        assert((std::ceil(std::log2(nCols)) - std::log2(nCols)) < 1e-10);
    }
    std::vector<std::complex<double>> vector_c(matrix.row(0).data(), matrix.row(0).data() + matrix.cols());
    const int power2 = std::ceil(std::log2(nCols));
    for (int i = 0; i < std::pow(2, power2) - nCols; i ++){
        vector_c.push_back(std::complex<double>(0.0, 0.0));
    }
    return vector_c;
}

std::vector<std::complex<double>> FFT1DMethod::getEven(std::vector<std::complex<double>>& vector_c) const{
    std::vector<std::complex<double>> even_vector_c;
    for (auto it = vector_c.begin(); it != vector_c.end(); it += 2){
        even_vector_c.push_back(*it);
    }
    return even_vector_c;
}

std::vector<std::complex<double>> FFT1DMethod::getOdd(std::vector<std::complex<double>>& vector_c) const{
    std::vector<std::complex<double>> odd_vector_c;
    for (auto it = vector_c.begin(); it != vector_c.end(); it += 2){
        odd_vector_c.push_back(*(it+1));
    }
    return odd_vector_c;
}

std::vector<std::complex<double>> FFT1DMethod::fftCompute(std::vector<std::complex<double>> vector_c) const{
    const int M = vector_c.size();
    if (M == 1){
        return vector_c;
    }
    
    const std::complex<double> omega = std::polar(1.0, -2*M_PI/M);

    std::vector<std::complex<double>> even_vector_c = this->getEven(vector_c);
    std::vector<std::complex<double>> odd_vector_c = this->getOdd(vector_c);

    std::vector<std::complex<double>> y_even = this->fftCompute(even_vector_c);
    std::vector<std::complex<double>> y_odd = this->fftCompute(odd_vector_c);

    std::vector<std::complex<double>> y(M, std::complex(0.0, 0.0));

    for (int k = 0; k < M/2; k ++){
        y[k] = y_even[k] + std::pow(omega, k) * y_odd[k];
        y[k + M/2] = y_even[k] - std::pow(omega, k) * y_odd[k];
    }

    return y;

}

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> FFT1DMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    assert(matrix.rows() == 1);
    const int N = matrix.cols();
    std::vector<std::complex<double>> expandedFFT = this->fftCompute(this->getCloserPowerTwo(matrix));
    const int M = expandedFFT.size();
    Eigen::Map<Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic>> mat(expandedFFT.data(), N);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> finalMat = mat;
    return finalMat;
}

////////////////////////////////////////////////

// Inv FFT Method

std::vector<std::complex<double>> Inv1DFFTMethod::getCloserPowerTwo(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    assert(matrix.rows() == 1);
    const int nCols = matrix.cols();
    std::vector<std::complex<double>> vector_c(matrix.row(0).data(), matrix.row(0).data() + matrix.cols());
    const int power2 = std::ceil(std::log2(nCols));
    for (int i = 0; i < std::pow(2, power2) - nCols; i ++){
        vector_c.push_back(std::complex<double>(0.0, 0.0));
    }
    return vector_c;
}

std::vector<std::complex<double>> Inv1DFFTMethod::getEven(std::vector<std::complex<double>>& vector_c) const{
    std::vector<std::complex<double>> even_vector_c;
    for (auto it = vector_c.begin(); it != vector_c.end(); it += 2){
        even_vector_c.push_back(*it);
    }
    return even_vector_c;
}

std::vector<std::complex<double>> Inv1DFFTMethod::getOdd(std::vector<std::complex<double>>& vector_c) const{
    std::vector<std::complex<double>> odd_vector_c;
    for (auto it = vector_c.begin(); it != vector_c.end(); it += 2){
        odd_vector_c.push_back(*(it+1));
    }
    return odd_vector_c;
}

std::vector<std::complex<double>> Inv1DFFTMethod::fftCompute(std::vector<std::complex<double>> vector_c) const{
    const int M = vector_c.size();
    if (M == 1){
        return vector_c;
    }
    
    const std::complex<double> omega = std::polar(1.0, 2*M_PI/M);

    std::vector<std::complex<double>> even_vector_c = this->getEven(vector_c);
    std::vector<std::complex<double>> odd_vector_c = this->getOdd(vector_c);

    std::vector<std::complex<double>> y_even = this->fftCompute(even_vector_c);
    std::vector<std::complex<double>> y_odd = this->fftCompute(odd_vector_c);

    std::vector<std::complex<double>> y(M, std::complex(0.0, 0.0));

    for (int k = 0; k < M/2; k ++){
        y[k] = y_even[k] + std::pow(omega, k) * y_odd[k];
        y[k + M/2] = y_even[k] - std::pow(omega, k) * y_odd[k];
    }

    return y;

}

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> Inv1DFFTMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    assert(matrix.rows() == 1);
    const int N = matrix.cols();
    std::vector<std::complex<double>> expandedFFT = this->fftCompute(this->getCloserPowerTwo(matrix));
    const int M = expandedFFT.size();
    Eigen::Map<Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic>> mat(expandedFFT.data(), N);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> finalMat = mat;
    finalMat = finalMat / N;
    return finalMat ;
}

////////////////////////////////////////////////

// BlueStein1D Method

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> BlueStein1DMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    assert(matrix.rows() == 1);
    const int N = matrix.cols();
    const int M = this->getNextPowerTwo(2*N - 1);   

    std::vector<std::complex<double>> a(M, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> b(M, std::complex<double>(0.0, 0.0));

    double angle = 0.0;
    for (int n = 0; n < N; n ++){
        angle = M_PI * (double)(n*n) / N;
        a[n] = matrix(0, n) * std::polar(1.0, -angle);
        b[n] = std::polar(1.0, angle);
    }

    for (int n = 1; n < N; n ++){
        b[M - n] = b[n];
    }

    std::vector<std::complex<double>> A = fftMethod_->fftCompute(a);
    std::vector<std::complex<double>> B = fftMethod_->fftCompute(b);
    std::vector<std::complex<double>> C(M, std::complex<double>(0.0, 0.0));

    for (int k = 0; k < M; k ++){
        C[k] = A[k] * B[k];
    }

    std::vector<std::complex<double>> c = invfftMethod_->fftCompute(C);

    std::vector<std::complex<double>> x(N, std::complex(0.0, 0.0));

    for (int k = 0; k < N; k ++){
        angle = M_PI * (k*k) / N;
        x[k] = (c[k] / std::polar((double)(M), 0.0))  * std::polar(1.0, -angle);
    }

    Eigen::Map<Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic>> mat(x.data(), N);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> finalMat = mat;
    return finalMat;
}

////////////////////////////////////////////////

// Inverse BlueStein1D Method

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> InvBlueStein1DMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    assert(matrix.rows() == 1);
    const int N = matrix.cols();
    const int M = this->getNextPowerTwo(2*N - 1);   

    std::vector<std::complex<double>> a(M, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> b(M, std::complex<double>(0.0, 0.0));

    double angle = 0.0;
    for (int n = 0; n < N; n ++){
        angle = M_PI * (double)(n*n) / N;
        a[n] = matrix(0, n) * std::polar(1.0, angle);
        b[n] = std::polar(1.0, -angle);
    }

    for (int n = 1; n < N; n ++){
        b[M - n] = b[n];
    }

    std::vector<std::complex<double>> A = fftMethod_->fftCompute(a);
    std::vector<std::complex<double>> B = fftMethod_->fftCompute(b);
    std::vector<std::complex<double>> C(M, std::complex<double>(0.0, 0.0));

    for (int k = 0; k < M; k ++){
        C[k] = A[k] * B[k];
    }

    std::vector<std::complex<double>> c = invfftMethod_->fftCompute(C);

    std::vector<std::complex<double>> x(N, std::complex(0.0, 0.0));

    for (int k = 0; k < N; k ++){
        angle = M_PI * (k*k) / N;
        x[k] = ((c[k] / std::polar((double)(M), 0.0))  * std::polar(1.0, angle) / (double)(N));
    }

    Eigen::Map<Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic>> mat(x.data(), N);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> finalMat = mat;
    return finalMat;
}

///////////////////////////////////////////////

// BlueStein Method

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> BlueSteinMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    const int nRows = matrix.rows();
    const int nCols = matrix.cols();

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> resultMatrix(nRows, nCols);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> currentRow;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> currentCol;

    for (int row = 0; row < nRows; row ++){
        currentRow = matrix.row(row).eval();
        currentRow = bluestein1d_->computeMatrixMethod(currentRow);
        resultMatrix.row(row) = currentRow;
    }

    

    for (int col = 0; col < nCols; col ++){
        currentCol = resultMatrix.col(col).transpose().eval();
        currentCol = bluestein1d_->computeMatrixMethod(currentCol);
        currentCol = currentCol.transpose();
        resultMatrix.col(col) = currentCol;
    }

    return resultMatrix;

}

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> LineOnlyBlueSteinMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    const int nRows = matrix.rows();
    const int nCols = matrix.cols();

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> resultMatrix(nRows, nCols);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> currentRow;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> currentCol;

    for (int row = 0; row < nRows; row ++){
        currentRow = matrix.row(row).eval();
        currentRow = bluestein1d_->computeMatrixMethod(currentRow);
        resultMatrix.row(row) = currentRow;
    }

    return resultMatrix;

}

///////////////////////////////////////////////

// Inverse BlueStein Method

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> InvBlueSteinMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    const int nRows = matrix.rows();
    const int nCols = matrix.cols();

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> resultMatrix(nRows, nCols);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> currentRow;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> currentCol;


    for (int row = 0; row < nRows; row ++){
        currentRow = matrix.row(row).eval();
        resultMatrix.row(row) = invbluestein1d_->computeMatrixMethod(currentRow);
    }

    for (int col = 0; col < nCols; col ++){
        currentCol = resultMatrix.col(col).eval().transpose();
        resultMatrix.col(col) = (invbluestein1d_->computeMatrixMethod(currentCol)).transpose();
    }

    return resultMatrix;

}

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> LineOnlyInvBlueSteinMethod::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const {
    const int nRows = matrix.rows();
    const int nCols = matrix.cols();

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> resultMatrix(nRows, nCols);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> currentRow;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> currentCol;

    for (int row = 0; row < nRows; row ++){
        currentRow = matrix.row(row).eval();
        currentRow = invbluestein1d_->computeMatrixMethod(currentRow);
        resultMatrix.row(row) = currentRow;
    }

    return resultMatrix;

}

//////////////////////////////////////////////

// Amplification and Attenuation Method

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> BandFiltering::getBandMatrix(int nRows, int nCols) const {
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> boolMatrix(nRows, nCols);
    for (int i = 0; i < nRows; i ++){
        for (int j = 0; j < nCols; j ++){
            if ((i >= lowerFractionX_ * (nRows-1) - 1)  && (i <= upperFractionX_ * (nRows-1) + 1) && (j >= lowerFractionY_ * (nCols-1) - 1) && (j <= upperFractionY_ * (nCols-1) + 1)){
                boolMatrix(i, j) = std::complex<double>(1.0, 0.0) * amplificationFactor_;
            } else {
                boolMatrix(i, j) = std::complex<double>(0.0, 0.0);
            }
        }
    }
    return boolMatrix;
}   

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> BandFiltering::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const{
    int nRows = matrix.rows();
    int nCols = matrix.cols();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> resultMatrix(nRows, nCols);
    resultMatrix = matrix.cwiseProduct(this->getBandMatrix(nRows, nCols));
    return resultMatrix;
}

//////////////////////////////////////////////


Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
ProbDensity::computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const 
{
    std::vector<std::complex<double>> hist(256, std::complex<double>(0.0, 0.0));
    int N = matrix.rows()*matrix.cols();
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            int val = static_cast<int>(std::round(std::real(matrix(i,j))));
            val = std::clamp(val, 0, 255); 
            hist[val] += 1.0;
        }
    }

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> out(2, 256);
    for (int k = 0; k < 256; ++k) {
        out(0, k) = k;      
        out(1, k) = hist[k]; 
    }

    return out;
}