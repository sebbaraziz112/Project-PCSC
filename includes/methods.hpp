/**
 * @file methods.h
 * @brief Déclaration des différentes méthodes de transformation et convolution.
 *
 * Contient les classes pour les méthodes de DFT, FFT, BlueStein, ainsi que des convolutions classiques et avancées.
 */

#ifndef METHODS
#define METHODS

#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <cassert>
#include <cmath>

/**
 * @class AbstractMethod
 * @brief Classe de base abstraite pour toutes les méthodes de transformation.
 *
 * Fournit des interfaces pour la conversion, le calcul et la manipulation des matrices complexes.
 */
class AbstractMethod
{
public: 
    /**
     * @brief Constructeur.
     * @param name Nom de la méthode.
     */
    AbstractMethod(std::string name): name_(name){};

    /**
     * @brief Destructeur.
     */
    ~AbstractMethod(){};

    /**
     * @brief Convertit un vecteur de matrices en matrices complexes.
     * @tparam T Type original des matrices.
     * @param matrices Vecteur de matrices à convertir.
     * @return Vecteur de matrices complexes.
     */
    template <typename T>
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> ComplexCast(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const;

    /**
     * @brief Convertit un vecteur de matrices complexes en un autre type.
     * @tparam O Type cible des matrices.
     * @param matrices Vecteur de matrices complexes.
     * @return Vecteur de matrices converties.
     */
    template <typename O>
    std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> OCast(std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const;

    /**
     * @brief Applique la méthode sur un vecteur de matrices.
     * @tparam T Type des matrices d'entrée.
     * @tparam O Type des matrices de sortie.
     * @param matrices Vecteur de matrices d'entrée.
     * @return Vecteur de matrices transformées.
     */
    template <typename T, typename O>
    std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> computeMethod(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const;

    /**
     * @brief Méthode abstraite : calcule la transformation sur une seule matrice.
     * @param matrix Matrice d'entrée.
     * @return Matrice transformée.
     */
    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const = 0;

    /**
     * @brief Retourne la prochaine puissance de deux supérieure ou égale à M.
     * @param M Nombre de départ.
     * @return Prochaine puissance de deux.
     */
    int getNextPowerTwo(int M) const;

    /**
     * @brief Retourne les facteurs premiers de N.
     * @param N Nombre entier.
     * @return Vecteur de facteurs premiers.
     */
    std::vector<int> getFactors(int N) const;

    /**
     * @brief Convertit un vecteur de matrices en vecteur de matrices selon les facteurs.
     * @tparam T Type des matrices.
     * @param vector_matrices Vecteur de matrices d'entrée.
     * @param factors Facteurs utilisés pour la conversion.
     * @return Vecteur de matrices converties.
     */
    template <typename T>
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> vectorToMatrix(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& vector_matrices, std::vector<int>& factors) const;

    /**
     * @brief Convertit un vecteur de matrices en vecteur linéarisé.
     * @tparam T Type des matrices.
     * @param matrices Vecteur de matrices à transformer.
     * @return Vecteur de matrices linéarisé.
     */
    template <typename T>
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> matrixToVector(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const;

    std::string name_; /**< Nom de la méthode */

protected: 
};

/**
 * @class IdentityMethod
 * @brief Méthode identité (retourne la matrice telle quelle).
 */
class IdentityMethod: public AbstractMethod
{
public: 
    IdentityMethod(): AbstractMethod("Identity Method"){};
    ~IdentityMethod(){};

    /**
     * @brief Applique la méthode identité sur la matrice.
     * @param matrix Matrice d'entrée.
     * @return La même matrice que celle d'entrée.
     */
    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

private: 
};

/**
 * @class NaiveDFTMethod
 * @brief Implémentation naïve de la DFT (Discrete Fourier Transform).
 */
class NaiveDFTMethod: public AbstractMethod
{
public: 
    NaiveDFTMethod(): AbstractMethod("DFT Method"){};
    ~NaiveDFTMethod(){};

    /**
     * @brief Calcule la DFT de la matrice.
     * @param matrix Matrice d'entrée.
     * @return Matrice transformée.
     */
    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

private: 
    /**
     * @brief Calcule un élément (u,v) de la DFT.
     * @param u Index ligne.
     * @param v Index colonne.
     * @param matrix Matrice d'entrée.
     * @return Élément de la DFT.
     */
    std::complex<double> getResultElement(int u, int v, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const;
};

/**
 * @class NaiveInvDFTMethod
 * @brief Implémentation naïve de la DFT inverse.
 */
class NaiveInvDFTMethod: public AbstractMethod
{
public: 
    NaiveInvDFTMethod(): AbstractMethod("Inverse DFT Method"){};
    ~NaiveInvDFTMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

private: 
    std::complex<double> getResultElement(int u, int v, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const;
};

/**
 * @class FFT1DMethod
 * @brief FFT 1D avec option de padding.
 */
class FFT1DMethod: public AbstractMethod
{
public: 
    /**
     * @brief Constructeur.
     * @param padd Indique si un padding est appliqué.
     */
    FFT1DMethod(bool padd): AbstractMethod("FFT Method"), padd_(padd){};
    ~FFT1DMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

    std::vector<std::complex<double>> getCloserPowerTwo(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const;
    std::vector<std::complex<double>> fftCompute(std::vector<std::complex<double>> vector_c) const;
    std::vector<std::complex<double>> getEven(std::vector<std::complex<double>>& vector_c) const;
    std::vector<std::complex<double>> getOdd(std::vector<std::complex<double>>& vector_c) const;

private:
    bool padd_; /**< Padding activé ou non */
};

/**
 * @class Inv1DFFTMethod
 * @brief FFT inverse 1D avec option de padding.
 */
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

/**
 * @class BlueStein1DMethod
 * @brief BlueStein 1D method pour les transformations FFT de taille non-puissance de deux.
 */
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

/**
 * @class InvBlueStein1DMethod
 * @brief BlueStein inverse 1D method.
 */
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

/**
 * @class BlueSteinMethod
 * @brief BlueStein méthode ND (1D et 2D).
 */
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

/**
 * @class InvBlueSteinMethod
 * @brief BlueStein inverse ND method.
 */
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


/**
 * @class LineOnlyBlueSteinMethod
 * @brief BlueStein method appliquée uniquement sur les lignes.
 */
class LineOnlyBlueSteinMethod: public AbstractMethod
{
public: 
    /**
     * @brief Constructeur.
     */
    LineOnlyBlueSteinMethod(): AbstractMethod("Line Only BlueStein Method"){
        bluestein1d_ = std::make_shared<BlueStein1DMethod>();
    }

    /**
     * @brief Destructeur.
     */
    ~LineOnlyBlueSteinMethod(){};

    /**
     * @brief Applique la méthode LineOnlyBlueStein sur une matrice.
     * @param matrix Matrice d'entrée.
     * @return Matrice transformée.
     */
    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

private:
    std::shared_ptr<BlueStein1DMethod> bluestein1d_; /**< Méthode BlueStein 1D utilisée pour les lignes */
};

/**
 * @class LineOnlyInvBlueSteinMethod
 * @brief BlueStein inverse appliquée uniquement sur les lignes.
 */
class LineOnlyInvBlueSteinMethod: public AbstractMethod
{
public: 
    LineOnlyInvBlueSteinMethod(): AbstractMethod("Line Only Inverse BlueStein Method"){
        invbluestein1d_ = std::make_shared<InvBlueStein1DMethod>();
    }

    ~LineOnlyInvBlueSteinMethod(){};

    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

private:
    std::shared_ptr<InvBlueStein1DMethod> invbluestein1d_; /**< Méthode inverse BlueStein 1D utilisée pour les lignes */
};

/**
 * @class Convolution
 * @brief Classe template pour la convolution 2D de matrices complexes.
 * @tparam DIM Dimension de la matrice de convolution (kernel).
 */
template <int DIM>
class Convolution: public AbstractMethod
{
public: 
    Convolution(): AbstractMethod("Convolution Method"){};
    ~Convolution(){};

    /**
     * @brief Applique la convolution sur une matrice.
     * @param matrix Matrice d'entrée.
     * @return Matrice convoluée.
     */
    virtual Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

    /**
     * @brief Définit un élément du kernel de convolution.
     * @param complexValue Valeur complexe.
     * @param i Index ligne.
     * @param j Index colonne.
     */
    void setKernel(std::complex<double> complexValue, int i, int j) {
        kernel_(i, j) = complexValue;
    }

    /**
     * @brief Normalise le kernel (somme des coefficients = 1).
     */
    void normalizeKernel(){
        kernel_ = kernel_ / kernel_.sum();
    }

    Eigen::Matrix<std::complex<double>, DIM, DIM> kernel_; /**< Kernel de convolution */

protected: 
    /**
     * @brief Calcule la valeur convoluée pour un élément spécifique.
     * @param matrix Matrice d'entrée.
     * @param i Index ligne.
     * @param j Index colonne.
     * @return Valeur convoluée.
     */
    std::complex<double> convolve(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix, int i, int j) const;
};

/**
 * @class MeanConv
 * @brief Convolution moyenne (filtre moyen) utilisant la classe Convolution.
 * @tparam DIM Dimension du kernel.
 */
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

/**
 * @class BlurrGaussianConv
 * @brief Convolution avec un filtre gaussien.
 * @tparam DIM Dimension du kernel.
 */
template <int DIM>
class BlurrGaussianConv: public Convolution<DIM>
{
public: 
    /**
     * @brief Constructeur.
     * @param sigmaX Écart type horizontal.
     * @param sigmaY Écart type vertical.
     */
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

    /**
     * @brief Calcule la valeur de la fonction gaussienne pour un point (x,y).
     * @param x Coordonnée x.
     * @param y Coordonnée y.
     * @param xc Centre x.
     * @param yc Centre y.
     * @return Valeur gaussienne.
     */
    double getGaussian(int x, int y, int xc, int yc){
        return (1.0 / (2 * M_PI * sigmaX_ * sigmaY_)) * std::exp((-std::pow((x - xc)/sigmaX_, 2) - std::pow((y - yc)/sigmaY_, 2))/2.0);
    }

    double sigmaX_; /**< Écart type horizontal */
    double sigmaY_; /**< Écart type vertical */
};

/**
 * @class SobelX
 * @brief Filtre Sobel horizontal (3x3) pour détection de contours.
 */
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


/**
 * @class SobelY
 * @brief Filtre Sobel vertical (3x3) pour la détection de contours.
 * 
 * La matrice de convolution applique le masque Sobel vertical :
 * \f[
 * \begin{bmatrix}
 *  1 & 2 & 1 \\
 *  0 & 0 & 0 \\
 * -1 & -2 & -1
 * \end{bmatrix}
 * \f]
 */
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

/**
 * @class LaplacianClassical
 * @brief Filtre de Laplacien classique (3x3) pour la détection des contours.
 * 
 * La matrice de convolution applique le masque :
 * \f[
 * \begin{bmatrix}
 * 0 & -1 & 0 \\
 * -1 & 4 & -1 \\
 * 0 & -1 & 0
 * \end{bmatrix}
 * \f]
 */
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

/**
 * @class Laplacian8Connected
 * @brief Filtre Laplacien 8-connecté (3x3) pour détection de contours plus précis.
 * 
 * La matrice de convolution applique le masque :
 * \f[
 * \begin{bmatrix}
 * -1 & -1 & -1 \\
 * -1 & 8 & -1 \\
 * -1 & -1 & -1
 * \end{bmatrix}
 * \f]
 */
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

/**
 * @class Laplacian4Connected
 * @brief Filtre Laplacien 4-connecté (3x3) pour détection de contours avec moins de diagonales.
 * 
 * La matrice de convolution applique le masque :
 * \f[
 * \begin{bmatrix}
 * 0 & 1 & 0 \\
 * 1 & -4 & 1 \\
 * 0 & 1 & 0
 * \end{bmatrix}
 * \f]
 */
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

/**
 * @class LaplacianOfGaussian
 * @brief Filtre Laplacien du Gaussien (LoG) (5x5) pour détection de contours avec pré-smoothing.
 * 
 * Ce filtre combine un lissage gaussien et un opérateur Laplacien pour détecter les contours.
 */
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
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> shift(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices, int rowPad, int colPad) const;

    template <typename T>
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> shiftCenter(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const{
        int nRows = matrices[0].rows();
        int nCols = matrices[0].cols();
        return this->shift<T>(matrices, nRows / 2, nCols / 2);
    }

    template <typename T>
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> shiftInverse(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const{
        int nRows = matrices[0].rows();
        int nCols = matrices[0].cols();
        int rowPad = - nRows / 2;
        int colPad = - nCols / 2;
        return this->shift<T>(matrices, rowPad, colPad);
    }
    
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
/**
 * @class LowPass
 * @brief Filtre passe-bas basé sur la classe BandFiltering.
 * 
 * Ce filtre conserve les basses fréquences et atténue les hautes fréquences.
 * Les fractions inférieures X et Y sont initialisées à 0.0 et le facteur
 * d'amplification à 1.0.
 */
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

/**
 * @class HighPass
 * @brief Filtre passe-haut basé sur la classe BandFiltering.
 * 
 * Ce filtre conserve les hautes fréquences et atténue les basses fréquences.
 * Les fractions supérieures X et Y sont initialisées à 0.0 et le facteur
 * d'amplification à 1.0.
 */
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

template <typename T>
std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> AbstractMethod::vectorToMatrix(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& vector_matrices, std::vector<int>& factors) const{
    const int nRows = vector_matrices[0].rows();
    const int nCols = vector_matrices[0].cols();
    const int nChann = vector_matrices.size();
    assert(nRows == 1);
    assert(factors[0] * factors[1] == nCols);
    
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> matrix_matrices;
    matrix_matrices.reserve(nChann);

    for (int chann = 0; chann < nChann; chann ++){
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat(factors[0], factors[1]);
        for (int i = 0; i < factors[0]; i ++){
            for (int j = 0; j < factors[1]; j ++){
                mat(i, j) = vector_matrices[chann](0, i * factors[1] + j);
            }
        }
        matrix_matrices.push_back(mat);
    }
    return matrix_matrices;
}
    
template <typename T>
std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> AbstractMethod::matrixToVector(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) const{
    const int nChann = matrices.size();
    const int nRows = matrices[0].rows();
    const int nCols = matrices[0].cols();
    const int N = nRows * nCols;

    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> vector_matrices;
    vector_matrices.reserve(nChann);
    for (int chann = 0; chann < nChann; chann ++){
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vec(1, N);
        for (int i = 0; i < nRows; i ++){
            for (int j = 0; j < nCols; j ++){
                vec(0, i * nCols + j) = matrices[chann](i, j);
            }
        }
        vector_matrices.push_back(vec);
    }
    return vector_matrices;
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
std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> BandFiltering::shift(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices, int rowPad, int colPad) const {
    const int nChann = matrices.size();
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> shiftedMatrices;
    shiftedMatrices.reserve(nChann);

    for (int chann = 0; chann < nChann; chann++) {
        const auto& mat = matrices[chann];
        const int nRows = mat.rows();
        const int nCols = mat.cols();

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> newMatrix(nRows, nCols);

        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                int newRow = (i + rowPad) % nRows;
                int newCol = (j + colPad) % nCols;

                if (newRow < 0) newRow += nRows;
                if (newCol < 0) newCol += nCols;

                newMatrix(newRow, newCol) = mat(i, j);
            }
        }

        

        shiftedMatrices.push_back(newMatrix);
    }

    return shiftedMatrices;
}


//////////////////////////////////////

// Probability density of intensity values 

class ProbDensity: public AbstractMethod
{

public:
    
    ProbDensity(): AbstractMethod("Probability Density Method"){};
    ~ProbDensity(){};

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> computeMatrixMethod(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& matrix) const override;

};


#endif