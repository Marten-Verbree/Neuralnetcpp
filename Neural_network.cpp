#include <cmath>
#include <initializer_list>
#include <iostream>
#include <list>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>     

template <typename T>
class Matrix
{
    public:

        Matrix()
        {
            rows = 0;
            cols = 0;
            data = nullptr;
            //std::cout << "Default construct" << std::endl; 
        }

        Matrix(int r, int c)
        {
            rows = r;
            cols = c;
            const int size = r * c;
            data = new T[size]{};
            //std::cout << "Size construct" << std::endl; 
        }

        ~Matrix()
        {
            rows = 0;
            cols = 0;
            delete[] data;
            //std::cout << "Deconstruct" << std::endl; 
        }

        Matrix(int r, int c, const std::initializer_list<T>& list)
        :Matrix(r, c)
            {
                int size = list.size();
                if (r * c != size){throw "List size differs from given dimensions.";}
                std::uninitialized_copy(list.begin(),list.end(),data);
                //std::cout << "List construct" << std::endl; 
            }

        Matrix(const Matrix& other)
        :Matrix(other.rows, other.cols)
        {
            for(auto i = 0; i < other.rows*other.cols; i++)
                data[i] = other.data[i];
            //std::cout << "Copy construct" << std::endl; 
        }

        Matrix (Matrix&& other)
        :Matrix(other.rows, other.cols)
        {
            delete[] data;
            data = other.data;
            other.rows = 0;
            other.cols = 0;
            other.data = nullptr;
            //std::cout << "Move construct" << std::endl; 
        }

        Matrix& operator=(const Matrix& other)
        {
            if (this != &other)
            {
                delete[] data;
                rows = other.rows;
                cols = other.cols;
                data = new T[rows*cols];
                for(auto i = 0; i < other.rows*other.cols; i++)
                    data[i] = other.data[i];
            }
            //std::cout << "Copy operator" << std::endl; 
            return *this;
        }

        Matrix& operator=(Matrix&& other)
        {
            if (this != &other)
            {
                delete[] data;
                rows = other.rows;
                cols = other.cols;
                data = other.data;
                other.rows = 0;
                other.cols = 0;
                other.data = nullptr;
            }
            //std::cout << "Move operator" << std::endl; 
            return *this;
        }

        T& operator[](const std::pair<int, int>& ij) 
        {
            if (std::get<0>(ij) >= rows || std::get<1>(ij) >= cols){throw "Index out of bounds.";}
            return data[std::get<1>(ij) + (std::get<0>(ij))*cols];
        }

        const T& operator[](const std::pair<int, int>& ij) const
        {
            if (std::get<0>(ij) >= rows || std::get<1>(ij) >= cols){throw "Index out of bounds.";}
            return data[std::get<1>(ij) + (std::get<0>(ij))*cols];
        }

        template<typename U>
        Matrix<typename std::common_type<T,U>::type> operator*(const U x) const 
        {
            Matrix<typename std::common_type<T,U>::type> A (rows, cols);
            for (auto i = 0; i < rows * cols; i++)
                A.data[i] = data[i]*x;
            return A;
        }

        template<typename U>
        Matrix<typename std::common_type<T,U>::type> operator*(const Matrix<U>& B) const 
        {
            if(cols != B.rows){throw "Matrices are incompatible for multiplication.";}
            Matrix<typename std::common_type<T,U>::type> result(rows,B.cols);
            for(auto i = 0; i < rows; i++)
                {
                    for(auto j = 0; j < B.cols; j++)
                        {
                            result.data[i*B.cols+j] = 0;
                            for(auto k = 0; k < cols; k++)
                                result.data[i*B.cols+j] += data[i*cols+k] * B.data[k*B.cols+j];
                        }
                }
            return result;
        }

        template<typename U>
        Matrix<typename std::common_type<T,U>::type> operator+(const Matrix<U>& B) const
        {
            if (B.rows == 1 && cols == B.cols)
            {
                Matrix<typename std::common_type<T,U>::type> result(rows, cols); 
                for(auto i = 0; i < rows; i++)
                    for(auto j = 0; j < cols; j++)
                        result.data[i*cols+j] = data[i*cols+j] + B.data[j];  
                return result;
            }
            else if(rows == B.rows && cols == B.cols)
            {
                Matrix<typename std::common_type<T,U>::type> result(rows, cols);
                for(auto i = 0; i < rows; i++)
                    for(auto j = 0; j < cols; j++)
                        result.data[i*cols+j] = data[i*cols+j] + B.data[i*cols+j];
                return result;
            }
            else{throw "Matrices are incompatible for addition.";}
            
        }

        template<typename U>
        Matrix<typename std::common_type<T,U>::type> operator-(const Matrix<U>& B) const
        {
            if (B.rows == 1 && cols == B.cols)
            {
                Matrix<typename std::common_type<T,U>::type> result(rows, cols); 
                for(auto i = 0; i < rows; i++)
                    for(auto j = 0; j < cols; j++)
                        result.data[i*cols+j] = data[i*cols+j] - B.data[j]; 
                return result;
            }
            else if(rows == B.rows && cols == B.cols)
            {
                Matrix<typename std::common_type<T,U>::type> result(rows, cols);
                for(auto i = 0; i < rows; i++)
                    for(auto j = 0; j < cols; j++)
                        result.data[i*cols+j] = data[i*cols+j] - B.data[i*cols+j];
                return result;
            }
            else{throw "Matrices are incompatible for subtraction.";}        
        }

        Matrix transpose() const 
        {
            Matrix<T> result(cols,rows);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    result.data[i + j*rows] = data[j + i *cols];
            return result;
        }

        int getRows() const 
        {
            return rows;
        }
        
        int getCols() const 
        {
            return cols;
        }
        
        int rows;
        int cols;
        T* data;
};

template<typename T>
class Layer
{
    virtual Matrix<T> forward(const Matrix<T>& x) = 0;
    virtual Matrix<T> backward(const Matrix<T>& dy) = 0;
};

template <typename T>
class Linear: public Layer<T>
{
    public:
        Linear(int in_features, int out_features, int n_samples, int seed) 
        {
            input = in_features;
            output = out_features;
            n = n_samples;
            s = seed;

            std::default_random_engine        generator(seed);
            std::normal_distribution<T>       distribution_normal(0.0, 1.0);
            std::uniform_real_distribution<T> distribution_uniform(0.0, 1.0);

            cache = Matrix<T>(n,input);
            bias = Matrix<T>(1,output);
            for (int i=0; i<bias.rows; ++i) 
            {
                for (int j=0; j<bias.cols; ++j) 
                {
                    bias.data[i*bias.cols+j] = distribution_uniform(generator);
                }
            }
            weights = Matrix<T>(input,output);
            for (int i=0; i<weights.rows; ++i) 
            {
                for (int j=0; j<weights.cols; ++j) 
                {
                    weights.data[i*weights.cols+j] = distribution_normal(generator);
                }
            }
            b_grad = Matrix<T>(1,output);
            w_grad = Matrix<T>(input,output);
        }
        
        ~Linear()
        {
            input = 0;
            output = 0;
            n = 0;
            s = 0;
        }

        virtual Matrix<T> forward(const Matrix<T>& x) override final 
        {
            Matrix<T> y = x * weights + bias;
            cache = x;
            return y;
        }

        virtual Matrix<T> backward(const Matrix<T>& dy) override final
        {
            Matrix<T> dx = dy * weights.transpose();
            w_grad = cache.transpose() * dy; 
            Matrix<T> one_vector(1,n);
            for(auto i = 0; i < n; i++)
               one_vector.data[i] = 1;
            b_grad = one_vector * dy;
            return dx;
        }

        void optimize(T learning_rate) 
        {
            weights = weights - w_grad * learning_rate;
            bias = bias - b_grad * learning_rate;
        }

        int input;
        int output;
        int n;
        int s;
        Matrix<T> cache;
        Matrix<T> bias;
        Matrix<T> weights;
        Matrix<T> b_grad;
        Matrix<T> w_grad;
};

template <typename T>
class ReLU: public Layer<T>
{
    public:
        ReLU(int in_features, int out_features, int n_samples) 
        {
            input = in_features;
            output = out_features;
            n = n_samples;

            cache = Matrix<T>(n,input); 
        }
        
        ~ReLU()
        {
            input = 0;
            output = 0;
            n = 0;
        }

        virtual Matrix<T> forward(const Matrix<T>& x) override final 
        {
            // We were told we can assume the input and output to be of equal dimensions.
            Matrix<T> y(x.rows,x.cols);
            T zero = 0;
            for(auto i = 0; i < x.rows * x.cols; i++)
                y.data[i] = std::max(zero,x.data[i]); 
            cache = x;
            return y;
        }

        virtual Matrix<T> backward(const Matrix<T>& dy) override final
        {
            Matrix<T> dx(dy.rows, dy.cols);
            for(auto i = 0; i < dy.rows*dy.cols; i++)
            {
                if(cache.data[i] > 0){dx.data[i] = dy.data[i];} //because dy * 1 = dy
                else{dx.data[i] = 0;}
            }    
            return dx;
        }

        int input;
        int output;
        int n;
        Matrix<T> cache;
};

template <typename T>
class Net 
{
    public:
        Net(int in_features, int hidden_dim, int out_features, int n_samples, int seed)
        :input_layer(in_features, hidden_dim, n_samples, seed), 
        hidden_layer(hidden_dim, hidden_dim, n_samples),
        output_layer(hidden_dim, out_features, n_samples, seed)
        {
            input = in_features;
            output = out_features;
            n = n_samples;
            s = seed;
            hidden = hidden_dim;    
        }

        ~Net()
        {
            input = 0;
            output = 0;
            n = 0;
            s = 0;
            hidden = 0;
        }

        Matrix<T> forward(const Matrix<T>& x) 
        {
            Matrix<T> input_to_hidden = input_layer.forward(x);
            Matrix<T> hidden_to_output = hidden_layer.forward(input_to_hidden);
            Matrix<T> y = output_layer.forward(hidden_to_output);
            return y;
        }

        Matrix<T> backward(const Matrix<T>& dy)
        {
            Matrix<T> output_to_hidden = output_layer.backward(dy);
            Matrix<T> hidden_to_input = hidden_layer.backward(output_to_hidden);
            Matrix<T> dx = input_layer.backward(hidden_to_input);
            return dx;
        }

        void optimize(T learning_rate) 
        {
            input_layer.optimize(learning_rate);
            output_layer.optimize(learning_rate);
        }

        int input;
        int output;
        int n;
        int s;
        int hidden;
        Linear<T> input_layer;
        ReLU<T> hidden_layer;
        Linear<T> output_layer; 
};

// Function to calculate the loss
template <typename T>
T MSEloss(const Matrix<T>& y_true, const Matrix<T>& y_pred) 
{
    if(y_true.cols != y_pred.cols || y_true.rows != y_pred.rows)
        {throw "Given output values have a different dimension from predicted values.";}
    int n = y_true.cols * y_true.rows;
    T result = 0;
    for(auto i = 0; i < n; i++)
        {
            result += std::pow(y_pred.data[i] - y_true.data[i],2);
        }
    result = result/float(n);
    return result;
}

// Function to calculate the gradients of the loss
template <typename T>
Matrix<T> MSEgrad(const Matrix<T>& y_true, const Matrix<T>& y_pred) 
{
    if(y_true.cols != y_pred.cols || y_true.rows != y_pred.rows)
        {throw "Given output values have a different dimension from predicted values.";}
    Matrix<T> grad (y_pred.rows,y_pred.cols);
    int n = y_pred.cols * y_pred.rows;
    for (int i = 0; i<n;i++)
    {
        grad.data[i] = (2*y_pred.data[i]-2*y_true.data[i]);
    }
    return grad;
}

// Calculate the argmax 
template <typename T>
Matrix<T> argmax(const Matrix<T>& y) 
{
    Matrix<T> result(1,y.rows);
    for(auto i = 0; i < y.rows; i++)
        {
            int index_max = 0;
            for(auto j = 1; j < y.cols; j++)
                {
                    if(y.data[i*y.cols+index_max] < y.data[i*y.cols+j]){index_max = j;}
                }
            result.data[i] = index_max;
        }
    return result; 
}

// Calculate the accuracy of the prediction, using the argmax
template <typename T>
T get_accuracy(const Matrix<T>& y_true, const Matrix<T>& y_pred)
{
    if(y_true.cols != y_pred.cols || y_true.rows != y_pred.rows)
        {throw "Given output values have a different dimension from predicted values.";}
    T result = 0;
    for(auto i = 0; i < y_true.rows; i++)
        {
            if(argmax(y_true).data[i] == argmax(y_pred).data[i]){result += 1;}
        }
    return result/(float)y_true.rows;
}

template <typename T>
void print_matrix(Matrix<T> M)
{
    std::cout << "Size: " << M.rows << "x" << M.cols << std::endl; 
    for (auto i = 0; i < M.rows; i++)
        {
            for (auto j = 0; j < M.cols; j++)
                std::cout << M.data[i*M.cols+j] << " ";
            std::cout << std::endl;
        }
    std::cout << std::endl;
}

int main(int argc, char* argv[])
{
    // Your training and testing of the Net class starts here
    float learning_rate = 0.0005;
    int optimizer_steps = 100;
    std::initializer_list<float> list = {0,0,0,1,1,0,1,1};
    std::initializer_list<float> list1 = {1,0,0,1,0,1,1,0};
    Matrix<float> x_xor(4,2,list);
    Matrix<float> y_xor(4,2,list1);
    Net<float> net(2,100,2,4,1);
    Matrix<float>loss_list(1,optimizer_steps);
    Matrix<float>accuracy_list(1,optimizer_steps);
    for(auto i = 0; i < optimizer_steps; i++)
        {
            Matrix<float> y_pred(net.forward(x_xor));
            loss_list.data[i] = MSEloss(y_xor,y_pred);
            Matrix<float> x_label = net.backward(MSEgrad(y_xor,y_pred));
            net.optimize(learning_rate);
            accuracy_list.data[i] = get_accuracy(y_xor,y_pred);
            if(i % 10 == 0)
                {
                    std::cout << "Predicted output at step " << i << ":" << std::endl;
                    print_matrix(y_pred);
                } 
        }
    std::cout << "Losses:" << std::endl;    
    print_matrix(loss_list);
    std::cout << "Accuracies:" << std::endl;
    print_matrix(accuracy_list);
    return 0;
}
