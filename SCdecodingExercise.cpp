/***************************
 SC decoding exercise

 Carlo Condo, PhD - 2017
 carlo.condo@mail.mcgill.ca
****************************/

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <limits>

using namespace std;

// functions list
void PolarCodeReadReliabilityFile(ifstream &file, vector<unsigned int>& reliability);
void PolarCodeEncoding(const vector<unsigned char>& u, vector<unsigned char>& x);
float sign(float x);
float f(float a, char b);
float g(float a, char b, float s);
//void SCdecoding(const vector<unsigned char>& alpha, vector<unsigned char>& beta, vector<unsigned char>& frozen);
void SCdecoding(const vector<float>& alpha, vector<unsigned char>& beta, vector<unsigned char>& frozen);

unsigned int n_leaf = 0;

int main()
{
	// code length
    const unsigned int N1 = 128;
    //const unsigned int N1 = 4;

    // number of information bits
    const unsigned int K1 = 16;
    //const unsigned int K1 = 120;

	// code rate
	double R1 = (double)K1 / (double)N1;

	// signal-to-noise ratio (SNR) values
    vector<double> SNR = {-1.2,-1.0,-0.8};
    //vector<double> SNR = {-10.2};

	// SNR standard deviation
	vector<double> stddev;
	for (unsigned int i = 0; i < SNR.size(); ++i) {
		stddev.push_back(sqrt(pow(10, -SNR[i] / 10)));
	}

	// polar code reliability file
	// sorted list of bit positions, from the most to the least reliable
	ifstream file1("n128_awgn_s0.5.pc");

    //output file
    ofstream file2("SCresults.txt");

	// maximum number of simulated frames
	const unsigned int maxIter = 10000;

	// minimum number of frame errors before the simulation is interrupted
	const unsigned int minError = 50;
	
	vector<unsigned int> reliability1;
	PolarCodeReadReliabilityFile(file1, reliability1);

	vector<unsigned char> frozen1(N1);
	vector<unsigned char> u1(N1);
	vector<unsigned char> x1(N1);
    //vector<unsigned char> y1(N1);
    vector<float> y1(N1);
	vector<unsigned char> u_hat1(N1);

    vector<double> BER(SNR.size());
    vector<double> FER(SNR.size());


	// information bits (all-zero codeword)
	for (unsigned int i = 0; i < K1; ++i) {
		frozen1[reliability1[i]] = 0;
		u1[reliability1[i]] = 0;
		}

	// frozen bits (0)
	for (unsigned int i = K1; i < N1; ++i) {
		frozen1[reliability1[i]] = 1;
		u1[reliability1[i]] = 0;
		}

	// for every SNR
	for (unsigned int SNR_cnt = 0; SNR_cnt < SNR.size(); ++SNR_cnt) {
		clock_t startTime = clock();
		unsigned int iter = 0;
        unsigned int n_err_FER = 0;
        unsigned int n_err_BER = 0;
        unsigned int n_err_tmp;
        //unsigned int n_leaf = 0;
	
        while ((iter < maxIter) || (n_err_FER < minError)) {
        //while ((iter < maxIter) && (n_err_FER < minError)) {

            n_leaf = 0;

			// encoding
			PolarCodeEncoding(u1, x1);

			// BPSK modulation and noise
			unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
			default_random_engine generator_n(seed);
			normal_distribution<double> distribution_n(0, 1);

			for (unsigned int i = 0; i < N1; ++i)
				y1[i] = 1 - 2 * x1[i] + stddev[SNR_cnt] * distribution_n(generator_n);

			// SC decoding
            vector<unsigned char> s1(N1);
            SCdecoding(y1,s1, frozen1);

			// Error rate calculation

            n_err_tmp = n_err_BER;  //BER of the previous frame
            for (unsigned int i=0; i<N1; i++)
            {
                if(s1[i] != x1[i])
                    n_err_BER++;    //updated BER
            }

            if (n_err_BER != n_err_tmp)
                n_err_FER++;        //updated FER

            iter++;
		} //iter loop

        BER[SNR_cnt] = (double)n_err_BER / (double)(iter*N1);
        FER[SNR_cnt] = (double)n_err_FER / (double)(iter);

        file2 << "SNR: \n";
        file2 << SNR[SNR_cnt];
        file2 << "\nBER: \n";
        file2 << BER[SNR_cnt];
        file2 << "\nFER: \n";
        file2 << FER[SNR_cnt];
        file2 << "\n\n";


	}	//SNR loop


    file1.close();
    file2.close();
	return 0;
}

//********************************************************************************************************************

void PolarCodeReadReliabilityFile(ifstream &file, vector<unsigned int>& reliability)
{
	string line;

	unsigned int value;

	// Read one line at a time into the variable line:
	if (file.is_open())
	{
		// first line contains code length
		getline(file, line);
		// second line ignored
		getline(file, line);
		// third line ignored
		getline(file, line);
		// forth line contains reliability values from the most to the least reliable bit
		getline(file, line);
		stringstream  lineStream(line);
		while (lineStream >> value)
		{
			// Add the integers from a line to a 1D array (vector)
			reliability.push_back(value);
		}
		file.close();
	}
	else cout << "Unable to open file"<<endl;
}


void PolarCodeEncoding(const vector<unsigned char>& u, vector<unsigned char>& x)
{
	unsigned int N = static_cast<unsigned int>(u.size());
	unsigned int log2N = 0;
	unsigned int N_tmp = N;
	while (N_tmp >>= 1) ++log2N;
	x = u;
	for (unsigned int i = 0; i < log2N; ++i) {
		for (unsigned int j = 0; j < (N >> (i + 1)); ++j) {
			for (unsigned int l = (j << (i + 1)); l < (((j << 1) + 1) << i); ++l) {
				x[l] ^= x[l + (1 << i)];
			}
		}
	}
}


float sign(float x)
{
    if ( x <0)
        return -1;
    else
        return 1;
}


char f(float a, float b)
{
    char sign_a = sign(a);
    char sign_b = sign(b);
    char abs_a = abs(a);
    char abs_b = abs(b);
    char min_ab = min(abs_a, abs_b);
    char f = sign_a*sign_b*min_ab;
    return f;
}


float g(float a, char b, float s)
{
    return (1-2*s)*a+b;
}


void SCdecoding(const vector<float>& alpha, vector<unsigned char>& beta, vector<unsigned char>& frozen)
{
    unsigned int N = static_cast<unsigned int>(alpha.size());
    unsigned int log2N = 0;
    unsigned int N_tmp = N;
    while (N_tmp >>= 1) ++log2N;

    if(log2N>=1)    //not leaf
    {
        vector<float> alpha_l(N/2);
        vector<float> alpha_r(N/2);
        vector<unsigned char> beta_l(N/2);
        vector<unsigned char> beta_r(N/2);

        //compute f and return alpha_l
        for (unsigned int i=0; i<N/2; i++)
            alpha_l[i] = f(alpha[i], alpha[i+pow(2,log2N-1)]);

        //SCdecoding of the left child
        SCdecoding(alpha_l, beta_l, frozen);

        //compute g and return alpha_r
        for (unsigned int i=0; i<N/2; i++)
            alpha_r[i] = g(alpha[i], alpha[i+pow(2,log2N-1)], beta_l[i]);

        //SCdecoding of the right child
        SCdecoding(alpha_r, beta_r, frozen);

        //compute the partial sum
        for (unsigned int i=0; i<N; i++)
        {
            if (i < N/2)
                beta[i] = beta_l[i] xor beta_r[i];
            else
                beta[i] = beta_r[i-pow(2,log2N-1)];
        }
    }
    else
    {
        if (frozen[n_leaf]==1)
            beta[0] = 0;
        else
        {
            if (alpha[0]>=0)
                beta[0]=0;
            else
                beta[0]=1;
        }
        n_leaf++;
    }
}




