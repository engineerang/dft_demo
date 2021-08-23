#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>

std::vector<std::complex<double> > dft(std::vector<std::complex<double> > X)
{
    // Determine number of samples. For simplicity, the number of input vs output samples
    // is the same i.e. our DFT is the length of the input signal.
    int N = X.size(); // N = Total number of samples from input singal X
    int K = N; // K = Total number of outputs of the DFT 

    // Allotcate memory for internals
    std::complex<double> intSum;

    // Allotcate memory for output
    std::vector<std::complex<double> > output;
    output.reserve(K);

    // Loop through each k (k = the index/bin of the DFT output in the frequency domain)
    for (int k=0; k<K; k++)
    {
        /* 
        Loop through each sample (n) from 0 to N-1 and plot a complex evaulation singal (w)
        with frequency, k. The output (intSum) will be a complex number that is a summation of 
        the product of the evaulation signal (w) and input signal (X) for each sample, n.
        */

        /* Note: The complex evaulation singal (w) is the signal we are comparing the input
        signal samples (X) against. The better the frequency (k) match between w & X the larger the
        resultant summation (intSum) will be. Hence a "response" in the index/bin (k) of the 
        DFT output in the frequency domain.
        */

        // reset the summation of n values for each k
        intSum = std::complex<double>(0, 0);

        for (int n=0; n<N; n++)
        {
            // Find amplitude value (y-axis) of a Cosine and Sine signal at sample n (x-axis)
            // with k cycles within 2*PI/N (total n points on x-axis normalised by N)
            // e.g. when k=1,n=0 : realPart = 1, imagePart = 0
            double realPart = cos(((2*M_PI)/N) *k * n);
            double imagPart = sin(((2*M_PI)/N) *k * n);
            // create complex signal, w, from the resulting amplitude values at sample n.
            // e.g. when k=1,n=0 : w = 1 - j0
            std::complex<double> w (realPart, -imagPart); // w = evaulation signal.

            //Note the "negative" sign from the DFT equation.

            // Multiple the complex values of our input signal, X, and the evaulation signal, w.
            // and add to summation of all n samples for value k.
            intSum += X[n] * w; // X[n] = the sequence of input samples from 0 to N-1.
        }

        // append the resultant summation (intSum) to the DFT output vector of size K (0 to K-1).
        output.push_back(intSum);
    }
    return output;
}

int main()
{
    // N = Total number of samples of signal (input signal).
    // Note - this is also your sampling frequency. e.g. 256 samples/second 
    // which will include frequencys upto ~ 120hz
    int N = 256;
    std::vector<std::complex<double> > signal; // Input signal
    signal.reserve(N);

    double sigK = 10.0; // Signal Frequency
    double sigAmp = 1.0; // Signal Amplitude
    double sigPhase = 0.0; // Signal Phase
    //double sigPhase = M_PI / 2.0; // Try this phase out - 90 degrees out of phase
    //double sigPhase = M_PI / 4.0; // Try this phase out - 45 degrees out of phase

    // Create complex Input Signal by computing each sample (currentSample) from 0 to N-1. 
    for (int x=0; x<N; ++x)
    {
        auto currentSample = std::complex<double>(sigAmp*cos((2*M_PI/static_cast<double>(N)) * sigK * static_cast<double>(x) + sigPhase), 0.0);

        // append each sample to the signal output vector of size N (0 to N-1).
        signal.push_back(currentSample);
    }

    // Compute the DFT
    std::vector<std::complex<double> > Fx = dft(signal);    

    // Display the real and imaginary components + magnitude - output is normalised
    std::cout << "****" << std::endl;
    std::cout << "DFT outputs.." << std::endl;
    std::cout << "\n" << "k (bin)\t" << std::setw(12) << "Real\t" << std::setw(12) << "Imag" << std::setw(12) << "Mag\t" << std::endl;
    for (int i=0; i<N/2; ++i)
    {
        std::cout << i << "\t" 
            << std::setw(12) << 2*(Fx[i].real() / static_cast<double>(N)) << "\t"
            << std::setw(12) << 2*(Fx[i].imag() / static_cast<double>(N)) << "\t"
            << std::setw(12) << sqrt(pow(2*(Fx[i].real() / static_cast<double>(N)),2) + pow(2*(Fx[i].imag() / static_cast<double>(N)),2))
            << std::endl; 
    }
}