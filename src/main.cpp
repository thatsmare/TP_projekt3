#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <set>
#include <math.h>
#include <vector>
#include <iostream>
#include <string> 
#include "AudioFile.h"
#include <opencv2/opencv.hpp>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#ifndef NOMINMAX
#define NOMINMAX
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int add(int i, int j) {
    return i + j;
}

int multiply(int i, int j) {
    AudioFile<double> audioFile;

    audioFile.load("C:/Users/marty/Documents/GitHub/TP_projekt3/AudioFile/examples/test-audio.wav");

    int sampleRate = audioFile.getSampleRate();
    std::cout << sampleRate << std::endl;

    return 0;
}


void signal_visualization() {
    AudioFile<double> audiosample;

    std::string file_location = "C:/Users/marty/Documents/GitHub/TP_projekt3/AudioFile/examples/test-audio.wav";

    std::cout << "Location of the file to visualize:" << std::endl;     //mozliwosc dobory sciezki do pliku
    std::cin >> file_location;

    audiosample.load(file_location);

    // Częstotliwość próbkowania
    int sampleRate = audiosample.getSampleRate();
    std::vector<std::vector<double>> samples = audiosample.samples;

    // Utworzenie osi czasu dla próbek
    size_t sample_number = samples[0].size();
    std::vector<double> time(sample_number);
    for (size_t i = 0; i < sample_number; ++i)
        time[i] = static_cast<double>(i) / sampleRate;

    using namespace matplot;
    figure();
    plot(time, samples[0]);
    title("Wizualizacja sygnału audio");
    xlabel("Czas");
    ylabel("Amplituda");
    show();

}

void signal_generate_sinusoidal() {
    double amplitude, frequency, phase, duration;
    int sampleRate;
    bool sine = true;
    char decision;

    std::cout << "Change to cos? (Y/N)" << std::endl;
    std::cin >> decision;
    if(decision == 'Y' || decision == 'y'){
        sine = false;
    }
    std::cout << "Amplitude [unit]: ";
    std::cin >> amplitude;
    std::cout << "Frequency [Hz]: ";
    std::cin >> frequency;
    std::cout << "Phase [degrees]: ";
    std::cin >> phase;
    std::cout << "Duration [s]: ";
    std::cin >> duration;
    std::cout << "Sample rate [Hz]: ";
    std::cin >> sampleRate;

    if(!sine){                   //adjusts to create cos
        phase += 90;
    }

    int samples_number = static_cast<int>(sampleRate * duration);

    // Utworzenie wektorów - oś czasu i sygnały
    std::vector<double> time(samples_number);
    std::vector<double> signal(samples_number);

    for (int i = 0; i < samples_number; ++i) {
        time[i] = i / static_cast<double>(sampleRate);
        signal[i] = amplitude * std::sin(2 * M_PI * frequency * time[i] + phase*M_PI/180.0); // aplitude x sin(2pi x f x time(i) + phase(changed to radians))
    }

    using namespace matplot;
    plot(time, signal);
    title(sine ? "Wygenerowany sygnał sinusoidalny" : "Wygenerowany sygnał cosinusoidalny");
    xlabel("Czas (s)");
    ylabel("Amplituda");
    show();
}

void signal_generate_square_wave() {
    double amplitude, period, duration, duty_cycle;
    int sampleRate;

    std::cout << "Amplitude [unit]: ";
    std::cin >> amplitude;
    std::cout << "Period [s]: ";
    std::cin >> period;
    std::cout << "Duty cycle [%]: ";
    std::cin >> duty_cycle;
    std::cout << "Duration [s]: ";
    std::cin >> duration;
    std::cout << "Sample rate [Hz]: ";
    std::cin >> sampleRate;

    duty_cycle /= 100;      //changes percentages to number

    int samples_number = static_cast<int>(sampleRate * duration);

    // Utworzenie wektorów - oś czasu i sygnały
    std::vector<double> time(samples_number);
    std::vector<double> signal(samples_number);

    for (int i = 0; i < samples_number; ++i) {
        time[i] = i / static_cast<double>(sampleRate);
        if (fmod(time[i], period) < duty_cycle * period) {
            signal[i] = amplitude;
        } else {
            signal[i] = 0.0;
        }
    }

    using namespace matplot;
    plot(time, signal);
    title("Wygenerowany sygnał prostokątny");
    xlabel("Czas (s)");
    ylabel("Amplituda");
    ylim({-0.2, amplitude + 0.2});
    show();
}

void saw_wave_generate(){
    double amplitude, frequency, duration, min_value, phase, sample_rate, fourier_approx;
    std::string unit;
    int harmonics = 30; //info from https://kconrad.math.uconn.edu/math1132s10/sawtooth.html#:~:text=Fairly%20general%2C%20even%20discontinuous%2C%20periodic,3x)%20%2B%20....&text=sin(x)%20-%201⁄,(6x)%20%2B%20...

    std::cout << "Unit (on y axis): " << std::endl;
    std::cin >> unit;
    std::cout << "Amplitude [unit]:" << std::endl;
    std::cin >> amplitude;
    std::cout << "Frequency [Hz]:" << std::endl;
    std::cin >> frequency;
     std::cout << "phase (degrees):" << std::endl;
    std::cin >> phase;
    std::cout << "Signal duration [s]: " << std::endl;
    std::cin >> duration;
    std::cout << "Minimal Value [unit]: " << std::endl;
    std::cin >> min_value;
    std::cout << "Sample Rate [Hz]: " << std::endl;
    std::cin >> sample_rate;

        phase *= M_PI/180;      //radians easier to count
        int all_measure_points = static_cast<int>(sample_rate*duration);

        std::vector<double> time(all_measure_points);           //time and signal vectors
        std::vector<double> signal(all_measure_points);     

        for(int t = 0; t < all_measure_points; t++ ){
            time[t] = t/sample_rate;
            fourier_approx =0;
            for(int i =1; i <= harmonics; i++){
            fourier_approx +=  pow(-1, i)*sin(2*M_PI*frequency*time[t]*i + phase)/i;
            }
            signal[t]= amplitude*(0.5 - (1/M_PI)*fourier_approx) + min_value; //from wikipedia + min_value
            if(signal[t] > min_value + amplitude){      //better accuracy
                signal[t] = min_value+amplitude;
            }
            if(signal[t] < min_value){
                signal[t] = min_value;
            }
        }


 using namespace matplot;
    plot(time, signal);
    title("Wygenerowany sygnał piłokształtny");
    xlabel("Czas (s)");
    ylabel("Amplituda [" + unit + "]");
    ylim({min_value - 1, amplitude + min_value + 1});
    show();

}

int twod_filter(){
    std::string file_path; //image location
    std::cout<< "What's the image's location?: " << std::endl;
    std::cin  >> file_path;
    
    file_path = "C:/Users/marty/Documents/projects/tp_projekt3/mayo.jpg";   //test

    cv::Mat image = cv::imread(file_path, cv::IMREAD_COLOR);   //loads img

      if (image.empty()) {                              //checks if loaded
        std::cerr << "Image not found" << std::endl;
        return -1;
    }
    int kernel_size=0;
    while(kernel_size % 2 != 1){
        std::cout<< "Choose the size of the kernel(an odd number): " << std::endl; //odd bc i need the middle 
        std::cin >> kernel_size;
    }

    std::vector<std::vector<int>> kernel(kernel_size, std::vector<int>(kernel_size)); //kernel vector (2D table) outer->row, inner->column

     for (int r = 0; r < kernel_size; ++r) {            //r-row
        for (int c = 0; c < kernel_size; ++c) {         //c-column
            kernel[r][c] = 1;
        }
    }

    int add_frame_size = kernel_size / 2;       //https://docs.opencv.org/4.8.0/dc/da3/tutorial_copyMakeBorder.html

   cv::Mat framed_image;

    cv::copyMakeBorder (image, framed_image, add_frame_size, add_frame_size, add_frame_size, add_frame_size, cv::BORDER_REPLICATE);
    
    // cv::imwrite("outpit_img.jpg", framed_image);  //saves extended img
    cv::Mat final_image = image.clone();

    //filtration aka multiplication of every element of the matrixes
    for(int img_rows = add_frame_size; img_rows < image.rows + add_frame_size; ++img_rows)              //https://www.youtube.com/watch?v=S4z-C-96xfU
        for(int img_cols = add_frame_size; img_cols < image.cols + add_frame_size; ++img_cols){
            int value[3] = {0,0,0};
            for(int kernel_y = 0; kernel_y < kernel_size; kernel_y++){
                for(int kernel_x = 0; kernel_x < kernel_size; kernel_x++){
                    cv::Vec3b pixel = framed_image.at<cv::Vec3b>(img_rows + kernel_y,img_rows + kernel_x);
                    value[0] += pixel[0] * kernel[kernel_y][kernel_x];
                    value[1] += pixel[1] * kernel[kernel_y][kernel_x];
                    value[2] += pixel[2] * kernel[kernel_y][kernel_x];
                }
            }
            final_image[img_rows - add_frame_size][img_cols - add_frame_size] = value;  

        }
    }


namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           multiply
           signal_visualization
           signal_generate_sinusoidal
           signal_generate_square_wave
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers
    )pbdoc");

    m.def("multiply", &multiply, R"pbdoc(
        Multiply two numbers
    )pbdoc");

    m.def("signal_visualization", &signal_visualization, R"pbdoc(
        Visualizes the signal with matplot library
    )pbdoc");

    m.def("signal_generate_sinusoidal", &signal_generate_sinusoidal, R"pbdoc(
        Generates sinusoidal (or cos) signal with matplot library
    )pbdoc");

    m.def("signal_generate_square_wave", &signal_generate_square_wave, R"pbdoc(
        Generates square signal with matplot library
    )pbdoc");

     m.def("saw_wave_generate", &saw_wave_generate, R"pbdoc(
        Generates sawtooth signal with matplot library
    )pbdoc");

     m.def("twod_filter", &twod_filter, R"pbdoc(
        1D filter function with opencv, puts an image through a filter, outputs the filtered image
    )pbdoc");



#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}

