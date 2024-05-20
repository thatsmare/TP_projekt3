#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <set>
#include <math.h>
#include <vector>
#include <iostream>
#include <string> 
#include "AudioFile.h"
#include <opencv2/opencv.hpp>
#include <algorithm>

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

    std::string file_location; // = "C:/Users/marty/Documents/GitHub/TP_projekt3/AudioFile/examples/test-audio.wav";

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

void signal_generate_sinusoidal(double amplitude, double frequency, double phase, double duration, int sampleRate) { //frequency [Hz], phase [degrees], duration [s], sampleRate [Hz]

    if (amplitude > 0 && frequency > 0 && duration > 0 && 1/sampleRate < duration) {
        
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
    title("Wygenerowany sygnał sinusoidalny");
    xlabel("Czas (s)");
    ylabel("Amplituda");
    ylim({-0.2 - amplitude, amplitude + 0.2});   
    show();
    }
    else std::cout<<"Function parameters are incorrect";
}

void signal_generate_square_wave(double amplitude, double period, double duration, double duty_cycle, int sampleRate) { //period [s], duration [s], duty_cycle [%], sampleRate [Hz]
    if (amplitude > 0 && period > 0 && duration > 0 && duty_cycle > 0 && duty_cycle < 100 && 1/sampleRate < duration) {
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
    else std::cout<<"Function parameters are incorrect";
}

void saw_wave_generate(){
    double amplitude, frequency, duration, min_value, phase, sample_rate, fourier_approx;
    std::string unit;
    int harmonics = 101; //info from https://kconrad.math.uconn.edu/math1132s10/sawtooth.html#:~:text=Fairly%20general%2C%20even%20discontinuous%2C%20periodic,3x)%20%2B%20....&text=sin(x)%20-%201⁄,(6x)%20%2B%20...

    std::cout << "Unit (on y axis): " << std::endl;
    std::cin >> unit;
    std::cout << "Amplitude [" << unit <<"]:" << std::endl;
    std::cin >> amplitude;
    std::cout << "Frequency [Hz]:" << std::endl;
    std::cin >> frequency;
     std::cout << "phase (degrees):" << std::endl;
    std::cin >> phase;
    std::cout << "Signal duration [s]: " << std::endl;
    std::cin >> duration;
    std::cout << "Minimal Value [" << unit <<"]: " << std::endl;
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
            fourier_approx +=  pow(-1, i + 1)*sin((2 *M_PI * frequency* time[t] + phase)*i)/i;
            }
            signal[t]= amplitude*(0.5 - (1/M_PI)*fourier_approx) + min_value; //from wikipedia + min_value https://en.wikipedia.org/wiki/Sawtooth_wave
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
  /*  std::string file_path; //image location
    std::cout<< "What's the image's location?: " << std::endl;
    std::cin  >> file_path;
    
    //file_path = "C:/Users/marty/Documents/projects/tp_projekt3/mayo.jpg";   //test

    cv::Mat image = cv::imread(file_path, cv::IMREAD_COLOR);   //loads img

      if (image.empty()) {                              //checks if loaded
        std::cerr << "Image not found" << std::endl;
        return -1;
    }
    */
    
  /*  int kernel_size=0;
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
                    cv::Vec3b pixel = framed_image.at<cv::Vec3b>(img_rows + kernel_y - add_frame_size, img_cols + kernel_x - add_frame_size);
                    value[0] += pixel[0] * kernel[kernel_y][kernel_x];              //different colors :)
                    value[1] += pixel[1] * kernel[kernel_y][kernel_x];
                    value[2] += pixel[2] * kernel[kernel_y][kernel_x];
                }
            }
           final_image.at<cv::Vec3b>(img_rows - add_frame_size, img_cols - add_frame_size) = cv::Vec3b(         //colors are from (0-255)
                std::min(std::max(value[0], 0), 255),
                std::min(std::max(value[1], 0), 255),               //bc colors are from 0 to 255, have to cap them
                std::min(std::max(value[2], 0), 255)
            );
        }

         if (!cv::imwrite("output_img.jpg", final_image)) {
        std::cerr << "Error saving the image" << std::endl;
        return -1;
        }

        cv::imshow("Filtered Image", final_image);
        cv::waitKey(0);
    */    
   return 0;
}

int bilinear_interpolation(){
 /*   std::string file_path; //image location
    std::cout<< "What's the image's location?: " << std::endl;
    std::cin  >> file_path;

      cv::Mat image = cv::imread(file_path, cv::IMREAD_COLOR);   //loads img

      if (image.empty()) {                              //checks if loaded
        std::cerr << "Image not found" << std::endl;
        return -1;
    }

    int old_width, old_height;  //old dimensions of the image
    int new_width, new_height;
    */
    return 0;
}

int oned_filtering(){
  AudioFile<double> audiosample;

   std::string file_location; // = "C:/Users/marty/Documents/GitHub/TP_projekt3/AudioFile/examples/test-audio.wav";

    std::cout << "Location of the file:" << std::endl;     //choose path 
    std::cin >> file_location;

    if(!audiosample.load(file_location)){
        std::cerr << "Audio could not be loaded" << std::endl;
        return 1;
    }

    if(audiosample.getNumChannels() !=1){
        std::cout << "1D means only one channel, the provided file contains more" << std::endl;
        return 1;
    }

    int numSamples = audiosample.getNumSamplesPerChannel();

    std::vector<double> audio(numSamples); 

    for (int i = 0; i < numSamples; i++)
    {
	    audio[i] = audiosample.samples[0][i]; //0 is for channels we have 1
    }
    
    int kernel_size;

    std::cout << "Kernel size: " << std::endl;
    while(!(std::cin >> kernel_size)){
        std::cout << "Input an integer, please" << std::endl;
    }

    std::vector<double> kernel(kernel_size,1);
    std::vector<double> output(numSamples);

    for(int i = 0; i< audio.size(); i++){               //https://www.youtube.com/watch?v=yd_j_zdLDWs
        int output_cell = 0;
        for(int j =0; j < kernel_size; j++){
            if(!((i+j) >audio.size())){
            output_cell += audio[i+j] * kernel[j];
            }
        }
        output[i] = output_cell;
    }


    AudioFile<double> post_filtering;                   //saving the audio
    post_filtering.setNumChannels(1);
    post_filtering.setSampleRate(audiosample.getSampleRate());          //everything about an audio file https://github.com/adamstark/AudioFile
    post_filtering.setBitDepth(audiosample.getBitDepth());
    post_filtering.setNumSamplesPerChannel(numSamples);

    for(int i = 0; i< numSamples; i++){
        post_filtering.samples[0][i] = output[i];
    }
    
    std::string post_filtering_loc = "post_filtering.wav";

     if (!post_filtering.save(post_filtering_loc)) {
        std::cerr << "Error saving the output audio file" << std::endl;
        return 1;
    }
    return 0;
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
        2d filtration
    )pbdoc");

    m.def("bilinear_interpolation", &bilinear_interpolation, R"pbdoc(
        Bilinear interpolation of an image
        )pbdoc");

    m.def("oned_filtering", &oned_filtering, R"pbdoc(
        1D filtering of a signal
        )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}

