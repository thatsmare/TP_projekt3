#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <set>
#include <cmath>
#include <vector>
#include "AudioFile.h"


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

int multiply(int i, int j){
    AudioFile<double> audioFile;

    audioFile.load ("C:/Users/marty/Documents/GitHub/TP_projekt3/AudioFile/examples/test-audio.wav");

    int sampleRate = audioFile.getSampleRate();
    std::cout << sampleRate;

    return 0;
}

void file_location_in(){

}

void signal_visualization(){
    AudioFile<double> audiosample;

    std::string file_location = C:/Users/marty/Documents/GitHub/TP_projekt3/AudioFile/examples/test-audio.wav;

    audiosample.load(file_location);

    //częstotliwość próbkowania
    int sampleRate = audiosample.getSampleRate();
    std::vector<std::vector<double>> sample = audiosample.samples;

    //utworzenie osi czasu dla próbek
    size_t sample_number = sample[0].size();
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

    return 0;
}

void signal_generate_sinusoidal()
{
    double amplitude, frequency, phase, duration;  
    int sampleRate;  

    std::cout<<"Amplitude: ";
    std::cin>>amplitude;
    std::cout<<endl<<"Frequency (Hz): ";
    std::cin>>fraquency;
    std::cout<<endl<<"Phase (in degrees): ";
    std::cin>>phase;
    std::cout<<endl<<"Duration (s): ";
    std::cin>>duration;
    std::cout<<endl<<"Sample rate (Hz): ";
    std::cin>>sampleRate;

    int samples_number = static_cast<int>(sampleRate * duration);
    
    //utworzenie wektorów - oś czasu i sygnały
    std::vector<double> time(samples_number);
    std::vector<double> signal(samples_number);

    for (int i = 0; i < samples_number; ++i) {
        time[i] = i / static_cast<double>(sampleRate);
        signal[i] = amplitude * std::sin(2 * M_PI * frequency * time[i] + phase);
    }

    using namespace matplot;
    plot(time, signal);
    title("Wygenerowany sygnał sinusoidalny");
    xlabel("Czas (s)");
    ylabel("Amplituda");
    show();

    return 0;
}

void signal_generate_square_wave()
{
    double amplitude, period, duration, duty_cycle;
    int sampleRate;

    std::cout<<"Amplitude: ";
    std::cin>>amplitude;
    std::cout<<endl<<"Period (s): ";
    std::cin>>fraquency;
    std::cout<<endl<<"Duty cycle (%): ";
    std::cin>>phase;
    std::cout<<endl<<"Duration (s): ";
    std::cin>>duration;
    std::cout<<endl<<"Sample rate (Hz): ";
    std::cin>>sampleRate;

    int samples_number = static_cast<int>(sampleRate * duration);
    
    //utworzenie wektorów - oś czasu i sygnały
    std::vector<double> time(samples_number);
    std::vector<double> signal(samples_number);

    for (int i = 0; i < samples_number; ++i) {
        time[i] = i / static_cast<double>(sampleRate);
        if (fmod(time[i], period) < dutyCycle * period) {
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
    ylim({-0.2, amplitude+0.2});
    show();

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
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("multiply", &multiply, R"pbdoc(
        Multiply two numbers
    )pbdoc");

     m.def("signal_visualization", &signal_visualization, R"pbdoc(
        Visualizes the signal with matplot library
    )pbdoc");

    m.def("signal_generate_sinusoidal", &signal_generate_sinusoidal,R"pbdoc(
        Generates sinusoidal signal with matplot library
    )pbdoc");

    m.def("signal_generate_square_wave", &signal_generate_square_wave,R"pbdoc(
        Generates square signal with matplot library
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
