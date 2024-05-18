#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <set>
#include <cmath>
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

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
