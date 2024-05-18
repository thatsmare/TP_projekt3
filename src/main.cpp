#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <set>
#include <cmath>
#include <vector>
#include <iostream>
#include "AudioFile.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#ifndef NOMINMAX
#define NOMINMAX
#endif

// Define M_PI if not defined
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

void file_location_in() {
    // Function implementation here
}

void signal_visualization() {
    AudioFile<double> audiosample;

    std::string file_location = "C:/Users/marty/Documents/GitHub/TP_projekt3/AudioFile/examples/test-audio.wav";

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

    std::cout << "Amplitude: ";
    std::cin >> amplitude;
    std::cout << "Frequency (Hz): ";
    std::cin >> frequency;
    std::cout << "Phase (in degrees): ";
    std::cin >> phase;
    std::cout << "Duration (s): ";
    std::cin >> duration;
    std::cout << "Sample rate (Hz): ";
    std::cin >> sampleRate;

    int samples_number = static_cast<int>(sampleRate * duration);

    // Utworzenie wektorów - oś czasu i sygnały
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
}

void signal_generate_square_wave() {
    double amplitude, period, duration, duty_cycle;
    int sampleRate;

    std::cout << "Amplitude: ";
    std::cin >> amplitude;
    std::cout << "Period (s): ";
    std::cin >> period;
    std::cout << "Duty cycle (%): ";
    std::cin >> duty_cycle;
    std::cout << "Duration (s): ";
    std::cin >> duration;
    std::cout << "Sample rate (Hz): ";
    std::cin >> sampleRate;

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
        Generates sinusoidal signal with matplot library
    )pbdoc");

    m.def("signal_generate_square_wave", &signal_generate_square_wave, R"pbdoc(
        Generates square signal with matplot library
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}

