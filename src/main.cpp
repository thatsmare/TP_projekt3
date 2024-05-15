#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <set>
#include <AudioFile.h>
#include <cmath>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

int multiply(int i, int j){
    using namespace matplot;

    std::set<std::vector<double>> Y = {
        {16, 5, 9, 4}, {2, 11, 7, 14}, {3, 10, 6, 15}, {13, 8, 12, 1}};
    plot(Y);

    show();
    return 0;
}

namespace audio
{
    void createAudiofile();
}

namespace audio
{
    void createAudioFile(){
        std::cout << "creating an audio file" << std::endl;

        AudioFile<float> aFile;
        aFile.setNumChannels (2);
        aFile.setNumSamplesPerChannel (44100);

           //---------------------------------------------------------------
        // 2. Create some variables to help us generate a sine wave
        
        const float sampleRate = 44100.f;
        const float frequencyInHz = 440.f;
        
        //---------------------------------------------------------------
        // 3. Write the samples to the AudioFile sample buffer
        
        for (int i = 0; i < aFile.getNumSamplesPerChannel(); i++)
        {
            for (int channel = 0; channel < aFile.getNumChannels(); channel++)
            {
                aFile.samples[channel][i] = sin ((static_cast<float> (i) / sampleRate) * frequencyInHz * 2.f * (float)M_PI);
            }
        }
        
        //---------------------------------------------------------------
        // 4. Save the AudioFile
        
        std::string filePath = "sine-wave.wav"; // change this to somewhere useful for you
        a.save ("sine-wave.wav", AudioFileFormat::Wave);


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
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("multiply", &multiply, R"pbdoc(
        Multiply two numbers
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
