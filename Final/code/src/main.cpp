#include "ppm.hpp"
#include "pt.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    for (int argNum = 1; argNum < argc; ++argNum) {
        cout << "Argument " << argNum << " is: " << argv[argNum] << endl;
    }

    string input_file, output_file;
    if (argc != 4) {
        cout << "Usage: ./bin/FINAL <input scene file> <output bmp file> <photon_num * 4>" << endl;
        return 1;
    } else {
        input_file = argv[1];
        output_file = argv[2];
    }
    int samples = samples = argc == 4 ? atoi(argv[3]) / 4 : 10000; // # samples
    SceneParser parser(input_file.c_str());
    parser.getGroup()->print();
#ifdef __AVX2__
    cout << "Using Intel AVX to accelerate BSpline curve..." << endl;
#endif
    // Image img = pathTrace(parser, samples);
    Image img = photonMap(parser, samples);
    img.SaveBMP(output_file.c_str());
    return 0;
}