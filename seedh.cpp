#include "AllSequences.h"

#include <boost/format.hpp>
#include <boost/program_options.hpp>

#define DEFAULT_HISTOGRAM_STEP "0.05"
#define DEFAULT_HISTOGRAM_LENGTH "2.00"

int main(int argc, char **argv) {
    string seedString;
    string fastqFn;
    string fastaFn;

    fastaDataT fasta;

    int directionF = 1;
    int directionB = 1;

    /* default parameters*/
    int overlappingAllowed = 0;
    float histogramStep = atof(DEFAULT_HISTOGRAM_STEP);
    float histogramLength = atof(DEFAULT_HISTOGRAM_LENGTH);

    /*
     * PARSE ARGUMENTS
     */
    try {
        namespace po = boost::program_options;

        po::options_description all("OPTIONS");

        all.add_options()
                //
                ("fasta,1", po::value<string>(&fastaFn)->required(), "FASTA reference")
                //
                ("reads,r", po::value<string>(&fastqFn)->required(), "reads")
                //
                ("seed,s", po::value<string>(&seedString)->required(), "seed")
                //
                //("allow-overlaps,o", "allow overlaps of seeds [disabled]")
                //
                ("forward,f", "search only in forward direction")
                //
                ("backward,b", "search only in backward direction")
                //
                //("histogram-step,d", po::value<float>(&histogramStep), "one step in histogram [" DEFAULT_HISTOGRAM_STEP "]")
                //
                //("histogram-length,l", po::value<float>(&histogramLength), "length of histogram [" DEFAULT_HISTOGRAM_LENGTH "]")
                ;

        po::variables_map vm;
        try {
            po::store(po::parse_command_line(argc, argv, all), vm); // can throw

            /*if (vm.count("allow-overlaps")) {
                overlappingAllowed = 1;
            }*/

            if (vm.count("forward")) {
                directionB = 0;
            }

            if (vm.count("backward")) {
                directionF = 0;
            }

            if (vm.count("forward") && vm.count("backward")) {
                directionB = 1;
                directionF = 1;
            }


            po::notify(vm);
        } catch (po::error &e) {
            cerr << "ERROR: " << e.what() << endl << endl;
            cerr << "Example: ./seed_meta -f bact.fa -r reads_fa1.bfast.fastq -s \"##-####-#-#-##-####-#-###\"" << endl;
            cerr << endl << all << endl;
            return EXIT_FAILURE;
        }

    } catch (exception &e) {
        cerr << "Unhandled Exception reached the top of main: " << e.what()
                << ", application will now exit" << std::endl;
        return EXIT_FAILURE;
    }

    Seed seed(seedString);

    /* 1. create indexes */

    cerr << endl << "CREATING INDEX" << endl;

    cerr << "Loading " << fastaFn << endl;
    loadFasta(fastaFn, fasta);
    cerr << "Creating index " << fastaFn << endl;
    AllSequences index(fasta, seed);

    /* 2. open fastq and iterate through all reads */

    cerr << endl << "PROCESSING FASTQ" << endl;

    try {
        cerr << "Processing " << fastqFn << endl;
        ifstream fastq(fastqFn);
        string line;
        string lineFiltered;

        string read_name="";

        int lineNb = 0;

        while (getline(fastq, line)) {
            switch (lineNb % 4) {
                case 1: {
                    /* line with actual sequence of the read */

                    int sharedOverlappingF=0;
                    int sharedNonoverlappingF=0;
                    int coveredPositionsF=0;

                    int sharedOverlappingB=0;
                    int sharedNonoverlappingB=0;
                    int coveredPositionsB=0;

                    /* a) cleaning */
                    cleanGenomicSequence(line);

                    /* b) seeding*/
                    string kmer(seed.weight, ' ');

                    {
                        /* b) i) forward */
                        string covered(line.size(), 0);
                        int lastNonOverlapping = -4242;

                        if (directionF) {
                            for (int i = 0; i <= (int) line.size() - (int) seed.seed.size(); i++) {
                                seed.apply(line, kmer, i);

                                if (index.isKmerPresent(kmer)) {
                                    sharedOverlappingF++;

                                    for (int j=0;j<seed.seed.size();j++){
                                        if (seed.seed[j]=='#'){
                                            covered[i+j]=1;
                                        }
                                    }

                                    if ((i - lastNonOverlapping) >= (int) seed.seed.size()){
                                        sharedNonoverlappingF++;
                                        lastNonOverlapping=i;
                                    }
                                }
                            }

                            for(int j=0;j<covered.size();j++){
                                if(covered[j]==1){
                                    coveredPositionsF++;
                                }
                            }
                        }
                    }

                    reverseComplement(line);

                    {
                        /* b) ii) backward */
                        string covered(line.size(), 0);
                        int lastNonOverlapping = -4242;

                        if (directionB) {
                            for (int i = 0; i <= (int) line.size() - (int) seed.seed.size(); i++) {
                                seed.apply(line, kmer, i);

                                if (index.isKmerPresent(kmer)) {
                                    sharedOverlappingB++;

                                    for (int j=0;j<seed.seed.size();j++){
                                        if (seed.seed[j]=='#'){
                                            covered[i+j]=1;
                                        }
                                    }

                                    if ((i - lastNonOverlapping) >= (int) seed.seed.size()){
                                        sharedNonoverlappingB++;
                                        lastNonOverlapping=i;
                                    }
                                }
                            }

                            for(int j=0;j<covered.size();j++){
                                if(covered[j]==1){
                                    coveredPositionsB++;
                                }
                            }
                        }
                    }

                    cout << read_name
                            << "\t" << max(sharedNonoverlappingF,sharedNonoverlappingB)
                            << "\t" << max(sharedOverlappingF,sharedOverlappingB)
                            << "\t" << max(coveredPositionsF,coveredPositionsB)
                            << endl;

                    break;
                }
                case 0:
                    /* read name */
                    assert(line[0] = '@');
                    read_name=line.substr(1,line.size()-1);

                    break;
                case 2:
                    /* qualities comment */
                    break;
                case 3:
                    /* base qualities */
                    break;
                default:
                    break;
            } // switch (lineNb % 4)

            lineNb++;
        } // while (getline(fastq, line))
        fastq.close();

    } catch (std::ifstream::failure &e) {
        std::cerr << "Exception opening/reading/closing file" << fastqFn
                << endl;
        exit(1);
    }

    return 0;
}
