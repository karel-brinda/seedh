#ifndef __UTILS_H__
#define __UTILS_H__

#include <algorithm>
#include <vector>
#include <list>
#include <utility>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <boost/algorithm/string/join.hpp>


using namespace std;

/* hashed nucleotides */
#define HASH_SIZE 8
/* size of the hash table */
//#define HASH_TABLE_SIZE ((int)(pow(4.0,HASH_SIZE)))
// TODO: solve it in more elegant way
#define HASH_TABLE_SIZE 4*4*4*4*4*4*4*4

#define TEST_SMALL_NUCL(__X__) ((__X__=='a') || (__X__=='c') || (__X__=='g') || (__X__=='t'))
#define TEST_CAP_NUCL(__X__) ((__X__=='A') || (__X__=='C') || (__X__=='G') || (__X__=='T'))

#define REPR_A 0
#define REPR_C 1
#define REPR_G 2
#define REPR_T 3


typedef int posT;
typedef pair<string, string> fastaPairT;
typedef vector<fastaPairT> fastaDataT;

inline void cleanGenomicSequence(string &seq) {
    int j = 0;
    for (int i = 0; i < (int) seq.size(); i++) {
        switch (seq[i]) {
            case 'A':
            case 'a':
                seq[j++] = 'a';
                break;
            case 'C':
            case 'c':
                seq[j++] = 'c';
                break;
            case 'G':
            case 'g':
                seq[j++] = 'g';
                break;
            case 'T':
            case 't':
                seq[j++] = 't';
                break;
            default:
                break;
        }
    }
    seq.resize(j);
}

inline void reverseComplement(string &seq) {
    int j = 0;
    for (int i = 0; i < (int) seq.size(); i++) {
        switch (seq[i]) {
            case 'A':
            case 'a':
                seq[i] = 't';
                break;
            case 'C':
            case 'c':
                seq[i] = 'g';
                break;
            case 'G':
            case 'g':
                seq[i] = 'c';
                break;
            case 'T':
            case 't':
                seq[i] = 'a';
                break;
            default:
                cerr << "unknown error" << endl;
                cerr << "file " << __FILE__ << ", line " << __LINE__ << endl;
                exit(1);
        }
    }
    reverse(seq.begin(), seq.end());
}

inline void loadFasta(const string &fastaFn, fastaDataT &fastaData) {
    fastaData.clear();
    try {
        ifstream fasta(fastaFn);
        string name("");
        list<string> seq;

        string line;
        fastaPairT fp;
        int charCount = 0;
        while (getline(fasta, line)) {
            if (line[0] == '>') {
                /* new seq*/

                /* save current seq */
                if (name != "") {
                    fp.first = name;
                    string *sec = new string();
                    fp.second = boost::algorithm::join(seq, "");
                    fastaData.push_back(fp);
                }

                seq.clear();

                /* extract seq name */
                int start = 1;
                int end;
                for (end = line.length() - 1;
                     line[end] == ' ' || line[end] == '\n' || line[end] == '\r' || line[end] == '\t'; end--) {
                }
                name = line.substr(start, end - start + 1);
            } else {
                cleanGenomicSequence(line);
                seq.push_back(string(line));
            }
        }
        fasta.close();

        fp.first = name;
        fp.second = boost::algorithm::join(seq, "");
        fastaData.push_back(fp);
    }

    catch (std::ifstream::failure &e) {
        cerr << "Exception opening/reading/closing file" << fastaFn << endl;
        cerr << "file " << __FILE__ << ", line " << __LINE__ << endl;
        exit(1);
    }

}

inline unsigned char nuclToId(const char &nucl) {
    //cerr << nucl;
    switch (nucl) {
        case 'a':
        case 'A':
            return REPR_A;
        case 'c':
        case 'C':
            return REPR_C;
        case 'g':
        case 'G':
            return REPR_G;
        case 't':
        case 'T':
            return REPR_T;
        default:
            cerr << "unknown nucleotide: '" << nucl << "' (number " << int(nucl) << ")" << endl;
            cerr << "file " << __FILE__ << ", line " << __LINE__ << endl;
            exit(1);
    }
}

inline char idToNucl(const unsigned char &nucl) {
    switch (nucl) {
        case REPR_A:
            return 'a';
        case REPR_C:
            return 'c';
        case REPR_G:
            return 'g';
        case REPR_T:
            return 't';
        default:
            cerr << "unknown error" << endl;
            cerr << "file " << __FILE__ << ", line " << __LINE__ << endl;
            exit(1);
    }
}

inline string getDehash(int hash) {
    assert(hash >= 0);
    assert(hash < HASH_TABLE_SIZE);
    string toRet(HASH_SIZE, nuclToId('a'));
    for (int i = 0; i < HASH_SIZE; i++) {
        toRet[HASH_SIZE - i - 1] = idToNucl(hash % 4);
        hash /= 4;
    }
    assert(hash == 0);
    return toRet;
}

inline unsigned short int getHash(const string &kmer) {
    int hash = 0;
    for (int i = 0; i < 8; i++) {
        hash *= 4;
        hash += nuclToId(kmer[i]);
    }
    return hash;
}

#endif
