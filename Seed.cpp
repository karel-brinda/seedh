#include "Seed.h"

Seed::Seed(string seed) :
        weight(0), seed(seed) {
    for (int i = 0; i < (int) seed.size(); i++) {
        assert(seed[i] == '#' || seed[i] == '-');
        if (seed[i] == '#') {
            weight++;
        }
    }
}

Seed::Seed(const Seed &seed) :
        weight(seed.weight), seed(seed.seed) {
}

void Seed::apply(const string &text, string &kmer, const int &applPos) const {
    assert((int) kmer.size() == weight);
    assert((0 <= applPos && applPos <= (int) text.size() - (int) seed.size()) || !(cerr << "applPos=" << applPos << endl));

    int wi = 0;
    for (int i = 0; wi < weight; i++) {
        assert(i < (int) seed.size());
        if (seed[i] == '#') {
            kmer[wi] = text[i + applPos];
            wi++;
        }
    }
}
