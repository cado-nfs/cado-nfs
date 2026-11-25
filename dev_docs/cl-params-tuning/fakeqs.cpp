// Synthetic SIQS relations for filtering simulations.
// Reads a sample of real siqs relations ("a,b:p,p,...:" hex primes), and
// writes R fake relations in the post-filter_galois format understood by
// purge ("a,b:[-]idx,..." hex, even indices).  Each fake relation is a
// sample relation whose primes >= T are replaced by a random prime ideal of
// about the same size (p' = p*exp(U(-eps,eps))).  Column numbering:
//   p <  T : exact rank of p among integers (2*p, sparse but harmless)
//   p >= T : T + (li(p)-li(T))/(2*shrink)   (half of the primes split)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <cstring>

static double li(double x) {
    double l = log(x);
    return x / l * (1 + 1 / l + 2 / (l * l) + 6 / (l * l * l) + 24 / (l * l * l * l));
}

int main(int argc, char **argv) {
    if (argc < 9) {
        fprintf(stderr, "usage: %s sample R lpb T eps shrink seed outprefix [nfiles]\n", argv[0]);
        return 1;
    }
    const char *sample = argv[1];
    uint64_t R = strtoull(argv[2], 0, 10);
    int lpb = atoi(argv[3]);
    double T = atof(argv[4]);
    double eps = atof(argv[5]);
    double shrink = atof(argv[6]);
    uint64_t seed = strtoull(argv[7], 0, 10);
    std::string prefix = argv[8];
    int nfiles = argc > 9 ? atoi(argv[9]) : 1;
    double L = ldexp(1.0, lpb);
    double liT = li(T);

    std::vector<std::vector<uint64_t>> rels;
    FILE *f = fopen(sample, "r");
    if (!f) { perror(sample); return 1; }
    char *line = NULL; size_t cap = 0;
    while (getline(&line, &cap, f) > 0) {
        if (line[0] == '#') continue;
        char *c = strchr(line, ':');
        if (!c) continue;
        std::vector<uint64_t> v;
        char *p = c + 1;
        while (*p && *p != ':' && *p != '\n') {
            char *e;
            uint64_t x = strtoull(p, &e, 16);
            if (e == p) break;
            v.push_back(x);
            p = e;
            if (*p == ',') p++;
        }
        if (!v.empty()) rels.push_back(v);
    }
    fclose(f);
    fprintf(stderr, "read %zu sample relations\n", rels.size());
    if (rels.empty()) return 1;

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(-eps, eps);
    std::vector<FILE *> out(nfiles);
    for (int i = 0; i < nfiles; i++) {
        std::string fn = prefix + "." + std::to_string(i);
        out[i] = fopen(fn.c_str(), "w");
        setvbuf(out[i], NULL, _IOFBF, 1 << 20);
    }
    std::vector<int64_t> idx;
    for (uint64_t r = 0; r < R; r++) {
        auto const &v = rels[rng() % rels.size()];
        idx.clear();
        for (size_t i = 0; i < v.size(); i++) {
            double p = (double)v[i];
            int64_t c;
            if (p < T) {
                c = (int64_t)v[i];
            } else {
                double q = p * exp(U(rng));
                if (q < T) q = T * T / q;       // reflect
                if (q >= L) q = L * L / q;
                if (q < T) q = T;
                c = (int64_t)T + (int64_t)((li(q) - liT) / (2 * shrink));
            }
            c *= 2;
            // random sign for jittered ideals; keep multiplicity
            if (p >= T && (rng() & 1)) c = -c;
            idx.push_back(c);
        }
        // combine entries on the same column: signed exponents add up
        std::sort(idx.begin(), idx.end(), [](int64_t a, int64_t b) {
            return llabs(a) < llabs(b) || (llabs(a) == llabs(b) && a < b);
        });
        FILE *o = out[r % nfiles];
        fprintf(o, "%lx,1:", (unsigned long)(r + 1));
        bool first = true;
        for (size_t i = 0; i < idx.size();) {
            size_t j = i; int e = 0;
            int64_t c = llabs(idx[i]);
            while (j < idx.size() && llabs(idx[j]) == c) { e += (idx[j] < 0 || (c == 0 && false)) ? -1 : 1; j++; }
            for (int k = 0; k < abs(e); k++) {
                if (!first) fputc(',', o);
                first = false;
                if (e < 0) fprintf(o, "-%lx", (unsigned long)c);
                else fprintf(o, "%lx", (unsigned long)c);
            }
            i = j;
        }
        fputc('\n', o);
    }
    for (auto o : out) fclose(o);
    return 0;
}
