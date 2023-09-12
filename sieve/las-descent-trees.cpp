#include "cado.h" // IWYU pragma: keep

#include <list>
#include <map>
#include <string>
// #include <iostream>
#include <sstream>
#include <cstring>      // strlcpy (openbsd) // IWYU pragma: keep

#include "las-descent-trees.hpp"
#include "portability.h"
#include "macros.h"


using namespace std;

int descent_tree::display_tree(FILE* o, tree * t, string const& prefix) {
    int res = 1;
    char comment[10] = {'\0'};
    if (!t->contender) {
        if (t->try_again) {
            snprintf(comment, sizeof(comment), " try%d", t->try_again);
        } else {
            size_t rc = strlcpy(comment, " ###", sizeof(comment));
            ASSERT_ALWAYS(rc < sizeof(comment));
        }
    }
    fprintf(o, "%s%s [%1.4f]%s\t\t%s\n",
            prefix.c_str(), t->label().c_str(), t->spent,
            comment,
            t->label.fullname().c_str());
    string new_prefix = prefix + "  ";
    typedef list<tree *>::iterator it_t;
    for(it_t i = t->children.begin() ; i != t->children.end() ; i++) {
        if (!display_tree(o, *i, new_prefix))
            res = 0;
    }
    if (!t->contender)
        res = 0;
    return res;
}

void descent_tree::display_last_tree(FILE * o) {
    bool xs = is_successful(forest.back());
    display_tree(o, forest.back(), xs ? "# ": "# FAILED ");
}

struct collected_stats {
    bool ok;
    double t;
    int d;
    int w;
    collected_stats(bool ok, double t, int d, int w) :
        ok(ok), t(t), d(d), w(w) {}
};

void descent_tree::display_all_trees(FILE * o)
{
    typedef list<tree *>::iterator it_t;
    int total = 0, good = 0;
    for(it_t i = forest.begin() ; i != forest.end() ; i++, total++) {
        bool xs = is_successful(*i);
        good += display_tree(o, *i, xs ? "# ": "# FAILED ");
    }
    fprintf(o, "# Success %d/%d (%1.2f%%)\n", good, total,
            100.0 * (double) good / total);
    list<collected_stats> foo;
    typedef map<tree_label, list<collected_stats> > stats_t;
    stats_t stats;
    typedef stats_t::iterator sit_t;
    for(it_t i = forest.begin() ; i != forest.end() ; i++, total++) {
        collected_stats w(is_successful(*i),
                (*i)->spent,
                tree_depth(*i),
                tree_weight(*i)
                );
        sit_t it = stats.find((*i)->label);
        if (it == stats.end()) {
            list<collected_stats> v;
            v.push_back(w);
            stats.insert(make_pair((*i)->label, v));
        } else {
            it->second.push_back(w);
        }
    }
#if 0
    /* (since we've dropped the dependence on siever_config, we no
     * longer do this)
     */
    /* We also use this to provide an updated stats file (whose
     * results are based on actual descents, and thus less
     * speculative) */
    
    list<string> new_hints;
    for(sit_t it = stats.begin() ; it != stats.end() ; it++) {
        /* Collect and print stats for this sq size */
        fprintf(o, "# Stats for %s\n", it->first().c_str());
        int n = it->second.size();
        int nok = 0;
        double t1 = 0;
        double t2 = 0;
        double tmin = 999999;
        double tmax = 0;
        int d1 = 0;
        int dmin = 0;
        int dmax = 0;
        int w1 = 0;
        int wmin = 0;
        int wmax = 0;
        typedef list<collected_stats>::iterator lit_t;
        for(lit_t i = it->second.begin() ; i != it->second.end() ; i++) {
            nok += i->ok;
            t1 += i->t;
            t2 += i->t * i->t;
            d1 += i->d;
            w1 += i->w;
            tmin = min(tmin, i->t);
            tmax = max(tmax, i->t);
            dmin = min(dmin, i->d);
            dmax = max(dmax, i->d);
            wmin = min(wmin, i->w);
            wmax = max(wmax, i->w);
        }
        fprintf(o, "#   success %.2f (%d/%d)\n", (double) nok / n, nok, n);
        fprintf(o, "#   depth avg %.1f, min %d, max %d\n",
                (double) d1/n, dmin, dmax);
        fprintf(o, "#   weight avg %.1f, min %d, max %d\n",
                (double) w1/n, wmin, wmax);
        double tt1 = (double) t1/n;
        double tt2 = (double) t2/n;
        fprintf(o, "#   time avg %.3f, sdev %.3f, min %.3f, max %.3f\n",
                tt1, sqrt(tt2-tt1*tt1), tmin, tmax);

        ostringstream os;
        os << it->first()
            << " " << tt1
            << " " << (double) nok / n 
            << " I=" << it->first.sc->logI;
        for(int i = 0 ; i < 2 ; i++) {
            os << " " << it->first.sc->sides[i]->lim
                << "," << it->first.sc->sides[i]->lpb
                << "," << it->first.sc->sides[i]->mfb
                << "," << it->first.sc->sides[i]->lambda;
        }
        new_hints.push_back(os.str());
    }
    fprintf(o, "# The following data is an _example_ which can be used to provide a hint file\n");
    for(list<string>::iterator i = new_hints.begin() ; i != new_hints.end() ; i++) {
        fprintf(o, "# %s\n", i->c_str());
    }
#endif
}

