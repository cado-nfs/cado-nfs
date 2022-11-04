#include "cado.h"
#include <istream>
#include <string>
#include "params.h"

int param_list_read(param_list_ptr pl, std::istream & is, bool stop_on_empty_line)
{
    int all_ok=1;
    for(std::string line ; std::getline(is, line, '\n') ; ) {
        if (line[0] == '#')
            continue;

        // remove possible comment at end of line.
        std::string::size_type p = line.find('#');
        if (p != std::string::npos)
            line.erase(p);

        for( ; !line.empty() && isspace(line[line.size()-1]) ; )
            line.erase(line.size()-1);

        for(p = 0 ; !line.empty() && isspace(line[p]) ; p++);
        if (p) line.erase(0, p);

        // empty ps are ignored (unless stop_on_empty_line is non-zero).
        if (line.empty()) {
            if (stop_on_empty_line)
                break;
            else
                continue;
        }

        // look for a left-hand-side. We grok anything that *BEGINS WITH
        // A DIGIT* as something that goes with the "NULL" token in the
        // pl dictionary. That looks like a pretty obscure hack, in fact.
        // Do we ever use it ?
        if (!(isalpha((int)(unsigned char)line[0]) || line[0] == '_' || line[0] == '-')) {
            param_list_add_key(pl, NULL, line.c_str(), PARAMETER_FROM_FILE);
            continue;
        }
        std::string::const_iterator q = line.begin();
        for( ; q != line.end() && (isalnum((int)(unsigned char)*q) || *q == '_' || *q == '-') ; ++q);
        if (q == line.begin()) {
            fprintf(stderr, "Parse error, no usable key for config line:\n%s\n",
                    line.c_str());
            all_ok=0;
            continue;
        }

        std::string key(line.cbegin(), q);

        /* Now we can match (whitespace+ | whitespace* separator whitespace*) data
         */
        for( ; q != line.end() && isspace((int)(unsigned char)*q) ; ++q);

        /* should we actually allow it, after all ? */
        if (q == line.end()) {
            fprintf(stderr, "Parse error, key with no value in config:\n%s\n",
                    line.c_str());
            continue;
        }

        /* match separator, which is one of : = := */
        if (*q == '=') {
            q++;
        } else if (*q == ':') {
            q++;
            if (*q == '=')
                q++;
        } else if (q == line.begin() + key.size()) {
            fprintf(stderr, "Parse error, no separator for config line:\n%s\n",
                    line.c_str());
            all_ok=0;
            continue;
        }
        for( ; *q && isspace((int)(unsigned char)*q) ; q++);

        std::string value(q, line.cend());

        param_list_add_key(pl, key.c_str(), value.c_str(), PARAMETER_FROM_FILE);
    }

    return all_ok;
}

