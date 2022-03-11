#include "cado.h" // IWYU pragma: keep
#include "macros.h"
#include <exception>
#include <sstream>
#include <limits.h>
#include "json.hpp"



struct json_parser {
    struct parse_error: public std::exception {
        const char * what() const noexcept override { return "parse error"; }
    };
    struct tokenizer_error: public std::exception {
        std::string message;
        tokenizer_error(std::istream& is) {
            std::string s;
            std::getline(is, s);
            std::ostringstream os;
            os << "tokenizer error before: " << s;
            message = os.str();
        }
        const char * what() const noexcept override { return message.c_str(); }
    };

private:
    /* https://tools.ietf.org/html/rfc7159 */
    enum expression_token {
        LEFT_BRACE,
        RIGHT_BRACE,
        LEFT_SQUARE_BRACKET,
        RIGHT_SQUARE_BRACKET,
        COLON,
        COMMA,
        /* escape U+0000 to U+001F, U+0022, U+005C, plus shortcuts:
         *
         *
              %x22 /          ; "    quotation mark  U+0022
              %x5C /          ; \    reverse solidus U+005C
              %x2F /          ; /    solidus         U+002F
              %x62 /          ; b    backspace       U+0008
              %x66 /          ; f    form feed       U+000C
              %x6E /          ; n    line feed       U+000A
              %x72 /          ; r    carriage return U+000D
              %x74 /          ; t    tab             U+0009
              %x75 4HEXDIG )  ; uXXXX                U+XXXX
        */
        STRING,
        NUMBER,
        JNULL, BOOL,
    };
    std::vector<expression_token> tokens;
    std::vector<std::string> strings;
    std::vector<double> numbers;
    std::vector<expression_token>::const_iterator ctok;
    std::vector<std::string>::const_iterator cstr;
    std::vector<double>::const_iterator cnum;
    inline bool dry() { return ctok == tokens.end(); }
    inline void next() { if (!dry()) ctok++; }
    inline bool test(expression_token s) { return !dry() && *ctok == s; }
    bool accept(expression_token s) {
        if (!test(s)) return false;
        next();
        return true;
    }
    int expect(expression_token s) {
        if (accept(s))
            return 1;
        throw parse_error();
        return 0;
    }
public:
    json parse_inner() {
        if (accept(LEFT_SQUARE_BRACKET)) {
            json_array * j;
            json jj(j = new json_array());
            bool trailing_comma = false;
            for(;!accept(RIGHT_SQUARE_BRACKET);) {
                j->data.push_back(parse_inner());
                trailing_comma = accept(COMMA);
            }
            if (trailing_comma) throw parse_error();
            return jj;
        } else if (accept(LEFT_BRACE)) {
            json_hash * j;
            json jj(j = new json_hash());
            bool trailing_comma = false;
            for(;!accept(RIGHT_BRACE);) {
                expect(STRING);
                std::string const& key(*cstr++);
                expect(COLON);
                j->data.emplace(key, parse_inner());
                trailing_comma = accept(COMMA);
            }
            if (trailing_comma) throw parse_error();
            return jj;
        } else if (accept(JNULL)) {
            return json(new json_null());
        } else if (accept(NUMBER)) {
            return json(new json_number(*cnum++));
        } else if (accept(BOOL)) {
            return json(new json_bool(*cnum++));
        } else if (accept(STRING)) {
            return json(new json_string(*cstr++));
        } else {
            throw parse_error();
        }
    }

    json parse() {
        json j = parse_inner();
        if (!dry())
            throw parse_error();
        return j;
    }


    bool tokenize(std::istream& is) {
        tokens.clear();
        strings.clear();
        numbers.clear();
        for( ; !is.eof() ; ) {
            int c;
            for(;;is.get()) {
                c = is.peek();
                if (is.eof()) break;
                if (!isspace(c)) break;
            }
            /* c is the next non-whitespace character */
            if (is.eof()) break;
            if (c == '"') {
                is.get(), c=is.peek();
                std::string s;
                for(;!is.eof();is.get(), c=is.peek()) {
                    if (c == '"') {
                        is.get(), c=is.peek();
                        break;
                    }
                    if (c == '\\') {
                        is.get(), c=is.peek();
                        if (c == '\\') { s += '\\'; }
                        else if (c == '/') { s += '/'; }
                        else if (c == 'n') { s += '\n'; }
                        else if (c == 't') { s += '\t'; }
                        else if (c == 'r') { s += '\r'; }
                        else if (c == 'b') { s += '\b'; }
                        else if (c == 'u') { throw tokenizer_error(is); }
                        else { s += c; /* any character may be escaped */ }
                    } else {
                        s += c;
                    }
                }
                tokens.push_back(STRING);
                strings.push_back(s);
            } else if (c == '[') { is.get(); tokens.push_back(LEFT_SQUARE_BRACKET);
            } else if (c == ']') { is.get(); tokens.push_back(RIGHT_SQUARE_BRACKET);
            } else if (c == '{') { is.get(); tokens.push_back(LEFT_BRACE);
            } else if (c == '}') { is.get(); tokens.push_back(RIGHT_BRACE);
            } else if (c == ',') { is.get(); tokens.push_back(COMMA);
            } else if (c == ':') { is.get(); tokens.push_back(COLON);
            } else if (isdigit(c) || c == '-' || c == '+') {
                int sign = 1;
                double d = 0;
                double m = 10;
                double f = 1;
                if (c == '-') { is.get(), c=is.peek(); sign = -1; }
                else if (c == '+') { is.get(), c=is.peek(); }
                
                /* we don't support fractionals very cautiously, we don't
                 * support exponents, etc. It's quite lame.
                 */
                if (!isdigit(c)) throw tokenizer_error(is);
                for(;!is.eof() && (isdigit(c) || (c == '.' && m == 10));is.get(), c=is.peek()) {
                    if (c == '.' && m == 10) { m = 1; f = 0.1; continue; }
                    d = d * m + (c - '0') * f;
                    if (f < 1) f = f*0.1;
                }
                d = d * sign;
                numbers.push_back(d);
                tokens.push_back(NUMBER);
            } else if (c == 'n') {
                /* the only thing we can do know is detect one of the
                 * values null, true, false
                 */
                is.get(), c=is.peek();
                if (c != 'u') { throw tokenizer_error(is); } is.get(), c=is.peek();
                if (c != 'l') { throw tokenizer_error(is); } is.get(), c=is.peek();
                if (c != 'l') { throw tokenizer_error(is); } is.get(), c=is.peek();
                tokens.push_back(JNULL);
            } else if (c == 't') {
                is.get(), c=is.peek();
                if (c != 'r') { throw tokenizer_error(is); } is.get(), c=is.peek();
                if (c != 'u') { throw tokenizer_error(is); } is.get(), c=is.peek();
                if (c != 'e') { throw tokenizer_error(is); } is.get(), c=is.peek();
                tokens.push_back(BOOL);
                numbers.push_back(1);
            } else if (c == 'f') {
                is.get(), c=is.peek();
                if (c != 'a') { throw tokenizer_error(is); } is.get(), c=is.peek();
                if (c != 'l') { throw tokenizer_error(is); } is.get(), c=is.peek();
                if (c != 's') { throw tokenizer_error(is); } is.get(), c=is.peek();
                if (c != 'e') { throw tokenizer_error(is); } is.get(), c=is.peek();
                tokens.push_back(BOOL);
                numbers.push_back(0);
            } else {
                throw tokenizer_error(is);
            }
            /* c is either a character, or EOF (but if it's EOF, then
             * is.eof() is set). The thing is that we insist on the code
             * pattern "is.get(), c=is.peek()", but that actually peeks
             * twice in a row in the general case, and static analyzers
             * don't like it. To fix that, we test c against INT_MAX, in
             * cases where it is known to always be in the [0..255]
             * range.
             */
            ASSERT_ALWAYS(is.eof() || c != INT_MAX);
        };
        ctok = tokens.begin();
        cnum = numbers.begin();
        cstr = strings.begin();
        return true;
    }
};

static std::string json_flatten_string(std::string const & s) {
    std::ostringstream ss;
    ss << "\"";
    for(auto c : s) {
        if (c == '"') { ss << "\\\""; }
        else if (c == '\\') { ss << "\\\\"; }
        else if (c == '\n') { ss << "\\n"; }
        else if (c == '\t') { ss << "\\t"; }
        else if (c == '\r') { ss << "\\r"; }
        else if (c == '\b') { ss << "\\b"; }
        else if ((int) c < 0x20) { throw std::runtime_error("implement me"); }
        else ss << c;
    }
    ss << "\"";
    return ss.str();
}

std::string json_string::flatten() const {
    return json_flatten_string(data);
}

std::string json_number::flatten() const {
    std::ostringstream ss;
    ss << data;
    return ss.str();
}

std::string json_bool::flatten() const {
    std::ostringstream ss;
    ss << (data ? "true" : "false");
    return ss.str();
}

std::string json_hash::flatten() const {
    std::ostringstream ss;
    ss << "{";
    int some = 0;
    for(auto const & x : data) {
        if (some++) ss << ",";
        ss << json_flatten_string(x.first) << ':' << x.second->flatten();
    }
    ss << "}";
    return ss.str();
}

std::string json_array::flatten() const {
    std::ostringstream ss;
    ss << "[";
    int some = 0;
    for(auto const & x : data) {
        if (some++) ss << ',';
        ss << x->flatten();
    }
    ss << "]";
    return ss.str();
}

std::istream& operator>>(std::istream& in, json & f)
{
    /* The tokenizer eats the whole stream */
    try {
        json_parser P;
        P.tokenize(in);
        f = P.parse();
    } catch (json_parser::tokenizer_error const & p) {
        in.setstate(std::ios_base::failbit);
        return in;
    } catch (json_parser::parse_error const & p) {
        in.setstate(std::ios_base::failbit);
        return in;
    }

    return in;
}

std::ostream& operator<<(std::ostream& o, json_base const & f) {
    return o << f.flatten();
}
std::ostream& operator<<(std::ostream& o, json const & f) {
    return o << f->flatten();
}
