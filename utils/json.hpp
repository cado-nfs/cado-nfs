#ifndef JSON_HPP_
#define JSON_HPP_

#include <stdexcept>
#include <vector>
#include <string>
#include <map>
#include <memory>

/*
 * This it _not_ the correct way to parse a json file. The correct way
 * would be to embed https://github.com/nlohmann/json ; but that means 1
 * megabyte to add, which is slightly large. So we'll do a q&d thing for
 * the moment, and think of the possibility of switching to something
 * better someday if the need arises.
 *
 * Regarding the super-verbose "using json_base::operator X" constructs,
 * see https://bbs.archlinux.org/viewtopic.php?pid=1830516#p1830516 and https://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings
 *
 */

struct json_base {
    struct error: public std::exception {
        std::string m;
        error(const char * s) : m(s) {}
        const char * what() const noexcept override { return m.c_str(); }
    };
    enum type_t { HASH, ARRAY, STRING, NUMBER, BOOL };
    type_t type;
    protected:
    json_base(type_t x) : type(x) {}
    public:
    inline type_t get_type() const { return type; }
    virtual std::string flatten() const = 0;
    virtual json_base & operator[](int) { throw error("index"); }
    virtual json_base & operator[](std::string const &) { throw error("key"); }
    virtual json_base & operator[](char const * a) { return operator[](std::string(a)); }
    virtual json_base const & operator[](int) const { throw error("index"); }
    virtual json_base const & operator[](std::string const &) const { throw error("key"); }
    virtual json_base const & operator[](char const * a) const { return operator[](std::string(a)); }
    /* conversion to string is not the same as flattening */
    virtual operator std::string() const { throw error("type"); }
    virtual operator double() const { throw error("type"); }
    virtual operator long() const { throw error("type"); }
    virtual size_t size() const { throw error("type"); } 
    virtual ~json_base() {}
    virtual json_base * clone() const { throw error("type"); }
};

struct json: public std::unique_ptr<json_base> {
    /* proxy the things that we like to have, and delegate that to the
     * subobject.
     */
    json() = default;
    json(json const & x) : std::unique_ptr<json_base>(x->clone()) {}
    json& operator=(json const & x) { reset(x->clone()); return *this; }
    json(json_base * p) : std::unique_ptr<json_base>(p) {}
    template<typename T> inline json_base const & operator[](T const &x) const { return (*get())[x]; }
    template<typename T> inline json_base & operator[](T const &x) { return (*get())[x]; }
};

struct json_hash : public json_base {
    std::map<std::string, json> data;
    json_hash() : json_base(HASH) {}
    virtual json_base & operator[](std::string const & x) override { return *data[x]; }
    std::pair<std::map<std::string, json>::iterator, bool> emplace(std::string const & x, json_base const & c) {
        return data.emplace(x, json(c.clone()));
    }
    virtual json_base const & operator[](std::string const & x) const override {
        return *data.find(x)->second;
    }
    virtual std::string flatten() const override;
    virtual ~json_hash() override {}
    virtual json_base * clone() const override { return new json_hash(*this); }

};

struct json_array : public json_base {
    std::vector<json> data;
    json_array() : json_base(ARRAY) {}
    virtual json_base & operator[](int x) override { return *data[x]; }
    virtual json_base const & operator[](int x) const override { return *data[x]; }
    virtual std::string flatten() const override;
    virtual ~json_array() override {}
    virtual size_t size() const override { return data.size(); }
    virtual json_base * clone() const override { return new json_array(*this); }
};

struct json_string : public json_base {
    std::string data;
    json_string() : json_base(STRING) {}
    json_string(std::string const & x) : json_base(STRING), data(x) {}
    virtual operator std::string() const override { return data; }
    virtual std::string flatten() const override;
    virtual ~json_string() override {}
    virtual json_base * clone() const override { return new json_string(*this); }
    private:
    using json_base::operator double;
    using json_base::operator long;
};

struct json_number : public json_base {
    double data;
    json_number() : json_base(NUMBER) {}
    json_number(double x) : json_base(NUMBER), data(x) {}
    virtual operator double() const override { return data; }
    virtual operator long() const override { return data; }
    virtual std::string flatten() const override;
    virtual ~json_number() override {}
    virtual json_base * clone() const override { return new json_number(*this); }
    private:
    using json_base::operator std::string;
};

struct json_bool : public json_base {
    bool data;
    json_bool() : json_base(BOOL) {}
    json_bool(double x) : json_base(BOOL), data(x) {}
    virtual operator long() const override { return data; }
    virtual std::string flatten() const override;
    virtual ~json_bool() override {}
    virtual json_base * clone() const override { return new json_bool(*this); }
    private:
    using json_base::operator double;
    using json_base::operator std::string;
};

struct json_null : public json_base {
    json_null() : json_base(STRING) {}
    virtual std::string flatten() const override { return "null"; }
    virtual ~json_null() override {}
    virtual json_base * clone() const override { return new json_null(); }
};

std::ostream& operator<<(std::ostream& o, json const & f);
std::istream& operator>>(std::istream& in, json & f);
std::ostream& operator<<(std::ostream& o, json_base const & f);

#endif	/* JSON_HPP_ */
