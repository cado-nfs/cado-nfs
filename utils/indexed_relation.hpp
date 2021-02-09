#ifndef INDEXED_RELATION_HPP_
#define INDEXED_RELATION_HPP_


/* indexed relations are those that appear after dup2. These are read by
 * purge or shrink_rels, for instance. They're also produced by fake_rels.
 *
 * Note that the purge tools do not use this c++ structure, they rely
 * rather on filter_io's earlyparsed_relation things.
 *
 * There are two possible underlying storage layouts for
 * indexed_relation:
 *
 * - either a simple vector of indices
 * - or a vector of vectors of indices, one vector for each side.
 *
 * Most methods are the same, and we prefer to use a bit of c++
 * templating so as to refactor some of the code.
 */

#include <vector>
#include "typedefs.h"
#include "relation.hpp"
#include "renumber.hpp"

struct indexed_relation_storage_base {
    // typedef void (*action_t)(unsigned int side, std::vector<index_t> &);
};

struct indexed_relation_normal_storage : public indexed_relation_storage_base {
    std::vector<index_t> data;
    protected:
    void set_nsides(unsigned int) {}
    unsigned int get_nsides() const { return 1; }
    std::vector<index_t> & operator[](unsigned int) { return data; }
    template<typename ctr>
    struct fake_tmpl {
        ctr b;
        fake_tmpl(ctr d) : b(d) {}
        ctr begin() { return b; }
        ctr end() { return b+1; }
    };
    typedef fake_tmpl<std::vector<index_t> *> fake;
    typedef fake_tmpl<std::vector<index_t> const *> const_fake;
    fake containers() { return fake(&data); }
    const_fake containers() const { return const_fake(&data); }
    static constexpr const bool separates_sides = false;
    protected:
    /* It doesn't make sense to expose a parsing function for the storage
     * that splits by side: we *can't* do such a thing !
     */
    int parse(relation_ab &, const char *line);

};

struct indexed_relation_byside_storage : public indexed_relation_storage_base {
    std::vector<std::vector<index_t>> sides;
    protected:
    void set_nsides(unsigned int n) {
        sides.assign(n, std::vector<index_t>());
    }
    unsigned int get_nsides() const { return sides.size(); }
    std::vector<index_t> & operator[](unsigned int side) { return sides[side]; }
    /* this is just for range-based for loops */
    std::vector<std::vector<index_t>> & containers() { return sides; }
    std::vector<std::vector<index_t>> const & containers() const { return sides; }
    static constexpr const bool separates_sides = true;
};

template<typename Storage>
struct indexed_relation_tmpl
    : public relation_ab
    , public Storage
{
    indexed_relation_tmpl() = default;
    indexed_relation_tmpl(indexed_relation_tmpl<Storage> const &) = default;
    indexed_relation_tmpl& operator=(indexed_relation_tmpl<Storage> const &) = default;
    indexed_relation_tmpl(indexed_relation_tmpl<Storage> &&) = default;
    indexed_relation_tmpl& operator=(indexed_relation_tmpl<Storage> &&) = default;

    indexed_relation_tmpl(relation_ab const & ab) : relation_ab(ab) {}
    indexed_relation_tmpl(relation const & rel, renumber_t const & R);

    void sort();
    void compress(bool dl);
    /* This is used by shrink_rels and fake_rels */
    void shrink(double shrink_factor, index_t shrink_threshold = 0);
    template<typename T>
    friend std::ostream& operator<<(std::ostream&, indexed_relation_tmpl<T> const&);
    template<typename T>
    friend std::istream& operator>>(std::istream&, indexed_relation_tmpl<T>&);
};

template<typename Storage>
std::ostream& operator<<(std::ostream& os, indexed_relation_tmpl<Storage> const & rel);

template<typename Storage>
std::istream& operator>>(std::istream& is, indexed_relation_tmpl<Storage>& rel);

/* forward-declare the template instantiation and the input and output
 * functions, which are actually defined in the cpp file
 */
extern template struct indexed_relation_tmpl<indexed_relation_normal_storage>;
extern template 
std::ostream& operator<<(std::ostream& os, indexed_relation_tmpl<indexed_relation_normal_storage> const & rel);
extern template
std::istream& operator>>(std::istream& is, indexed_relation_tmpl<indexed_relation_normal_storage>& rel);
typedef indexed_relation_tmpl<indexed_relation_normal_storage> indexed_relation;

/* forward-declare the template instantiation and the output
 * function, which are actually defined in the cpp file
 */
extern template struct indexed_relation_tmpl<indexed_relation_byside_storage>;
extern template 
std::ostream& operator<<(std::ostream& os, indexed_relation_tmpl<indexed_relation_byside_storage> const & rel);
typedef indexed_relation_tmpl<indexed_relation_byside_storage> indexed_relation_byside;

#endif	/* INDEXED_RELATION_HPP_ */
