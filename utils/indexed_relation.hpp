#ifndef CADO_INDEXED_RELATION_HPP
#define CADO_INDEXED_RELATION_HPP


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
 *   NOTE: In that case, we need to know what are the active sides as
 *   well.
 *
 * Most methods are the same, and we prefer to use a bit of c++
 * templating so as to refactor some of the code.
 */

#include <cstddef>

#include <array>
#include <istream>
#include <ostream>
#include <vector>

#include "typedefs.h"
#include "relation.hpp"
#include "renumber.hpp"

struct indexed_relation_storage_base {
    // using &) = void (*action_t)(int side, std::vector<index_t>;
    std::array<int, 2> active_sides;
    protected:
    void set_active_sides(std::array<int, 2> const & aa) {
        active_sides = aa;
    }
    std::array<int, 2> get_active_sides() const { return active_sides; }
};

struct indexed_relation_normal_storage : public indexed_relation_storage_base {
    std::vector<index_t> data;
    protected:
    std::vector<index_t> & operator[](unsigned int) { return data; }
    template<typename ctr>
    struct fake_tmpl {
        ctr b;
        fake_tmpl(ctr d) : b(d) {}
        ctr begin() { return b; }
        ctr end() { return b+1; }
    };
    using fake = fake_tmpl<std::vector<index_t> *>;
    using const_fake = fake_tmpl<std::vector<index_t> const *>;
    fake containers() { return { &data }; }
    const_fake containers() const { return { &data }; }
    static constexpr const bool separates_sides = false;
    protected:
    /* It doesn't make sense to expose a parsing function for the storage
     * that splits by side: we *can't* do such a thing !
     *
     * XXX what does that mean after all?
     */
    int parse(relation_ab &, const char *line);

};

struct indexed_relation_byside_storage : public indexed_relation_storage_base {
    using sides_t = std::array<std::vector<index_t>, 2>;
    sides_t sides;
    protected:
    std::vector<index_t> & operator[](size_t side) { return sides[side]; }
    /* this is just for range-based for loops */
    sides_t & containers() { return sides; }
    sides_t const & containers() const { return sides; }
    static constexpr const bool separates_sides = true;
};

template<typename Storage>
struct indexed_relation_tmpl
    : public relation_ab
    , public Storage
{
    indexed_relation_tmpl() = default;
    ~indexed_relation_tmpl() = default;
    indexed_relation_tmpl(indexed_relation_tmpl<Storage> const &) = default;
    indexed_relation_tmpl& operator=(indexed_relation_tmpl<Storage> const &) = default;
    indexed_relation_tmpl(indexed_relation_tmpl<Storage> &&) = default;
    indexed_relation_tmpl& operator=(indexed_relation_tmpl<Storage> &&) = default;

    explicit indexed_relation_tmpl(relation_ab const & ab) : relation_ab(ab) {}
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
using indexed_relation = indexed_relation_tmpl<indexed_relation_normal_storage>;

/* forward-declare the template instantiation and the output
 * function, which are actually defined in the cpp file
 */
extern template struct indexed_relation_tmpl<indexed_relation_byside_storage>;
extern template 
std::ostream& operator<<(std::ostream& os, indexed_relation_tmpl<indexed_relation_byside_storage> const & rel);
using indexed_relation_byside = indexed_relation_tmpl<indexed_relation_byside_storage>;

#endif	/* CADO_INDEXED_RELATION_HPP */
