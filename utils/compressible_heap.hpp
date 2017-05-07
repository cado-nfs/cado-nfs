#ifndef COMPRESSIBLE_HEAP_HPP_
#define COMPRESSIBLE_HEAP_HPP_

#include <vector>
#include <cstddef>
#include <cassert>

/* This is functionally similar to a subset of vector<vector<T>>.
 *
 * We may access an item by its index:
 *                      H[i] returns a pair T*, Size
 *
 * We may also do       H.push_back(T * p, Size n)
 *                      which adds a "row" [p..p+n) to the heap.
 *                      
 * And                  H.kill(i)
 *                      which kills the i-th inserted row.
 *
 * Or                   H.shrink_value(i, n)
 *                      which modifies the stored length of the i-th
 *                      inserted row.
 *
 * The two latter operations might trigger compression of the underlying
 * data, if it so happens that its size has shrunk by more than 10% since
 * the last compression operation.
 *
 */
template<typename T, typename Size = size_t, int batch_size = 1024>
struct compressible_heap
{
    typedef std::pair<T *, Size> value_type;
    typedef std::pair<T *, Size> & reference;
    typedef std::pair<const T *, Size> const& const_reference;
    typedef std::pair<T *, Size> * pointer;
    typedef std::pair<const T *, Size> const * const_pointer;
public:
    class chunk {
        std::vector<T> data;
        value_type items[batch_size];
        size_t _size, killed;
        public:
        size_t allocated_bytes() const { return data.capacity()*sizeof(T) + batch_size * sizeof(value_type)+ sizeof(data); }
        size_t size() const { return _size; }
        /* full can also mean "full of emptiness"... */
        bool full() const { return _size == batch_size; }
        reference operator[](size_t i) { return items[i]; }
        const_reference operator[](size_t i) const { return ((const_pointer)items)[i]; }
        void compress() {
            if (killed <= data.size() / 10) return;
            /* compress the chunk in place, and update the pointers */
            size_t top = 0;
            ptrdiff_t temp[batch_size];
            for(int i = 0 ; i < batch_size ; i++) {
                temp[i] = -1;
                if (items[i].first == NULL) continue;
                const T * row = items[i].first;
                const T * endrow = row + items[i].second;
                /* It is guaranteed that this copy() won't go past the
                 * end of the data vector */
                std::copy(row, endrow, data.begin() + top);
                temp[i] = top;
                top += items[i].second;
            }
            data.resize(top);
            data.shrink_to_fit();
            for(int i = 0 ; i < batch_size ; i++)
                if (temp[i] >= 0) items[i].first = &data.front() + temp[i];
            killed = 0;
        }
        size_t kill(size_t i) {
            if (!items[i].first) return 0;
            items[i].first = NULL;
            killed += items[i].second;
            Size ret = items[i].second;
            // items[i].second = 0; // useless anyway.
            compress();
            return ret;
        }
        void shrink_value_unlocked(size_t i, size_t k) {
            assert(k <= items[i].second);
            killed += items[i].second - k;
            items[i].second = k;
        }
        void shrink_value(size_t i, size_t k) {
            shrink_value_unlocked(i, k);
            compress();
        }
        T * prepare_push_back(size_t n) {
            assert(_size < batch_size);
            items[_size].second = n;
            if (data.size() + n > data.capacity()) {
                /* we're going to invalidate the pointers... */
                ptrdiff_t temp[batch_size];
                for(size_t i = 0 ; i < _size ; i++) {
                    temp[i]=-1;
                    if (items[i].first == NULL) continue;
                    temp[i] = items[i].first - &data.front();
                }
                temp[_size++] = data.size();
                data.insert(data.end(), n, T());
                for(size_t i = 0 ; i < _size ; i++)
                    if (temp[i] >= 0) items[i].first = &data.front() + temp[i];
            } else {
                items[_size++].first = &data.front() + data.size();
                data.insert(data.end(), n, T());
            }
            return items[_size-1].first;
        }
        T * push_back(const T * p, size_t n) {
            return push_back(p, p+n);
        }
        template<typename RandomAccessIterator>
        T * push_back(RandomAccessIterator p, RandomAccessIterator q) {
            assert(_size < batch_size);
            assert(q >= p);
            items[_size].second = q - p;
            if (data.size() + (size_t) (q - p) > data.capacity()) {
                /* we're going to invalidate the pointers... */
                ptrdiff_t temp[batch_size];
                for(size_t i = 0 ; i < _size ; i++) {
                    temp[i]=-1;
                    if (items[i].first == NULL) continue;
                    temp[i] = items[i].first - &data.front();
                }
                temp[_size++] = data.size();
                data.insert(data.end(), p, q);
                for(size_t i = 0 ; i < _size ; i++)
                    if (temp[i] >= 0) items[i].first = &data.front() + temp[i];
            } else {
                /* If data.capacity() is zero, then we'll end up storing
                 * the null pointer. In turn, this might be
                 * misinterpreted as meaning that we have a killed row in
                 * that position.
                 */
                if (data.capacity() == 0)
                    data.reserve(64);
                items[_size++].first = &data.front() + data.size();
                data.insert(data.end(), p, q);
            }
            return items[_size-1].first;
        }
    };
    std::vector<chunk> chunks;
public:
    class iterator {/*{{{*/
        compressible_heap& H;
    public:
        size_t i;
    private:
        iterator(compressible_heap& H, size_t i) : H(H), i(i) {}
        iterator nextvalid() {
            for(; i < H.size() && !H[i].first ; ++i);
            return *this;
        }
        public:
        iterator() {}
        iterator operator++(int) {
            iterator v = *this;
            ++i;
            nextvalid();
            return v;
        }
        iterator operator++() {
            ++i;
            return nextvalid();
        }
        reference operator*() { return H[i]; }
        pointer operator->() { return &H[i]; }
        bool operator==(iterator const& o) const { return &H == &o.H && i == o.i; }
        bool operator!=(iterator const& o) const { return !(i == o.i); }
        size_t index() const { return i; }
        friend class compressible_heap;
    };
    iterator begin() { return iterator(*this, 0).nextvalid(); }
    iterator end() { return iterator(*this, size()); }
/*}}}*/
    class const_iterator {/*{{{*/
        compressible_heap const & H;
        size_t i;
        const_iterator(compressible_heap const & H, size_t i) : H(H), i(i) {}
        const_iterator nextvalid() {
            for(; i < H.size() && !H[i].first ; ++i);
            return *this;
        }
        public:
        const_iterator() {}
        const_iterator operator++(int) {
            iterator v = *this;
            size_t n = H.size();
            ++i;
            nextvalid();
            return v;
        }
        const_iterator operator++() {
            size_t n = H.size();
            ++i;
            nextvalid();
            return *this;
        }
        const_reference operator*() { return H[i]; }
        const_pointer operator->() { return &H[i]; }
        bool operator==(iterator const& o) const { return &H == &o.H && i == o.i; }
        bool operator!=(iterator const& o) const { return !(i == o.i); }
        size_t index() const { return i; }
        friend class compressible_heap;
    };
    const_iterator begin() const { return const_iterator(*this, 0).nextvalid(); }
    const_iterator end() const { return const_iterator(*this, size()); }
/*}}}*/
    T * push_back(std::vector<T> const & row) { return push_back(&row.front(), row.size()); }
    T * prepare_push_back(size_t n) {
        if (chunks.empty() || chunks.back().full())
            chunks.push_back(chunk());
        return chunks.back().prepare_push_back(n);
    }
     T * push_back(const T * p, size_t n) {
        if (chunks.empty() || chunks.back().full())
            chunks.push_back(chunk());
        return chunks.back().push_back(p, n);
    }
    template<typename RandomAccessIterator>
    T * push_back(RandomAccessIterator p, RandomAccessIterator q) {
        if (chunks.empty() || chunks.back().full())
            chunks.push_back(chunk());
        return chunks.back().push_back(p, q);
    }
    reference operator[](size_t i) { return chunks[i / batch_size][i % batch_size]; }
    const_reference operator[](size_t i) const { return chunks[i / batch_size][i % batch_size]; }
    size_t kill(size_t i) { return chunks[i / batch_size].kill(i % batch_size); }
    void shrink_value_unlocked(size_t i, size_t k) {
        chunks[i / batch_size].shrink_value_unlocked(i % batch_size, k);
    }
    void shrink_value_unlocked(iterator const& it, size_t k) { shrink_value_unlocked(it.i, k); }
    void shrink_value(size_t i, size_t k) {
        chunks[i / batch_size].shrink_value(i % batch_size, k);
    }
    void shrink_value(iterator const& it, size_t k) { shrink_value(it.i, k); }
    bool empty() const { return chunks.empty(); }
    size_t size() const {
        if (empty()) return 0;
        return (chunks.size()-1)*batch_size + chunks.back().size();
    }
    void compress() { for(auto & c: chunks) c.compress(); }
    size_t allocated_bytes() const {
        size_t res=sizeof(*this);
        for(auto const& c: chunks) res += c.allocated_bytes();
        return res;
    }
};


#endif	/* COMPRESSIBLE_HEAP_HPP_ */
