__attribute__((noinline))
long frob(const long *t, long *out, int stride, int n) {
    long x = 0;
    for (int i = 0; i < n; i++, t += stride) {
        x -= *t;
        out[i] = *t * 2;
    }
    return x;
}

int main(void) {
    long data[2] = { -10, -20 };
    long out[2];
    long x = frob(data, out, 1, 2);
    return (x == 30) ? 0 : 1;
}
