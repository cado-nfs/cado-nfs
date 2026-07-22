
get_property(c90_21505_chrono GLOBAL PROPERTY c90_21505_chrono)
cado_define_test(c90_21505_render_chronogram
    PROGRAM
    ${CADO_NFS_SOURCE_DIR}/scripts/chronograms/render.py ${c90_21505_chrono} -o
    ${c90_21505_chrono}.html
    FIXTURES_REQUIRED c90_21505.chrono
)
