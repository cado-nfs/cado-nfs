cado_define_test(test_mod_ul.c LIBRARIES utils tests ARGUMENTS 1000)
cado_define_test(test_modredc_ul.c LIBRARIES utils tests ARGUMENTS 1000)
# These two seem to time out easily in coverage mode
cado_define_test(test_modredc_15ul.c LIBRARIES utils tests ARGUMENTS 1000
    TIMEOUT 60
)
cado_define_test(test_modredc_2ul2.c LIBRARIES utils tests ARGUMENTS 1000
    TIMEOUT 120
)
cado_define_test(test_modul_poly.c LIBRARIES utils tests)
cado_define_test(test_modredcul_batchinv.c LIBRARIES utils tests)
cado_define_test(test_modredcul_batchinv_ul.c LIBRARIES utils tests)
cado_define_test(test_modredc2ul2_Q_to_Fp.c LIBRARIES utils tests)
