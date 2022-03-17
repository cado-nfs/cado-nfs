# This is a buggy strategy file
define blah
    [S0:PP1-27, 150, 4000]
    S0: ECM-B12(12345),200,2000
define foo
    S0: ECM-B12,200,2000
# foo after blah should fail, because a default parameter choice cannot
# work after a non-default one.
r0=10,r1=20
    use blah
    use foo
