# This is an example strategy file. The longer one is F9.st in the
# directory above. This is deliberately tortured in order to exert some corner cases of the parser. All of this should parse successfully.

r0=10,r1=20]
    [S0:PM1,200,2000

define 0a [S0:ECM-B12,50,5000]
define 0b [S0:ECM-B12,50,5000]
define 0c S0: ECM-B12, 50, 5000

define 0A [S0:PP1-65,150,4500]  use 0a use 0a use 0a
define 0B S0:PP1-65,150,4500    use 0a use 0b
define 0C S0: PP1-65, 150, 4500 use 0c use 0b

[r0=1, r1=45] : (p = 1.000000, t = 13.149463)
    use 0A
[r0=1, r1=46] # comment at end of line
    use 0B
[r0=1, r1=47]
    use 0C
[r0=2, r1=47]
        S0:ECM-B12,200,2000
        S0:ECM-B12(12345),200,2000
        S0:ECM-B12(6789),200,2000
        S0:ECM-M12,200,2000
        S0:ECM-M12(12345),200,2000
        S0:ECM-M16,200,2000
        S0:ECM-TM12(12345),200,2000
