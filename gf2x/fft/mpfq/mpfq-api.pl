# This is a trimmed down version of the mpfq api. We select only what
# gf2x-fft needs.
'#ID'	=> 'MPFQ Field API (trimmed down for gf2x-fft)',

# This file has the form of a parseable perl input file, but clearly it
# needs not be so. Parsing is done in Mpfq/conf.pm


# The first line is magical, and specifies the types that have to be defined.
'#TYPES'        => [ qw/field dst_field src_field/ ],

'#TYPES'        => [ qw/elt dst_elt src_elt/ ],

restrict('+URE'),
'#TYPES'        => [ qw/elt_ur dst_elt_ur src_elt_ur/ ],
restrict('-URE'),

'#COMMON_ARGS'  => [ ],

'#COMMON_ARGS'  => [ 'dst_field' ],


###############################################
hdr('Elementary assignment functions'),

'set'		=> 'dst_elt src_elt',
'set_zero'	=> 'dst_elt',

############################
hdr('Comparison functions'),

'is_zero'	=> 'int <- src_elt',
###############################################
hdr('Arithmetic operations on elements'),

'add'		=> 'dst_elt src_elt src_elt',
'mul'		=> 'dst_elt src_elt src_elt',
'sqr'		=> 'dst_elt src_elt',

###############################################
restrict('+URE'),
hdr('Operations involving unreduced elements'),

'elt_ur_set'	=> 'dst_elt_ur src_elt_ur',
'elt_ur_set_elt'	=> 'dst_elt_ur src_elt',
'elt_ur_set_zero'	=> 'dst_elt_ur',
'elt_ur_add'	=> 'dst_elt_ur src_elt_ur src_elt_ur',

'mul_ur'	=> 'dst_elt_ur src_elt src_elt',
'sqr_ur'	=> 'dst_elt_ur src_elt',

'reduce'	=> [ 'dst_elt dst_elt_ur',
	"reduces the dst_elt_ur operand, store the reduced operand in the dst_elt
	operand. Note that the unreduced operand is clobbered." ],

restrict('-URE'),

###############################################
#hdr('Vector functions'),
'#TYPES'        => [ qw/vec dst_vec src_vec/ ],
'vec_set_zero'          => [ 'dst_vec uint', "zeroes out a vector" ],

###############################################
# we only export types (mpfq bug)
#hdr('Polynomials'),
'#TYPES'        => [ qw/poly dst_poly src_poly/ ],


# It is normal for this file to end with a comma.
