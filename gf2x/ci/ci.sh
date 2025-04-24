
# Most of the common steps of the build are specified here. The
# lower-level plumbing is done by the `ci-backend` project.

export project_name=gf2x 
: ${CI_PROJECT_URL=https://gitlab.inria.fr/gf2x/gf2x} 
export build_system=autotools 
 
# needs_bc=1 
# needs_python=1 
needs_perl=1 
# needs_optional_hwloc=1 
# needs_optional_ecm=1 
# needs_optional_fmt=1 
needs_gmp=1 

project_package_selection() { : ; }

