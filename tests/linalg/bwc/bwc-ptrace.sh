#!/usr/bin/env bash

old_setx=$(shopt -po xtrace || :)

set -e

exec < /dev/null

if [ "$CADO_DEBUG" ] ; then set -x ; fi

# This script is intended *for testing only*. It's used for, e.g.,
# coverage tests. This even goes with dumping all intermediary data to
# magma-parseable format. So don't try to use this script for real
# examples. Real examples require the user to craft his own bwc.pl
# command line with care (command-lines as created by this tool might be
# a source of inspiration, though).

pass_bwcpl_args=()

while [ $# -gt 0 ] ; do
    a="$1"
    shift
    if [ "$a" = "--" ] ; then
        break
    else
        if [[ $a =~ ^(prime|m|n|wdir|interval|seed|nullspace|mpi|thr|simd)= ]] ; then
            # There are quite a few parameters that we parse here and
            # pass to bwcpl anyway, so let's not put them in
            # pass_bwcpl_args.
            # (in fact, this is the case of most parameters)
            :
        # (try to remove both the exception here and the offending
        # exception later on)
        # elif [[ $a =~ ^(tolerate_failure|stop_at_step|keep_rolling_checkpoints|checkpoint_precious|skip_online_checks|interleaving) ]] ; then
        #     # basically the same story here, except that there seem to
        #     # be two ways these arguments get passed to bwc.pl : we
        #     # forcibly add them to ${common[@]} later on in this file,
        #     # but why do we do that? Is it only to provide for the case
        #     # where these arguments are exported via the shell
        #     # environment?
        #     :
        elif [[ $a =~ ^(bindir|mats|pre_wipe|random_matrix_size|random_matrix_minkernel|script_steps|nrhs|sage|magma|wordsize) ]] ; then
            # and there are even parameters that only make sense here.
            :
        else
            pass_bwcpl_args+=("$a")
        fi
        eval "$a"
    fi
done


# various configuration variables. environment can be used to override them
: ${scriptpath=$0}
: ${m=8}
: ${n=4}
: ${prime=4148386731260605647525186547488842396461625774241327567978137}
: ${mpi:=1x1}
: ${lingen_mpi:=1x1}
: ${thr:=2x2}
# Set the "matrix" variable in order to work on a real matrix.
: ${matrix=}
# Set the "bindir" variable to use pre-built binaries (must point to the
# directories with the bwc binaries)
: ${bindir=}

# This is related to the random matrix which gets automatically created
# if no matrix was given on the cmdline.
: ${random_matrix_size=1000}
# there is no random_matrix_coeffs_per_row ; see random_matrix.c
: ${random_matrix_maxcoeff=10}
: ${random_matrix_minkernel=10}
: ${mats=$HOME/Local/mats}
: ${pre_wipe=}
: ${seed=$RANDOM}
: ${balancing_options=reorder=columns}
: ${script_steps=wipecheck,matrix,bwc.pl/:complete,bwccheck,magma}

pass_bwcpl_args+=("seed=$seed")

: ${wordsize=64}

# XXX note that $wdir is wiped out by this script !
: ${wdir=/tmp/bwcp}

# By default we don't enable magma. Just say magma=magma (or
# magma=/path/to/magma) to get (very expensive) magma checks.
: ${magma=}

# We also have sagemath testing. For the moment it's quite limited, but
# eventually we expect that it will catch up with what we do in magma.
: ${sage=}

# This has to be provided as auxiliary data for the GF(p) case (or we'll
# do without any rhs whatsoever).
: ${rhs=}

# For the rest of this script, we'll prefer the variable name "rhsfile"
rhsfile="$rhs"
: ${nullspace=right}
: ${interval=50}
# : ${mm_impl=basicp}

nullspace_lowercase=$(echo $nullspace | tr A-Z a-z)

if ! type -p seq >/dev/null ; then
    seq() {
        first="$1"
        shift
        if [ "$2" ] ; then
            incr="$1"
            shift
        else
            incr=1
        fi
        last="$1"
        while [ "$first" -le "$last" ] ; do
            echo "$first"
            let first=first+incr
        done
    }
fi

# inject the variables that were provided by guess_mpi_configs
if [ "$mpi" ] ; then
    eval "$exporter_mpirun"
    eval "$exporter_mpi_extra_args"
    pass_bwcpl_args+=(mpi_extra_args="${mpi_extra_args[*]}")
fi

usage() {
    cat >&2 <<-EOF
    Usage: $scriptpath [param1=value param2=value2 ...]
    This runs a complete block Wiedemann algorithm, given the parameters
    specified by various environment variables (the command line can also
    be used to set those variables).
EOF
    exit 1
}

argument_checking() {
    case "$nullspace" in
        left|right) ;;
        LEFT|RIGHT) ;;
        *) echo "\$nullspace must be left or right" >&2; usage;;
    esac
    if [ "$prime" != 2 ] ; then
        if [ "$matrix" ] && [ "$rhsfile" ] && [ "$nrhs" ] ; then
            # inhomogeneous mod p
            :
        elif [ "$matrix" ] && ! [ "$rhsfile" ] && ! [ "$nrhs" ] ; then
            # homogeneous mod p (bogus as of 755e9c5).
            :
        elif ! [ "$matrix" ] && ! [ "$rhsfile" ] ; then
            if ! [ "$random_matrix_size" ] || ! [ "$random_matrix_maxcoeff" ] ; then
                echo "Please give random matrix dimensions" >&2
                exit 1
            fi
            if ! [ "$nrhs" ] && ! [ "$random_matrix_minkernel" ] ; then
                echo "Please give random matrix dimensions" >&2
                exit 1
            fi
        else
            echo "Unsupported combination of arguments for specifying system" >&2
            echo "Detected arguments:" >&2
            for a in matrix rhsfile nrhs random_matrix_size random_matrix_maxcoeff ; do
                if [ "${!a}" ] ; then echo "$a=${!a}" >&2 ; fi
            done
            echo "Supported combinations:">&2
            echo "matrix rhsfile nrhs" >&2
            echo "matrix" >&2
            echo "random_matrix_size random_matrix_maxcoeff random_matrix_minkernel" >&2
            echo "random_matrix_size random_matrix_maxcoeff nrhs" >&2
            usage
        fi
    fi
}

derived_variables() {
    if ! [ -d $mats ] ; then mats=$wdir; fi
    if [ "$prime" = 2 ] ; then
        splitwidth=64
    else
        splitwidth=1
    fi
    : ${simd=$splitwidth}
}

prepare_wdir() {
    if [[ $script_steps =~ wipecheck ]] ; then
        if [ -d $wdir ] ; then
            if [ "$pre_wipe" ] ; then
                rm -rf $wdir 2>/dev/null
            else
                echo "Won't wipe $wdir unless \$pre_wipe is set" >&2
                exit 1
            fi
        fi
        mkdir $wdir
    elif [[ $script_steps =~ keepdir ]] ; then
        if ! [ -d $wdir ] ; then
            echo "cannot work with keepdir: $wdir does not exist" >&2
            exit 1
        fi
    else
        echo "Want either wipecheck or keepdir in script_steps" >&2
        exit 1
    fi
}


create_test_matrix_if_needed() {
    if [ "$matrix" ] ; then
        if ! [ -e "$matrix" ] ; then
            echo "matrix=$matrix inaccessible" >&2
            usage
        fi
        # Get absolute path.
        matrix=$(readlink /proc/self/fd/99 99< $matrix)
        if [ -e "$rhsfile" ] ; then
            rhsfile=$(readlink /proc/self/fd/99 99< $rhsfile)
        fi
        return
    fi

    # It's better to look for a kernel which is not trivial. Thus
    # specifying --kright for random generation is a good move prior to
    # running this script for nullspace=right
    kside="--k${nullspace_lowercase}"

    # We only care the binary matrix, really. Nevertheless, the random
    # matrix is created as text, and later transformed to binary format.
    # We also create the auxiliary .rw and .cw files.

    # We're not setting the density, as there is an automatic setting in
    # random_matrix for that (not mandatory though, since
    # nrows,ncols,density, may be specified in full).

    rmargs=()
    # defaults, some of the subcases below tweak that.
    nrows=`echo $random_matrix_size | cut -d, -f1`
    ncols=`echo $random_matrix_size | cut -d, -f2`   # = nrows if no comma
    outer_nrows=`echo $random_matrix_size | cut -d, -f1`
    outer_ncols=`echo $random_matrix_size | cut -d, -f2`   # = nrows if no comma
    escaped_size=$(echo $random_matrix_size | tr , _)
    if [ "$prime" = 2 ] ; then
        basename=$mats/t${escaped_size}
        matrix="$basename.matrix.bin"
        rmargs+=(--k${nullspace_lowercase} ${random_matrix_minkernel})
        # ncols=
    elif ! [ "$nrhs" ] ; then
        basename=$mats/t${escaped_size}p
        matrix="$basename.matrix.bin"
        rmargs+=(--k${nullspace_lowercase} ${random_matrix_minkernel})
        rmargs+=(-c ${random_matrix_maxcoeff})
        rmargs+=(-Z)
    else
        # This is an experimental mode. In the DLP context, we have a
        # matrix with N rows, N-r ideal columns, and r Schirokauer maps.   
        # let's say random_matrix_size is N and nrhs is r. We'll generate
        # a matrix with N rows and N-r columns, and later pad the column
        # width data with r zeroes.
        ncols=$((nrows-nrhs))
        basename=$mats/t${escaped_size}p+${nrhs}
        # We want something proportional to log(N)^2 as a density.
        # There's no undebatable support for this heuristic, it's just an
        # arbitrary choice. We'll do that without relying on external
        # tools for computing the log: just plain stupid shell. That's
        # gonna be good enough. The ratio below is taken as 2, which fits
        # well the relevant density ranges for the matrices we wish to
        # consider.
        density=$((${#random_matrix_size}*${#random_matrix_size}*2))
        if [ "$density" -lt 12 ] ; then density=12; fi
        rmargs+=(-d $density -Z)
        matrix="$basename.matrix.bin"
        rhsfile="$basename.rhs.txt"
        rmargs+=(-c ${random_matrix_maxcoeff})
        rmargs+=(rhs="$nrhs,$prime,$rhsfile,${nullspace_lowercase}")
    fi
    rmargs=($nrows $ncols -s $seed "${rmargs[@]}" --freq --binary --output "$matrix")
    rwfile=${matrix%%bin}rw.bin
    cwfile=${matrix%%bin}cw.bin
    ncols=$outer_ncols
    nrows=$outer_nrows
    if [[ $script_steps =~ matrix ]] ; then
        ${bindir}/random_matrix "${rmargs[@]}"
        data_ncols=$((`wc -c < $cwfile` / 4))
        data_nrows=$((`wc -c < $rwfile` / 4))
        if [ "$data_ncols" -lt "$ncols" ] ; then
            if [ "$prime" = 2 ] || ! [ "$nrhs" ] ; then
                echo "padding $cwfile with $((ncols-data_ncols)) zero columns"
                dd if=/dev/zero bs=4 count=$((ncols-data_ncols)) >> $cwfile
            fi
        fi
        if [ "$data_nrows" -lt "$nrows" ] ; then
            echo "padding $cwfile with $((nrows-data_nrows)) zero rows"
            dd if=/dev/zero bs=4 count=$((nrows-data_nrows)) >> $rwfile
        fi
    fi
}

# This is only useful in the situation where an input matrix has been
# provided by the user (which may be ither in text or binary format). In
# this case, we must make sure that we have the .cw and .rw files too.
create_auxiliary_weight_files() {
    if ! [[ $script_steps =~ matrix ]] ; then
        return
    fi
    if [ "$prime" != 2 ] ; then withcoeffs=--withcoeffs ; fi
    case "$matrix" in
        *.txt)
            echo "please convert the matrix to binary first" >&2
            exit 1
            ;;
        *.bin)
            rwfile=${matrix%%bin}rw.bin
            cwfile=${matrix%%bin}cw.bin
            if [ "$rwfile" -nt "$matrix" ] && [ "$cwfile" -nt "$matrix" ] ; then
                echo "Taking existing $rwfile, $cwfile as accompanying $matrix"
            else
                $bindir/mf_scan2  $withcoeffs --mfile $matrix --freq
            fi
            ;;
    esac
    ncols=$((`wc -c < $cwfile` / 4))
    nrows=$((`wc -c < $rwfile` / 4))
}

prepare_common_arguments() {
    # This sets the common arguments for all bwc binaries
    common=(
        matrix=$matrix
        balancing_options=$balancing_options
        mpi=$mpi
        thr=$thr
        m=$m
        n=$n
        wdir=$wdir
        prime=$prime
        nullspace=$nullspace
        interval=$interval
        simd=$simd
    )
    if [ "$mm_impl" ] ; then common+=(mm_impl=$mm_impl) ; fi
    if [ "$prime" != 2 ] ; then common+=(lingen_mpi=$lingen_mpi) ; fi
}

argument_checking
derived_variables
prepare_wdir

if ! [ "$matrix" ] ; then
    # This also sets rwfile cwfile nrows ncols
    create_test_matrix_if_needed
else
    # This also sets rwfile cwfile nrows ncols
    create_auxiliary_weight_files
fi
for f in matrix cwfile rwfile ; do
    if ! [ -f "${!f}" ] ; then echo "Missing file $f=${!f}" >&2 ; exit 1 ; fi
done

# if [ "$magma" ] ; then interval=1 ; fi

prepare_common_arguments

if [ "$rhsfile" ] ; then
    common+=(rhs=$rhsfile)
fi

# for v in tolerate_failure stop_at_step keep_rolling_checkpoints checkpoint_precious skip_online_checks interleaving ; do
#     if [ "${!v}" ] ; then common+=("$v=${!v}") ; fi
# done

if ! [[ $script_steps =~ magma ]] ; then
    magma=
fi

if [ "$magma" ] || [ "$sage" ] ; then
    common+=(save_submatrices=1)
fi

if [ "$magma" ] ; then
    echo "### Enabling magma checking ###"
    if [ -d "$magma" ] ; then
        if [ -x "$magma/magma" ] ; then
            magma="$magma/magma"
        else
            echo "If \$magma is a directory, we want to find \$magma/magma !" >&2
            exit 1
        fi
    fi
fi

if [ "$sage" ] ; then
    echo "### Enabling SageMath checking ###"
fi

if ! [ "$magma" ] && ! [ "$sage" ] ; then
    rc=0
    if [[ $script_steps =~ bwc\.pl/([a-z:/]*) ]] ; then
        IFS=/ read -a bwcpl_steps <<< "${BASH_REMATCH[1]}"
        for s in "${bwcpl_steps[@]}" ; do
            if $bindir/bwc.pl "$s" "${common[@]}" "${pass_bwcpl_args[@]}" ; then
                :
            else
                rc=$?
                break
            fi
        done
    fi
    if [[ $script_steps =~ bwccheck ]] ; then
        $bindir/bwc.pl :mpirun_single -- $bindir/bwccheck prime=$prime m=$m n=$n -- $wdir/[ACVFS]* > $wdir/bwccheck.log
        grep NOK $wdir/bwccheck.log || :
    fi
    exit $rc
else
    set +e
    rc=0
    if [[ $script_steps =~ bwc\.pl/([a-z:/]*) ]] ; then
        IFS=/ read -a bwcpl_steps <<< "${BASH_REMATCH[1]}"
        for s in "${bwcpl_steps}" ; do
            if $bindir/bwc.pl "$s" "${common[@]}" "${pass_bwcpl_args[@]}" ; then
                :
            else
                rc=$?
                break
            fi
        done
    fi
    if [[ $script_steps =~ bwccheck ]] ; then
        $bindir/bwc.pl :mpirun_single -- $bindir/bwccheck prime=$prime m=$m n=$n -- $wdir/[ACVFS]* > $wdir/bwccheck.log
        grep NOK $wdir/bwccheck.log
    fi
    eval $old_setx
    if [[ $script_steps =~ bwc\.pl([a-z:/]*) ]] ; then
        if [ "$rc" = 0 ] ; then
            echo " ========== SUCCESS ! bwc.pl returned true ========== "
            echo " ========== SUCCESS ! bwc.pl returned true ========== "
            echo " ========== SUCCESS ! bwc.pl returned true ========== "
        else
            echo " ########## FAILURE ! bwc.pl returned false ########## "
            echo " ########## FAILURE ! bwc.pl returned false ########## "
            echo " ########## FAILURE ! bwc.pl returned false ########## "
        fi
    fi
fi

eval $old_setx

magma_sage_check_parameters() { # {{{
    # Required parameters below this point:

    # wdir
    # matrix
    # m
    # n
    # prime
    # interval
    # splitwidth
    # cmd (pay attention to dirname $0 above !)
    # rwfile
    # cwfile
    # nullspace

    for v in wdir matrix m n prime interval splitwidth cmd rwfile cwfile nullspace ; do
        if ! [ "${!v}" ] ; then echo "Missing parameter \$$v" >&2 ; fi
    done

    # nrows and ncols are also needed, but can be deduced from rwfile and
    # cwfile easily.
    : ${ncols=$((`wc -c < $cwfile` / 4))}
    : ${nrows=$((`wc -c < $rwfile` / 4))}


    mdir=$wdir

    if ! [[ "$mpi,$thr" =~ ^([0-9]*)x([0-9]*),([0-9]*)x([0-9]*)$ ]] ; then
        echo "bad format for mpi=$mpi and thr=$thr" >&2
        exit 1
    fi

    Nh=$((${BASH_REMATCH[1]}*${BASH_REMATCH[3]}))
    Nv=$((${BASH_REMATCH[2]}*${BASH_REMATCH[4]}))

    bfile="`basename $matrix .bin`.${Nh}x${Nv}/`basename $matrix .bin`.${Nh}x${Nv}.bin"

    # This is for the **unbalanced** matrix !!

    if [ "$prime" = 2 ] ; then
        magmaprintmode=vector
    else
        magmaprintmode=spvector64
    fi

    : ${nrhs:=0}
}
# }}}
magma_print_main_parameters() { # {{{
    echo "m:=$m;n:=$n;interval:=$interval;"
    echo "nrhs:=$nrhs;"
} # }}}

magma_save_matrix() { # {{{
    echo "Saving matrix to magma format"
    $cmd weights < $rwfile > $mdir/rw.m
    $cmd weights < $cwfile > $mdir/cw.m
    if [ "$prime" = 2 ] ; then
        (echo "nc_orig:=$ncols;nr_orig:=$nrows;" ; $cmd bmatrix < $matrix) > $mdir/t.m
    else
        $cmd bpmatrix_${nrows}_${ncols} < $matrix > $mdir/t.m
    fi
    if ! [ -f "$wdir/$bfile" ] ; then
        echo "no balancing file found $wdir/$bfile" >&2
        exit $rc
    fi
    $cmd balancing < "$wdir/$bfile" > $mdir/b.m
    checksum=$(perl -ne '/checksum (\w+)/ && print "$1\n";' < $mdir/b.m)
    echo "ncpad:=nc;nrpad:=nr;nc:=nc_orig;nr:=nr_orig;" >> $mdir/b.m
    if [ "$prime" = 2 ] ; then
        echo "VS:=KMatrixSpace(GF(2), 64, nr_orig);"
    else
        echo "VS:=VectorSpace(GF(p), nr_orig);" 
    fi > $mdir/vectorspace.m

    placemats() {
        if [ "${nullspace_lowercase}" = left ] ; then
            transpose_if_left="Transpose"
        else
            transpose_if_left=""
        fi
        cat <<-EOF
            p:=$prime;
            nullspace:="${nullspace_lowercase}";
            xtr:=func<x|$transpose_if_left(x)>;
            M:=Matrix(GF(p),Matrix (M));
            nr:=Nrows(M);
            nc:=Ncols(M);
            nh:=$Nh;
            nv:=$Nv;
            ALIGNMENT_ON_ALL_BWC_VECTORS:=64;
            MINIMUM_ITEMS_IN_BWC_CHUNKS:=4;
            chunk:=ALIGNMENT_ON_ALL_BWC_VECTORS div MINIMUM_ITEMS_IN_BWC_CHUNKS;
            nr:=nh*nv*(chunk*Ceiling(x/chunk)) where x is Ceiling(Maximum(nr, nc)/(nh*nv));
            nc:=nr;
            x:=Matrix(GF(p),nr,nc,[]);InsertBlock(~x,M,1,1);M:=x;
            nrp:=nv*(chunk*Ceiling(x/chunk)) where x is Ceiling (nr/(nh*nv));
            ncp:=nh*(chunk*Ceiling(x/chunk)) where x is Ceiling (nc/(nh*nv));
            Mt:=Matrix(GF(p),nh*nrp,nv*ncp,[]);
EOF
        for i in `seq 0 $((Nh-1))` ; do
            cat <<-EOF
                nr$i:=nrp; /* nr div nh + ($i lt nr mod nh select 1 else 0); */
                snr$i:=$i*nrp; /* $i*(nr div nh) + Min($i, nr mod nh); */
EOF
            for j in `seq 0 $((Nv-1))` ; do
                if [ "$prime" = 2 ] ; then
                    $cmd bmatrix < $wdir/${bfile%%.bin}.h$i.v$j.bin > $mdir/t$i$j.m
                else
                    $cmd bpmatrix < $wdir/${bfile%%.bin}.h$i.v$j.bin > $mdir/t$i$j.m
                fi
                cat <<-EOF
                    nc$j:=ncp; /*  div nv + ($j lt nc mod nv select 1 else 0); */
                    snc$j:=$j*ncp; /* (nc div nv) + Min($j, nc mod nv); */
                    load "$mdir/t$i$j.m";
                    M$i$j:=Matrix(GF(p),Matrix(var));
                    x:=RMatrixSpace(GF(p),nr$i,nc$j)!0;
                    InsertBlock(~x,$transpose_if_left(M$i$j),1,1);
                    M$i$j:=x;
                    InsertBlock(~Mt,M$i$j,1+snr$i,1+snc$j);
EOF
            done
        done
        echo "mlist:=["
        for i in `seq 0 $((Nh-1))` ; do
            if [ "$i" != 0 ] ; then echo "," ; fi
            echo -n "["
            for j in `seq 0 $((Nv-1))` ; do
                if [ "$j" != 0 ] ; then echo -n ", " ; fi
                echo -n "M$i$j";
            done
            echo -n "]"
        done
        echo "];"
    }

    placemats > $mdir/placemats.m
}
# }}}

magma_save_all_vectors() { # {{{
    echo "Saving vectors to magma format"
    $cmd x $wdir/X > $mdir/x.m

    # prints to stdout magma code defining all vectors V whose iteration
    # number is equal to the first specified parameter.
    # The second parameter tells which is the name of the magma variable
    # where this is put.
    print_all_v() {
        echo "$2:=[VS|];"
        for j in `seq 0 $splitwidth $((n-1))` ; do
            let j1=j+splitwidth
            if ! [ -f $wdir/V${j}-${j1}.$1 ] ; then
                break;
            fi
            $cmd $magmaprintmode < $wdir/V${j}-${j1}.$1 > $mdir/V${j}-${j1}.$1.m
            cat <<-EOF
            load "$mdir/V${j}-${j1}.${1}.m";
            Append(~$2, g(var));
EOF
        done
    }

    # When magma debug has been activated, we normally have the first 10
    # iterates of each sequence.
    echo "VV:=[];" > $mdir/VV.m
    for i in {0..10} ; do
        ii=$((i*interval))
        print_all_v $ii V_$ii > $mdir/V.$ii.m
        if ! [ -f $mdir/V0-$splitwidth.$ii.m ] ; then
            break;
        fi
        cat >> $mdir/VV.m <<-EOF
        load "$mdir/V.$ii.m";
        Append(~VV, V_$ii);
EOF
    done
}
# }}}
magma_save_check_vectors() { # {{{
    echo "Saving check vectors to magma format"
    if [ -f "$wdir/Cv0-$splitwidth.0" ] ; then
        $cmd $magmaprintmode < $wdir/Cv0-$splitwidth.0 > $mdir/Cv0.m
        $cmd $magmaprintmode < $wdir/Cv0-$splitwidth.$interval > $mdir/Cvi.m
    else
        echo "var:=0;" > $mdir/Cv0.m
        echo "var:=0;" > $mdir/Cvi.m
    fi
    if [ -f "$wdir/Cd0-$splitwidth.$interval" ] ; then
        $cmd $magmaprintmode < $wdir/Cd0-$splitwidth.$interval > $mdir/Cdi.m
    else
        echo "var:=0;" > $mdir/Cdi.m
    fi
    if [ -f "$wdir/Cr0-$splitwidth.0-$splitwidth" ] ; then
        $cmd $magmaprintmode < $wdir/Cr0-$splitwidth.0-$splitwidth > $mdir/Cr.m
    else
        echo "var:=0;" > $mdir/Cr.m
    fi
    if [ -f "$wdir/Ct0-$splitwidth.0-$m" ] ; then
        $cmd $magmaprintmode < $wdir/Ct0-$splitwidth.0-$m > $mdir/Ct.m
    else
        echo "var:=0;" > $mdir/Ct.m
    fi
}
# }}}

magma_save_krylov_sequence() { # {{{
    echo "Saving krylov sequence to magma format"
    # in case $(basename $wdir) matches the format string, we want to also
    # check for being a regular file.
    afile=(`find $wdir -name 'A*[0-9]' -a -type f`)
    if [ "${#afile[@]}" != 1 ] ; then
        echo "########## FAILURE ! Krylov sequence not finished #############"
        find $wdir -name 'A*[0-9]' -a -type f | xargs -r -n 1 echo "### Found file "
    else
        afile=$(basename "${afile[0]}")
        $cmd $magmaprintmode < $wdir/$afile > $mdir/A.m
    fi
}
    # }}}
magma_print_all_Ffiles() { # {{{
#    $cmd spvector64 < $wdir/$afile.gen > $mdir/F.m
#    if [ -f $wdir/$afile.gen.rhs ] ; then
#        $cmd spvector64 < $wdir/$afile.gen.rhs > $mdir/rhscoeffs.m
#    else
#        # We must create an empty file, because otherwise magma won't accept
#        # "load" being conditional...
#        echo > $mdir/rhscoeffs.m
#    fi
    # Ffiles=(`ls $wdir | perl -ne '/^F.sols\d+-\d+.\d+-\d+$/ && print;' | sort -n`)
    echo "Saving linear generator to magma format" >&${text_output_fd}
    echo "Fchunks:=[];";
    echo "Rchunks:=[];";
    # Solutions correspond to *columns* of F, so let us be consistent
    # and use that everywhere, including in the format of the data we
    # pass to magma for verification.
    for j in `seq 0 $splitwidth $((n-1))` ; do
        let j1=j+splitwidth
        echo "Fj:=[];"
        echo "Rj:=[];"
        for s in `seq 0 $splitwidth $((n-1))` ; do
            let s1=s+splitwidth
            basename=F.sols${s}-${s1}.${j}-${j1}
            if [ -f "$wdir/$basename" ] ; then
            $cmd $magmaprintmode < $wdir/$basename > $mdir/$basename.m
            cat <<-EOF
                load "$mdir/$basename.m";
                Append(~Fj, g(var));
EOF
            if [ "$nrhs" ] ; then
                if [ "$j" -lt "$nrhs" ] ; then
                    basename=F.sols${s}-${s1}.${j}-${j1}.rhs
                    $cmd $magmaprintmode < $wdir/$basename > $mdir/$basename.m
                    cat <<-EOF
                        load "$mdir/$basename.m";
                        Append(~Rj, g(var));
EOF
                fi
            fi
            fi
        done
        echo "Append(~Fchunks, Fj);"
        if [ "$j" -lt "$nrhs" ] ; then
            echo "Append(~Rchunks, Rj);"
        fi
    done
}
# }}}
magma_save_mksol_data() { # {{{
    echo "Saving mksol data to magma format" 

    echo "vars:=[];" > $mdir/S.m
    for i in `seq 0 $splitwidth $((n-1))` ; do
        k=0;
        echo "solvars:=[];" >> $mdir/S.m
        while true ; do
            k0=$k
            k1=$((k+interval))
            let i1=i+splitwidth
            binfile="S.sols${i}-${i1}.${k0}-${k1}"
            magmafile="S${i}.${k0}-${k1}.m"
            if ! [ -f "$wdir/$binfile" ] ; then
                break
            fi
            $cmd $magmaprintmode < "$wdir/$binfile" > "$mdir/$magmafile"
            echo "load \"$mdir/$magmafile\"; Append(~solvars, var);" >> "$mdir/S.m"
            let k=k1
        done
        if [ $k = 0 ] ; then
            break
        else
            echo "Append(~vars, solvars);" >>  $mdir/S.m
        fi
    done
}
# }}}
magma_save_gather_data() { # {{{
    echo "Saving gather data to magma format"
    echo "vars:=[];" > $mdir/K.m
    (cd $wdir ; find  -name K.sols\*[0-9]) | while read f ; do
        $cmd $magmaprintmode < $wdir/$f > $mdir/$f.m
        echo "load \"$mdir/$f.m\"; Append(~vars, var);" >> $mdir/K.m
    done
}
# }}}
magma_run_script() { # {{{
    echo "Running magma verification script"

    s="`dirname $(readlink -f $0)`/bwc-ptrace.m"
    cd $wdir
    # magma does not exit with a useful return code, so we have to grep its
    # output.
    $magma -b $s < /dev/null | tee magma.out | grep -v '^Loading'
    if egrep -q "(Runtime error|Assertion failed)" magma.out ; then
        exit 1
    fi
}
# }}}

if [ "$magma" ] ; then
    echo "Converting files to magma format"

    cmd=`dirname $0`/convert_magma.pl

    exec {text_output_fd}>&1
    magma_sage_check_parameters
    magma_print_main_parameters > $mdir/mn.m
    magma_save_matrix
    magma_save_all_vectors
    magma_save_check_vectors
    magma_save_krylov_sequence
    magma_print_all_Ffiles > $mdir/Fchunks.m
    magma_save_mksol_data
    magma_save_gather_data

    #if [ "$prime" = 2 ] ; then
    #    convert_and_display() {
    #        text="$1"
    #        pattern="$2"
    #        terse="$3"
    #        n=($(find $wdir -maxdepth 1 -name "$pattern"'*' -a \! -name '*.m' -printf '%f\n'))
    #        echo -n "Saving $text to magma format ; ${#n[@]} files"
    #        if ! [ "$terse" ] && [ "${#n[@]}" -lt 10 ] && [ "${#n[@]}" -gt 0 ] ; then echo -n ": ${n[@]}" ; fi
    #        echo
    #        for f in "${n[@]}" ; do $cmd vector < $wdir/$f > $mdir/`basename $f`.m ; done
    #    }
    #
    #    convert_and_display "check vectors" "H"
    #    # convert_and_display "krylov sequence" "[AY]"
    #    # convert_and_display "linear generator" "F"
    #    # convert_and_display "mksol sequence" "S"
    #    convert_and_display "gather data" "[WK]"
    #
    #else
    #    :
    #fi

    check_script_diagnostic_fd=1
    if [ "$FORCE_BWC_EXTERNAL_CHECKS_OUTPUT_ON_FD3" ] && (exec 1>&3) 2>&- ; then
        check_script_diagnostic_fd=3
    fi
    if [ "$CADO_DEBUG" ] ; then set -x ; fi
    magma_run_script >&${check_script_diagnostic_fd}
    eval $old_setx
fi

if [ "$sage" ] ; then
    cmd=/bin/true
    magma_sage_check_parameters
    cd `dirname "$0"`/../..
    export PYTHONUNBUFFERED=true
    # we want cado_sage to be accessible. Note that as we customarily run
    # sage within a docker container, it's best if this PYTHONPATH is a
    # relative subdirectory.
    export PYTHONPATH=sagemath
    check_script_diagnostic_fd=1
    if [ "$FORCE_BWC_EXTERNAL_CHECKS_OUTPUT_ON_FD3" ] && (exec 1>&3) 2>&- ; then
        check_script_diagnostic_fd=3
    fi
    sage_args=( m=$m n=$n p=$prime
                wdir=$wdir matrix=$matrix
                nh=$Nh nv=$Nv
                nullspace=$nullspace
    )
    if [ "$wordsize" != 64 ] ; then
        sage_args+=(wordsize=$wordsize)
    fi
    if [ "$rhsfile" ] ; then
        sage_args+=(rhs=$rhsfile)
    fi
    if [ "$CADO_DEBUG" ] ; then set -x ; fi
    set -eo pipefail
    "$sage" linalg/bwc/bwc.sage "${sage_args[@]}" >&${check_script_diagnostic_fd}
    eval $old_setx
fi
