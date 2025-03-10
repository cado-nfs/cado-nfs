# Montgomery square root.

## required inputs:

 * A file `BLAH.purged` that contains the (`npurged`) $(a,b)$ pairs. This
   file typically also contains the relations, but we don't need them.
   This is used in order to compute the initial logarithmic embeddings.
   (producer: `purge`)

 * A `BLAH.sparse.bin` matrix of (`small_nrows`) rows and (`small_ncols`)
   columns, which is the one that was fed to the linear algebra step.
   (producer: `replay-dl`)

 * A `BLAH.index` matrix of (`small_nrows`) rows and (`npurged`) cols,
   which indicates the contents of each relation-set as it was produced
   by `merge-dl`. (producer: `replay-dl`)

 * A file `BLAH.ideals` that gives the correspondence between the column
   number in the matrix and the ideal number that is found in the
   relations. (note that `replay-dl`, which produces this file, sorts its
   output by decreasing column weight, hence the odd-looking contents
   of this file). (producer: `replay-dl`)

 * A file `BLAH.renumber` file that can be used to interpret the ideal
   numbers in relation into something that is mathematically meaningful.
   (producer: `freerel`)

 * One or (often) several files `BLAH.ker` which are each a single vector
   of the left kernel of the matrix, thus a list of `small_nrows`
   integers mod $\ell$.  (producer: `bwc.pl`, under names such as
   `K.sols0-1.0.txt`)


# Creating a data set

To create a data set for this program, do as follows.
In all steps below, we assume that the following env variables are set:
 - `$N` is our modulus.
 - `$e` is our exponent.
 - `$wdir` is the directory with the computation files.
 - `$name` is the short prefix of all computation files.
 - `$build_tree` is the path to the cado-nfs binaries.

Examples:

```
N=123412341234618726123412341234611
e=29 ; # we don't need it until the linear algebra.
wdir=/tmp/blah
name=c30
eval `make show` ; # build_tree=/tmp/cado
```

### start a factoring job with the requested modulus.

This will, in particular, choose a polynomial, and obviously collect
relations.  It's about all there is that we will keep. In fact, it's even
possible to run the computation with `tasks.filter.run=False`

```
./cado-nfs.py $N tasks.filter.run=False --wdir $wdir
```

### Redo the freerel step

We do so in order to obtain a look-up table for ideals that includes
_all_ ideals, even those that we usually skip in the factoring context.
The important bit in the freerel invocation is the `-dl` flag. The first
line in the example below just picks the lpb parameters values from the
parameter files.

```
eval $(perl -ne '/(lpb[01])\s*=\s*(\d+)/ && print "$1=$2\n";' $wdir/$name.parameters_snapshot.*)
$build_tree/sieve/freerel -poly $wdir/$name.poly -renumber $wdir/$name.renumber-dl.gz -lpb0 $lpb0 -lpb1 $lpb1 -out $wdir/$name.freerel-dl.gz -dl
```

[FIXME: is it really necessary? The renumber file is almost identical,
except for one comment line (but maybe it does get parsed). And the
freerel file does indeed show a shifted numbering, but are we going to
use the free relations anyway?

For the record, it's not necessary for the programs, but humans do like
to read the correspondence between ideal indices and their definition.
It's given by the following command:

```
$build_tree/misc/debug_renumber -poly $wdir/$name.poly  -lpbs $lpb0,$lpb1 -build
```

### The `dup1` and `dup2` steps.

In this simple experiment, `dup1` only greps the non-comment lines in the
relation files, but we keep it in the workflow for consistency.

```
RELATION_FILES=(`find $wdir/$name.upload -name "$name.[0-9]*"`)
mkdir -p $wdir/$name.dup1/0 ; $build_tree/filter/dup1 -prefix $name.dup1 -out $wdir/$name.dup1 -n 0 "${RELATION_FILES[@]}"
```

Then, `dup2` reads these relations and renumbers them:

```
nrels=$(zcat $wdir/$name.dup1/0/$name.dup1.0000.gz | wc -l)
$build_tree/filter/dup2 -poly $wdir/$name.poly -nrels $nrels -renumber $wdir/$name.renumber-dl.gz -dl $wdir/$name.dup1/0/$name.dup1.0000.gz
```

Again, we have some debug tools. Consistency of the relations can be
checked as follows:

```
zcat $wdir/$name.dup1/0/$name.dup1.0000.gz  | head | $build_tree/misc/explain_indexed_relation -renumber $wdir/$name.renumber-dl.gz -poly $wdir/$name.poly -dl -relations /dev/stdin > /tmp/ck.sage
sage /tmp/ck.sage

$build_tree/misc/explain_indexed_relation -renumber $wdir/$name.renumber-dl.gz -poly $wdir/$name.poly -dl -python -all -skip-ideal-checks < /dev/null > $wdir/$name.debug-ideals.py
$build_tree/misc/explain_indexed_relation -renumber $wdir/$name.renumber-dl.gz -poly $wdir/$name.poly -dl -raw -all -skip-ideal-checks < /dev/null > $wdir/$name.debug-ideals.txt
```

### redo the `purge` / `merge-dl` / `replay-dl` steps

There's a stupid catch, which is that we must tell the purge program
ahead of time what is the size of the renumber table. This limitation
shouldn't exist, really, but for the time being we have to work around
it. Any program that reads the renumber table (purge itself doesn't!) can
get this information easily.

Note that for `merge-dl` and `replay-dl`, the dl version is necessary
because it's the one that keeps correct track of the exponents.

```
nideals=$($build_tree/misc/explain_indexed_relation -renumber $wdir/$name.renumber-dl.gz -poly $wdir/$name.poly -dl < /dev/null | perl -ne '/INFO.*size = (\d+)/ && print "$1\n";')

nrels_nodup=$(zcat $wdir/$name.dup1/0/$name.dup1.0000.gz | wc -l)
nfreerels=$(zcat $wdir/$name.freerel-dl.gz | wc -l)
nrels_input=$((nrels_nodup+nfreerels))

min_index=$((nideals/40))
if [ $min_index -lt 500 ] ; then min_index=500; fi

$build_tree/filter/purge -out $wdir/$name.purged.gz -nrels $nrels_input -keep 160 -col-min-index $min_index -col-max-index $nideals -t 8 -required_excess 0.0 $wdir/$name.dup1/0/$name.dup1.0000.gz $wdir/$name.freerel-dl.gz

$build_tree/filter/merge-dl -mat $wdir/$name.purged.gz -out $wdir/$name.history.gz -target_density 170.0 -t 8

$build_tree/filter/replay-dl -purged $wdir/$name.purged.gz -his $wdir/$name.history.gz -index $wdir/$name.index.gz -ideals $wdir/$name.ideals.gz -out $wdir/$name.sparse.bin
```

### Do the linear algebra step

This works for a baby example. A larger example probably requires some
tweaking, but it seems at least that we have something that sorta makes
sense, which is pretty much a good sign.

The blocking factor n must be at least as large as the number of
solutions that are requested, which in turn must be at least as large as
the unit rank.

Here we choose exponent `e=29`. This works for the example I tried, but
it might not work for yours. Most probably it would work for a larger
exponent such as 65537 with overwhelming probability. (but for testing,
small e is preferred, obviously)

```
e=29; ./build/coffee/linalg/bwc/bwc.pl :complete matrix=/tmp/blah/c30.sparse.bin prime=$e wdir=/tmp/blah/c30.bwc.$e nullspace=LEFT m=8 n=4 solutions=0-3
```

### Compute the Schirokauer maps

The sagemath code is capable of computing the Schirokauer maps, but as a
matter of fact the cado-nfs C programs can also do this. I'm just not
sure that the latter correctly deals with the case where the exponent
(called ell in the Schirokauer maps contexts, e in the e-th root
computation context) is within the factor base (it appears to be the
case, but more thorough checking would be needed). For sure the sagemath
code has provision to deal with this case. Anyway, to get the Schirokauer
maps from the C code, here is the command line.

Further analysis tells that if the Schirokauer maps that `sm_append`
computes use pieces that involve an ideal of degree 1, then there is a
chance that this ideal is encountered in the factorizations of
`$a-b*\alpha$`, and the output is clearly rubbish in that case. And the
computation fails, so we must fix this.

```
$build_tree/filter/sm_append -ell $e -in $wdir/$name.purged.gz -poly $wdir/$name.poly | gzip >  $wdir/$name.purged_withsm${e}.gz
```

# Debugging

To keep track of what the linear algebra step did, see `montgomery.sage`

