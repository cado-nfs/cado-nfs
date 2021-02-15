#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

disksize_watch() {
    read -d -r perl_code <<'EOF'
#!/usr/bin/perl

my $inv_avail_ratio = 2;        # first warning at 50%
while (1) {
    open F, "df -k . |";
    my $a;
    my $df;
    while (defined($_=<F>)) {
        next unless /^\S+\s+(\d+)\s+(\d+)\s+/;
        $a = 1/(1-$2/$1);
        $df = $_;
    }
    close F;
    if (defined($a) && $a > $inv_avail_ratio) {
        print STDERR "Filesystem is %.2f%% full\n", 100.0*(1-1/$a);
        $inv_avail_ratio *= 2;
    }
    sleep(1);
}

EOF
    perl -e "$perl_code" &
}

disksize_watch &
watch_process=$!
trap "kill $watch_process" EXIT ERR

"$(dirname $0)"/01-conf.sh
"$(dirname $0)"/02-build1.sh
"$(dirname $0)"/02-build2.sh
"$(dirname $0)"/03-check.sh

kill $watch_process
