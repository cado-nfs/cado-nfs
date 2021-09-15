# This must be sourced (because of traps), and we're bash. 
start_disksize_watch() {
    read -d -r perl_code <<'EOF'
#!/usr/bin/perl

my $inv_avail_ratio = 2;        # first warning at 50%
while (1) {
    open F, "df -k $ARGV[0] |";
    my $a;
    my $df;
    while (defined($_=<F>)) {
        next unless /^\S+\s+(\d+)\s+(\d+)\s+/;
        $a = 1/(1-$2/$1);
        $df = $_;
    }
    close F;
    if (defined($a) && $a > $inv_avail_ratio) {
        printf STDERR "Filesystem is %.2f%% full\n", 100.0*(1-1/$a);
        $inv_avail_ratio *= 2;
    }
    sleep(1);
}

EOF
    perl -e "$perl_code" ${1:-.} &
    watch_process=$!
    major_message "Started asynchronous disk usage watchdog process (pid $watch_process)"
    trap_add "kill $watch_process || :" EXIT ERR
}

start_disksize_watch .
