# This is only meant to be *sourced* from a shell

# Some helper shell functions are provided here. These may be handy when
# one wants to do some maintenance tasks with the cado-nfs workunit
# database while a large computation is running.

# WARNING: While the functions below were designed to perform rather
# safe operations, you must be aware of the fact that you're tinkering
# with a database. In essence, this is a dangerous operation.


# The code below assumes that:
#     - your cado-nfs server is running with a mysql back-end (making
#       that work with sqlite3 would be possible -- maybe a few syntax
#       tweaks).
#
#     - you have put the necessary database authentication stuff
#       in ~/.db_auth, namely:
#             DB_HOST=__put_the_hostname_of_your_mysql_server_here__
#             DB_USER=__put_your_mysql_username_here__
#             DB_PASSWORD=__put_your_mysql_password_here__
#             DB_NAME=__put_your_mysql_database_name_here__
#
#     - you have sourced the file ~/.db_auth and the resubmit-tools.sh
#       script, e.g. with the shell command:
#             . ~/.db_auth
#             . ~/resubmit-tools.sh
# 

do_mysql() {
    mysql -h $DB_HOST --user=$DB_USER --password=$DB_PASSWORD ${DB_NAME} -B -e "$1"
}
sql_list_allfiles() {
    table1="(select wurowid from workunits where status=5) x"
    do_mysql "select path from $table1 inner join files using (wurowid);"
}
dispatch() {
    perl -e '$a{$_}=1 for (@ARGV); die unless @ARGV; while (defined($_=<STDIN>)) { /(?:rsa200\.|rsa200_sieving_)(\d{3})/ && do { next unless $a{$1}; my $fh = $h{$1}; if(!defined($fh)) { open $fh, ">/tmp/good$1.txt"; $h{$1}=$fh; } print $fh $_; };} ' "$@"
}
extract_goodfiles() {
    sql_list_allfiles | dispatch "$@"
}


gather_wu_history_in_database() {
    . ~/.db_auth
    read -d \~ -r query <<'EOF'
    drop table wu_byid;
    create table wu_byid (wuid_base varchar(512), wurowid integer);
    create index wu_byid_baseid_index on wu_byid(wuid_base);
    insert into wu_byid 
        select case p when 0 then wuid else substring(wuid,1,p-1) end,wurowid from
        (select position('#' in wuid) as p, wuid, wurowid from workunits) x
        ;
~
EOF
    do_mysql "$query"
}

put_dead_wus_in_a_table() {
    read -d \~ -r query <<'EOF'
    drop table deadid;
    create table deadid (path varchar(512) PRIMARY KEY);
    load data local infile '/tmp/miss.txt' into table deadid;
~
EOF
    do_mysql "$query"
}
list_full_history_of_dead_wus() {
    read -d \~ -r query <<'EOF'
    select wurowid,wuid_base,wuid,timeassigned,timeresult,assignedclient,status
        from workunits inner join
        (select wurowid,wuid_base from wu_byid inner join
            (select path as wuid_base from deadid) y
            using (wuid_base)
        ) x
        using (wurowid)
        order by wuid_base,timecreated;
~
EOF
    do_mysql "$query"
}
list_terminal_state_of_dead_wus() {
    read -d \~ -r query <<'EOF'
    drop table temp;
    create table temp (wurowid integer primary key);
    insert into temp
        select max(wurowid)
        from (
            select wuid_base,wurowid
            from workunits inner join
            (select wuid_base,wurowid from wu_byid inner join deadid
            on wuid_base = path) tt
            using (wurowid)
        ) x
        group by wuid_base ;
    select wurowid,wuid,timeassigned,timeresult,assignedclient,status from workunits inner join temp using (wurowid) order by wuid;
~
EOF
    do_mysql "$query"
}

reschedule_dead_wus_with_status6()
{
    # We also need to "un-receive" these WUs.
    #
    # Use case: a WU wasn't resubmitted beause attempt>=maxwuerror and status=6
    #
    # XXX: this doesn't properly update wu_range_received
    read -d \~ -r query <<'EOF'
    update workunits inner join temp using (wurowid) set status=2 where status=6;
    update sieving set value=value-row_count() where kkey='wu_received';
~
EOF
    do_mysql "$query"
}
reschedule_dead_wus_with_status7()
{
    # Exactly as above. The use case is different. A priori those
    # correspond to WUs that "will time out soon" (because we know the
    # job's already dead), or that timed out and didn't respawn because
    # of maxresubmit.
    read -d \~ -r query <<'EOF'
    update workunits inner join temp using (wurowid) set status=2 where status=7;
    update sieving set value=value-row_count() where kkey='wu_received';
~
EOF
    do_mysql "$query"
}
reschedule_dead_wus_with_status1()
{
    # Those have never been received, so we should just schedule them for
    # timeout.
    read -d \~ -r query <<'EOF'
    update workunits inner join temp using (wurowid) set status=2 where status=1;
~
EOF
    do_mysql "$query"
}

typical_flow() {
    extract_goodfiles "$@"
    gather_wu_history_in_database &
    wait
    put_dead_wus_in_a_table
    list_terminal_state_of_dead_wus
    # maybe followed by reschedule_dead_wus_with_status6
}

