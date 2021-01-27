#!/bin/bash

# This is only used by .docker-script.sh


uid=$(stat --printf '%u' /host)
gid=$(stat --printf '%g' /host)
groupadd -g $gid hostgroup
useradd -s /bin/bash -m -g $gid -u $uid hostuser

## touch /tmp/trampoline.sh
## chmod 755 /tmp/trampoline.sh
echo "hostuser	ALL=(ALL:ALL) NOPASSWD:ALL" > /etc/sudoers.d/hostuser

cd /host

## if ! [ -f build ] && ! [ -e build ] ; then
    ## echo "export force_build_tree=/tmp/b" >> /tmp/trampoline.sh
## fi

# if [ -f /etc/fedora-release ] ; then
    # echo "cd /host" >> /tmp/trampoline.sh
    # echo "exec bash -i" >> /tmp/trampoline.sh
# elif [ -f /etc/debian_version ] ; then
    # # echo "cd /host" >> /tmp/trampoline.sh
    # echo "exec bash -i" >> /tmp/trampoline.sh
    # echo "exec bash" >> /tmp/trampoline.sh
# fi

# if [ -f /etc/fedora-release ] ; then
# set -x
# exec su -P -l hostuser -c /tmp/trampoline.sh
# else
# set -x
# # exec su -l hostuser -c /tmp/trampoline.sh
# fi
# 
# cat <<EOF
# # NOTE: we're not sure how to give you a shell that Just Works (tm).
# EOF
# if [ -s /tmp/trampoline.sh ] ; then
# cat <<EOF
# # You might need to run the commands from /tmp/trampoline.sh:
# `cat /tmp/trampoline.sh`
# #
# EOF
# fi

cat > /etc/profile.d/99-cado-nfs-build-tree.sh <<'EOF'
echo
export force_build_tree=/tmp/b 
echo "# NOTE: cado-nfs build tree has just been set to $force_build_tree"
echo "# NOTE: cado-nfs source tree is under /host"
EOF

exec su - hostuser
