#!/usr/bin/env bash

cd /host

cat > $HOME/.bashrc <<'EOF'
export force_build_tree=/tmp/b 
echo "# NOTE: cado-nfs build tree has just been set to $force_build_tree"
echo "# NOTE: cado-nfs source tree is under /host"
EOF

# Do we want hooks in .bashrc or .bash_profile ? alpine linux doesn't
# seem to honour .bashrc for some reason.

cat >> $HOME/.bash_profile <<EOF
cd /host
EOF

if [ "$BUILD_NAME" ] ; then
cat >> $HOME/.bash_profile <<EOF
BUILD_NAME="$BUILD_NAME"
export BUILD_NAME
DISPLAY_CONFIG=1
export DISPLAY_CONFIG
. ci/000-functions.sh
. ci/001-environment.sh
. ci/005-build-environment.sh
set +e
EOF
fi

exec /usr/bin/env bash --login
