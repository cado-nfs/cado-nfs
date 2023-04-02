# This file must be sourced.
#
# We're /bin/sh, not bash.
#
# We **must** be silent, because we're read by 00-dockerfile.sh. Note
# though that 00-dockerfile.sh sets HUSH_STDOUT, so that ECHO_E is
# actually a no-op (i.e. ECHO_E or major_message are ok to use safely,
# just don't use plain echo)

export CLICOLOR_FORCE=1

# Note that our set of scripts reacts on BUILD_NAME, and the logic for
# this is here. BUILD_NAME must follow the regexp below.
#
# (note that we **CANNOT** do regexp matching in this script because
# we're /bin/sh, not bash !)
#
# (container for )?(coverage tests on )?((expensive )?checks)?( on (osx|alpine|(debian|fedora|centos|freebsd)[0-9]*) system)?( with ((32[- ]bit )?gcc|clang|icc))?
#
# with one exception which is "merge coverage tests"
#
# a downside is that in the gitlab page, this makes many pipeline steps
# with similar names. We could change to abbreviated names (e.g.
# "coverage tests on " would be V, "container for " would be L, "checks"
# and "expensive checks" would be C and XC, and so on. But it would be
# really cryptic. to have, e.g. LC/debian10+gcc ; wouldn't it ?

if [ "$DISPLAY_CONFIG" ] ; then
    display_config() { major_message "$@" ; }
else
    display_config() { : ; }
fi

if type -p hostname > /dev/null 2>&1 ; then
    HOSTNAME=$(hostname)
else
    HOSTNAME="[[placeholder]]"
fi
case "$HOSTNAME" in
    # some of our very slow machines have so little ram that clearly, we
    # must not tax them too much.
    genepi|calva|pine64) export NCPUS_FAKE=1;;
    *) : ;;
esac

#### set COMMIT_SHORT_SHA to CI_COMMIT_SHORT_SHA or GITHUB_SHA
if [ "$CI_COMMIT_SHORT_SHA" ] ; then
    COMMIT_SHORT_SHA="$CI_COMMIT_SHORT_SHA"
elif [ "$GITHUB_SHA" ] ; then
    COMMIT_SHORT_SHA="$GITHUB_SHA"
elif [ -d .git ] && type -p git > /dev/null 2>&1 ; then
    COMMIT_SHORT_SHA="$(git rev-parse --short HEAD)"
    display_config "Setting COMMIT_SHORT_SHA=\"$COMMIT_SHORT_SHA\""
fi
display_config "COMMIT_SHORT_SHA=$COMMIT_SHORT_SHA"

#### set JOB_ID to either CI_JOB_ID or GITHUB_RUN_ID
if [ "$CI_JOB_ID" ] ; then
    JOB_ID="$CI_JOB_ID"
elif [ "$GITHUB_RUN_ID" ] ; then
    JOB_ID="$GITHUB_RUN_ID"
else
    JOB_ID=0
    display_config "Setting JOB_ID=\"$JOB_ID\""
fi
display_config "JOB_ID=$JOB_ID"

### set REPOSITORY to either $CI_PROJECT_NAMESPACE/$CI_PROJECT_NAME or GITHUB_REPOSITORY

if [ "$CI_PROJECT_NAMESPACE" ] && [ "$CI_PROJECT_NAME" ] ; then
    REPOSITORY="$CI_PROJECT_NAMESPACE/$CI_PROJECT_NAME"
elif [ "$GITHUB_REPOSITORY" ] ; then
    REPOSITORY="$GITHUB_REPOSITORY"
else
    # no default
    REPOSITORY=
fi
display_config "REPOSITORY=$REPOSITORY"
    
if [ -x /opt/homebrew/bin/brew ] ; then
    eval `/opt/homebrew/bin/brew shellenv`
fi

MAKE=make
if type -p gmake > /dev/null 2>&1 ; then
    MAKE=gmake
fi

export CC CXX
export CFLAGS CXXFLAGS
export MAKE
export ENABLE_SHARED
