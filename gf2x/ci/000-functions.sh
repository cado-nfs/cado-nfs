# This file must be sourced.
#
# We're /bin/sh, not bash.
#
# We **must** be silent, because we're read by 00-dockerfile.sh. Note
# though that 00-dockerfile.sh sets HUSH_STDOUT, so that ECHO_E is
# actually a no-op (i.e. ECHO_E is ok to use safely, just don't use plain
# echo)

set -e

is_debian() { [ -f /etc/debian_version ] ; }
is_fedora() { [ -f /etc/fedora-release ] ; }
is_alpine() { [ -f /etc/alpine-release ] ; }
is_opensuse() { type -p zypper >/dev/null 2>&1 ; }
is_osx() { case "`uname -s`" in Darwin) true;; *)false;; esac; }
is_freebsd() {
    if [ -f /var/run/os-release ] ; then
        . "/var/run/os-release"
        [ "$NAME" = "FreeBSD" ]
    else
        false
    fi
}

CSI_RED="\e[01;31m"
CSI_BLUE="\e[01;34m"
CSI_RESET="\e[00;39m\e[m"
CSI_KILLLINE="\e[0K"

ECHO_E=echo
if [ "$BASH_VERSION" ] ; then
    ECHO_E="echo -e"
    if is_osx ; then
        # bash3 on osx does not like \e
        CSI_RED="[01;31m"
        CSI_BLUE="[01;34m"
        CSI_RESET="[00;39m[m"
        CSI_KILLLINE="[0K"
    fi
elif [ -f /proc/$$/exe ] && [ `readlink /proc/$$/exe` = /bin/busybox ] ; then
    ECHO_E="echo -e"
elif is_freebsd ; then
    # /bin/sh on freebsd does grok echo -e, but we have seemingly no way
    # to check /bin/sh's version. Presumably it's attached to the system
    # as a whole...
    ECHO_E="echo -e"
fi

if [ "$HUSH_STDOUT" ] ; then
    ECHO_E=:
fi

major_message()
{
    $ECHO_E "${CSI_BLUE}$*${CSI_RESET}"
}

pushed_sections=""

# Usage: enter_section [internal name] [message]
#
# The message is optional and defaults to the internal name
enter_section() {
    internal_name="$1"
    shift
    message="$*"
    : ${message:="$internal_name"}
    set -- "$internal_name" $pushed_sections
    pushed_sections="$*"
    $ECHO_E "section_start:`date +%s`:$internal_name\r${CSI_KILLLINE}${CSI_BLUE}$message${CSI_RESET}"
}

# Usage: leave_section [internal name] [message]
#
# Both arguments are optional. The message defaults to nothing. The
# internal name defaults to the last pushed section. If an inconsistency
# is detected, error out.
leave_section() {
    if ! [ "$pushed_sections" ] ; then
        echo "script error, no section stack !" >&2
        exit 1
    fi
    set -- $pushed_sections
    current_section="$1"
    shift
    pushed_sections="$*"
    $ECHO_E "section_end:`date +%s`:$current_section\r${CSI_KILLLINE}${CSI_BLUE}$2${CSI_RESET}"
}

# succeed if **ALL** of the listed tools exist
check_mandatory_tools() {
    fail=
    for tool in "$@" ; do
        if ! type -p "$tool" > /dev/null ; then
            echo "Missing tool: $tool" >&2
            fail=1
        else
            echo "ok - $tool"
        fi
    done
    if [ "$fail" ] ; then
        $ECHO_E "${CSI_RED}Fix these missing tools on the runner host, and try again${CSI_RESET}" >&2
        exit 1
    fi
}

# succeed if any of the provided files exists
check_mandatory_file() {
    for file in "$@" ; do
        if [ -f "$file" ] ; then
            echo "ok - $file"
            return
        fi
    done
    $ECHO_E "${CSI_RED}Fix these missing files on the runner host, and try again${CSI_RESET}" >&2
    exit 1
}

# succeed anyway, but report if none of the provided files exists
check_optional_file() {
    for file in "$@" ; do
        if [ -f "$file" ] ; then
            echo "ok - $file"
            return
        fi
    done
    $ECHO_E "${CSI_RED}Some optional files could not be found. This is not a fatal error${CSI_RESET}" >&2
}

check_optional_nonzero_output_shell() {
    out="$(eval "$@" 2>/dev/null || :)"
    if [ "$out" ] ; then
        echo "ok - $* [$out]"
        return
    fi
    echo "shell test failed: $*" >&2
    $ECHO_E "${CSI_RED}Some optional files could not be found. This is not a fatal error${CSI_RESET}" >&2
}

check_mandatory_nonzero_output_shell() {
    out="$(eval "$@" 2>/dev/null || :)"
    if [ "$out" ] ; then
        echo "ok - $* [$out]"
        return
    fi
    echo "shell test failed: $*" >&2
    $ECHO_E "${CSI_RED}Fix these missing tests on the runner host, and try again${CSI_RESET}" >&2
    exit 1
}

