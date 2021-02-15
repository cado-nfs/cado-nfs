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
elif [ -f /proc/$$/exe ] && [ `readlink /proc/$$/exe` = /bin/busybox ] ; then
    ECHO_E="echo -e"
elif is_freebsd ; then
    # /bin/sh on freebsd does grok echo -e, but we have seemingly no way
    # to check /bin/sh's version. Presumably it's attached to the system
    # as a whole...
    ECHO_E="echo -e"
elif is_osx ; then
    # bash3 on osx does not like \e
    CSI_RED="[01;31m"
    CSI_BLUE="[01;34m"
    CSI_RESET="[00;39m\e[m"
    CSI_KILLLINE="[0K"
    ECHO_E="echo -e"
fi

if [ "$HUSH_STDOUT" ] ; then
    ECHO_E=:
fi

major_message()
{
    $ECHO_E "${CSI_BLUE}$*${CSI_RESET}"
}

# Usage: enter_section [internal name] [message]
#
# The message is optional and defaults to the internal name
enter_section() {
    internal_name="$1"
    shift
    message="$*"
    : ${message:="$internal_name"}
    current_section="$1"
    $ECHO_E "section_start:`date +%s`:$internal_name\r${CSI_KILLLINE}${CSI_BLUE}$message${CSI_RESET}"
}

# Usage: leave_section [internal name] [message]
#
# Both arguments are optional. The message defaults to nothing. The
# internal name defaults to the last pushed section. If an inconsistency
# is detected, error out.
leave_section() {
    if ! [ "$current_section" ] ; then
        echo "script error, no section stack !" >&2
        exit 1
    fi
    if [ "$1" ] && [ "$current_section" != "$1" ] ; then
        echo "script error, last pushed section is $current_section, not $1 !" >&2
        exit 1
    fi
    $ECHO_E "section_end:`date +%s`:$1\r${CSI_KILLLINE}${CSI_BLUE}$2${CSI_RESET}"
    unset current_section
}

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

check_mandatory_files() {
    fail=
    for file in "$@" ; do
        if ! [ -f "$file" ] ; then
            echo "Missing file: $file" >&2
            fail=1
        else
            echo "ok - $file"
        fi
    done
    if [ "$fail" ] ; then
        $ECHO_E "${CSI_RED}Fix these missing files on the runner host, and try again${CSI_RESET}" >&2
        exit 1
    fi
}

check_optional_files() {
    fail=
    for file in "$@" ; do
        if ! [ -f "$file" ] ; then
            echo "Optional file not found: $file" >&2
            fail=1
        else
            echo "ok - $file"
        fi
    done
    if [ "$fail" ] ; then
        $ECHO_E "${CSI_RED}Some optional files could not be found. This is not a fatal error${CSI_RESET}" >&2
    fi
}

check_optional_nonzero_output_shell() {
    fail=
    for cmd in "$@" ; do
        if ! eval "$cmd" | grep -q . ; then
            echo "shell test failed: $cmd" >&2
            fail=1
        else
            echo "ok - $cmd"
        fi
    done
    if [ "$fail" ] ; then
        $ECHO_E "${CSI_RED}Some optional files could not be found. This is not a fatal error${CSI_RESET}" >&2
    fi
}

check_mandatory_nonzero_output_shell() {
    fail=
    for command in "$@" ; do
        if ! eval "$command" | grep -q . ; then
            echo "shell test failed: $command" >&2
            fail=1
        else
            echo "ok - $command"
        fi
    done
    if [ "$fail" ] ; then
        $ECHO_E "${CSI_RED}Fix these missing tests on the runner host, and try again${CSI_RESET}" >&2
    fi
}

