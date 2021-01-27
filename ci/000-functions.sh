# This file must be sourced.

section_stack=()

set -e

CSI_RED="\e[01;31m"
CSI_BLUE="\e[01;34m"
CSI_RESET="\e[00;39m\e[m"
CSI_KILLLINE="\e[0K"

# Usage: enter_section [internal name] [message]
#
# The message is optional and defaults to the internal name
enter_section() {
    internal_name="$1"
    message="$2"
    : ${message:="$1"}
    section_stack+=("$1")
    echo -e "section_start:`date +%s`:$internal_name\r${CSI_KILLLINE}${CSI_BLUE}$message${CSI_RESET}"
}

# Usage: leave_section [internal name] [message]
#
# Both arguments are optional. The message defaults to nothing. The
# internal name defaults to the last pushed section. If an inconsistency
# is detected, error out.
leave_section() {
    if ! [ "${#section_stack[@]}" -gt 0 ] ; then
        echo "script error, no section stack !" >&2
        exit 1
    fi
    last_index=$((${#section_stack[@]}-1))
    last_section="${section_stack[$last_index]}"
    section_stack="${section_stack[@]::$last_index}"

    if [ "$internal_name" ] && [ "$last_section" != "$internal_name" ] ; then
        echo "script error, last pushed section is $last_section, not $internal_name !" >&2
        exit 1
    fi
    internal_name="$last_section"
    message="$2"
    echo -e "section_end:`date +%s`:$internal_name\r${CSI_KILLLINE}${CSI_BLUE}$message${CSI_RESET}"
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
        echo -e "${CSI_RED}Fix these missing tools on the runner host, and try again${CSI_RESET}" >&2
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
        echo -e "${CSI_RED}Fix these missing files on the runner host, and try again${CSI_RESET}" >&2
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
        echo -e "${CSI_RED}Some optional files could not be found. This is not a fatal error${CSI_RESET}" >&2
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
        echo -e "${CSI_RED}Some optional files could not be found. This is not a fatal error${CSI_RESET}" >&2
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
        echo -e "${CSI_RED}Fix these missing tests on the runner host, and try again${CSI_RESET}" >&2
    fi
}

