#!/usr/bin/env python3

import argparse
import re
import json
import subprocess

headers_c_language = \
    set("cassert cctype cerrno cfloat cinttypes climits"
        " cmath csignal cstdarg cstdbool cstddef cstdint"
        " cstdio cstdlib cstring ctime complex.h".split())

# complex.h is tricky because ccomplex was removed from c++10

headers_c_language_not_converted = \
    set([c[1:] + ".h" for c in headers_c_language])

headers_cxx_language = \
    set("algorithm array atomic condition_variable deque exception fstream"
        " functional initializer_list iomanip ios iosfwd iostream istream"
        " iterator limits list map memory mutex new ostream queue set"
        " sstream stack stdexcept streambuf string thread tuple type_traits"
        " typeinfo unordered_map unordered_set utility vector".split())

# posix, compiler, and platform headers are usually few and can be presented
# grouped.
headers_posix = \
    set("dirent.h dlfcn.h strings.h fcntl.h pthread.h regex.h sys/mman.h"
        " sys/endian.h"
        " sys/resource.h sys/socket.h sys/stat.h sys/statvfs.h sys/syscall.h"
        " sys/time.h sys/types.h sys/utsname.h sys/wait.h unistd.h"
        " sys/user.h semaphore.h".split())

headers_compiler = \
    set("arm_neon.h emmintrin.h immintrin.h mmintrin.h nmmintrin.h"
        " pmmintrin.h smmintrin.h tmmintrin.h xmmintrin.h"
        " x86intrin.h"
        " omp.h".split())

headers_platform = \
    set("mach/mach.h execinfo.h cxxabi.h features.h"
        " linux/binfmts.h linux/limits.h".split())

headers_user_libraries = \
    set("fmt/format.h gmp.h hwloc.h hwloc/bitmap.h mpi.h".split())

angle_headers_book = dict(
    c_language=headers_c_language,
    c_language_not_converted=headers_c_language_not_converted,
    cxx_language=headers_cxx_language,
    compiler=headers_compiler,
    platform=headers_platform,
    posix=headers_posix,
    user_libraries=headers_user_libraries)


def recognize_angle_header(line):
    regexp = r"(?:// )?#include\s<([a-z0-9_/\.]+)>( *(\/\/.*|\/\*.*\*\/))?"
    if m := re.match(regexp, line):
        hdr = m.group(1)
        # comment = m.group(3)
        # print(f"Found <> header {hdr}")
        for k, s in angle_headers_book.items():
            if hdr in s:
                return [k, hdr, line]
        return
        # raise NotImplementedError(f"Cannot recognize {hdr}")


def recognize_user_header(line):
    regexp = r'(?:// )?#include\s"([a-zA-Z0-9_/\.-]+)"( *(\/\/.*|\/\*.*\*\/))?'
    if m := re.match(regexp, line):
        hdr = m.group(1)
        # comment = m.group(3)
        # print(f"Found \"\" header {hdr}")
        return ['user', hdr, line]


def detect_header_macro_protect(line, header_macros):
    m = re.match(r'#(?:ifndef|define)\s([A-Z0-9_]+)', line)
    if not m:
        return False
    return m.group(1) in header_macros


class parser_base(object):
    def __init__(self):
        pass

    def empty(self, T):
        n = 0
        while n < len(T):
            if T[n][0] in ['empty']:
                n += 1
            else:
                break
        return n

    def conditional(self, T, *sub_parsers):
        n = 0
        parsed = []
        if n >= len(T) or T[n][0] != 'if':
            return
        parsed.append(T[n])
        n += 1
        while n < len(T):
            n0 = n
            for s in sub_parsers:
                if m := getattr(self, s)(T[n:], False):
                    n += m[0]
                    parsed += m[1]
            if n < len(T) and T[n][0] == 'elif/else':
                parsed.append(T[n])
                n += 1
            if n == n0:
                break
        if n >= len(T) or T[n][0] != 'fi':
            return
        parsed.append(T[n])
        n += 1
        return n, parsed


class parser_common(parser_base):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def prologue(self, T, strip_trailing_empty_lines=True):
        n = 0
        parsed = []
        while n < len(T):
            if T[n][0] in ['empty', 'comment', 'header-protect']:
                parsed.append(T[n])
                n += 1
            else:
                break
        if strip_trailing_empty_lines:
            n += self.empty(T[n:])
        if not n:
            return
        return n, parsed

    def pragma_err_def(self, T, strip_trailing_empty_lines=True):
        n = 0
        parsed = []
        while n < len(T):
            if T[n][0] in ['empty', 'comment', 'pragma-error', 'define']:
                parsed.append(T[n])
                n += 1
            elif T[n][0] in ['empty']:
                n += 1
            elif m := self.conditional(T[n:], 'pragma_err_def'):
                parsed.append(['conditional', m[1]])
                n += m[0]
            else:
                break
        if strip_trailing_empty_lines:
            n += self.empty(T[n:])
        if not n:
            return
        return n, parsed

    def user_early(self, T, strip_trailing_empty_lines=True):
        n = 0
        parsed = []
        while n < len(T):
            if T[n][0] in ['empty', 'comment']:
                parsed.append(T[n])
                n += 1
            else:
                break
        while n < len(T):
            if T[n][0] in ['user-early']:
                parsed.append(T[n])
                n += 1
            elif T[n][0] == 'user' and re.search(r'config.h(?:pp)?', T[n][1]):
                # things such as las-config.h, for instance.
                T[n][0] = 'user-early'
                parsed.append(T[n])
                n += 1
            elif T[n][0:2] == ['user', 'macros.h']:
                # It's a hack, but it can happen if we want some specific
                # macros to decide on conditionals.
                T[n][0] = 'user-early'
                parsed.append(T[n])
                n += 1
            else:
                break
        if strip_trailing_empty_lines:
            n += self.empty(T[n:])
        if not n:
            return
        return n, parsed

    def c_language(self, T, strip_trailing_empty_lines=True):
        n = 0
        parsed = []
        while n < len(T):
            if T[n][0] in ['empty', 'comment']:
                parsed.append(T[n])
                n += 1
            else:
                break
        while n < len(T):
            if T[n][0] in ['c_language']:
                parsed.append(T[n])
                n += 1
            elif T[n][0] in ['c_language_not_converted']:
                parsed.append(T[n])
                n += 1
            elif T[n][0] in ['empty']:
                n += 1
            elif m := self.conditional(T[n:], 'c_language'):
                parsed.append(['conditional', m[1]])
                n += m[0]
            else:
                break
        if strip_trailing_empty_lines:
            n += self.empty(T[n:])
        if not n:
            return
        return n, parsed

    def cxx_language(self, T, strip_trailing_empty_lines=True):
        n = 0
        parsed = []
        while n < len(T):
            if T[n][0] in ['empty', 'comment']:
                parsed.append(T[n])
                n += 1
            else:
                break
        while n < len(T):
            if T[n][0] in ['cxx_language']:
                parsed.append(T[n])
                n += 1
            elif T[n][0] in ['empty']:
                n += 1
            elif m := self.conditional(T[n:], 'cxx_language'):
                parsed.append(['conditional', m[1]])
                n += m[0]
            else:
                break
        if strip_trailing_empty_lines:
            n += self.empty(T[n:])
        if not n:
            return
        return n, parsed

    def platform_etc(self, T, strip_trailing_empty_lines=True):
        n = 0
        parsed = []
        while n < len(T):
            if T[n][0] in ['empty', 'comment']:
                parsed.append(T[n])
                n += 1
            else:
                break
        while n < len(T):
            if T[n][0] in ['compiler', 'platform', 'posix']:
                parsed.append(T[n])
                n += 1
            elif T[n][0] in ['empty']:
                n += 1
            elif m := self.conditional(T[n:],
                                       'c_language', 'cxx_language',
                                       'platform_etc'):
                parsed.append(['conditional', m[1]])
                n += m[0]
            else:
                break
        if strip_trailing_empty_lines:
            n += self.empty(T[n:])
        if not n:
            return
        return n, parsed

    def user_libraries(self, T, strip_trailing_empty_lines=True):
        n = 0
        parsed = []
        while n < len(T):
            if T[n][0] in ['empty', 'comment']:
                parsed.append(T[n])
                n += 1
            else:
                break
        while n < len(T):
            if T[n][0] in ['user_libraries']:
                parsed.append(T[n])
                n += 1
            elif T[n][0] in ['empty']:
                n += 1
            elif m := self.conditional(T[n:], 'user_libraries'):
                parsed.append(['conditional', m[1]])
                n += m[0]
            else:
                break
        if not n:
            return
        if strip_trailing_empty_lines:
            n += self.empty(T[n:])
        return n, parsed

    def user(self, T, strip_trailing_empty_lines=True):
        n = 0
        parsed = []
        while n < len(T):
            if T[n][0] in ['empty', 'comment']:
                parsed.append(T[n])
                n += 1
            else:
                break
        while n < len(T):
            if T[n][0] in ['user']:
                parsed.append(T[n])
                n += 1
            elif T[n][0] in ['empty']:
                n += 1
            elif m := self.conditional(T[n:], 'user'):
                parsed.append(['conditional', m[1]])
                n += m[0]
            else:
                break
        if not n:
            return
        if strip_trailing_empty_lines:
            n += self.empty(T[n:])
        return n, parsed

    def unsorted(self, T, strip_trailing_empty_lines=True):
        n = 0
        parsed = []
        while n < len(T):
            if T[n][0] in ['empty', 'comment']:
                parsed.append(T[n])
                n += 1
            else:
                break
        while n < len(T):
            if T[n][0] in ['user', 'c_language', 'cxx_language',
                           'c_language_not_converted', 'user_libraries',
                           'platform', 'compiler', 'posix', 'user']:
                parsed.append(T[n])
                n += 1
            elif T[n][0] in ['empty']:
                n += 1
            elif m := self.conditional(T[n:], 'unsorted'):
                parsed.append(['conditional', m[1]])
                n += m[0]
            else:
                break
        if not n:
            return
        if strip_trailing_empty_lines:
            n += self.empty(T[n:])
        return n, parsed


class parser_header(parser_common):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __call__(self, T):
        n = 0
        parsed = []

        for phase in ['prologue',
                      'user_early',
                      'pragma_err_def',
                      'c_language', 'cxx_language',
                      'platform_etc',
                      'user_libraries',
                      'user', 'unsorted']:
            if m := getattr(self, phase)(T[n:]):
                parsed.append([phase, m[1]])
                if args.verbose:
                    print(f"After {phase} n={n}")
                n += m[0]
        if not n:
            return
        return n, parsed


class parser_compilation_unit(parser_common):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __call__(self, T):
        n = 0
        parsed = []

        for phase in ['prologue',
                      'user_early',
                      'pragma_err_def',
                      'c_language', 'cxx_language',
                      'platform_etc',
                      'user_libraries',
                      'user', 'unsorted',
                      'pragma_err_def',
                      ]:
            if m := getattr(self, phase)(T[n:]):
                parsed.append([phase, m[1]])
                if args.verbose:
                    print(f"After {phase} n={n}")
                n += m[0]
        if not n:
            return
        return n, parsed


def scan_file(args, filename):
    in_prologue = True
    in_comment = False
    structure = []
    is_header = re.search(r'\.h(pp)?$', filename)
    is_compilation_unit = re.search(r'\.c(pp)?$', filename)
    in_condition = 0
    pragma_magic_stop = False
    if not (is_header or is_compilation_unit):
        raise RuntimeError(f"What is {filename}")
    header_macros = []
    if is_header:
        c = 0
        uu = re.sub(r"[^a-zA-Z0-9]", "_", filename[c:]).upper()
        header_macros += [uu, uu + "_"]
        while (d := filename[c:].find('/')) >= 0:
            c += d + 1
            uu = re.sub(r"[^a-zA-Z0-9]", "_", filename[c:]).upper()
            header_macros += [uu, uu + "_"]
        header_macros += ["CADO_" + x for x in header_macros]
    if args.verbose:
        print(header_macros)
    with open(filename) as f:
        for i, line in enumerate(f):
            line = line.rstrip()
            if not in_prologue:
                if not pragma_magic_stop and re.search(r"^#include", line):
                    print(f"{filename}:{i+1}:"
                          " ABORT on unknown header:"
                          " " + line)
                    return
                structure.append(['rest', line])
            elif line == '// scan-headers: skip':
                print(f"{filename}:{i+1}: SKIP")
                return
            elif line == '// scan-headers: stop here':
                pragma_magic_stop = True
                structure.append(['rest', line])
                in_prologue = False
            elif line == '// scan-headers:':
                raise NotImplementedError("Unknown pragma: " + line)
            elif not in_comment and re.match(r'/\*.*\*/', line):
                structure.append(['comment', line])
            elif not in_comment and re.search(r'^/\*', line):
                structure.append(['comment', line])
                in_comment = True
            elif in_comment:
                if re.search(r'\*/$', line):
                    in_comment = False
                structure.append(['comment', line])
            elif not line:
                structure.append(['empty', line])
            elif detect_header_macro_protect(line, header_macros):
                structure.append(['header-protect', line])
            elif re.match(r'#define', line):
                structure.append(['define', line])
            elif re.match(r'#(pragma|error)', line):
                structure.append(['pragma-error', line])
            elif re.match(r'#if', line):
                structure.append(['if', line])
                in_condition += 1
            elif in_condition and re.match(r'#el', line):
                structure.append(['elif/else', line])
            elif in_condition and re.match(r'#endif', line):
                structure.append(['fi', line])
                in_condition -= 1
            elif m := recognize_angle_header(line):
                structure.append(m)
            elif m := recognize_user_header(line):
                if is_compilation_unit and m[1] == 'cado.h':
                    m[0] = 'user-early'
                elif is_header and m[1] == 'cado_config.h':
                    m[0] = 'user-early'
                structure.append(m)
            elif not in_comment and re.search(r'^//', line):
                structure.append(['comment', line])
            else:
                i = 1
                while i < len(structure):
                    if structure[-i][0] in ['comment',
                                            'empty',
                                            'if',
                                            'define'
                                            ]:
                        structure[-i][0] = 'rest'
                        i += 1
                    else:
                        break
                # print(f"Rest of {filename} begins with {line}")
                structure.append(['rest', line])
                in_prologue = False

    if args.verbose:
        print(filename, ' '.join([x[0] for x in structure if x[0] != 'rest']))

    if is_header:
        if m := parser_header()(structure):
            n, parsed = m
            if args.verbose:
                print(json.dumps(parsed, indent=4))
                print(n)
            if n < len(structure):
                if structure[n][0] != 'rest':
                    print(f"{filename}:{n+1}: stopped at", structure[n])
        else:
            print(f"{filename}:1: incomplete parse")

    elif is_compilation_unit:
        if m := parser_compilation_unit()(structure):
            n, parsed = m
            if args.verbose:
                print(json.dumps(parsed, indent=4))
                print(n)
            if n < len(structure):
                if structure[n][0] != 'rest':
                    print(f"{filename}:{n+1}: stopped at", structure[n])
        else:
            print(f"{filename}:1: incomplete parse")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="This script scans the headers of a C source file and"
                    " rearranges them")

    parser.add_argument("file", nargs='*')
    parser.add_argument("-v", "--verbose", action='store_true', default=False)
    parser.add_argument("--all", action='store_true', default=False)

    args = parser.parse_args()

    L = args.file
    if args.all:
        S = subprocess.Popen("git ls-files *.[ch] *.[ch]pp".split(),
                             stdout=subprocess.PIPE)
        L = []
        excluded = r'\b(gf2x|flint-fft|fmt|out-of-core|misc|config)/'
        for b in S.stdout:
            b = b.decode().strip()
            if not re.search(excluded, b):
                L.append(b)

    for f in L:
        scan_file(args, f)
