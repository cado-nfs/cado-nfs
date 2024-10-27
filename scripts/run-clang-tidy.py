#!/usr/bin/env python3

import json
import os
import re
import tempfile

def keep_entry(e):
    if re.search(r'/flint-fft/transform_interface\.c$', e['file']):
        return True
    if re.search(r'/(flint-fft|embedding)/', e['file']):
        return False
    if not re.search(r'\.[ch](pp)?$', e['file']):
        return False
    return True



if __name__ == '__main__':
    D = json.load(open('compile_commands.json'))
    D = [ e for e in D if keep_entry(e) ]
    with tempfile.TemporaryDirectory() as tmpdirname:
        json.dump(D, open(f'{tmpdirname}/compile_commands.json', 'w'))
        os.system(f"run-clang-tidy -p {tmpdirname} -j4 | tee clang-tidy.log")


