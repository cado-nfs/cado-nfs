"""This backwards compatible code will eventually go away, but we
still need it for some outdated OSes"""

import io
import os
import subprocess
from queue import Queue, Empty
from threading import Thread


class important_file(object):
    def __init__(self, outfile, call_that, temporary_is_reusable=False):
        self.child = None
        print("command line:\n" + " ".join(call_that))
        self.outfile = outfile
        if temporary_is_reusable:
            self.outfile_tmp = outfile
        else:
            self.outfile_tmp = outfile + ".tmp"
        if os.path.exists(outfile):
            print("reusing file %s" % outfile)
            self.reader = open(outfile, 'r')
            self.writer = None
        else:
            print("running program, saving output to %s" % outfile)
            self.child = subprocess.Popen(call_that, stdout=subprocess.PIPE)
            self.reader = io.TextIOWrapper(self.child.stdout, 'utf-8')
            self.writer = open(self.outfile_tmp, 'w')

    def streams(self):
        return self.reader, self.writer

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.reader)
        if self.writer is not None:
            self.writer.write(line)
            self.writer.flush()
        return line

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self.child is not None:
            # self.writer.close()
            print("Waiting for child to finish")
            # self.child.kill()
            # self.reader.close()   # do I need to put it before ? broken pipe

            for line in self.reader:
                self.writer.write(line)
                self.writer.flush()
            # self.reader.close()
            # self.writer.close()
            # I'm not sure that there's still a point in babysitting the
            # child output as we do above. Maybe a simple communicate()
            # as we do below is all we need. Furthermore, communicate()
            # is the thing that we have to do in order to properly catch
            # the return code.
            self.child.communicate()
            if self.outfile != self.outfile_tmp:
                os.rename(self.outfile_tmp, self.outfile)
            if self.child.returncode != 0:
                raise RuntimeError("Child process failed with return code " +
                                   f" {self.child.returncode} ;" +
                                   " failed command line was:\n" +
                                   " ".join(self.child.args))
            print("ok, done")
        else:
            self.reader.close()


def monitor_important_files(sources,
                            consumer,
                            consumer_args,
                            temporary_is_reusable=False):
    if len(sources) == 1:
        filename, args = sources[0]
        with important_file(filename, args) as f:
            for line in f:
                consumer(*consumer_args, 0, line)
    else:
        q = Queue()

        def enqueue_output(i, out, q):
            for line in out:
                q.put((i, line))
            # q.close()

        threads = [Thread(target=enqueue_output, args=(i, process, q))
                   for (i, process) in enumerate(sources)]

        for t in threads:
            t.daemon = True
            t.start()

        while True:
            try:
                i, line = q.get_nowait()
            except Empty:
                if all([not t.is_alive() for t in threads]):
                    if q.empty():
                        break
                    else:
                        continue
            else:
                consumer(*consumer_args, i, line)
