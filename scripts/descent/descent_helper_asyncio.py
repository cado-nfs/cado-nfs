import asyncio
import os
import sys


class important_files_async_loop():
    def __init__(self, calls,
                 consumer,
                 consumer_args=(),
                 temporary_is_reusable=False):
        self.consumer = consumer
        self.consumer_args = consumer_args
        self.calls = calls
        self.temporary_is_reusable = temporary_is_reusable
        self.q = None
        self.lk = None
        self.consumer_done = None

    async def producer_task(self, i, outfile, args) -> None:
        if os.path.exists(outfile):
            print(f"reusing file {outfile}")
            with open(outfile, "r") as f:
                for line in f.readlines():
                    async with self.lk:
                        if self.consumer_done.is_set():
                            break
                        await self.q.put((i, line.strip()))
        else:
            if self.temporary_is_reusable:
                outfile_tmp = outfile
            else:
                outfile_tmp = outfile + ".tmp"

            errfile = open(outfile + ".stderr", "wb")
            errfile.write((' '.join(args) + "\n").encode())
            errfile.flush()

            print(f"running program, saving output to {outfile}")
            print("calling:", *args, file=sys.stderr)

            with open(outfile_tmp, 'w') as save:
                proc = await asyncio.create_subprocess_exec(
                    *args,
                    stderr=errfile,
                    stdout=asyncio.subprocess.PIPE)
                while True:
                    line = await proc.stdout.readline()
                    if not line:
                        break
                    line = line.decode().strip()
                    async with self.lk:
                        if self.consumer_done.is_set():
                            proc.kill()
                            break
                        await self.q.put((i, line))
                        print(line, file=save)
                await proc.wait()

            if not self.temporary_is_reusable:
                os.rename(outfile_tmp, outfile)

    async def consumer_task(self) -> None:
        while True:
            i, t = await self.q.get()
            self.q.task_done()
            if self.consumer(*self.consumer_args, i, t):
                break
        async with self.lk:
            self.consumer_done.set()
        # make sure we drain the queue
        while not self.q.empty():
            await self.q.get()
            self.q.task_done()

    async def monitor(self):
        self.q = asyncio.Queue()
        self.consumer_done = asyncio.Event()
        self.lk = asyncio.Lock()
        sources = [asyncio.create_task(self.producer_task(i, *c))
                   for i, c in enumerate(self.calls)]
        sinks = [asyncio.create_task(self.consumer_task())]
        await asyncio.gather(*sources)
        await self.q.join()
        for c in sinks:
            c.cancel()
        return


def monitor_important_files(*args, **kwargs):
    """
    This function calls one or several subprograms (listed in the
    "calls" array) and calls the consumer function on each output line,
    with both the caller-provided consumer_args, the index of the
    subprocess that produced the line, and of course the line itself.
    Processing terminates as soon as the consumer function returns a True
    value, or when all programs exit (thus the consumer function may not
    be called at all).

    Each entry in the "calls" array is a pair (filename, tuple of
    arguments). The output of each call is saved to the filename.
    Depending on whether temporary_is_reusable is True or False, this
    filename is written directly, or via a temporary file.

    The limit argument can be used to limit the number of data lines that
    are produced. limit=0 means unlimited
    """
    F = important_files_async_loop(*args, **kwargs)
    asyncio.run(F.monitor())
