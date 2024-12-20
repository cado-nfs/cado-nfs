#!/usr/bin/env python3

import signal
from cadofactor import patterns

# A signal handler receives as parameters only the signal number and the
# stack frame of the function that was interruped. This means we can't
# register instance methods, as there is no way to pass the instance
# reference. This module defines a module-scope function which is used as
# the signal handler and registers it with signal.signal(). The signal
# handler calls a module- scope Observable instance, to which other
# objects can subscribe, and calls any previously registered signal
# handler. Since modules are loaded at most once, this is kind of a poor
# man's Singleton. The Observable thus acts as a kind of message relay
# that can forward signal notifications to bound methods. The
# notification sent to subscribed Observers has the signal number as the
# message.

signal_usr1_relay = patterns.Observable()


def signal_usr1_handler(signum, frame):
    signal_usr1_relay.notifyObservers(signum)
    if callable(old_signal_usr1_handler):
        old_signal_usr1_handler(signum, frame)


old_signal_usr1_handler = signal.signal(signal.SIGUSR1, signal_usr1_handler)
