import abc

"""
This module defines some well-known design patterns, mostly taken from the
Gang-of-Four book.
"""


class Observable(object):
    """
    Implements the subject part of the Observer pattern.

    An Observer can register with the subject, the subject then can notify
    the Observer of state changes by calling the Observer's updateObserver()
    method.
    """
    def __init__(self, *args, **kwargs):
        self.__observers = []
        super().__init__(*args, **kwargs)

    def subscribeObserver(self, observer):
        assert isinstance(observer, Observer)
        if observer not in self.__observers:
            self.__observers.append(observer)

    def unsubscribeObserver(self, observer):
        try:
            self.__observers.remove(observer)
        except ValueError:
            pass

    def notifyObservers(self, message):
        was_received = False
        for observer in self.__observers:
            if observer.updateObserver(message):
                was_received = True
        return was_received


class Observer(object):
    """
    Defines the interface of the client part of the Observer pattern.

    This is an abstract base class; the subclass must implement
    updateObserver.
    """
    @abc.abstractmethod
    def updateObserver(self, message):
        # Subclasses should implement this
        return


class Mediator(object, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def relay_notification(self, notification):
        pass

    @abc.abstractmethod
    def answer_request(self, request):
        pass


class Colleague(object):
    def __init__(self, *args, mediator, **kwargs):
        super().__init__(*args, **kwargs)
        self.__mediator = mediator

    def send_notification(self, notification):
        return self.__mediator.relay_notification(notification)

    def send_request(self, request):
        return self.__mediator.answer_request(request)


class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args,
                                                                 **kwargs)
        return cls._instances[cls]
