class BwcParameters(object):
    NULLSPACE_LEFT = 0
    NULLSPACE_RIGHT = 1

    def __init__(self, *args, **kwargs):
        try:
            self.m = kwargs['m']
            self.n = kwargs['n']
            self.p = kwargs['p']
        except KeyError:
            raise KeyError("Arguments m n p to BwcParameters are required")
        self.wordsize = kwargs.get('wordsize', 64)
        if self.p == 2:
            self.splitwidth = 64
            self.nullspace = self.NULLSPACE_LEFT
        else:
            self.splitwidth = 1
            self.nullspace = self.NULLSPACE_RIGHT

        word = 2**self.wordsize
        self.p_words = len((self.p**self.splitwidth-1).digits(word))
        self.p_bytes = self.p_words * (self.wordsize // 8)
        n = kwargs.get('nullspace', None)
        if n == 'left':
            self.nullspace = self.NULLSPACE_LEFT
        elif n == 'right':
            self.nullspace = self.NULLSPACE_RIGHT
        else:
            # default behavior is to let the prime decide.
            pass

    def is_nullspace_left(self):
        return self.nullspace == self.NULLSPACE_LEFT

    def is_nullspace_right(self):
        # of course it's just for convenience matters that we have both.
        return self.nullspace == self.NULLSPACE_RIGHT
