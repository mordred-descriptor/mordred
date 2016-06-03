class MordredException(Exception):
    critical = False


class FragmentError(MordredException, ValueError):
    def __init__(self, cxt, desc):
        self.cxt = cxt
        self.desc = desc

    def __str__(self):
        return '{!r} require just 1 fragmented molecule(molecule: {}).'.format(self.desc, self.cxt)
