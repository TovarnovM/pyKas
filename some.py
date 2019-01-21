class Class4Tst(object):
    def __init__(self, *args, **kwargs):
        self.n = args[0]

#some changes 4 tst git token2

if __name__ == "__main__":
    c = Class4Tst(10)
    print(c.n)