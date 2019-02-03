cdef int foo1(int a):
    print(a)
    return 2*a

def foo(a):
    return foo1(a)