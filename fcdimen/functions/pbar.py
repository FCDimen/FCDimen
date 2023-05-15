import sys
def progressbar(it, prefix="", size=60, out=sys.stdout):
    count = len(it)
    def show(j):
        x = int(size*j/count)
        y = int((j/count) * 100)
        #print("{}[{}{}] {}/{}".format(prefix, "="*x, "."*(size-x), j, count), end='\r', file=out, flush=True)
        print("{}[{}{}] {}/100".format(prefix, "="*x, "."*(size-x), y), end='\r', file=out, flush=True)
        #print("{}[{}{}]".format(prefix, "="*x, "."*(size-x)), end='\r', file=out, flush=True)
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=out)
